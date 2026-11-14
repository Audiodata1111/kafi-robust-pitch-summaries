# =============================================================================
# MONTE CARLO ROBUSTNESS ASSESSMENT FOR KAFI PITCH PARAMETERS
# =============================================================================
#
# Authors : Prena Singh, Soubhik Chakraborty and Pushkal Kumar
#           Birla Institute of Technology Mesra, Ranchi, India
# Licence : MIT (see LICENSE)
# Part of : Robust Summaries of Tonic-Normalised Pitch Events, with an
#           Application to Performance Variation in a Raga Kafi Bandish
#
# Purpose
# -------
# Evaluate how strongly conventional and robust P1--P3 change when a controlled
# proportion of extracted note events is given a simulated octave error.
#
# Scientific scope
# ----------------
# This study validates the stability of the event-to-parameter stage. It does
# not validate the absolute F0 accuracy of pYIN, which would require manually
# annotated or otherwise known ground-truth F0 contours.
#
# Each measure is compared with its own uncontaminated value. Conventional and
# robust parameters are not treated as estimators of one identical target.
#
# Input
# -----
# The script searches for:
#   output/all_note_events.csv
# and then:
#   all_note_events.csv
#
# Required input columns:
#   recording_id, Event_ID, Duration, Pitch_Semitones
#
# Primary simulation design
# -------------------------
# * All available recordings are used, including R.
# * Contamination rates: 1%, 2.5%, 5%, and 10% of note events.
# * At least one event is contaminated at every positive rate.
# * The number contaminated is ceiling(rate * number of events).
# * Selected events are sampled uniformly without replacement.
# * Optional duration-weighted design: EVENT_SELECTION <- "duration" draws
#   events with probability proportional to duration until their combined
#   share of total event duration reaches the contamination rate. This tests
#   the duration-share breakdown of robust P1 and P2 (Proposition 2). Its
#   outputs and checkpoints go to separate *_duration folders.
# * Each selected pitch receives either -12 or +12 semitones with equal
#   probability; event order and event durations remain unchanged.
# * Monte Carlo repetitions: 1000 per rate and recording.
# * Recording-specific checkpoints make interrupted runs resumable.
#
# Required packages
# -----------------
# install.packages(c("data.table", "ggplot2"))
# =============================================================================


# -----------------------------------------------------------------------------
# 1. SETTINGS
# -----------------------------------------------------------------------------

SIMULATION_VERSION <- "1.0"
SIMULATION_SEED <- 20260829L
N_REPLICATES <- 1000L
CONTAMINATION_RATES <- c(0.01, 0.025, 0.05, 0.10)
OCTAVE_SHIFT_SEMITONES <- 12
RELATIVE_CHANGE_WARNING_PERCENT <- 10

ROBUST_PITCH_LOWER_PROB <- 0.05
ROBUST_PITCH_UPPER_PROB <- 0.95
ROBUST_TRANSITION_RETAINED_PROB <- 0.90

# All recordings, including R, are used for methodological robustness.
INCLUDE_REFERENCE <- TRUE
REFERENCE_ID <- "R"

# Reuse completed recording-level simulations unless inputs or settings change.
FORCE_REPROCESS <- FALSE
REPROCESS_IDS <- character()

INPUT_CANDIDATES <- c(
  file.path("output", "all_note_events.csv"),
  "all_note_events.csv"
)

# "uniform" runs the primary simulation reported in the article;
# "duration" runs the duration-weighted design.
EVENT_SELECTION <- "uniform"

if (!(EVENT_SELECTION %in% c("uniform", "duration"))) {
  stop("EVENT_SELECTION must be 'uniform' or 'duration'.", call. = FALSE)
}

# The duration-weighted design writes to its own folders.
RUN_SUFFIX <- if (EVENT_SELECTION == "duration") "_duration" else ""
OUTPUT_DIR <- file.path("output", paste0("robustness_simulation", RUN_SUFFIX))
CACHE_DIR <- file.path("kafi_cache", paste0("robustness_simulation", RUN_SUFFIX))


# -----------------------------------------------------------------------------
# 2. PACKAGE AND INPUT CHECKS
# -----------------------------------------------------------------------------

required_packages <- c("data.table", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    "\nRun: install.packages(c(",
    paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

if (N_REPLICATES < 1L) {
  stop("N_REPLICATES must be at least 1.", call. = FALSE)
}

if (any(!is.finite(CONTAMINATION_RATES)) ||
    any(CONTAMINATION_RATES <= 0 | CONTAMINATION_RATES >= 1)) {
  stop("Every contamination rate must lie strictly between 0 and 1.",
       call. = FALSE)
}

if (!is.finite(OCTAVE_SHIFT_SEMITONES) || OCTAVE_SHIFT_SEMITONES <= 0) {
  stop("OCTAVE_SHIFT_SEMITONES must be positive.", call. = FALSE)
}

existing_inputs <- INPUT_CANDIDATES[file.exists(INPUT_CANDIDATES)]

if (length(existing_inputs) == 0L) {
  stop(
    "Could not find all_note_events.csv. Place it either in the project ",
    "directory or in the output directory.",
    call. = FALSE
  )
}

INPUT_FILE <- existing_inputs[1L]
INPUT_MD5 <- unname(tools::md5sum(INPUT_FILE))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

message("Reading note events from: ", normalizePath(INPUT_FILE))


# -----------------------------------------------------------------------------
# 3. SMALL UTILITY FUNCTIONS
# -----------------------------------------------------------------------------

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]", "_", as.character(x))
}

save_rds_safely <- function(object, path) {
  temporary_path <- paste0(path, ".tmp")
  saveRDS(object, temporary_path, compress = TRUE)

  copied <- file.copy(temporary_path, path, overwrite = TRUE)
  unlink(temporary_path)

  if (!isTRUE(copied) || !file.exists(path)) {
    stop("Could not finalize checkpoint: ", path, call. = FALSE)
  }
}

weighted_quantile <- function(x, weights, probabilities) {
  keep <- is.finite(x) & is.finite(weights) & weights > 0
  x <- x[keep]
  weights <- weights[keep]

  if (length(x) == 0L) {
    return(rep(NA_real_, length(probabilities)))
  }

  if (length(x) == 1L) {
    return(rep(x, length(probabilities)))
  }

  ordering <- order(x)
  x <- x[ordering]
  weights <- weights[ordering]

  midpoint_probability <- (
    cumsum(weights) - 0.5 * weights
  ) / sum(weights)

  stats::approx(
    x = midpoint_probability,
    y = x,
    xout = probabilities,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )$y
}

calculate_pitch_parameters <- function(pitch, duration) {
  keep <- is.finite(pitch) & is.finite(duration) & duration > 0
  pitch <- pitch[keep]
  duration <- duration[keep]

  if (length(pitch) < 2L) {
    stop("At least two valid events are required for P3.", call. = FALSE)
  }

  conventional_p1 <- (
    sum(pitch * duration) / sum(duration)
  ) - min(pitch)
  conventional_p2 <- max(pitch) - min(pitch)

  transitions <- abs(diff(pitch))
  conventional_p3 <- mean(transitions)

  pitch_quantiles <- weighted_quantile(
    pitch,
    duration,
    c(ROBUST_PITCH_LOWER_PROB, ROBUST_PITCH_UPPER_PROB)
  )

  pitch_q05 <- pitch_quantiles[1L]
  pitch_q95 <- pitch_quantiles[2L]
  winsorised_pitch <- pmin(pmax(pitch, pitch_q05), pitch_q95)

  robust_p1 <- (
    sum(winsorised_pitch * duration) / sum(duration)
  ) - pitch_q05
  robust_p2 <- pitch_q95 - pitch_q05

  retained_transition_count <- max(
    1L,
    as.integer(ceiling(
      ROBUST_TRANSITION_RETAINED_PROB * length(transitions)
    ))
  )

  robust_p3 <- mean(
    sort(transitions)[seq_len(retained_transition_count)]
  )

  c(
    P1_Conventional = conventional_p1,
    P1_Robust = robust_p1,
    P2_Conventional = conventional_p2,
    P2_Robust = robust_p2,
    P3_Conventional = conventional_p3,
    P3_Robust = robust_p3
  )
}

metric_information <- data.table::data.table(
  Metric = c(
    "P1_Conventional", "P1_Robust",
    "P2_Conventional", "P2_Robust",
    "P3_Conventional", "P3_Robust"
  ),
  Parameter = rep(c("P1", "P2", "P3"), each = 2L),
  Method = rep(c("Conventional", "Robust"), 3L)
)

checkpoint_path <- function(recording_id) {
  file.path(
    CACHE_DIR,
    paste0(safe_filename(recording_id), "_robustness_simulation.rds")
  )
}

make_signature <- function(recording_id, n_events) {
  recording_seed <- as.integer(
    SIMULATION_SEED +
      sum(utf8ToInt(as.character(recording_id)) *
            seq_along(utf8ToInt(as.character(recording_id))))
  )

  signature <- list(
    simulation_version = SIMULATION_VERSION,
    input_md5 = INPUT_MD5,
    recording_id = as.character(recording_id),
    n_events = as.integer(n_events),
    global_seed = as.integer(SIMULATION_SEED),
    recording_seed = recording_seed,
    include_reference = as.logical(INCLUDE_REFERENCE),
    n_replicates = as.integer(N_REPLICATES),
    contamination_rates = as.numeric(CONTAMINATION_RATES),
    octave_shift_semitones = as.numeric(OCTAVE_SHIFT_SEMITONES),
    robust_pitch_lower_probability = as.numeric(ROBUST_PITCH_LOWER_PROB),
    robust_pitch_upper_probability = as.numeric(ROBUST_PITCH_UPPER_PROB),
    robust_transition_retained_probability = as.numeric(
      ROBUST_TRANSITION_RETAINED_PROB
    )
  )

  # The selection mode enters the signature only for the duration-weighted
  # design, so checkpoints of the uniform design remain valid.
  if (EVENT_SELECTION != "uniform") {
    signature$event_selection <- EVENT_SELECTION
  }

  signature
}

# Draw events with probability proportional to duration, without
# replacement, until their combined duration share reaches the target rate.
select_events_by_duration <- function(duration, target_share) {
  draw_order <- sample.int(length(duration), prob = duration)
  cumulative_share <- cumsum(duration[draw_order]) / sum(duration)
  draw_order[seq_len(which(cumulative_share >= target_share)[1L])]
}

read_valid_checkpoint <- function(recording_id, n_events) {
  path <- checkpoint_path(recording_id)

  if (!file.exists(path)) {
    return(NULL)
  }

  checkpoint <- tryCatch(readRDS(path), error = function(e) NULL)

  if (is.null(checkpoint) || !is.list(checkpoint)) {
    return(NULL)
  }

  if (!all(c("signature", "baseline", "replicates") %in% names(checkpoint))) {
    return(NULL)
  }

  if (!identical(
    checkpoint$signature,
    make_signature(recording_id, n_events)
  )) {
    return(NULL)
  }

  checkpoint
}


# -----------------------------------------------------------------------------
# 4. READ AND VALIDATE THE EVENT DATA
# -----------------------------------------------------------------------------

events <- data.table::fread(INPUT_FILE)

required_columns <- c(
  "recording_id", "Event_ID", "Duration", "Pitch_Semitones"
)
missing_columns <- setdiff(required_columns, names(events))

if (length(missing_columns) > 0L) {
  stop(
    "The event file is missing: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

events[, recording_id := as.character(recording_id)]
events[, Event_ID := as.integer(Event_ID)]
events[, Duration := as.numeric(Duration)]
events[, Pitch_Semitones := as.numeric(Pitch_Semitones)]

if (anyNA(events[, .(
  recording_id, Event_ID, Duration, Pitch_Semitones
)])) {
  stop("Required event columns contain missing values.", call. = FALSE)
}

if (any(!is.finite(events$Duration) | events$Duration <= 0)) {
  stop("All event durations must be finite and positive.", call. = FALSE)
}

if (any(!is.finite(events$Pitch_Semitones))) {
  stop("All event pitches must be finite.", call. = FALSE)
}

if (events[, anyDuplicated(paste(recording_id, Event_ID))] > 0L) {
  stop("Event_ID values must be unique within each recording.", call. = FALSE)
}

if (!INCLUDE_REFERENCE) {
  events <- events[recording_id != REFERENCE_ID]
}

if (nrow(events) == 0L) {
  stop("No events remain after applying the inclusion settings.", call. = FALSE)
}

event_counts <- events[, .(N_Events = .N), by = recording_id]

if (any(event_counts$N_Events < 2L)) {
  stop(
    "Every included recording must contain at least two note events.",
    call. = FALSE
  )
}

participant_ids <- sort(setdiff(unique(events$recording_id), REFERENCE_ID))
recording_ids <- c(
  if (REFERENCE_ID %in% events$recording_id) REFERENCE_ID else character(),
  participant_ids
)

data.table::setorder(events, recording_id, Event_ID)


# -----------------------------------------------------------------------------
# 5. RECORDING-LEVEL MONTE CARLO SIMULATION
# -----------------------------------------------------------------------------

status <- data.table::data.table(
  recording_id = recording_ids,
  status = "pending",
  message = "",
  checkpoint = vapply(recording_ids, checkpoint_path, character(1))
)

status_path <- file.path(CACHE_DIR, "simulation_status.csv")
data.table::fwrite(status, status_path)

for (recording_index in seq_along(recording_ids)) {
  current_id <- recording_ids[recording_index]
  current_events <- events[recording_id == current_id]
  data.table::setorder(current_events, Event_ID)

  n_events <- nrow(current_events)
  pitch <- current_events$Pitch_Semitones
  duration <- current_events$Duration

  use_cache <- !FORCE_REPROCESS &&
    !(current_id %in% REPROCESS_IDS)

  checkpoint <- if (use_cache) {
    read_valid_checkpoint(current_id, n_events)
  } else {
    NULL
  }

  if (!is.null(checkpoint)) {
    message("Using cached simulation for ", current_id, " ...")
    status[recording_id == current_id, `:=`(
      status = "cached",
      message = "Valid recording-level simulation checkpoint reused."
    )]
    data.table::fwrite(status, status_path)
    next
  }

  message("Simulating octave contamination for ", current_id, " ...")
  status[recording_id == current_id, `:=`(
    status = "processing",
    message = "Monte Carlo simulation in progress."
  )]
  data.table::fwrite(status, status_path)

  simulation_result <- tryCatch(
    {
      baseline <- calculate_pitch_parameters(pitch, duration)

      # The recording-specific seed makes results reproducible even when a
      # previous run is resumed from checkpoints.
      set.seed(make_signature(current_id, n_events)$recording_seed)

      rate_results <- vector("list", length(CONTAMINATION_RATES))

      for (rate_index in seq_along(CONTAMINATION_RATES)) {
        contamination_rate <- CONTAMINATION_RATES[rate_index]
        contaminated_count <- max(
          1L,
          as.integer(ceiling(contamination_rate * n_events))
        )

        result_matrix <- matrix(
          NA_real_,
          nrow = N_REPLICATES,
          ncol = length(baseline),
          dimnames = list(NULL, names(baseline))
        )

        # Record, per replicate, how many events were shifted and what share
        # of total duration they carried.
        selected_count <- integer(N_REPLICATES)
        selected_duration_share <- numeric(N_REPLICATES)

        for (replicate_id in seq_len(N_REPLICATES)) {
          contaminated_pitch <- pitch
          selected_events <- if (EVENT_SELECTION == "duration") {
            select_events_by_duration(duration, contamination_rate)
          } else {
            sample.int(
              n_events,
              size = contaminated_count,
              replace = FALSE
            )
          }
          selected_count[replicate_id] <- length(selected_events)
          selected_duration_share[replicate_id] <-
            sum(duration[selected_events]) / sum(duration)

          shift_direction <- sample(
            c(-1, 1),
            size = length(selected_events),
            replace = TRUE
          )

          contaminated_pitch[selected_events] <-
            contaminated_pitch[selected_events] +
            shift_direction * OCTAVE_SHIFT_SEMITONES

          result_matrix[replicate_id, ] <- calculate_pitch_parameters(
            contaminated_pitch,
            duration
          )
        }

        rate_results[[rate_index]] <- data.table::as.data.table(result_matrix)
        # Event counts and duration shares are stored per replicate; in the
        # uniform design the count is constant.
        rate_results[[rate_index]][, `:=`(
          recording_id = current_id,
          contamination_rate = contamination_rate,
          actual_contamination_rate = selected_count / n_events,
          contaminated_events = selected_count,
          duration_contamination_share = selected_duration_share,
          replicate = seq_len(N_REPLICATES)
        )]

        data.table::setcolorder(
          rate_results[[rate_index]],
          c(
            "recording_id", "contamination_rate",
            "actual_contamination_rate", "contaminated_events",
            "duration_contamination_share",
            "replicate", names(baseline)
          )
        )
      }

      list(
        signature = make_signature(current_id, n_events),
        baseline = baseline,
        replicates = data.table::rbindlist(rate_results, use.names = TRUE),
        completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
      )
    },
    error = function(e) e
  )

  if (inherits(simulation_result, "error")) {
    status[recording_id == current_id, `:=`(
      status = "failed",
      message = conditionMessage(simulation_result)
    )]
    data.table::fwrite(status, status_path)
    stop(
      "Robustness simulation failed for ", current_id, ": ",
      conditionMessage(simulation_result),
      "\nCompleted recording checkpoints were preserved.",
      call. = FALSE
    )
  }

  save_rds_safely(
    simulation_result,
    checkpoint_path(current_id)
  )

  status[recording_id == current_id, `:=`(
    status = "processed",
    message = "Simulation checkpoint created successfully."
  )]
  data.table::fwrite(status, status_path)
}


# -----------------------------------------------------------------------------
# 6. COMBINE CHECKPOINTS AND CALCULATE CHANGES FROM CLEAN BASELINES
# -----------------------------------------------------------------------------

checkpoints <- lapply(recording_ids, function(current_id) {
  current_n <- event_counts[recording_id == current_id, N_Events]
  checkpoint <- read_valid_checkpoint(current_id, current_n)

  if (is.null(checkpoint)) {
    stop("A valid final checkpoint is missing for ", current_id,
         call. = FALSE)
  }

  checkpoint
})

baseline_table <- data.table::rbindlist(
  lapply(seq_along(recording_ids), function(index) {
    values <- checkpoints[[index]]$baseline
    data.table::data.table(
      recording_id = recording_ids[index],
      Metric = names(values),
      Baseline_Value = as.numeric(values)
    )
  })
)

baseline_table <- merge(
  baseline_table,
  metric_information,
  by = "Metric",
  all.x = TRUE,
  sort = FALSE
)

replicate_wide <- data.table::rbindlist(
  lapply(checkpoints, function(x) x$replicates),
  use.names = TRUE,
  fill = TRUE
)

# Checkpoints from earlier versions of this script lack the duration share.
if (!("duration_contamination_share" %in% names(replicate_wide))) {
  replicate_wide[, duration_contamination_share := NA_real_]
}

replicate_long <- data.table::melt(
  replicate_wide,
  id.vars = c(
    "recording_id", "contamination_rate",
    "actual_contamination_rate", "contaminated_events",
    "duration_contamination_share", "replicate"
  ),
  measure.vars = metric_information$Metric,
  variable.name = "Metric",
  value.name = "Contaminated_Value"
)

replicate_long <- merge(
  replicate_long,
  baseline_table,
  by = c("recording_id", "Metric"),
  all.x = TRUE,
  sort = FALSE
)

replicate_long <- merge(
  replicate_long,
  event_counts,
  by = "recording_id",
  all.x = TRUE,
  sort = FALSE
)

replicate_long[, Signed_Change := Contaminated_Value - Baseline_Value]
replicate_long[, Absolute_Change := abs(Signed_Change)]
replicate_long[, Relative_Change_Percent := 100 * Absolute_Change /
                 pmax(abs(Baseline_Value), sqrt(.Machine$double.eps))]


# -----------------------------------------------------------------------------
# 7. RECORDING-LEVEL AND OVERALL ROBUSTNESS SUMMARIES
# -----------------------------------------------------------------------------

recording_summary <- replicate_long[, .(
  N_Replicates = .N,
  N_Events = unique(N_Events),
  # Medians, because the duration-weighted design varies these per
  # replicate; in the uniform design they are constant.
  Actual_Contamination_Rate = stats::median(actual_contamination_rate),
  Contaminated_Events = stats::median(contaminated_events),
  Median_Duration_Contamination_Share = stats::median(
    duration_contamination_share
  ),
  Baseline_Value = unique(Baseline_Value),
  Median_Signed_Change = stats::median(Signed_Change),
  Median_Absolute_Change = stats::median(Absolute_Change),
  IQR_Absolute_Change = stats::IQR(Absolute_Change),
  Q95_Absolute_Change = unname(stats::quantile(
    Absolute_Change,
    probs = 0.95,
    names = FALSE,
    type = 7
  )),
  RMSE_Change = sqrt(mean(Signed_Change^2)),
  Median_Relative_Change_Percent = stats::median(
    Relative_Change_Percent
  ),
  IQR_Relative_Change_Percent = stats::IQR(Relative_Change_Percent),
  Proportion_Relative_Change_GT_10pct = mean(
    Relative_Change_Percent > RELATIVE_CHANGE_WARNING_PERCENT
  )
), by = .(
  recording_id,
  contamination_rate,
  Metric,
  Parameter,
  Method
)]

overall_summary <- recording_summary[, .(
  N_Recordings = .N,
  Median_Recording_Median_Absolute_Change = stats::median(
    Median_Absolute_Change
  ),
  Q1_Recording_Median_Absolute_Change = unname(stats::quantile(
    Median_Absolute_Change, 0.25, names = FALSE
  )),
  Q3_Recording_Median_Absolute_Change = unname(stats::quantile(
    Median_Absolute_Change, 0.75, names = FALSE
  )),
  Median_Recording_RMSE_Change = stats::median(RMSE_Change),
  Median_Recording_Median_Relative_Change_Percent = stats::median(
    Median_Relative_Change_Percent
  ),
  Q1_Recording_Median_Relative_Change_Percent = unname(stats::quantile(
    Median_Relative_Change_Percent, 0.25, names = FALSE
  )),
  Q3_Recording_Median_Relative_Change_Percent = unname(stats::quantile(
    Median_Relative_Change_Percent, 0.75, names = FALSE
  )),
  Mean_Proportion_Relative_Change_GT_10pct = mean(
    Proportion_Relative_Change_GT_10pct
  )
), by = .(
  contamination_rate,
  Metric,
  Parameter,
  Method
)]

# Structural validation: one summary row is required for every
# recording x contamination rate x metric combination.
expected_recording_rows <- length(recording_ids) *
  length(CONTAMINATION_RATES) * nrow(metric_information)

if (nrow(recording_summary) != expected_recording_rows) {
  stop(
    paste0(
      "Recording-level robustness summary has ",
      nrow(recording_summary), " rows; expected ",
      expected_recording_rows, "."
    ),
    call. = FALSE
  )
}

if (recording_summary[, anyDuplicated(
  paste(recording_id, contamination_rate, Metric, sep = "|")
)] > 0L) {
  stop(
    "Recording-level robustness summary contains duplicate keys.",
    call. = FALSE
  )
}

expected_overall_rows <- length(CONTAMINATION_RATES) *
  nrow(metric_information)

if (nrow(overall_summary) != expected_overall_rows ||
    any(overall_summary$N_Recordings != length(recording_ids))) {
  stop(
    paste0(
      "Overall robustness summary failed structural validation; ",
      "each row must represent all ", length(recording_ids),
      " recordings."
    ),
    call. = FALSE
  )
}

data.table::setorder(
  recording_summary,
  Parameter,
  Method,
  contamination_rate,
  recording_id
)
data.table::setorder(
  overall_summary,
  Parameter,
  Method,
  contamination_rate
)

gain_wide <- data.table::dcast(
  overall_summary,
  Parameter + contamination_rate ~ Method,
  value.var = "Median_Recording_Median_Relative_Change_Percent"
)

gain_wide[, Robustness_Gain_Ratio := data.table::fifelse(
  is.finite(Robust) & Robust > 0,
  Conventional / Robust,
  NA_real_
)]
data.table::setorder(gain_wide, Parameter, contamination_rate)


# -----------------------------------------------------------------------------
# 8. WRITE NUMERICAL OUTPUTS
# -----------------------------------------------------------------------------

data.table::fwrite(
  baseline_table,
  file.path(OUTPUT_DIR, "baseline_pitch_parameters.csv")
)

data.table::fwrite(
  recording_summary,
  file.path(OUTPUT_DIR, "recording_level_robustness_summary.csv")
)

data.table::fwrite(
  overall_summary,
  file.path(OUTPUT_DIR, "overall_robustness_summary.csv")
)

data.table::fwrite(
  gain_wide,
  file.path(OUTPUT_DIR, "robustness_gain.csv")
)

saveRDS(
  replicate_wide,
  file.path(OUTPUT_DIR, "simulation_replicate_results.rds"),
  compress = TRUE
)

simulation_settings <- data.table::data.table(
  Setting = c(
    "Simulation version",
    "Input file",
    "Input MD5",
    "Number of recordings",
    "Include reference",
    "Monte Carlo repetitions per rate and recording",
    "Contamination rates",
    "Pitch shift (semitones)",
    "Event selection",
    "Shift direction",
    "Robust pitch lower quantile",
    "Robust pitch upper quantile",
    "Robust transition retained proportion",
    "Relative-change warning threshold (percent)",
    "Random seed"
  ),
  Value = c(
    SIMULATION_VERSION,
    INPUT_FILE,
    INPUT_MD5,
    as.character(length(recording_ids)),
    as.character(INCLUDE_REFERENCE),
    as.character(N_REPLICATES),
    paste(CONTAMINATION_RATES, collapse = ", "),
    as.character(OCTAVE_SHIFT_SEMITONES),
    # Report the selection mode actually used.
    if (EVENT_SELECTION == "duration") {
      "Probability proportional to duration, without replacement"
    } else {
      "Uniform without replacement"
    },
    paste0(
      "-", OCTAVE_SHIFT_SEMITONES, " or +", OCTAVE_SHIFT_SEMITONES,
      " semitones with equal probability"
    ),
    as.character(ROBUST_PITCH_LOWER_PROB),
    as.character(ROBUST_PITCH_UPPER_PROB),
    as.character(ROBUST_TRANSITION_RETAINED_PROB),
    as.character(RELATIVE_CHANGE_WARNING_PERCENT),
    as.character(SIMULATION_SEED)
  )
)

data.table::fwrite(
  simulation_settings,
  file.path(OUTPUT_DIR, "simulation_settings.csv")
)

writeLines(
  capture.output(utils::sessionInfo()),
  file.path(OUTPUT_DIR, "session_info.txt")
)


# -----------------------------------------------------------------------------
# 9. PUBLICATION FIGURE
# -----------------------------------------------------------------------------

zero_rows <- unique(
  overall_summary[, .(Parameter, Method)]
)
zero_rows[, `:=`(
  contamination_rate = 0,
  Median_Recording_Median_Relative_Change_Percent = 0,
  Q1_Recording_Median_Relative_Change_Percent = 0,
  Q3_Recording_Median_Relative_Change_Percent = 0
)]

plot_data <- data.table::rbindlist(
  list(
    zero_rows,
    overall_summary[, .(
      Parameter,
      Method,
      contamination_rate,
      Median_Recording_Median_Relative_Change_Percent,
      Q1_Recording_Median_Relative_Change_Percent,
      Q3_Recording_Median_Relative_Change_Percent
    )]
  ),
  use.names = TRUE,
  fill = TRUE
)

plot_data[, Contamination_Percent := 100 * contamination_rate]
plot_data[, Parameter := factor(Parameter, levels = c("P1", "P2", "P3"))]
plot_data[, Method := factor(
  Method,
  levels = c("Conventional", "Robust")
)]

robustness_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = Contamination_Percent,
    y = Median_Recording_Median_Relative_Change_Percent,
    colour = Method,
    fill = Method,
    group = Method
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = Q1_Recording_Median_Relative_Change_Percent,
      ymax = Q3_Recording_Median_Relative_Change_Percent
    ),
    alpha = 0.14,
    colour = NA
  ) +
  ggplot2::geom_line(linewidth = 1.05) +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::facet_wrap(~Parameter, nrow = 1) +
  ggplot2::scale_colour_manual(
    values = c("Conventional" = "#C0392B", "Robust" = "#2166AC")
  ) +
  ggplot2::scale_fill_manual(
    values = c("Conventional" = "#C0392B", "Robust" = "#2166AC")
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(0, 1, 2.5, 5, 10),
    labels = c("0", "1", "2.5", "5", "10")
  ) +
  ggplot2::labs(
    title = "Robustness to simulated octave-tracking errors",
    subtitle = paste0(
      "Lines: median of recording-specific Monte Carlo summaries; ",
      "bands: interquartile range"
    ),
    # In the duration-weighted design the axis is a share of duration.
    x = if (EVENT_SELECTION == "duration") {
      "Contaminated share of event duration (%)"
    } else {
      "Contaminated note events (%)"
    },
    y = "Median relative change from baseline (%)",
    colour = "Parameter form",
    fill = "Parameter form"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(
      fill = "#E9EEF5",
      colour = "#9AA7B2",
      linewidth = 0.5
    ),
    strip.text = ggplot2::element_text(face = "bold", size = 12),
    axis.title = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold", size = 15),
    plot.subtitle = ggplot2::element_text(size = 10.5, colour = "#3F4C5A"),
    legend.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "robustness_simulation_relative_change.png"
  ),
  plot = robustness_plot,
  width = 12,
  height = 5.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "robustness_simulation_relative_change.pdf"
  ),
  plot = robustness_plot,
  width = 12,
  height = 5.2,
  units = "in",
  device = "pdf",
  bg = "white"
)


# -----------------------------------------------------------------------------
# 10. COMPLETION MESSAGE
# -----------------------------------------------------------------------------

status[status == "processed" | status == "cached", status := "complete"]
data.table::fwrite(status, status_path)

message("Robustness simulation completed successfully.")
message("Recordings included: ", length(recording_ids))
message("Replicates per rate and recording: ", N_REPLICATES)
message("Contamination rates: ", paste(100 * CONTAMINATION_RATES, collapse = ", "), "%")
message("Results are available in: ", normalizePath(OUTPUT_DIR))
