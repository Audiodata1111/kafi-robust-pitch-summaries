# =============================================================================
# KAFI INFERENTIAL ANALYSES
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
# Inferential analyses built on the outputs of Kafi_Swara_Analysis.R. Run
# that script first, then run this one from the same project directory.
#
#   Article item                        Output file
#   ----------------------------------  ---------------------------------------
#   Table S1 (lengths, sampling rates)  recording_metadata.csv
#   Section 4.3, Table S15              centrality_typicality.csv
#   Table 3, Table S16                  swara_tests.csv
#   Table S17                           confound_spearman.csv
#   Table S18                           breakdown_diagnostics.csv
#   Section 4.4, Table 4                markov_homogeneity.csv
#   Figure 2                            swara_share_vs_reference.png / .pdf
#   Optional source-group check         source_group_tests.csv,
#                                       swara_tests_direct_only.csv
#
# Inputs (searched in output/ first, then the project directory)
# --------------------------------------------------------------
#   swara_distributions.csv  recording_id, Note, Count_Probability,
#                            Duration_Probability
#   all_note_events.csv      recording_id, Event_ID, Duration
#   all_transitions.csv      recording_id, From, To, Count
#   P1_to_P5.csv             recording_id, n_events, P2_melodic_span,
#                            P5_mean_note_duration
#   quality_summary.csv      recording_id, total_frames, raw_voiced_fraction,
#                            frame_duration_sec
#   recordings.csv           recording_id, file_name, sa_hz   (project dir)
#   kafi_cache/f0/*_pyin_f0.rds  optional; supplies each file's sampling rate
#                            (otherwise read from the WAV headers in audio/)
#   recording_sources.csv    optional; recording_id, source, where source is
#                            original, direct or demucs. If absent, a template
#                            is written and the source-group check is skipped.
#
# Notes
# -----
# * The reference R is a fixed benchmark in the per-swara tests; the 33
#   participant recordings are the observational units.
# * All permutation tests use a fixed seed (ANALYSIS_SEED) for reproducibility.
#
# Required packages
# -----------------
# install.packages(c("data.table", "ggplot2"))
# =============================================================================


# -----------------------------------------------------------------------------
# 1. SETTINGS
# -----------------------------------------------------------------------------

REFERENCE_ID <- "R"
SWARA_ORDER <- c("Sa", "Re", "g", "Ma", "Pa", "Dha", "n")

N_PERMUTATIONS <- 2000L
ANALYSIS_SEED <- 20261022L

# Must match the settings used for the robust parameters and the simulation.
ROBUST_TRANSITION_RETAINED_PROB <- 0.90
CONTAMINATION_RATES <- c(0.01, 0.025, 0.05, 0.10)

ALPHA <- 0.05

CONFIG_FILE <- "recordings.csv"
SOURCES_FILE <- "recording_sources.csv"
SOURCE_LEVELS <- c("original", "direct", "demucs")
F0_CACHE_DIR <- file.path("kafi_cache", "f0")
OUTPUT_DIR <- file.path("output", "inference_analyses")


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

find_input <- function(file_name) {
  candidates <- c(file.path("output", file_name), file_name)
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop("Could not find ", file_name, " in output/ or the project directory.",
         call. = FALSE)
  }
  existing[1L]
}

check_columns <- function(data, columns, label) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(label, " is missing: ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

distributions <- data.table::fread(find_input("swara_distributions.csv"))
events <- data.table::fread(find_input("all_note_events.csv"))
transitions_raw <- data.table::fread(find_input("all_transitions.csv"))
parameters <- data.table::fread(find_input("P1_to_P5.csv"))
quality <- data.table::fread(find_input("quality_summary.csv"))

check_columns(distributions, c("recording_id", "Note", "Count_Probability",
                               "Duration_Probability"), "swara_distributions.csv")
check_columns(events, c("recording_id", "Event_ID", "Duration"),
              "all_note_events.csv")
check_columns(transitions_raw, c("recording_id", "From", "To", "Count"),
              "all_transitions.csv")
check_columns(parameters, c("recording_id", "n_events", "P2_melodic_span",
                            "P5_mean_note_duration"), "P1_to_P5.csv")
check_columns(quality, c("recording_id", "total_frames", "raw_voiced_fraction",
                         "frame_duration_sec"), "quality_summary.csv")

distributions[, recording_id := as.character(recording_id)]
events[, recording_id := as.character(recording_id)]
transitions_raw[, recording_id := as.character(recording_id)]
parameters[, recording_id := as.character(recording_id)]
quality[, recording_id := as.character(recording_id)]

recording_ids <- c(
  REFERENCE_ID,
  sort(setdiff(unique(distributions$recording_id), REFERENCE_ID))
)
participant_ids <- setdiff(recording_ids, REFERENCE_ID)

if (!(REFERENCE_ID %in% distributions$recording_id)) {
  stop("Reference recording R is missing from swara_distributions.csv.",
       call. = FALSE)
}

message("Recordings found: ", length(recording_ids),
        " (", length(participant_ids), " participants)")


# -----------------------------------------------------------------------------
# 3. SMALL UTILITY FUNCTIONS
# -----------------------------------------------------------------------------

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_-]", "_", x)
}

bhattacharyya_distance <- function(p, q) {
  coefficient <- sum(sqrt(p * q))
  coefficient <- min(max(coefficient, 0), 1)
  if (coefficient <= 0) return(Inf)
  -log(coefficient)
}

# Recording x swara probability matrix, rows in recording_ids order.
probability_matrix <- function(column) {
  wide <- data.table::dcast(
    distributions[Note %in% SWARA_ORDER],
    recording_id ~ Note,
    value.var = column,
    fill = 0
  )
  for (s in setdiff(SWARA_ORDER, names(wide))) wide[, (s) := 0]
  m <- as.matrix(wide[, SWARA_ORDER, with = FALSE])
  rownames(m) <- wide$recording_id
  m <- m[recording_ids, , drop = FALSE]
  m / rowSums(m)
}

# Matched-pairs rank-biserial correlation (Kerby 2014).
rank_biserial <- function(x) {
  x <- x[x != 0]
  if (length(x) == 0L) return(NA_real_)
  r <- rank(abs(x))
  (sum(r[x > 0]) - sum(r[x < 0])) / sum(r)
}


# -----------------------------------------------------------------------------
# 4. RECORDING METADATA FOR TABLE S1
# -----------------------------------------------------------------------------
# Recording length is total frames x frame step. The sampling rate is stored
# in each raw pYIN checkpoint by Kafi_Swara_Analysis.R.

metadata <- quality[recording_id %in% recording_ids, .(
  recording_id,
  total_frames,
  frame_duration_sec,
  length_sec = total_frames * frame_duration_sec
)]

metadata[, sample_rate_hz := vapply(recording_id, function(id) {
  path <- file.path(F0_CACHE_DIR, paste0(safe_filename(id), "_pyin_f0.rds"))
  if (!file.exists(path)) return(NA_real_)
  checkpoint <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(checkpoint) || is.null(checkpoint$sample_rate_hz)) {
    return(NA_real_)
  }
  as.numeric(checkpoint$sample_rate_hz)
}, numeric(1))]

if (file.exists(CONFIG_FILE)) {
  config <- data.table::fread(CONFIG_FILE)
  config[, recording_id := as.character(recording_id)]
  metadata <- merge(
    metadata,
    config[, .(recording_id, file_name, sa_hz)],
    by = "recording_id",
    all.x = TRUE
  )

  # Fallback: read the sampling rate from the WAV header when the pYIN cache
  # does not provide it. Walks the RIFF chunks to the "fmt " chunk, whose
  # bytes 5-8 hold the sampling rate (little-endian). No package needed.
  read_wav_sample_rate <- function(path) {
    if (!file.exists(path)) return(NA_real_)
    con <- file(path, "rb")
    on.exit(close(con))
    riff <- readChar(con, 4L, useBytes = TRUE)
    readBin(con, "integer", n = 1L, size = 4L, endian = "little")
    wave <- readChar(con, 4L, useBytes = TRUE)
    if (!identical(riff, "RIFF") || !identical(wave, "WAVE")) return(NA_real_)
    repeat {
      chunk_id <- readChar(con, 4L, useBytes = TRUE)
      if (length(chunk_id) == 0L || nchar(chunk_id, type = "bytes") < 4L) {
        return(NA_real_)
      }
      chunk_size <- readBin(con, "integer", n = 1L, size = 4L, endian = "little")
      if (identical(chunk_id, "fmt ")) {
        readBin(con, "integer", n = 1L, size = 2L, endian = "little")
        readBin(con, "integer", n = 1L, size = 2L, endian = "little")
        rate <- readBin(con, "integer", n = 1L, size = 4L, endian = "little")
        return(as.numeric(rate))
      }
      seek(con, chunk_size + (chunk_size %% 2L), origin = "current")
    }
  }

  AUDIO_DIR <- "audio"
  metadata[is.na(sample_rate_hz) & !is.na(file_name),
           sample_rate_hz := vapply(file.path(AUDIO_DIR, file_name),
                                    read_wav_sample_rate, numeric(1))]
}

if (file.exists(SOURCES_FILE)) {
  sources <- data.table::fread(SOURCES_FILE)
  check_columns(sources, c("recording_id", "source"), SOURCES_FILE)
  sources[, recording_id := as.character(recording_id)]
  sources[, source := tolower(trimws(as.character(source)))]

  bad_sources <- setdiff(unique(sources$source), SOURCE_LEVELS)
  if (length(bad_sources) > 0L) {
    stop("Unexpected source labels in ", SOURCES_FILE, ": ",
         paste(bad_sources, collapse = ", "),
         ". Use: ", paste(SOURCE_LEVELS, collapse = ", "), call. = FALSE)
  }
  missing_ids <- setdiff(recording_ids, sources$recording_id)
  if (length(missing_ids) > 0L) {
    stop(SOURCES_FILE, " has no source for: ",
         paste(missing_ids, collapse = ", "), call. = FALSE)
  }

  metadata <- merge(metadata, sources[, .(recording_id, source)],
                    by = "recording_id", all.x = TRUE)
  message("Recording sources: ",
          paste(names(table(sources$source)), table(sources$source),
                sep = " = ", collapse = ", "))
} else {
  sources <- NULL
  data.table::fwrite(
    data.table::data.table(
      recording_id = recording_ids,
      source = ifelse(recording_ids == REFERENCE_ID, "original", "")
    ),
    "recording_sources_TEMPLATE.csv"
  )
  warning(
    SOURCES_FILE, " not found. A template was written to ",
    "recording_sources_TEMPLATE.csv; fill in direct or demucs for each ",
    "participant, save it as ", SOURCES_FILE, ", and rerun. ",
    "Source-group analyses are skipped."
  )
}

metadata[, Order := match(recording_id, recording_ids)]
data.table::setorder(metadata, Order)
metadata[, Order := NULL]

if (anyNA(metadata$sample_rate_hz)) {
  warning("Sampling rate unavailable for: ",
          paste(metadata[is.na(sample_rate_hz), recording_id], collapse = ", "),
          ". Check kafi_cache/f0/.")
}

data.table::fwrite(metadata, file.path(OUTPUT_DIR, "recording_metadata.csv"))


# -----------------------------------------------------------------------------
# 5. IS THE REFERENCE TYPICAL? (CENTRALITY RANK TEST, TABLE S15)
# -----------------------------------------------------------------------------
# c_i = mean Bhattacharyya distance from recording i to the other recordings.
# Under exchangeability of R with the participants, rank(c_R) is uniform, so
# p = #{i : c_i >= c_R} / N is an exact one-sided p-value for "R is less
# central than a typical performance".

centrality_for <- function(column, label) {
  P <- probability_matrix(column)
  n <- nrow(P)
  D <- matrix(0, n, n, dimnames = list(rownames(P), rownames(P)))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j) D[i, j] <- bhattacharyya_distance(P[i, ], P[j, ])
    }
  }
  centrality <- rowSums(D) / (n - 1L)
  data.table::data.table(
    weighting = label,
    recording_id = rownames(P),
    centrality = centrality,
    rank = data.table::frank(centrality, ties.method = "min"),
    p_typicality = NA_real_
  )
}

centrality_table <- data.table::rbindlist(list(
  centrality_for("Duration_Probability", "duration"),
  centrality_for("Count_Probability", "count")
))

centrality_table[, p_typicality := mean(
  centrality >= centrality[recording_id == REFERENCE_ID]
), by = weighting]
centrality_table[recording_id != REFERENCE_ID, p_typicality := NA_real_]
data.table::setorder(centrality_table, weighting, rank)

data.table::fwrite(centrality_table,
                   file.path(OUTPUT_DIR, "centrality_typicality.csv"))

typicality_summary <- centrality_table[recording_id == REFERENCE_ID,
                                       .(weighting, rank, centrality, p_typicality)]
message("Reference centrality rank (of ", length(recording_ids), "):")
print(typicality_summary)


# -----------------------------------------------------------------------------
# 6. PER-SWARA TESTS AGAINST THE REFERENCE (TABLE 3 AND TABLE S16)
# -----------------------------------------------------------------------------
# H0: the median participant share equals R's share. Sign test (primary),
# Wilcoxon signed-rank (sensitivity; continuity correction off), Holm over
# the seven swaras within each weighting.

swara_tests_for <- function(column, label, ids = participant_ids) {
  P <- probability_matrix(column)
  reference_share <- P[REFERENCE_ID, ]
  participant_shares <- P[ids, , drop = FALSE]

  result <- data.table::rbindlist(lapply(SWARA_ORDER, function(s) {
    shares <- participant_shares[, s]
    x <- shares - reference_share[s]
    n_above <- sum(x > 0)
    n_below <- sum(x < 0)
    sign_p <- if (n_above + n_below > 0L) {
      stats::binom.test(n_above, n_above + n_below, p = 0.5)$p.value
    } else {
      NA_real_
    }
    wilcoxon_p <- suppressWarnings(
      stats::wilcox.test(x, mu = 0, correct = FALSE)$p.value
    )
    data.table::data.table(
      weighting = label,
      swara = s,
      reference_share = reference_share[s],
      participant_median = stats::median(shares),
      participant_q1 = unname(stats::quantile(shares, 0.25)),
      participant_q3 = unname(stats::quantile(shares, 0.75)),
      median_difference = stats::median(x),
      n_above = n_above,
      n_below = n_below,
      n_tied = sum(x == 0),
      sign_p = sign_p,
      wilcoxon_p = wilcoxon_p,
      rank_biserial = rank_biserial(x)
    )
  }))

  result[, sign_p_holm := stats::p.adjust(sign_p, method = "holm")]
  result[, wilcoxon_p_holm := stats::p.adjust(wilcoxon_p, method = "holm")]
  result
}

swara_tests <- data.table::rbindlist(list(
  swara_tests_for("Duration_Probability", "duration"),
  swara_tests_for("Count_Probability", "count")
))

data.table::fwrite(swara_tests, file.path(OUTPUT_DIR, "swara_tests.csv"))

message("Per-swara tests (duration weighting):")
print(swara_tests[weighting == "duration", .(
  swara, reference_share, participant_median, median_difference,
  n_above, n_below, sign_p_holm, wilcoxon_p_holm, rank_biserial
)], digits = 3)


# -----------------------------------------------------------------------------
# 7. CONFOUND CHECKS (TABLE S17)
# -----------------------------------------------------------------------------

duration_matrix <- probability_matrix("Duration_Probability")
distance_to_reference <- vapply(participant_ids, function(id) {
  bhattacharyya_distance(duration_matrix[REFERENCE_ID, ], duration_matrix[id, ])
}, numeric(1))

characteristics <- merge(
  metadata[recording_id %in% participant_ids, .(recording_id, length_sec)],
  quality[, .(recording_id, voiced_percent = 100 * raw_voiced_fraction)],
  by = "recording_id"
)
characteristics <- merge(
  characteristics,
  parameters[, .(recording_id, n_events, P5_mean_note_duration,
                 P2_melodic_span)],
  by = "recording_id"
)
characteristics[, `:=`(
  sa_share = duration_matrix[recording_id, "Sa"],
  dha_share = duration_matrix[recording_id, "Dha"],
  distance_to_R = distance_to_reference[recording_id]
)]

predictors <- c(
  "Recording length (s)" = "length_sec",
  "Number of events" = "n_events",
  "Voiced frames (%)" = "voiced_percent",
  "Mean event duration P5" = "P5_mean_note_duration",
  "Conventional span P2" = "P2_melodic_span"
)
outcomes <- c("Sa share" = "sa_share", "Dha share" = "dha_share",
              "Distance to R" = "distance_to_R")

confounds <- data.table::rbindlist(lapply(names(predictors), function(pl) {
  data.table::rbindlist(lapply(names(outcomes), function(ol) {
    test <- suppressWarnings(stats::cor.test(
      characteristics[[predictors[[pl]]]],
      characteristics[[outcomes[[ol]]]],
      method = "spearman",
      exact = FALSE
    ))
    data.table::data.table(characteristic = pl, outcome = ol,
                           rho = unname(test$estimate), p_value = test$p.value)
  }))
}))

data.table::fwrite(confounds, file.path(OUTPUT_DIR, "confound_spearman.csv"))


# -----------------------------------------------------------------------------
# 7A. SOURCE-SEPARATION CHECK (SECTION 4.3)
# -----------------------------------------------------------------------------
# R and the directly supplied recordings were analysed as received; media
# recordings were analysed as Demucs v4 vocal stems. Two checks:
#   (a) Wilcoxon rank-sum, separated vs direct participants, with the
#       Hodges-Lehmann shift (demucs minus direct) and its 95% interval;
#   (b) per-swara sign tests against R within the direct recordings only.

if (!is.null(sources)) {
  direct_ids <- intersect(participant_ids, sources[source == "direct", recording_id])
  demucs_ids <- intersect(participant_ids, sources[source == "demucs", recording_id])

  message("Source groups: direct = ", length(direct_ids),
          ", demucs = ", length(demucs_ids))

  group_data <- data.table::data.table(
    recording_id = participant_ids,
    sa_share = duration_matrix[participant_ids, "Sa"],
    dha_share = duration_matrix[participant_ids, "Dha"],
    re_share = duration_matrix[participant_ids, "Re"],
    distance_to_R = distance_to_reference[participant_ids]
  )
  group_data[, source := sources$source[match(recording_id, sources$recording_id)]]

  source_outcomes <- c("Sa share" = "sa_share", "Re share" = "re_share",
                       "Dha share" = "dha_share", "Distance to R" = "distance_to_R")

  source_group_tests <- data.table::rbindlist(lapply(names(source_outcomes), function(ol) {
    y_demucs <- group_data[source == "demucs"][[source_outcomes[[ol]]]]
    y_direct <- group_data[source == "direct"][[source_outcomes[[ol]]]]
    if (length(y_demucs) < 2L || length(y_direct) < 2L) {
      return(data.table::data.table(outcome = ol, n_demucs = length(y_demucs),
                                    n_direct = length(y_direct)))
    }
    test <- suppressWarnings(stats::wilcox.test(
      y_demucs, y_direct, conf.int = TRUE, exact = FALSE, correct = FALSE
    ))
    data.table::data.table(
      outcome = ol,
      n_demucs = length(y_demucs),
      n_direct = length(y_direct),
      median_demucs = stats::median(y_demucs),
      median_direct = stats::median(y_direct),
      reference_value = if (ol == "Distance to R") 0 else
        duration_matrix[REFERENCE_ID, c(sa_share = "Sa", re_share = "Re",
                                        dha_share = "Dha")[[source_outcomes[[ol]]]]],
      hodges_lehmann_shift = unname(test$estimate),
      shift_ci_lower = test$conf.int[1],
      shift_ci_upper = test$conf.int[2],
      rank_sum_p = test$p.value
    )
  }), fill = TRUE)

  data.table::fwrite(source_group_tests,
                     file.path(OUTPUT_DIR, "source_group_tests.csv"))
  message("Separated vs direct participants (duration weighting):")
  print(source_group_tests, digits = 3)

  if (length(direct_ids) >= 2L) {
    direct_only <- data.table::rbindlist(list(
      swara_tests_for("Duration_Probability", "duration", ids = direct_ids),
      swara_tests_for("Count_Probability", "count", ids = direct_ids)
    ))
    direct_only[, subgroup := "direct"]
    data.table::fwrite(direct_only,
                       file.path(OUTPUT_DIR, "swara_tests_direct_only.csv"))
    message("Direct recordings only vs R (duration weighting):")
    print(direct_only[weighting == "duration", .(
      swara, reference_share, participant_median, n_above, n_below,
      sign_p, sign_p_holm
    )], digits = 3)
  }
}


# -----------------------------------------------------------------------------
# 8. BREAKDOWN DIAGNOSTICS (TABLE S18)
# -----------------------------------------------------------------------------
# w_max is computed exactly from event durations. t is the number of trimmed
# differences in robust P3; the worst-case condition of Proposition 3 fails
# when 2q > t, with q = ceiling(rate x n) shifted events.

breakdown <- events[recording_id %in% recording_ids, {
  n <- .N
  n_differences <- n - 1L
  retained <- max(1L, as.integer(ceiling(
    ROBUST_TRANSITION_RETAINED_PROB * n_differences
  )))
  trimmed <- n_differences - retained
  w <- Duration / sum(Duration)
  out <- list(
    n_events = n,
    w_max = max(w),
    trimmed_differences = trimmed
  )
  for (rate in CONTAMINATION_RATES) {
    q <- max(1L, as.integer(ceiling(rate * n)))
    out[[sprintf("q_%s", format(100 * rate))]] <- q
    out[[sprintf("trim_exhausted_%s", format(100 * rate))]] <- 2L * q > trimmed
  }
  out
}, by = recording_id]

breakdown[, Order := match(recording_id, recording_ids)]
data.table::setorder(breakdown, Order)
breakdown[, Order := NULL]

data.table::fwrite(breakdown, file.path(OUTPUT_DIR, "breakdown_diagnostics.csv"))

message("Recordings whose longest event exceeds 5% of duration: ",
        breakdown[w_max > 0.05, .N], " of ", nrow(breakdown),
        " (median w_max = ", round(stats::median(breakdown$w_max), 3), ")")


# -----------------------------------------------------------------------------
# 9. MARKOV-CHAIN HOMOGENEITY (SECTION 4.4)
# -----------------------------------------------------------------------------
# Anderson-Goodman likelihood-ratio statistic, calibrated by permutation:
# successors are permuted across recordings within each originating swara,
# which keeps every recording's row totals fixed.

# Aggregate defensively in case a From--To pair appears more than once, then
# expand to one row per observed transition.
transition_counts <- transitions_raw[
  recording_id %in% recording_ids,
  .(Count = as.integer(round(sum(Count)))),
  by = .(recording_id, From, To)
][Count > 0]

transition_long <- transition_counts[
  rep(seq_len(.N), Count),
  .(recording_id, From, To)
]

g2_statistic <- function(long) {
  d <- long[, .(N = .N), by = .(recording_id, From, To)]
  d[, recording_row_total := sum(N), by = .(recording_id, From)]
  d[, pooled_cell := sum(N), by = .(From, To)]
  d[, pooled_row_total := sum(N), by = From]
  2 * sum(d$N * log((d$N / d$recording_row_total) /
                      (d$pooled_cell / d$pooled_row_total)))
}

permutation_test <- function(long, label, n_permutations = N_PERMUTATIONS) {
  observed <- g2_statistic(long)
  null_values <- numeric(n_permutations)
  permuted <- data.table::copy(long)
  for (b in seq_len(n_permutations)) {
    permuted[, To := To[sample.int(.N)], by = From]
    null_values[b] <- g2_statistic(permuted)
  }
  data.table::data.table(
    test = label,
    n_transitions = nrow(long),
    G2 = observed,
    permutation_p = (1 + sum(null_values >= observed)) / (n_permutations + 1),
    n_permutations = n_permutations
  )
}

group_reference <- function(long) {
  grouped <- data.table::copy(long)
  grouped[, recording_id := data.table::fifelse(
    recording_id == REFERENCE_ID, "R", "Participants"
  )]
  grouped
}

set.seed(ANALYSIS_SEED)
without_self <- transition_long[From != To]

markov_results <- data.table::rbindlist(list(
  permutation_test(transition_long, "All recordings, with self-transitions"),
  permutation_test(without_self, "All recordings, without self-transitions"),
  permutation_test(group_reference(transition_long),
                   "R vs participants, with self-transitions"),
  permutation_test(group_reference(without_self),
                   "R vs participants, without self-transitions")
))

data.table::fwrite(markov_results, file.path(OUTPUT_DIR, "markov_homogeneity.csv"))
message("Markov homogeneity tests:")
print(markov_results, digits = 4)


# -----------------------------------------------------------------------------
# 10. FIGURE 2: DURATION-WEIGHTED SWARA SHARES AGAINST THE REFERENCE
# -----------------------------------------------------------------------------

plot_data <- swara_tests[weighting == "duration"]
plot_data[, swara := factor(swara, levels = rev(SWARA_ORDER))]
plot_data[, significant := sign_p_holm < ALPHA]
plot_data[, above_label := sprintf("%d/%d", n_above, length(participant_ids))]

label_x <- 1.08 * max(plot_data$participant_q3, plot_data$reference_share)

swara_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(y = swara)) +
  ggplot2::geom_segment(
    ggplot2::aes(x = participant_q1, xend = participant_q3, yend = swara,
                 colour = significant),
    linewidth = 1.4, lineend = "round"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = participant_median, colour = significant),
    size = 2.8
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = reference_share),
    shape = 124, size = 6, colour = "#333333"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = label_x, label = above_label),
    hjust = 0, size = 3.4, colour = "#333333"
  ) +
  ggplot2::annotate("text", x = label_x, y = length(SWARA_ORDER) + 0.7,
                    label = "Participants\nabove R", hjust = 0, size = 3,
                    colour = "#555555", lineheight = 0.9) +
  ggplot2::scale_colour_manual(
    values = c(`TRUE` = "#1F5FBF", `FALSE` = "#9A9A9A"),
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, label_x * 1.18),
    breaks = seq(0, 0.5, by = 0.1),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    x = "Share of extracted event duration",
    # The explanation of marks lives in the LaTeX caption of Figure 2.
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    plot.caption = ggplot2::element_text(hjust = 0, colour = "#555555"),
    plot.margin = ggplot2::margin(18, 10, 6, 6)
  )

ggplot2::ggsave(
  file.path(OUTPUT_DIR, "swara_share_vs_reference.png"),
  swara_plot, width = 6.5, height = 3.6, units = "in", dpi = 600, bg = "white"
)
ggplot2::ggsave(
  file.path(OUTPUT_DIR, "swara_share_vs_reference.pdf"),
  swara_plot, width = 6.5, height = 3.6, units = "in", device = "pdf",
  bg = "white"
)


# -----------------------------------------------------------------------------
# 11. SETTINGS, SESSION INFORMATION AND COMPLETION MESSAGE
# -----------------------------------------------------------------------------

data.table::fwrite(
  data.table::data.table(
    Setting = c("Reference ID", "Recordings", "Permutations", "Seed",
                "Robust transition retained proportion",
                "Contamination rates", "Significance level"),
    Value = c(REFERENCE_ID, length(recording_ids), N_PERMUTATIONS,
              ANALYSIS_SEED, ROBUST_TRANSITION_RETAINED_PROB,
              paste(CONTAMINATION_RATES, collapse = ", "), ALPHA)
  ),
  file.path(OUTPUT_DIR, "analysis_settings.csv")
)

writeLines(
  capture.output(utils::sessionInfo()),
  file.path(OUTPUT_DIR, "session_info.txt")
)

message("Inferential analyses completed successfully.")
message("Results are available in: ", normalizePath(OUTPUT_DIR))
