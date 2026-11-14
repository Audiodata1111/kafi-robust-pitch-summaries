# =============================================================================
# KAFI SWARA-TRANSITION MATRICES AND HEAT MAP
# =============================================================================
#
# Authors : Prena Singh, Soubhik Chakraborty and Pushkal Kumar
#           Birla Institute of Technology Mesra, Ranchi, India
# Licence : MIT (see LICENSE)
# Part of : Robust Summaries of Tonic-Normalised Pitch Events, with an
#           Application to Performance Variation in a Raga Kafi Bandish
#
# Input
# -----
# The script searches for:
#   output/all_transitions.csv
# and then:
#   all_transitions.csv
#
# Required input columns:
#   recording_id, From, To, Count
#
# The existing Probability column, when present, is not trusted blindly;
# probabilities are recomputed from the transition counts.
#
# Outputs
# -------
# output/transition_summary/
#   all_transition_matrices_complete.csv
#   reference_transition_matrix.csv
#   participant_mean_transition_matrix.csv
#   dominant_transitions_reference.csv
#   dominant_transitions_participants.csv
#   transition_heatmap_R_vs_participants.png
#   transition_heatmap_R_vs_participants.pdf
#
# Required packages
# -----------------
# install.packages(c("data.table", "ggplot2"))
# =============================================================================


# -----------------------------------------------------------------------------
# 1. SETTINGS
# -----------------------------------------------------------------------------

REFERENCE_ID <- "R"
EXPECTED_PARTICIPANTS <- 33L
TOP_N_TRANSITIONS <- 10L

SWARA_ORDER <- c("Sa", "Re", "g", "Ma", "Pa", "Dha", "n")

INPUT_CANDIDATES <- c(
  file.path("output", "all_transitions.csv"),
  "all_transitions.csv"
)

OUTPUT_DIR <- file.path("output", "transition_summary")


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

existing_inputs <- INPUT_CANDIDATES[file.exists(INPUT_CANDIDATES)]

if (length(existing_inputs) == 0L) {
  stop(
    "Could not find all_transitions.csv. Place it either in the project ",
    "directory or in the output directory.",
    call. = FALSE
  )
}

INPUT_FILE <- existing_inputs[1L]
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Reading transition data from: ", normalizePath(INPUT_FILE))


# -----------------------------------------------------------------------------
# 3. READ AND VALIDATE THE OBSERVED TRANSITIONS
# -----------------------------------------------------------------------------

transitions_raw <- data.table::fread(INPUT_FILE)

required_columns <- c("recording_id", "From", "To", "Count")
missing_columns <- setdiff(required_columns, names(transitions_raw))

if (length(missing_columns) > 0L) {
  stop(
    "The transition file is missing: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

transitions_raw[, recording_id := as.character(recording_id)]
transitions_raw[, From := as.character(From)]
transitions_raw[, To := as.character(To)]
transitions_raw[, Count := as.numeric(Count)]

if (anyNA(transitions_raw[, .(recording_id, From, To, Count)])) {
  stop("The required transition columns contain missing values.", call. = FALSE)
}

if (any(!is.finite(transitions_raw$Count) | transitions_raw$Count < 0)) {
  stop("Transition counts must be finite and non-negative.", call. = FALSE)
}

unexpected_swaras <- setdiff(
  union(unique(transitions_raw$From), unique(transitions_raw$To)),
  SWARA_ORDER
)

if (length(unexpected_swaras) > 0L) {
  stop(
    "Unexpected swara labels were found: ",
    paste(unexpected_swaras, collapse = ", "),
    call. = FALSE
  )
}

recording_ids <- unique(transitions_raw$recording_id)

if (!(REFERENCE_ID %in% recording_ids)) {
  stop("Reference recording R is not present in the transition file.",
       call. = FALSE)
}

participant_ids <- setdiff(recording_ids, REFERENCE_ID)

if (length(participant_ids) != EXPECTED_PARTICIPANTS) {
  warning(
    "Expected ", EXPECTED_PARTICIPANTS, " participant recordings, but found ",
    length(participant_ids), ". The available participant recordings will be used."
  )
}


# -----------------------------------------------------------------------------
# 4. COMPLETE EVERY RECORDING TO A 7 x 7 MATRIX
# -----------------------------------------------------------------------------

# Aggregate defensively in case a From--To pair appears more than once.
transition_counts <- transitions_raw[, .(
  Count = sum(Count)
), by = .(recording_id, From, To)]

complete_grid <- data.table::CJ(
  recording_id = recording_ids,
  From = SWARA_ORDER,
  To = SWARA_ORDER,
  unique = TRUE
)

complete_transitions <- merge(
  complete_grid,
  transition_counts,
  by = c("recording_id", "From", "To"),
  all.x = TRUE,
  sort = FALSE
)

complete_transitions[is.na(Count), Count := 0]
complete_transitions[, Row_Total := sum(Count),
                     by = .(recording_id, From)]
complete_transitions[, Probability := data.table::fifelse(
  Row_Total > 0,
  Count / Row_Total,
  NA_real_
)]

complete_transitions[, From_Order := match(From, SWARA_ORDER)]
complete_transitions[, To_Order := match(To, SWARA_ORDER)]
data.table::setorder(
  complete_transitions,
  recording_id,
  From_Order,
  To_Order
)

# Numerical validation: every estimable row must sum to one.
row_checks <- complete_transitions[Row_Total > 0, .(
  Probability_Sum = sum(Probability)
), by = .(recording_id, From)]

bad_rows <- row_checks[abs(Probability_Sum - 1) > 1e-10]

if (nrow(bad_rows) > 0L) {
  stop("At least one transition-probability row does not sum to one.",
       call. = FALSE)
}


# -----------------------------------------------------------------------------
# 5. REFERENCE AND EQUAL-RECORDING-WEIGHT PARTICIPANT MATRICES
# -----------------------------------------------------------------------------

reference_matrix <- data.table::copy(
  complete_transitions[recording_id == REFERENCE_ID]
)

# For each originating swara, average only across participant recordings in
# which that swara has at least one outgoing transition. Because an absent row
# is NA for all seven destination cells, every retained recording contributes
# one complete conditional-probability row with equal weight.
participant_mean_matrix <- complete_transitions[
  recording_id %in% participant_ids,
  .(
    Mean_Probability = if (all(is.na(Probability))) {
      NA_real_
    } else {
      mean(Probability, na.rm = TRUE)
    },
    Median_Probability = if (all(is.na(Probability))) {
      NA_real_
    } else {
      stats::median(Probability, na.rm = TRUE)
    },
    Total_Count = sum(Count),
    N_Recordings_With_Origin = sum(Row_Total > 0),
    N_Recordings_With_Transition = sum(Count > 0)
  ),
  by = .(From, To)
]

participant_mean_matrix[, Transition_Prevalence := data.table::fifelse(
  N_Recordings_With_Origin > 0,
  N_Recordings_With_Transition / N_Recordings_With_Origin,
  NA_real_
)]
participant_mean_matrix[, From_Order := match(From, SWARA_ORDER)]
participant_mean_matrix[, To_Order := match(To, SWARA_ORDER)]
data.table::setorder(participant_mean_matrix, From_Order, To_Order)

participant_row_checks <- participant_mean_matrix[
  !is.na(Mean_Probability),
  .(Probability_Sum = sum(Mean_Probability)),
  by = From
]

if (any(abs(participant_row_checks$Probability_Sum - 1) > 1e-10)) {
  stop("At least one mean participant transition row does not sum to one.",
       call. = FALSE)
}


# -----------------------------------------------------------------------------
# 6. DOMINANT-TRANSITION TABLES
# -----------------------------------------------------------------------------

reference_dominant <- data.table::copy(reference_matrix[Count > 0])
# Rank by conditional probability, as for the participant ranking below.
# Count and Row_Total are written out so that sparse rows (for example
# Ma -> g, 5 of 10) can be identified.
data.table::setorder(reference_dominant, -Probability, -Count, From_Order,
                     To_Order)
reference_dominant <- head(reference_dominant, TOP_N_TRANSITIONS)
reference_dominant[, Rank := seq_len(.N)]
reference_dominant <- reference_dominant[, .(
  Rank, From, To, Count, Row_Total, Probability
)]

participant_dominant <- data.table::copy(
  participant_mean_matrix[N_Recordings_With_Transition > 0]
)
data.table::setorder(
  participant_dominant,
  -Mean_Probability,
  -Transition_Prevalence,
  -Total_Count,
  From_Order,
  To_Order
)
participant_dominant <- head(participant_dominant, TOP_N_TRANSITIONS)
participant_dominant[, Rank := seq_len(.N)]
participant_dominant <- participant_dominant[, .(
  Rank,
  From,
  To,
  Mean_Probability,
  Median_Probability,
  Total_Count,
  N_Recordings_With_Origin,
  N_Recordings_With_Transition,
  Transition_Prevalence
)]


# -----------------------------------------------------------------------------
# 7. WRITE NUMERICAL OUTPUTS
# -----------------------------------------------------------------------------

data.table::fwrite(
  complete_transitions[, .(
    recording_id, From, To, Count, Row_Total, Probability
  )],
  file.path(OUTPUT_DIR, "all_transition_matrices_complete.csv")
)

reference_wide <- data.table::dcast(
  reference_matrix,
  From ~ To,
  value.var = "Probability"
)
data.table::setcolorder(reference_wide, c("From", SWARA_ORDER))
reference_wide[, From_Order := match(From, SWARA_ORDER)]
data.table::setorder(reference_wide, From_Order)
reference_wide[, From_Order := NULL]

participant_wide <- data.table::dcast(
  participant_mean_matrix,
  From ~ To,
  value.var = "Mean_Probability"
)
data.table::setcolorder(participant_wide, c("From", SWARA_ORDER))
participant_wide[, From_Order := match(From, SWARA_ORDER)]
data.table::setorder(participant_wide, From_Order)
participant_wide[, From_Order := NULL]

data.table::fwrite(
  reference_wide,
  file.path(OUTPUT_DIR, "reference_transition_matrix.csv")
)

data.table::fwrite(
  participant_wide,
  file.path(OUTPUT_DIR, "participant_mean_transition_matrix.csv")
)

data.table::fwrite(
  reference_dominant,
  file.path(OUTPUT_DIR, "dominant_transitions_reference.csv")
)

data.table::fwrite(
  participant_dominant,
  file.path(OUTPUT_DIR, "dominant_transitions_participants.csv")
)


# -----------------------------------------------------------------------------
# 8. TWO-PANEL PUBLICATION HEATMAP
# -----------------------------------------------------------------------------

plot_reference <- reference_matrix[, .(
  Panel = "(a) Reference R",
  From,
  To,
  Probability
)]

plot_participants <- participant_mean_matrix[, .(
  Panel = "(b) Mean of A01--A33",
  From,
  To,
  Probability = Mean_Probability
)]

plot_data <- data.table::rbindlist(
  list(plot_reference, plot_participants),
  use.names = TRUE
)

plot_data[, Panel := factor(
  Panel,
  levels = c("(a) Reference R", "(b) Mean of A01--A33")
)]
plot_data[, To := factor(To, levels = SWARA_ORDER)]
plot_data[, From := factor(From, levels = rev(SWARA_ORDER))]
plot_data[, Label := ifelse(
  is.na(Probability),
  "--",
  sprintf("%.2f", Probability)
)]
plot_data[, Label_Colour := ifelse(
  !is.na(Probability) & Probability >= 0.50,
  "white",
  "#17202A"
)]

transition_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = To, y = From, fill = Probability)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = Label, colour = Label_Colour),
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_fill_gradientn(
    colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    na.value = "#EEEEEE",
    name = "Conditional\nprobability"
  ) +
  ggplot2::scale_x_discrete(drop = FALSE) +
  ggplot2::scale_y_discrete(drop = FALSE) +
  ggplot2::facet_wrap(~Panel, nrow = 1) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Sequential swara-transition probabilities",
    subtitle = "Rows show the current swara; columns show the succeeding swara",
    x = "Next swara",
    y = "Current swara"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.spacing = grid::unit(1.1, "lines"),
    strip.background = ggplot2::element_rect(
      fill = "#E9EEF5",
      colour = "#9AA7B2",
      linewidth = 0.5
    ),
    strip.text = ggplot2::element_text(face = "bold", size = 12),
    axis.text = ggplot2::element_text(colour = "#17202A", size = 10.5),
    axis.title = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold", size = 15),
    plot.subtitle = ggplot2::element_text(size = 10.5, colour = "#3F4C5A"),
    legend.title = ggplot2::element_text(face = "bold"),
    legend.position = "right"
  )

ggplot2::ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "transition_heatmap_R_vs_participants.png"
  ),
  plot = transition_plot,
  width = 12,
  height = 5.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "transition_heatmap_R_vs_participants.pdf"
  ),
  plot = transition_plot,
  width = 12,
  height = 5.8,
  units = "in",
  device = "pdf",
  bg = "white"
)


# -----------------------------------------------------------------------------
# 9. COMPLETION MESSAGE
# -----------------------------------------------------------------------------

message("Transition analysis completed successfully.")
message("Recordings found: ", length(recording_ids))
message("Participant recordings found: ", length(participant_ids))
message("Results are available in: ", normalizePath(OUTPUT_DIR))

