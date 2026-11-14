# =============================================================================
# RUN THE COMPLETE KAFI ANALYSIS
# =============================================================================
#
# Authors : Prena Singh, Soubhik Chakraborty and Pushkal Kumar
# Licence : MIT (see LICENSE)
#
# Runs the four analysis scripts in order, each in a fresh R process, from the
# project directory. Every script caches its work, so a rerun after an
# interruption resumes where it stopped.
#
# Usage: Rscript run_all.R
# =============================================================================

scripts <- c(
  "Kafi_Swara_Analysis.R",
  "Kafi_Transition_Heatmap_Analysis.R",
  "Kafi_Robustness_Simulation.R",
  "Kafi_Inference_Analyses.R"
)

missing_scripts <- scripts[!file.exists(scripts)]
if (length(missing_scripts) > 0L) {
  stop("Run this file from the project directory. Missing: ",
       paste(missing_scripts, collapse = ", "), call. = FALSE)
}

rscript <- file.path(R.home("bin"), "Rscript")

for (script in scripts) {
  message("\n==== Running ", script, " ====")
  started <- Sys.time()
  status <- system2(rscript, shQuote(script))
  if (!identical(status, 0L)) {
    stop(script, " stopped with exit status ", status,
         ". Fix the problem and rerun; completed steps are cached.",
         call. = FALSE)
  }
  message("Finished ", script, " in ",
          format(round(difftime(Sys.time(), started, units = "mins"), 1)))
}

message("\nAll analyses completed. Results are in output/.")
