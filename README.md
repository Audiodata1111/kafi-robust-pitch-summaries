# Robust summaries of tonic-normalised pitch events: Raga Kafi case study

R code for the article:

> Singh, P., Chakraborty, S. and Kumar, P. *Robust Summaries of Tonic-Normalised Pitch Events, with an Application to Performance Variation in a Raga Kafi Bandish.* (Submitted.)

The code turns vocal recordings into tonic-normalised note events, computes conventional and robust pitch summaries, studies their sensitivity to octave-tracking errors by Monte Carlo simulation, and carries out the inferential analyses reported in the article.

## Contents

| File | Purpose |
| --- | --- |
| `Kafi_Swara_Analysis.R` | Main pipeline: pYIN pitch tracking, tonic normalisation, swara labelling, temporal cleaning, note events, parameters P1-P5 (conventional and robust), swara distributions and Bhattacharyya distances |
| `Kafi_Transition_Heatmap_Analysis.R` | First-order swara-transition matrices, dominant transitions and the transition heat map |
| `Kafi_Robustness_Simulation.R` | Monte Carlo study of conventional and robust P1-P3 under simulated octave errors |
| `Kafi_Inference_Analyses.R` | Typicality rank test, per-swara sign and Wilcoxon tests with Holm adjustment, confound checks, breakdown diagnostics, Markov homogeneity tests and Figure 2 |
| `run_all.R` | Runs the four scripts in order |
| `recordings.csv` | Recording identifiers, file names and tonic (Sa) of each recording |
| `LICENSE` | MIT licence |

## Requirements

- R (version 4.2 or later recommended) with the packages `data.table`, `ggplot2` and `reticulate`:

  ```r
  install.packages(c("data.table", "ggplot2", "reticulate"))
  ```

- Python with `librosa` 0.11.0. This is installed automatically by `reticulate` the first time pitch tracking runs; that first run needs an internet connection. Later runs work offline.

## Folder layout

```
project/
  Kafi_Swara_Analysis.R
  Kafi_Transition_Heatmap_Analysis.R
  Kafi_Robustness_Simulation.R
  Kafi_Inference_Analyses.R
  run_all.R
  recordings.csv
  audio/              reference.wav, A01.wav, ..., A33.wav
```

The folders `output/` and `kafi_cache/` are created when the scripts run. Run everything from the project folder (an RStudio project, or `Rscript` from that folder).

## Data

The corpus contains 34 vocal recordings of one bandish in Raga Kafi: a reference performance (R) and 33 participant performances (A01-A33).

- The reference and nine participant recordings were supplied directly by the performers and are unaccompanied.
- The other 24 participant recordings were taken from official online media. Their vocal part was separated from the accompaniment with Demucs v4 (Rouard, Massa and Défossez, 2023) before analysis.

**Audio files are not included in this repository.** <!-- Replace with the access route for the audio: the data repository link, the original source links for the media recordings, or "available from the authors on request". -->

Place the WAV files in `audio/` with the names given in `recordings.csv`.

Each tonic in `recordings.csv` is the median pYIN F0 over a stable region in which the performer sustains Sa, identified on the pitch contour with expert guidance.

## Running the analysis

Run all steps:

```bash
Rscript run_all.R
```

or run the scripts one at a time, in this order:

1. `Kafi_Swara_Analysis.R` (must run first; all other scripts use its outputs)
2. `Kafi_Transition_Heatmap_Analysis.R`
3. `Kafi_Robustness_Simulation.R`
4. `Kafi_Inference_Analyses.R`

Steps 2 to 4 are independent of one another.

**Run time.** Pitch tracking is the slowest step. It runs once per recording and is cached in `kafi_cache/f0/`. The simulation (34 recordings, 4 contamination levels, 1,000 replicates each) is the next slowest; it is checkpointed per recording, so an interrupted run resumes where it stopped.

## Where each result comes from

| Article item | Script | Output |
| --- | --- | --- |
| Table 2, Tables S2-S7 (extraction quality, parameters, descriptive statistics, Shapiro-Wilk) | `Kafi_Swara_Analysis.R` | `output/quality_summary.csv`, `output/P1_to_P5.csv`, `output/descriptive_statistics.csv`, `output/shapiro_wilk_results.csv` |
| Figure 1, Tables S8-S9 (Bhattacharyya distances, swara distributions) | `Kafi_Swara_Analysis.R` | `output/bhattacharyya_distances.csv`, `output/swara_distributions.csv` |
| Figure 3, Tables S10-S12 (transitions) | `Kafi_Transition_Heatmap_Analysis.R` | `output/transition_summary/` |
| Figure 4, Tables S13-S14 (robustness) | `Kafi_Robustness_Simulation.R` | `output/robustness_simulation/` |
| Table 3, Figure 2, Tables S15-S18 | `Kafi_Inference_Analyses.R` | `output/inference_analyses/` |
| Table 4 (Markov homogeneity) | `Kafi_Inference_Analyses.R` | `output/inference_analyses/markov_homogeneity.csv` |

<!-- Check the table and figure numbers against the final version of the article. -->

## Main settings

| Setting | Value | Script |
| --- | --- | --- |
| pYIN search range | 70-700 Hz | `Kafi_Swara_Analysis.R` (`recordings.csv`) |
| Frame length, frame step | 2048 samples, about 5 ms | `Kafi_Swara_Analysis.R` |
| Swara template | Just intonation; tolerance 45 cents | `Kafi_Swara_Analysis.R` |
| Gap bridging, minimum event length | 25 ms, 50 ms | `Kafi_Swara_Analysis.R` |
| Robust quantiles, trimming | 5% and 95%; 10% of largest steps | `Kafi_Swara_Analysis.R`, `Kafi_Robustness_Simulation.R` |
| Simulation | 1,000 replicates; 1, 2.5, 5, 10% contamination; ±12 semitones; seed 20260829 | `Kafi_Robustness_Simulation.R` |
| Permutation tests | 2,000 permutations; seed 20261022 | `Kafi_Inference_Analyses.R` |

Each output folder contains the settings used and `session_info.txt` with package versions.

**Optional analyses.**

- Set `EVENT_SELECTION <- "duration"` in `Kafi_Robustness_Simulation.R` for the duration-weighted contamination design (results go to `output/robustness_simulation_duration/`).
- Set `TUNING_SYSTEM <- "equal"` in `Kafi_Swara_Analysis.R` to repeat the analysis with an equal-tempered template.

## Troubleshooting

| Message | Cause | Fix |
| --- | --- | --- |
| "Supply a positive sa_hz" | Missing tonic in `recordings.csv` | Enter the tonic |
| "The following audio files do not exist" | Wrong file name or folder | Check `audio/` and `recordings.csv` |
| "pYIN extraction failed" | Python or librosa problem | Read `kafi_cache/logs/<ID>_pyin.log` |
| "Stage 2 cannot begin" | A checkpoint is missing or out of date | Set `RUN_STAGE <- "both"` in `Kafi_Swara_Analysis.R` and rerun |

## Citation

If you use this code, please cite the article above. <!-- Add the DOI after publication, and the Zenodo DOI of this repository. -->

## Licence

The code is released under the MIT licence; see `LICENSE`. The licence covers the code only, not the audio recordings.

## Contact

Pushkal Kumar, Department of Quantitative Economics and Data Science, Birla Institute of Technology Mesra, Ranchi, India. <!-- Add a contact email if you wish. -->
