# Calibration pipeline of FINNmodelBWI

Self-contained R scripts (house rules: one script per step, `PARAMS` block,
`Rscript <script>` from this repository's root, no scientific parameters
through the environment). `run_all.R` runs the whole chain in order:

| step | script | what | runtime |
|---|---|---|---|
| 1 | `run_all.R` 1a-1d | downloads the raw layers (BWI ACCESS2003 release, Thuenen inventory-archive API, DWD grids, WorldClim) into the bwi-dataprep checkout; `PROVENANCE.md` has sources, licences, dates and sha256 | hours (API) |
| 2 | `run_all.R` 2a-2l, `prep/02_*.R`, `prep/04_*.R` | bwi-dataprep (pinned commit) 00, 01, 02 (all tiers, literal copy), 03, 04 (literal copy, first pass), 16, 04 (second pass, with drought), 40, 70, 91, 71 -> `two-interval-v8`; checked against `two-interval-v8.sha256` | 1-2 h |
| 3 | `01_stage_data_and_weights.R` | copies the prepared sets (and, from a FINN-bwi study checkout, the fitted weights) into `FINN-zoo-data/bwi/`, verifies sha256 | 1 min |
| 4 | `07_expected_preparation.R` | writes / compares `expected_preparation.csv` | 1 min |
| 5 | `02_fit_stage1.R` | stage 1: per-process prefit (L-BFGS) | minutes |
| 6 | `03_fit_stage2.R` | stage 2: inverse calibration of growth, mortality, competition, Adam 1000 epochs, `weight_scale = "raw"` | ~6 h (1 thread) |
| 7 | `04_fit_stage3.R` | stage 3: regeneration + Beverton-Holt cap from the epoch-700 checkpoint, `weight_scale = "raw"` | ~30 min |
| 8 | `06_extract_parameters.R` | `parameters.pt` (named list of tensors, the release asset) + reproduction check of `build_model()` | 1 min |
| 9 | `05_validation.R` | holdout metrics and figures (`validation/`) from the rebuilt model | 1 min |
| 10 | `09_example_sites.R` | the four example stands shipped in `inst/extdata` (CC BY-SA 4.0 with attribution; 1 km cells rounded, no ids) | seconds |
| - | `08_quickstart_check.R` | the vignette's quickstart as a plain script | seconds |

Outputs that are committed: `PREPARATION.md`, `PROVENANCE.md`,
`provenance_sha256.txt`, `two-interval-v8.sha256`, `expected_preparation.csv`,
`parameters_summary.csv`, `parameters_extraction.txt`, `validation/` (metrics
csv, `report.md`, figures). Data, fits and `.pt` files live outside git in
`../FINN-zoo-data/bwi/` (sibling of the repository's parent, or
`../FINN-zoo-data/bwi/` next to it).

FINN: 0.2.0 (main `0848b60`), `.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))`
in every script. The released v0.1.0 parameters come from the same recipe run
with the 0.2.0 release candidate `6705a50` in the FINN-bwi study (stage 1 is
bit-identical under both; `weight_scale = "raw"` makes stages 2 and 3 score the
same objective under 0848b60).
