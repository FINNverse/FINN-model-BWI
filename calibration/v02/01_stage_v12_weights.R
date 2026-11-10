#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINNmodelBWI | calibration/v02/01_stage_v12_weights.R
##
## WHAT   Stages the fitted weights of the study's v12 chain (the same recipe
##        as v0.1.0, refitted on two-interval-v12 under FINN 0.2.0 main 0848b60
##        with fit(weight_scale = "raw")) into the zoo weights folder of the
##        v0.2.0 release, with sha256 sums. After it, run
##        calibration/06_extract_parameters.R (parameters.pt + reproduction
##        check) and calibration/05_validation.R with weights.json's tag set to
##        v0.2.0.
## WHY    v0.2 for BWI takes its parameters from the study's v12 refit (running
##        on the laptop and on Max's node); a separate zoo fit would duplicate
##        that compute. The candidate is the user's choice; the default below is
##        the chain with the same structure as v0.1.0 (linear env responses,
##        single drought index drought2, fitted recruitment cap, stage 2 epoch
##        700). The two other v12 chains are listed for the case the user picks
##        one of them (both need build_model() changes: the two-term drought
##        chain has env columns drought + drought_lag1 instead of drought2; the
##        DNN-climate chain has a network env response).
## INPUT  SOURCE_DIR/fits/<stage 1, 2, 3 folders> of the FINN-bwi study checkout (read-only)
## OUTPUT DATA_DIR/weights/v0.2.0/{model.pt, model_stage1.pt, *.csv, *.md, SHA256SUMS, README.md}
## RUN    Rscript calibration/v02/01_stage_v12_weights.R   (from the FINNmodelBWI repo root, seconds; stops while the chain is unfinished)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))
DATA_DIR   <- Filter(dir.exists, c(file.path("..", "..", "FINN-zoo-data", "bwi"), file.path("..", "FINN-zoo-data", "bwi")))[1]
SOURCE_DIR <- Filter(dir.exists, c(file.path("..", "..", "FINN-bwi"), file.path("..", "FINN-bwi")))[1]
if (is.na(SOURCE_DIR)) stop("FINN-bwi checkout not found next to this repository")

PARAMS <- list(
  release    = "v0.2.0",                                  # the next semver tag: same inputs and structure, refit on corrected data = minor bump
  ## the v12 chain with v0.1.0's structure (default candidate; the user decides)
  stage1_dir = "fits/2026-11-10-A-stage1-v12-drought2-mgmt-masked-finn0848b60",
  stage2_dir = "fits/2026-11-10-A-stage2-v12-drought2-ushape-lambda0.1-mortgrowth-masked-finn0848b60-raw",
  stage3_dir = "fits/2026-11-10-A-stage3-v12-drought2-ushape-lambda0.1-mortgrowth-masked-regcap-ep700-finn0848b60-raw",
  ## the other two v12 chains of the study (not usable with build_model() as written):
  ##   two drought terms: fits/2026-11-10-A-stage1-v12-mgmt-masked-finn0848b60 ->
  ##     fits/2026-11-10-A-stage2-v12-ushape-lambda0.1-mortgrowth-masked-finn0848b60-raw ->
  ##     fits/2026-11-10-A-stage3-v12-ushape-lambda0.1-mortgrowth-masked-regcap-ep700-finn0848b60-raw
  ##   DNN climate + linear management: fits/2026-11-10-A-stage2-v12-drought2-ushape-dnnclim-linmgmt-mortgrowth-masked-finn0848b60-raw ->
  ##     fits/2026-11-10-A-stage3-v12-drought2-ushape-dnnclim-linmgmt-mortgrowth-masked-regcap-ep700-finn0848b60-raw
  weight_files = c(
    "model.pt"                   = "stage3_dir/model.pt",
    "loss_history_stage3.csv"    = "stage3_dir/loss_history.csv",
    "regeneration_vs_stage1.csv" = "stage3_dir/regeneration_vs_stage1.csv",
    "fit_summary_stage3.md"      = "stage3_dir/fit_summary.md",
    "species_parameters.csv"     = "stage2_dir/species_parameters.csv",
    "mortality_size_curve.csv"   = "stage2_dir/mortality_size_curve.csv",
    "l2_deviation.csv"           = "stage2_dir/l2_deviation.csv",
    "loss_history_stage2.csv"    = "stage2_dir/loss_history.csv",
    "fit_summary_stage2.md"      = "stage2_dir/fit_summary.md",
    "model_stage1.pt"            = "stage1_dir/model_stage1.pt",
    "stage1_coefficients.csv"    = "stage1_dir/stage1_coefficients.csv",
    "dispersion_stage1.csv"      = "stage1_dir/dispersion_stage1.csv",
    "fit_summary_stage1.md"      = "stage1_dir/fit_summary.md")
)
missing <- Filter(function(nm) { s <- strsplit(PARAMS$weight_files[[nm]], "/", fixed = TRUE)[[1]]; !file.exists(file.path(SOURCE_DIR, PARAMS[[s[1]]], s[2])) }, names(PARAMS$weight_files))
if (length(missing)) stop("the v12 chain is not finished yet; missing: ", paste(missing, collapse = ", "))
dst <- file.path(DATA_DIR, "weights", PARAMS$release)
dir.create(dst, recursive = TRUE, showWarnings = FALSE)
for (nm in names(PARAMS$weight_files)) {
  s <- strsplit(PARAMS$weight_files[[nm]], "/", fixed = TRUE)[[1]]
  stopifnot(file.copy(file.path(SOURCE_DIR, PARAMS[[s[1]]], s[2]), file.path(dst, nm), overwrite = TRUE))
}
sums <- data.table(file = names(PARAMS$weight_files))[
  , sha256 := sapply(file.path(dst, file), digest::digest, algo = "sha256", file = TRUE)][, bytes := file.size(file.path(dst, file))]
writeLines(sprintf("%s  %s", sums$sha256, sums$file), file.path(dst, "SHA256SUMS"))
writeLines(c(sprintf("# %s", PARAMS$release), "",
             sprintf("Staged %s from the FINN-bwi fit folders `%s`, `%s`, `%s` (two-interval-v12, FINN 0.2.0 main 0848b60, weight_scale raw).",
                     format(Sys.time(), "%Y-%m-%d %H:%M"), PARAMS$stage1_dir, PARAMS$stage2_dir, PARAMS$stage3_dir), "",
             "Next: set `tag` to this release in inst/model/weights.json, run calibration/06_extract_parameters.R and calibration/05_validation.R.", "",
             "| file | bytes | sha256 |", "|---|---|---|", sprintf("| %s | %d | `%s` |", sums$file, sums$bytes, sums$sha256)), file.path(dst, "README.md"))
print(sums[, .(file, bytes, sha256 = paste0(substr(sha256, 1, 12), "..."))])
message("done: ", dst)
