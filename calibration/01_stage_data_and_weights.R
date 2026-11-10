#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINNmodelBWI | calibration/01_stage_data_and_weights.R
##
## WHAT   Stages the prepared BWI data and the fitted weights of the zoo model
##        into the zoo data folder (outside git) and verifies them against the
##        SHA-256 manifests, so that every other script of this model folder
##        reads from ONE documented location.
##
##        The zoo model is a documented copy of the FINN-bwi study pipeline
##        (branch inventory-simulator). Its prepared data are the versioned set
##        `two-interval-v8` (PREPARATION.md inside the version folder), written
##        from the open BWI base data by the public preparation route
##        `github.com/ykaeber/bwi-dataprep` (scripts 00-04, 16, 91 and
##        71_prepare_two_interval_v8.R, branch port-two-interval-v8). This
##        script does NOT rebuild the data: the raw download (~280 MB Access
##        databases + the Thuenen inventory-archive API + DWD grids + WorldClim)
##        and the preparation take hours and are documented in that repository
##        and in the staged PREPARATION.md. See the model card, section Data.
##
## WHY    The zoo convention keeps raw and prepared data and *.pt outside git in
##        FINN-zoo-data/bwi/ (raw/, prepared/<version>/<set>/, weights/). The
##        model card cites sha256 sums of the weights; this script computes and
##        writes them.
##
## INPUT  SOURCE_DIR = a checkout of FINN-bwi with
##          outputs/prepared/two-interval-v8/{train,holdout_blocked,pilot}/*.csv
##          outputs/prepared/two-interval-v8.sha256     (46-file manifest)
##          fits/<stage 1, 2, 3 folders of the candidate>  (see PARAMS)
##        or an unpacked copy of the same files from the data deposit.
## OUTPUT DATA_DIR/prepared/two-interval-v8/{PREPARATION.md, two-interval-v8.sha256,
##          split_blocks.csv, <set>/*.csv}
##        DATA_DIR/weights/<release>/{model.pt, model_stage1.pt, *.csv, *.md, SHA256SUMS}
##        DATA_DIR/raw/README.md  (where the raw data come from; nothing downloaded here)
## RUN    Rscript calibration/01_stage_data_and_weights.R   (from the FINNmodelBWI repo root, ~1 min)
##        Alternative without a FINN-bwi checkout: run calibration/run_all.R, which
##        rebuilds two-interval-v8 from the public sources and refits (hours).
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))

## zoo data folder: the sibling of the FINN-zoo checkout (worktree layout first,
## then a plain sibling checkout); no home paths
DATA_DIR <- Filter(dir.exists, c(file.path("..", "..", "FINN-zoo-data", "bwi"),
                                 file.path("..", "FINN-zoo-data", "bwi")))[1]
if (is.na(DATA_DIR)) stop("FINN-zoo-data/bwi not found next to the FINN-zoo checkout; create it first")
## the FINN-bwi study checkout (or the unpacked data deposit) next to the zoo
SOURCE_DIR <- Filter(dir.exists, c(file.path("..", "..", "FINN-bwi"), file.path("..", "FINN-bwi")))[1]
if (is.na(SOURCE_DIR)) stop("FINN-bwi (or the unpacked data deposit) not found next to the FINN-zoo checkout")

PARAMS <- list(
  prep_version = "two-interval-v12",                    # v8 without the spurious env row year 0 (PREPARATION.md): 10 env rows = 10 annual transitions per decade
  sets         = c("train", "holdout_blocked", "pilot", "national"),   # train = fitted sites, holdout_blocked = spatially blocked holdout, pilot = 2 k site-periods for smoke tests, national = every site-period (expected_preparation.csv, drought-by-tier experiments)
  release      = "v0.1.0",                              # the model repo's semver release tag (= weights folder name = GitHub release)
  ## the candidate (pending the final choice): linear, single drought index drought2, fitted recruitment cap, data v8
  stage1_dir   = "fits/2026-11-08-A-stage1-v8-drought2-mgmt-masked",
  stage2_dir   = "fits/2026-11-08-A-stage2-v8-drought2-ushape-lambda0.1-mortgrowth-masked",
  stage3_dir   = "fits/2026-11-08-A-stage3-v8-drought2-ushape-lambda0.1-mortgrowth-masked-regcap-ep700",
  weight_files = c(                                     # what the zoo needs: the simulation model + the tables the card and the fit scripts read
    "model.pt"                   = "stage3_dir/model.pt",                    # stage-2 demography + stage-3 regeneration: THE model
    "loss_history_stage3.csv"    = "stage3_dir/loss_history.csv",
    "regeneration_vs_stage1.csv" = "stage3_dir/regeneration_vs_stage1.csv",  # regeneration light, env slopes and fitted cap K per species
    "fit_summary_stage3.md"      = "stage3_dir/fit_summary.md",
    "species_parameters.csv"     = "stage2_dir/species_parameters.csv",      # competition, growth, mortality species parameters (stage 2)
    "mortality_size_curve.csv"   = "stage2_dir/mortality_size_curve.csv",    # U-shape curvature per species
    "l2_deviation.csv"           = "stage2_dir/l2_deviation.csv",            # env weights, stage 2 minus stage 1
    "loss_history_stage2.csv"    = "stage2_dir/loss_history.csv",
    "fit_summary_stage2.md"      = "stage2_dir/fit_summary.md",
    "model_stage1.pt"            = "stage1_dir/model_stage1.pt",             # the prefit (L2 anchor); needed only to refit stages 2 and 3
    "stage1_coefficients.csv"    = "stage1_dir/stage1_coefficients.csv",
    "dispersion_stage1.csv"      = "stage1_dir/dispersion_stage1.csv",
    "fit_summary_stage1.md"      = "stage1_dir/fit_summary.md")
)

## ---- 1. prepared data ---------------------------------------------------------
src_prep <- file.path(SOURCE_DIR, "outputs", "prepared", PARAMS$prep_version)
dst_prep <- file.path(DATA_DIR, "prepared", PARAMS$prep_version)
## shasum format: "<sha256>  ./<path>" (two spaces)
man_lines <- readLines(file.path(src_prep, "..", paste0(PARAMS$prep_version, ".sha256")))
manifest  <- data.table(sha256 = sub("^([0-9a-f]{64}) .*$", "\\1", man_lines),
                        path   = sub("^[0-9a-f]{64} +\\./", "", man_lines))
dir.create(dst_prep, recursive = TRUE, showWarnings = FALSE)
for (f in c("PREPARATION.md", "split_blocks.csv")) stopifnot(file.copy(file.path(src_prep, f), dst_prep, overwrite = TRUE))
stopifnot(file.copy(file.path(src_prep, "..", paste0(PARAMS$prep_version, ".sha256")), dst_prep, overwrite = TRUE))
for (s in PARAMS$sets) {
  dir.create(file.path(dst_prep, s), showWarnings = FALSE)
  files <- list.files(file.path(src_prep, s), pattern = "\\.csv$")
  stopifnot(length(files) > 0, all(file.copy(file.path(src_prep, s, files), file.path(dst_prep, s), overwrite = TRUE)))
  message(sprintf("staged %s/%s: %d files", PARAMS$prep_version, s, length(files)))
}
## verify every staged csv against the manifest
staged <- manifest[path %like% paste0("^(", paste(PARAMS$sets, collapse = "|"), ")/")]
staged[, sha_here := unname(sapply(file.path(dst_prep, path), function(p) if (file.exists(p)) digest::digest(p, algo = "sha256", file = TRUE) else NA_character_))]
bad <- staged[is.na(sha_here) | sha_here != sha256]
if (nrow(bad)) { print(bad); stop("prepared data differ from the manifest") }
message(sprintf("manifest check: %d of %d files verified (sets %s)", nrow(staged), nrow(manifest), paste(PARAMS$sets, collapse = ", ")))

## ---- 2. weights -----------------------------------------------------------------
dst_w <- file.path(DATA_DIR, "weights", PARAMS$release)
dir.create(dst_w, recursive = TRUE, showWarnings = FALSE)
for (nm in names(PARAMS$weight_files)) {
  spec <- strsplit(PARAMS$weight_files[[nm]], "/", fixed = TRUE)[[1]]
  src  <- file.path(SOURCE_DIR, PARAMS[[spec[1]]], spec[2])
  if (!file.exists(src)) stop("missing: ", src)
  stopifnot(file.copy(src, file.path(dst_w, nm), overwrite = TRUE))
}
sums <- data.table(file = names(PARAMS$weight_files))[
  , sha256 := sapply(file.path(dst_w, file), digest::digest, algo = "sha256", file = TRUE)][
  , bytes := file.size(file.path(dst_w, file))]
writeLines(sprintf("%s  %s", sums$sha256, sums$file), file.path(dst_w, "SHA256SUMS"))
writeLines(c(
  sprintf("# %s", PARAMS$release), "",
  sprintf("Staged %s from `%s` (FINN-bwi fit folders `%s`, `%s`, `%s`).", format(Sys.time(), "%Y-%m-%d %H:%M"),
          basename(normalizePath(SOURCE_DIR)), PARAMS$stage1_dir, PARAMS$stage2_dir, PARAMS$stage3_dir), "",
  "Candidate, pending the final choice of the study's model (see the model card). `model.pt` is the",
  "simulation model (stage-2 demography + stage-3 regeneration); `model_stage1.pt` and the stage-1 csv",
  "are needed only to refit stages 2 and 3 (`03_fit_stage2.R`, `04_fit_stage3.R`).", "",
  "| file | bytes | sha256 |", "|---|---|---|", sprintf("| %s | %d | `%s` |", sums$file, sums$bytes, sums$sha256)
), file.path(dst_w, "README.md"))
print(sums[, .(file, bytes, sha256 = paste0(substr(sha256, 1, 12), "..."))])

## ---- 3. where the raw data come from (nothing is downloaded here) -----------------
writeLines(c(
  "# raw/ (not staged; documentation only)", "",
  "The zoo model reads the versioned prepared data in `../prepared/two-interval-v8/`. The raw",
  "sources, their licences and the download commands are documented in the public preparation",
  "repository https://github.com/ykaeber/bwi-dataprep (`data/README.md`) and in the model card:", "",
  "- BWI 2002/2012 base data, ACCESS2003 release: https://bwi.info/Download/de/BWI-Basisdaten/ACCESS2003/",
  "  (bwi20150320_alle_daten2002.zip, ...2012.zip, ...2002-2012.zip, ...metadaten.zip, ShapeFile_Tnr_INSPIRE_Poly1000.zip;",
  "  anonymised: tract positions as 1 km INSPIRE cells only). Cite as the file Zitiervorschlag-BWI-Basisdaten-ACCESS2003-(2015).txt says.",
  "- BWI 2022 and the 2017 carbon inventory: Thuenen inventory-archive REST API (https://ci.thuenen.de/rest/v1/,",
  "  documentation https://thuenen-forest-ecosystems.github.io/TFM-Documentation/), filtered queries only; bulk use via the SQL dump.",
  "- Soil moisture (drought index): DWD Climate Data Center, 1 km monthly grids, https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/ (CC BY 4.0).",
  "- Climate normals: WorldClim 2.1, 30 arc-seconds, BIO1 and BIO12 (https://www.worldclim.org; academic / non-commercial use, no redistribution of the rasters).", "",
  "Never redistribute the field-crew name column (TruppFuehrer) or any other person-related column of the raw tables."
), file.path(DATA_DIR, "raw", "README.md"))
message("done")
