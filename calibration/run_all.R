#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINNmodelBWI | calibration/run_all.R
##
## WHAT   The executable recipe: raw public sources -> prepared data
##        two-interval-v8 -> three fit stages -> parameter release asset ->
##        holdout validation -> example sites, in a fixed order with every
##        repeat spelled out. Each step is one self-contained script; this
##        file only sequences them, checks the pinned inputs and compares the
##        prepared data against the committed expectations.
## WHY    Model-repository requirement: a rerun elsewhere must need no
##        knowledge beyond this file. Nothing scientific is passed through
##        environment variables; the two preparation scripts of bwi-dataprep
##        that read switches from the environment are replaced by the literal
##        copies in calibration/prep/ (all management tiers; management as a
##        covariate; no site subsample). The only environment variables set
##        here are technical (thread count, the WorldClim folder, the data
##        folder of bwi-dataprep).
## NEEDS  - R >= 4.1 with FINN 0.2.0 (main 0848b60) in ~/R/lib-finn-0848b60, torch,
##          data.table, digest, sf, terra, FNN, ggplot2, patchwork; mdbtools on PATH
##          (`mdb-export`, for the one-time Access import); curl; bash.
##        - A checkout of github.com/ykaeber/bwi-dataprep at the pinned commit,
##          next to this repository (../../bwi-dataprep or ../bwi-dataprep).
##        - The Thuenen inventory-archive anonymous API key in
##          <bwi-dataprep>/data/raw/bwi4/.apikey (public, see data/README.md there).
##        - WorldClim 2.1 30 arc-second bio and elevation rasters (10.4 GB zip +
##          elevation; WorldClim's licence allows academic use but not
##          redistribution, so they are downloaded here, not shipped).
##        - Disk: ~3 GB raw + ~1.5 GB prepared and fits. Time: download hours
##          (API rate limits), preparation ~1-2 h, stage 1 minutes, stage 2
##          ~6 h single-threaded (1000 epochs x ~22 s), stage 3 ~30 min,
##          validation minutes.
## RUN    Rscript calibration/run_all.R   (from the FINNmodelBWI repo root)
##        Steps already done are skipped when their outputs exist (delete an
##        output to redo a step).
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))

PARAMS <- list(
  dataprep_commit = "3dc3c059f190acda92421b06fbe595cd2d1f02e6",   # ykaeber/bwi-dataprep, branch port-two-interval-v8: scripts 91 + 71 reproduce the study's v8 tables byte for byte
  dataprep_url    = "https://github.com/ykaeber/bwi-dataprep",
  access2003_base = "https://bwi.info/Download/de/BWI-Basisdaten/ACCESS2003",
  access2003_files = c("bwi20150320_alle_daten2002.zip", "bwi20150320_alle_daten2012.zip", "bwi20150320_alle_daten2002-2012.zip",
                       "bwi20150320_alle_metadaten.zip", "ShapeFile_Tnr_INSPIRE_Poly1000.zip",
                       "Zitiervorschlag-BWI-Basisdaten-ACCESS2003-(2015).txt", "_Hinweis_BWI-DatenUnterACCESS2003.pdf"),
  worldclim_base  = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base",   # wc2.1_30s_bio.zip (10.4 GB), wc2.1_30s_elev.zip
  threads         = "1",                 # technical: single-threaded torch, reproducible timing
  prep_version    = "two-interval-v8",
  expected_sha    = file.path("calibration", "two-interval-v8.sha256"),      # the study's manifest of the 46 prepared files
  expected_table  = file.path("calibration", "expected_preparation.csv"),    # summary table a rerun is diffed against (07)
  provenance_sha  = file.path("calibration", "provenance_sha256.txt")        # sha256 of every raw file as pulled for the release
)
Sys.setenv(OMP_NUM_THREADS = PARAMS$threads, MKL_NUM_THREADS = PARAMS$threads)   # technical only

## ---- locations -----------------------------------------------------------------
DATA_DIR <- Filter(dir.exists, c(file.path("..", "..", "FINN-zoo-data", "bwi"), file.path("..", "FINN-zoo-data", "bwi")))[1]
if (is.na(DATA_DIR)) { DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "bwi"); dir.create(DATA_DIR, recursive = TRUE) }
PREP <- Filter(dir.exists, c(file.path("..", "..", "bwi-dataprep"), file.path("..", "bwi-dataprep")))[1]
if (is.na(PREP)) stop(sprintf("bwi-dataprep not found next to this repository: git clone %s ../../bwi-dataprep && git -C ../../bwi-dataprep checkout %s",
                              PARAMS$dataprep_url, PARAMS$dataprep_commit))
head_sha <- system2("git", c("-C", PREP, "rev-parse", "HEAD"), stdout = TRUE)
if (!identical(head_sha, PARAMS$dataprep_commit)) stop("bwi-dataprep is at ", head_sha, ", pinned: ", PARAMS$dataprep_commit)
WORLDCLIM <- Filter(dir.exists, c(file.path(PREP, "..", "worldclim"), file.path("..", "..", "worldclim")))[1]
if (is.na(WORLDCLIM)) { WORLDCLIM <- file.path(PREP, "..", "worldclim"); dir.create(WORLDCLIM) }
Sys.setenv(FINN_WORLDCLIM_DIR = normalizePath(WORLDCLIM))                      # technical: where bwi-dataprep finds the rasters
raw <- function(...) file.path(PREP, "data", "raw", ...)
here <- normalizePath(".")

step <- function(label, outputs, expr) {
  if (length(outputs) && all(file.exists(outputs))) { message(sprintf("[skip] %-48s outputs exist", label)); return(invisible()) }
  message(sprintf("[run ] %-48s %s", label, format(Sys.time(), "%H:%M")))
  t0 <- Sys.time(); force(expr)
  if (length(outputs) && !all(file.exists(outputs))) stop(label, ": expected output missing: ", paste(outputs[!file.exists(outputs)], collapse = ", "))
  message(sprintf("[done] %-48s %.1f min", label, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
rscript <- function(script, wd = here) {           # one self-contained R script, run from its project root
  status <- system2("Rscript", shQuote(normalizePath(script, mustWork = TRUE)), wait = TRUE, env = character(),
                    stdout = "", stderr = "", invisible = FALSE)
  if (!identical(wd, here)) NULL
  if (status != 0) stop("failed: ", script)
}
in_dir <- function(wd, f) { old <- setwd(wd); on.exit(setwd(old), add = TRUE); f() }
sha_ok <- function(path, label) {                  # compare with the release's provenance list; mismatch = warn (live sources change)
  ref <- fread(PARAMS$provenance_sha, header = FALSE, sep = "\n")[[1]]
  ref <- ref[!grepl("^#", ref)]
  hit <- ref[endsWith(ref, label)]
  if (!length(hit)) return(invisible())
  now <- digest::digest(path, algo = "sha256", file = TRUE)
  if (!startsWith(hit[1], now)) warning(sprintf("%s differs from the release's copy (sha256 %s vs %s)", label, substr(now, 1, 12), substr(hit[1], 1, 12)), call. = FALSE)
}

## ---- 1. raw layers (public sources; PROVENANCE.md) --------------------------------------
step("1a ACCESS2003 release (bwi.info, ~280 MB)", raw(PARAMS$access2003_files), {
  dir.create(raw(), recursive = TRUE, showWarnings = FALSE)
  for (f in PARAMS$access2003_files) if (!file.exists(raw(f)))
    download.file(file.path(PARAMS$access2003_base, utils::URLencode(f)), raw(f), mode = "wb", quiet = TRUE)
  for (f in PARAMS$access2003_files) sha_ok(raw(f), f)
  for (z in c("bwi20150320_alle_daten2002", "bwi20150320_alle_daten2012", "bwi20150320_alle_daten2002-2012"))
    unzip(raw(paste0(z, ".zip")), exdir = raw(z))
  unzip(raw("bwi20150320_alle_metadaten.zip"), exdir = raw("metadaten"))
  unzip(raw("ShapeFile_Tnr_INSPIRE_Poly1000.zip"), exdir = raw("inspire"))
})
step("1b Thuenen inventory archive via API (BWI 2022, CI 2017)", raw("bwi4", c("plot.csv", "tree.csv", "regeneration.csv", "lookup_tree_species.csv")), {
  if (!file.exists(raw("bwi4", ".apikey"))) stop("put Thuenen's public anonymous API key into ", raw("bwi4", ".apikey"), " (see bwi-dataprep/data/README.md)")
  in_dir(PREP, function() { stopifnot(system2("bash", "scripts/fetch_bwi4.sh") == 0, system2("bash", "scripts/fetch_bwi4_tree.sh") == 0) })
  for (f in list.files(raw("bwi4"), pattern = "\\.csv$")) sha_ok(raw("bwi4", f), file.path("bwi4", f))
})
step("1c DWD monthly 1 km grids 2000-2022 (~240 MB, 1104 files)", raw("dwd_drought", "grids_germany_monthly_soil_moist_202212.asc.gz"), {
  in_dir(PREP, function() stopifnot(system2("bash", c("scripts/fetch_dwd_drought.sh", "2000", "2022")) == 0))
})
step("1d WorldClim 2.1 30s bio + elevation (10.4 GB + 0.3 GB)", file.path(WORLDCLIM, c("wc2.1_30s_bio/wc2.1_30s_bio_1.tif", "wc2.1_30s_bio/wc2.1_30s_bio_12.tif", "wc2.1_30s_elev.tif")), {
  for (z in c("wc2.1_30s_bio.zip", "wc2.1_30s_elev.zip")) if (!file.exists(file.path(WORLDCLIM, z)))
    download.file(file.path(PARAMS$worldclim_base, z), file.path(WORLDCLIM, z), mode = "wb", quiet = TRUE)
  unzip(file.path(WORLDCLIM, "wc2.1_30s_bio.zip"), exdir = file.path(WORLDCLIM, "wc2.1_30s_bio"))
  unzip(file.path(WORLDCLIM, "wc2.1_30s_elev.zip"), exdir = WORLDCLIM)
})

## ---- 2. preparation (bwi-dataprep at the pinned commit; run from ITS root) -------------------
## Known traps, all handled by the order below:
##  - scripts/02 and scripts/04 read their switches from environment variables in
##    bwi-dataprep; the copies in calibration/prep/ carry them as literal PARAMS
##    (all four management tiers M0-M3; management as an env covariate; every site).
##  - scripts/04 runs TWICE: once before scripts/16 (it writes outputs/finn/species_dt.csv,
##    which scripts/40 and 91 need for the species lumping, and outputs/finn without a
##    drought column), and again after scripts/16 so that outputs/finn carries the
##    per-year drought anomaly that scripts/70 reads.
##  - Field dates (`Datum_feld` in b2_tab_w / b3_tab_w) come out of mdbtools as
##    "MM/DD/YY HH:MM:SS"; R/bwi_trees.R parses the two-digit year from that
##    format. A different mdbtools version (ISO dates) silently gives NA years and
##    nominal 10-year periods everywhere: the check after step 2b catches it.
prep_script <- function(s) rscript(file.path(PREP, "scripts", s), wd = PREP)
local_prep  <- function(s) rscript(file.path("calibration", "prep", s), wd = PREP)
run_prep <- function(label, outputs, scripts, local = FALSE) step(label, file.path(PREP, outputs), in_dir(PREP, function()
  for (s in scripts) if (local) rscript(normalizePath(file.path(here, "calibration", "prep", s))) else rscript(file.path("scripts", s))))
run_prep("2a 00 Access -> csv (mdbtools)", "data/bwi_csv/b3_baeume.csv", "00_import_bwi.R")
run_prep("2b 01 corner metadata", "outputs/corner_metadata.csv", "01_prepare_metadata.R")
step("2b' field-date check (Datum_feld parsing)", character(0), {
  cm <- fread(file.path(PREP, "outputs", "corner_metadata.csv"))
  pl <- cm$period_length[is.finite(cm$period_length)]
  if (median(pl) != 10 || mean(pl == 10) > 0.95)   # the realised interval spreads over 8-12 years; all-10 means the dates were not parsed
    stop("period_length degenerate (", round(100 * mean(pl == 10)), " % exactly 10): Datum_feld not parsed, check the mdbtools date format")
})
run_prep("2c 02 calibration set, ALL tiers (literal copy)", "outputs/calibration_set_corners.csv", "02_select_calibration_set.R", local = TRUE)
run_prep("2d 03 WorldClim + elevation + INSPIRE coordinates", "outputs/corner_env.csv", "03_attach_environment.R")
run_prep("2e 04 FINN inputs, first pass (species list; no drought yet)", "outputs/finn/species_dt.csv", "04_prepare_finn_inputs.R", local = TRUE)
run_prep("2f 16 per-year drought anomaly (DWD soil moisture)", "outputs/site_drought_annual.csv", "16_drought_index.R")
step("2g 04 FINN inputs, SECOND pass (with drought)", character(0), in_dir(PREP, function() {
  rscript(normalizePath(file.path(here, "calibration", "prep", "04_prepare_finn_inputs.R")))
  stopifnot("drought" %in% names(fread(file.path("outputs", "finn", "env_dt.csv"), nrows = 1)))
}))
run_prep("2h 40 BWI3->BWI4 tree panel (archive schema)", "outputs/new_schema/panel.csv", "40_import_new_schema.R")
run_prep("2i 70 two-interval-v3 (its split_blocks.csv is reused by v8)", "outputs/prepared/two-interval-v3/split_blocks.csv", "70_prepare_two_interval.R")
run_prep("2j 91 decade-2 panel with the fates fixed", "outputs/new_schema_chained/panel.csv", "91_import_new_schema_chained.R")
run_prep("2k 71 two-interval-v8", file.path("outputs", "prepared", PARAMS$prep_version, "PREPARATION.md"), "71_prepare_two_interval_v8.R")
step("2l v8 against the study's manifest", character(0), {
  man <- readLines(PARAMS$expected_sha)
  ref <- data.table(sha = sub("^([0-9a-f]{64}) .*$", "\\1", man), path = sub("^[0-9a-f]{64} +\\./", "", man))
  ref[, now := sapply(file.path(PREP, "outputs", "prepared", PARAMS$prep_version, path), function(p) if (file.exists(p)) digest::digest(p, algo = "sha256", file = TRUE) else NA_character_)]
  bad <- ref[is.na(now) | now != sha]
  message(sprintf("%d of %d prepared files identical to the study's", nrow(ref) - nrow(bad), nrow(ref)))
  if (nrow(bad)) message("differ: ", paste(bad$path, collapse = ", "), " (split_blocks.csv differs by design, see bwi-dataprep README)")
})

## ---- 3. stage, summarise, fit, release, validate (this repository's scripts) -------------------
step("3a 01 stage prepared data (+ weights if a study checkout is present)", file.path(DATA_DIR, "prepared", PARAMS$prep_version, "train", "obs_dt.csv"),
     rscript(file.path("calibration", "01_stage_data_and_weights.R")))
step("3b 07 expected_preparation.csv (compares with the committed table)", character(0), rscript(file.path("calibration", "07_expected_preparation.R")))
step("3c 02 stage 1 prefit (minutes)", file.path(DATA_DIR, "fits", "stage1", "model_stage1.pt"), rscript(file.path("calibration", "02_fit_stage1.R")))
step("3d 03 stage 2 inverse calibration (~6 h single-threaded)", file.path(DATA_DIR, "fits", "stage2", "model.pt"), rscript(file.path("calibration", "03_fit_stage2.R")))
step("3e 04 stage 3 regeneration + cap (~30 min)", file.path(DATA_DIR, "fits", "stage3", "model.pt"), rscript(file.path("calibration", "04_fit_stage3.R")))
step("3f 06 parameter list = release asset", file.path(DATA_DIR, "weights", "v0.1.0", "parameters.pt"), rscript(file.path("calibration", "06_extract_parameters.R")))
step("3g 05 holdout validation", file.path("calibration", "validation", "metrics.csv"), rscript(file.path("calibration", "05_validation.R")))
step("3h 09 real example stands (local, not shipped)", file.path(DATA_DIR, "example_sites", "example_env.csv"), rscript(file.path("calibration", "09_example_sites.R")))
message("all steps done ", format(Sys.time(), "%Y-%m-%d %H:%M"))
