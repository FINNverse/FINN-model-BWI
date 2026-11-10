#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINNmodelBWI | calibration/07_expected_preparation.R
##
## WHAT   Summarises the prepared data two-interval-v8 (the output of the last
##        preparation step of run_all.R) into one committed table,
##        calibration/expected_preparation.csv: site-periods, observations,
##        tallied trees, sums of basal area and stems, recruits, removals and
##        management tiers, per set and decade and per species. A rerun of the
##        preparation elsewhere is diffed against this table; the sha256
##        manifest (two-interval-v8.sha256) is the byte-level check.
## WHY    Executable-recipe requirement of the model repositories.
## INPUT  DATA_DIR/prepared/two-interval-v8/{national,train,holdout_blocked,pilot}/*.csv
## OUTPUT calibration/expected_preparation.csv (long format: set, interval, species, measure, value)
##        When the same file exists, the script compares and reports the largest relative difference.
## RUN    Rscript calibration/07_expected_preparation.R   (from the FINNmodelBWI repo root, ~1 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))
DATA_DIR <- Filter(dir.exists, c(file.path("..", "..", "FINN-zoo-data", "bwi"), file.path("..", "FINN-zoo-data", "bwi")))[1]
PARAMS <- list(
  prep_version = "two-interval-v12",      # v8 without env row year 0; only env_rows and the drought2 moments differ from v8
  sets         = c("national", "train", "holdout_blocked", "pilot"),
  period       = 10L,
  out          = file.path("calibration", "expected_preparation.csv"),
  rel_tol      = 1e-6                     # relative difference above which a rerun is reported as different
)
summarise_set <- function(s) {
  d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, s)
  obs  <- fread(file.path(d, "obs_dt.csv")); env <- fread(file.path(d, "env_dt.csv")); envu <- fread(file.path(d, "env_unscaled_dt.csv"))
  init <- fread(file.path(d, "init_cohort_dt.csv")); tree <- fread(file.path(d, "tree_obs_dt.csv"))
  sid  <- fread(file.path(d, "siteID_dt.csv")); sp <- fread(file.path(d, "species_dt.csv"))
  iv <- unique(sid[, .(siteID, interval)])
  obs <- merge(obs, iv, by = "siteID"); tree <- merge(tree, iv, by = "siteID"); init <- merge(init, iv, by = "siteID")
  mg <- merge(envu[year == 1L, .(siteID, management)], iv, by = "siteID")   # management is constant over the decade
  per_int <- rbindlist(list(
    obs[year == PARAMS$period, .(measure = c("site_periods", "obs_rows", "sum_ba_end_m2_per_patch", "sum_trees_end", "sum_recruits_per_ha_period",
                                             "sum_n_at_risk", "sum_n_died", "mean_growth"),
                                 value = c(uniqueN(siteID), .N, sum(ba), sum(trees), sum(pmax(reg, 0), na.rm = TRUE), sum(n_at_risk), sum(n_died),
                                           mean(growth, na.rm = TRUE))), by = interval],
    tree[, .(measure = c("tallied_trees", "tallied_trees_died", "sum_tallied_stems_per_patch", "mean_tallied_dbh_cm"),
             value = c(.N, sum(died, na.rm = TRUE), sum(trees), mean(dbh))), by = interval],
    init[living == TRUE, .(measure = c("init_cohort_rows", "sum_init_stems_per_patch", "sum_init_ba_m2_per_patch"),
                           value = c(.N, sum(trees), sum(trees * pi * (dbh / 200)^2))), by = interval],
    mg[, .(measure = paste0("sites_tier_", 0:3), value = sapply(0:3, function(t) sum(management == t))), by = interval],
    merge(env, iv, by = "siteID")[, .(measure = c("env_rows", "mean_drought2_scaled", "sd_drought2_scaled"),
                                      value = c(.N, mean(drought2), sd(drought2))), by = interval]
  ))[, species := 0L]
  per_sp <- merge(obs[year == PARAMS$period, .(measure = c("sum_ba_end_m2_per_patch", "sum_trees_end", "sum_recruits_per_ha_period", "sum_n_at_risk", "sum_n_died"),
                                                 value = c(sum(ba), sum(trees), sum(pmax(reg, 0), na.rm = TRUE), sum(n_at_risk), sum(n_died))), by = .(interval, species)],
                  sp, by = "species")[, species_name := NULL]
  rbind(per_int, per_sp)[, set := s][]
}
EXP <- rbindlist(lapply(PARAMS$sets, summarise_set))
setcolorder(EXP, c("set", "interval", "species", "measure", "value"))
setorder(EXP, set, interval, species, measure)
EXP[, value := signif(value, 10)]
if (file.exists(PARAMS$out)) {
  OLD <- fread(PARAMS$out)
  CMP <- merge(OLD, EXP, by = c("set", "interval", "species", "measure"), all = TRUE, suffixes = c("_expected", "_now"))
  CMP[, rel := abs(value_now - value_expected) / pmax(abs(value_expected), 1e-12)]
  bad <- CMP[is.na(rel) | rel > PARAMS$rel_tol]
  message(sprintf("compared %d rows with the committed table: %d differ (rel > %g)", nrow(CMP), nrow(bad), PARAMS$rel_tol))
  if (nrow(bad)) print(bad[order(-rel)][1:min(20, nrow(bad))])
}
fwrite(EXP, PARAMS$out)
message(sprintf("wrote %s: %d rows; national site-periods %s", PARAMS$out, nrow(EXP),
                paste(EXP[set == "national" & measure == "site_periods", sprintf("decade %d: %d", interval, as.integer(value))], collapse = ", ")))
