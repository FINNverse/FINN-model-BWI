#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINNmodelBWI | calibration/v02/audit_residuals.R
##
## WHAT   Evidence tables for AUDIT_v02.md: the holdout residuals of the released
##        parameters (rebuilt model, v12 layout) broken down by decade,
##        management tier, start basal-area class and species: predicted /
##        observed removals, growth and recruits, change metrics (delta BA,
##        delta stems: predicted vs observed, against the persistence null) and
##        the stage-2 / stage-3 loss trends of the released fit.
## WHY    v0.2 audit (V02_BRIEF.md): residuals by species, size, interval and
##        management; change metrics next to states; convergence.
## INPUT  DATA_DIR/weights/<tag>/{parameters.pt, loss_history_stage2.csv, loss_history_stage3.csv}
##        DATA_DIR/prepared/two-interval-v12/holdout_blocked/*
## OUTPUT calibration/v02/audit/{residuals_by_decade_tier.csv, residuals_by_ba_class.csv,
##          residuals_by_species_decade.csv, change_metrics.csv, loss_trend.csv, fig_residuals.png}
## RUN    OMP_NUM_THREADS=1 Rscript calibration/v02/audit_residuals.R   (from the FINNmodelBWI repo root, ~1 min, one thread)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0 = main 0848b60
suppressPackageStartupMessages({ library(FINN); library(torch); library(data.table); library(ggplot2) })
torch_set_num_threads(1L)
source("R/definition.R")

DATA_DIR <- Filter(dir.exists, c(file.path("..", "..", "FINN-zoo-data", "bwi"), file.path("..", "FINN-zoo-data", "bwi")))[1]
OUT <- file.path("calibration", "v02", "audit"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PARAMS <- list(
  release       = jsonlite::read_json(file.path("inst", "model", "weights.json"))$variants$linear$tag,   # the parameters under audit
  prep_version  = "two-interval-v12",
  holdout_set   = "holdout_blocked",
  env_vars      = c("mat_c", "map_mm", "elev_m", "management", "drought2"),
  period        = 10L, patch_size_ha = 0.1, dbh_binsize = 2, sim_batchsize = 1000L,
  theta_sim_k   = 3.5,
  ba_breaks     = c(0, 10, 20, 30, 40, Inf),   # start basal area classes, m2/ha
  trend_window  = 100L,                        # epochs for the loss trend (last vs previous window)
  seed          = 1L
)
P <- PARAMS$period
W_DIR <- file.path(DATA_DIR, "weights", PARAMS$release)
model <- build_model("linear", torch_load(file.path(W_DIR, "parameters.pt")), theta_sim_factor = PARAMS$theta_sim_k)
d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$holdout_set)
obs <- fread(file.path(d, "obs_dt.csv")); env <- fread(file.path(d, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
envu <- fread(file.path(d, "env_unscaled_dt.csv")); init <- fread(file.path(d, "init_cohort_dt.csv"))
sid <- unique(fread(file.path(d, "siteID_dt.csv"))[, .(siteID, interval)]); species <- fread(file.path(d, "species_dt.csv"))
n_species <- nrow(species)

## ---- 1. simulate the holdout decade (v12: 10 env rows = 10 transitions) ------------------
FINN.seed(PARAMS$seed)
S <- as.data.table(simulateForest(model, env = env, init_cohort = makeInitCohorts(init, dbh_binsize = PARAMS$dbh_binsize, Nspecies = n_species),
                                  patch_size = PARAMS$patch_size_ha, batchsize = PARAMS$sim_batchsize, device = "cpu")$wide$site)
for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
st <- S[year == P, .(siteID, species, ba_m = ba, trees_m = trees)]
rt <- S[year <= P, .(growth_m = mean(growth, na.rm = TRUE), mort_m = mean(mort, na.rm = TRUE), reg_m = sum(r_mean_ha, na.rm = TRUE)), by = .(siteID, species)]
M <- merge(merge(obs[year == P], st, by = c("siteID", "species")), rt, by = c("siteID", "species"))
ini <- init[living == TRUE, .(ba_0 = sum(trees * pi * (dbh / 200)^2), trees_0 = sum(trees)), by = .(siteID, species)]
M <- merge(M, ini, by = c("siteID", "species"), all.x = TRUE)[is.na(ba_0), `:=`(ba_0 = 0, trees_0 = 0)]
M <- merge(M, sid, by = "siteID")
M <- merge(M, envu[year == 1L, .(siteID, tier = as.integer(management))], by = "siteID")
M <- merge(M, species, by = "species")
M[, `:=`(removed = fifelse(n_at_risk > 0 & is.finite(mort), n_at_risk * (1 - (1 - pmin(pmax(mort, 0), 1))^P), 0),
         removed_m = fifelse(n_at_risk > 0, n_at_risk * (1 - (1 - pmin(pmax(mort_m, 0), 1))^P), 0))]
site0 <- M[, .(ba0_ha = sum(ba_0) / PARAMS$patch_size_ha), by = siteID][, ba_class := cut(ba0_ha, PARAMS$ba_breaks, right = FALSE)]
M <- merge(M, site0[, .(siteID, ba_class)], by = "siteID")

r2 <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)
summ <- function(x, by) x[, .(
  n_site_periods = uniqueN(siteID), n_trees_at_risk = sum(n_at_risk),
  removals_pred_over_obs = sum(removed_m) / sum(removed),
  growth_pred_over_obs = sum(growth_m[is.finite(growth) & growth_n > 0]) / sum(growth[is.finite(growth) & growth_n > 0]),
  recruits_pred_over_obs = sum(reg_m) / sum(pmax(reg, 0), na.rm = TRUE),
  ba_end_pred_over_obs = sum(ba_m) / sum(ba), stems_end_pred_over_obs = sum(trees_m) / sum(trees),
  ba_r2_species = r2(ba, ba_m), ba_r2_persistence = r2(ba, ba_0), stems_r2_species = r2(trees, trees_m), stems_r2_persistence = r2(trees, trees_0)
), by = by]
fwrite(summ(M, c("interval", "tier"))[order(interval, tier)], file.path(OUT, "residuals_by_decade_tier.csv"))
fwrite(summ(M, c("ba_class"))[order(ba_class)], file.path(OUT, "residuals_by_ba_class.csv"))
fwrite(summ(M, c("species_name", "interval"))[order(species_name, interval)], file.path(OUT, "residuals_by_species_decade.csv"))

## ---- 2. change metrics: delta BA and delta stems per site-period (total) -----------------
TOT <- M[, .(ba = sum(ba), ba_m = sum(ba_m), ba_0 = sum(ba_0), trees = sum(trees), trees_m = sum(trees_m), trees_0 = sum(trees_0)), by = .(siteID, interval, tier)]
TOT[, `:=`(dba_obs = ba - ba_0, dba_pred = ba_m - ba_0, dst_obs = trees - trees_0, dst_pred = trees_m - trees_0)]
CH <- rbind(
  TOT[, .(group = "all", n = .N, dba_obs_mean = mean(dba_obs) / PARAMS$patch_size_ha, dba_pred_mean = mean(dba_pred) / PARAMS$patch_size_ha,
          dba_r2 = r2(dba_obs, dba_pred), dba_r_pearson = cor(dba_obs, dba_pred),
          dst_obs_mean = mean(dst_obs) / PARAMS$patch_size_ha, dst_pred_mean = mean(dst_pred) / PARAMS$patch_size_ha,
          dst_r2 = r2(dst_obs, dst_pred), dst_r_pearson = cor(dst_obs, dst_pred))],
  TOT[, .(group = paste0("decade ", interval, ", tier ", tier), n = .N, dba_obs_mean = mean(dba_obs) / PARAMS$patch_size_ha, dba_pred_mean = mean(dba_pred) / PARAMS$patch_size_ha,
          dba_r2 = r2(dba_obs, dba_pred), dba_r_pearson = cor(dba_obs, dba_pred),
          dst_obs_mean = mean(dst_obs) / PARAMS$patch_size_ha, dst_pred_mean = mean(dst_pred) / PARAMS$patch_size_ha,
          dst_r2 = r2(dst_obs, dst_pred), dst_r_pearson = cor(dst_obs, dst_pred)), by = .(interval, tier)][order(interval, tier)][, .(group, n, dba_obs_mean, dba_pred_mean, dba_r2, dba_r_pearson, dst_obs_mean, dst_pred_mean, dst_r2, dst_r_pearson)])
fwrite(CH, file.path(OUT, "change_metrics.csv"))

## ---- 3. loss trends of the released fit -----------------------------------------------------
h2 <- fread(file.path(W_DIR, "loss_history_stage2.csv")); h3 <- fread(file.path(W_DIR, "loss_history_stage3.csv"))
trend <- function(x, w) { n <- length(x); last <- mean(x[(n - w + 1):n]); prev <- mean(x[(n - 2 * w + 1):(n - w)]); c(last = last, previous = prev, rel_change = last / prev - 1) }
LT <- rbind(
  data.table(stage = "stage 2", loss = "growth", epochs = nrow(h2), t(trend(h2$growth, PARAMS$trend_window))),
  data.table(stage = "stage 2", loss = "mortality", epochs = nrow(h2), t(trend(h2$mort, PARAMS$trend_window))),
  data.table(stage = "stage 2 (to epoch 700, the checkpoint used)", loss = "growth", epochs = 700L, t(trend(h2$growth[1:700], PARAMS$trend_window))),
  data.table(stage = "stage 2 (to epoch 700, the checkpoint used)", loss = "mortality", epochs = 700L, t(trend(h2$mort[1:700], PARAMS$trend_window))),
  data.table(stage = "stage 3", loss = "recruitment", epochs = nrow(h3), t(trend(h3$reg, 25L))))
fwrite(LT, file.path(OUT, "loss_trend.csv"))

## ---- 4. figure: pred/obs ratios by decade and tier -------------------------------------------
R <- melt(summ(M, c("interval", "tier"))[, .(interval, tier, removals = removals_pred_over_obs, growth = growth_pred_over_obs,
                                             recruits = recruits_pred_over_obs, `basal area (end)` = ba_end_pred_over_obs, `stems (end)` = stems_end_pred_over_obs)],
          id.vars = c("interval", "tier"), variable.name = "quantity", value.name = "ratio")
g <- ggplot(R, aes(factor(tier), ratio, fill = factor(interval))) + geom_col(position = "dodge") + geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ quantity, nrow = 1) + scale_fill_brewer(palette = "Set2", name = "decade") + theme_minimal(base_size = 10) +
  labs(x = "management tier (recorded)", y = "predicted / observed", title = sprintf("Holdout, released parameters %s, v12 layout", PARAMS$release))
ggsave(file.path(OUT, "fig_residuals.png"), g, width = 11, height = 3.5, dpi = 120, bg = "white")
print(summ(M, c("interval", "tier"))[order(interval, tier), .(interval, tier, n_site_periods, removals_pred_over_obs, growth_pred_over_obs, recruits_pred_over_obs, ba_end_pred_over_obs, stems_end_pred_over_obs)])
print(CH[, .(group, n, dba_obs_mean, dba_pred_mean, dba_r2, dst_obs_mean, dst_pred_mean, dst_r2)])
print(LT)
message("done")
