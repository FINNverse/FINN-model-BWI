#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINNmodelBWI | calibration/05_validation.R
##
## WHAT   Holdout validation of the released parameters, read-only: the model is
##        rebuilt from parameters.pt through R/definition.R::build_model() (what
##        users get from load_model()), and
##        the spatially blocked holdout (`holdout_blocked`, 8,596 site-periods
##        in 11 of 53 100-km blocks never used in any fit) is simulated over
##        its decade from the observed initial stand and scored the way the
##        study's fit() scores it; a bare-ground succession at five sites
##        across the climate space checks the long-term behaviour. Writes the
##        metrics tables and the figures the model card shows.
## WHY    Zoo convention: every card reports holdout R2 / RMSE per response
##        (species level AND stand total), the mortality AUC, an observed vs
##        predicted figure, a bare-ground succession figure and the loss
##        curves, produced by a script that loads the released weights. The
##        numbers here are the ones in metrics.csv; the study's own, larger
##        validation (train + holdout, drought experiments per management
##        tier, Harz) is in the FINN-bwi fit folder and summarised in the card.
## INPUT  DATA_DIR/weights/v0.1.0/{parameters.pt, loss_history_stage2.csv, loss_history_stage3.csv}
##        DATA_DIR/prepared/two-interval-v8/holdout_blocked/{obs_dt, env_dt, env_unscaled_dt, env_scales,
##          init_cohort_dt, species_dt}.csv
## OUTPUT calibration/validation/metrics.csv            holdout metrics, one row per response x level
##        calibration/validation/metrics_species.csv    per species: BA R2, stems R2, growth r, mortality AUC, biases
##        calibration/validation/bare_ground_equilibrium.csv
##        calibration/validation/report.md              the tables above as text + provenance
##        calibration/validation/figures/{loss_curves, holdout_obs_pred, holdout_totals, mortality_roc, bare_ground}.png
##        (metrics.csv and the figures are copied into the FINN-zoo model folder BWI-DE_staged-linear/ for the card)
## RUN    OMP_NUM_THREADS=1 Rscript calibration/05_validation.R   (from the FINNmodelBWI repo root, ~3 min, one thread)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0 = main 0848b60
suppressPackageStartupMessages({ library(FINN); library(torch); library(data.table); library(ggplot2) })
torch_set_num_threads(1L)
source("R/definition.R")                 # build_model(): the model definition this package ships

DATA_DIR  <- Filter(dir.exists, c(file.path("..", "..", "FINN-zoo-data", "bwi"), file.path("..", "FINN-zoo-data", "bwi")))[1]
MODEL_DIR <- file.path("calibration", "validation")
FIG_DIR   <- file.path(MODEL_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  release        = "v0.1.0",                # weights folder (= the model repo's release tag)
  prep_version   = "two-interval-v8",
  holdout_set    = "holdout_blocked",       # spatially blocked, never fitted
  env_vars       = c("mat_c", "map_mm", "elev_m", "management", "drought2"),   # the model's env columns
  period         = 10L,                     # every observation is a decade
  patch_size_ha  = 0.1,                     # the fitted patch area
  dbh_binsize    = 2,                       # angle-count virtual trees -> cohorts, as the fits
  sim_batchsize  = 1000L,                   # sites per simulateForest batch
  theta_sim_k    = 3.5,                     # recruit-draw dispersion = 3.5 x the fitted circle-count theta (study's validation, scripts/79)
  bare_locations = 5L,                      # sites at the 10..90 % quantiles of mean annual temperature
  bare_years     = 1000L, bare_patches = 15L,
  bare_natural   = c(management = 0, drought2 = 0),   # natural scale: no regular harvest, soil moisture at its long-term mean
  n_top_species  = 6L,                      # coloured in the succession figure; the rest grey
  seed           = 1L
)
P <- PARAMS$period

## ---- 1. weights and data -------------------------------------------------------
W_DIR <- file.path(DATA_DIR, "weights", PARAMS$release)
model <- build_model("linear", torch_load(file.path(W_DIR, "parameters.pt")), theta_sim_factor = PARAMS$theta_sim_k)
sha <- digest::digest(file.path(W_DIR, "parameters.pt"), algo = "sha256", file = TRUE)
d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$holdout_set)
obs    <- fread(file.path(d, "obs_dt.csv"))
env    <- fread(file.path(d, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
envu   <- fread(file.path(d, "env_unscaled_dt.csv"))
init   <- fread(file.path(d, "init_cohort_dt.csv"))
species <- fread(file.path(d, "species_dt.csv"))
scales  <- fread(file.path(d, "env_scales.csv"))
n_species <- nrow(species)
zs <- function(v, x) (x - scales[variable == v, mean]) / scales[variable == v, sd]
stopifnot(all(obs$period_length == P), nrow(env) == uniqueN(obs$siteID) * (P + 1L))
message(sprintf("holdout: %d site-periods, %d species, parameters.pt sha256 %s...", uniqueN(obs$siteID), n_species, substr(sha, 1, 12)))

## ---- 2. simulate the holdout decade from the observed stand ---------------------------
## each observation (year 10) is scored from the simulation as fit() scores it:
## state at sim year 11, growth and mortality as the mean annual rate over sim
## years 2..11, recruitment as their sum (env row k drives sim year k + 1)
FINN.seed(PARAMS$seed)
t0 <- Sys.time()
sim <- simulateForest(model, env = env, init_cohort = makeInitCohorts(init, dbh_binsize = PARAMS$dbh_binsize, Nspecies = n_species),
                      patch_size = PARAMS$patch_size_ha, batchsize = PARAMS$sim_batchsize, device = "cpu")
S <- as.data.table(sim$wide$site)
for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
message(sprintf("simulated in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
st <- S[year == P + 1L, .(siteID, species, dbh_m = dbh, ba_m = ba, trees_m = trees)]
rt <- S[year >= 2L & year <= P + 1L, .(growth_m = mean(growth, na.rm = TRUE), mort_m = mean(mort, na.rm = TRUE),
                                        reg_m = sum(r_mean_ha, na.rm = TRUE)), by = .(siteID, species)]
M <- merge(merge(obs[year == P], st, by = c("siteID", "species")), rt, by = c("siteID", "species"))
M <- merge(M, species, by = "species")
M[, mort_period_m := 1 - (1 - pmin(pmax(mort_m, 0), 1))^P]
## persistence reference: the decade-start stand as the prediction
ini <- init[living == TRUE, .(ba_0 = sum(trees * pi * (dbh / 200)^2), trees_0 = sum(trees)), by = .(siteID, species)]
M <- merge(M, ini, by = c("siteID", "species"), all.x = TRUE)
M[is.na(ba_0), `:=`(ba_0 = 0, trees_0 = 0)]

## ---- 3. metrics ------------------------------------------------------------------------
metrics <- function(o, p) {
  ok <- is.finite(o) & is.finite(p); o <- o[ok]; p <- p[ok]
  data.table(n = length(o), rmse = sqrt(mean((o - p)^2)), r2 = 1 - sum((o - p)^2) / sum((o - mean(o))^2),
             pearson_r = suppressWarnings(cor(o, p)), bias = mean(p - o), obs_mean = mean(o), pred_over_obs = sum(p) / sum(o))
}
## per response, species level (one row per site-period x species, as the loss scores them)
VARS <- list(ba = list("basal area (m2 / 0.1 ha)", quote(TRUE), "ba", "ba_m"),
             trees = list("stems (per 0.1 ha)", quote(TRUE), "trees", "trees_m"),
             dbh = list("mean dbh (cm)", quote(trees > 0 & trees_m > 0), "dbh", "dbh_m"),
             growth = list("growth (rel. increment / yr)", quote(is.finite(growth) & growth_n > 0), "growth", "growth_m"),
             reg = list("recruitment (per ha per decade)", quote(is.finite(reg)), "reg", "reg_m"))
MET <- rbindlist(lapply(names(VARS), function(v) { x <- VARS[[v]]; r <- M[eval(x[[2]])]
  cbind(response = v, label = x[[1]], level = "species", metrics(r[[x[[3]]]], r[[x[[4]]]])) }))
## stand totals (site-period sums over species)
TOT <- M[, .(ba = sum(ba), ba_m = sum(ba_m), trees = sum(trees), trees_m = sum(trees_m), ba_0 = sum(ba_0), trees_0 = sum(trees_0),
             reg = sum(reg, na.rm = TRUE), reg_m = sum(reg_m)), by = siteID]
MET <- rbind(MET,
             cbind(response = "ba", label = VARS$ba[[1]], level = "total", metrics(TOT$ba, TOT$ba_m)),
             cbind(response = "trees", label = VARS$trees[[1]], level = "total", metrics(TOT$trees, TOT$trees_m)),
             cbind(response = "reg", label = VARS$reg[[1]], level = "total", metrics(TOT$reg, TOT$reg_m)),
             cbind(response = "ba", label = VARS$ba[[1]], level = "species, persistence reference", metrics(M$ba, M$ba_0)),
             cbind(response = "trees", label = VARS$trees[[1]], level = "species, persistence reference", metrics(M$trees, M$trees_0)),
             cbind(response = "ba", label = VARS$ba[[1]], level = "total, persistence reference", metrics(TOT$ba, TOT$ba_0)),
             cbind(response = "trees", label = VARS$trees[[1]], level = "total, persistence reference", metrics(TOT$trees, TOT$trees_0)))
## mortality: the target is TOTAL removal (natural death + harvest; obs_dt$mort is
## the annualised removal rate of the tallied trees, n_died counts natural deaths
## only). Observed removals per row = n_at_risk * (1 - (1 - mort)^10), as the
## study's validation scores them. Tree-weighted AUC (each row `removed`
## positives and n_at_risk - removed negatives at one predicted decade removal
## probability; ties count half) and calibration (predicted / observed removals).
auc_w <- function(p, pos, neg) {
  a <- data.table(p, pos, neg)[, .(pos = sum(pos), neg = sum(neg)), by = p][order(p)]
  a[, below := cumsum(neg) - neg]
  sum(a$pos * (a$below + 0.5 * a$neg)) / (sum(a$pos) * sum(a$neg))
}
M[, removed := fifelse(n_at_risk > 0 & is.finite(mort), round(n_at_risk * (1 - (1 - pmin(pmax(mort, 0), 1))^P)), 0)]
MR <- M[n_at_risk > 0 & is.finite(mort) & is.finite(mort_period_m)][, `:=`(pos = removed, neg = n_at_risk - removed)]
auc <- auc_w(MR$mort_period_m, MR$pos, MR$neg)
calib <- sum(MR$mort_period_m * MR$n_at_risk) / sum(MR$removed)
MET <- rbind(MET, data.table(response = "mortality", label = "decade removal probability per tree", level = "tree-weighted",
                             n = sum(MR$n_at_risk), rmse = NA_real_, r2 = NA_real_, pearson_r = NA_real_, bias = NA_real_,
                             obs_mean = sum(MR$removed) / sum(MR$n_at_risk), pred_over_obs = calib))
MET[, `:=`(auc = fifelse(response == "mortality", auc, NA_real_), set = PARAMS$holdout_set)]
fwrite(MET, file.path(MODEL_DIR, "metrics.csv"))
## per species
SPM <- M[, {
  g <- is.finite(growth) & growth_n > 0; r <- n_at_risk > 0 & is.finite(mort)
  .(n_site_periods = .N, n_trees_at_risk = sum(n_at_risk),
    ba_r2 = 1 - sum((ba - ba_m)^2) / sum((ba - mean(ba))^2), ba_r2_persistence = 1 - sum((ba - ba_0)^2) / sum((ba - mean(ba))^2),
    stems_r2 = 1 - sum((trees - trees_m)^2) / sum((trees - mean(trees))^2),
    growth_r = suppressWarnings(cor(growth[g], growth_m[g])), growth_pred_over_obs = sum(growth_m[g]) / sum(growth[g]),
    mort_auc = if (sum(removed[r]) > 0 && sum(n_at_risk[r] - removed[r]) > 0) auc_w(mort_period_m[r], removed[r], n_at_risk[r] - removed[r]) else NA_real_,
    mort_pred_over_obs = sum(mort_period_m[r] * n_at_risk[r]) / sum(removed[r]),
    stems_pred_over_obs = sum(trees_m) / sum(trees), recruits_pred_over_obs = sum(reg_m) / sum(reg, na.rm = TRUE))
}, by = .(species, species_name)][order(-n_trees_at_risk)]
fwrite(SPM, file.path(MODEL_DIR, "metrics_species.csv"))

## ---- 4. succession from bare ground -------------------------------------------------
site_env <- merge(env[year == 0L, .(siteID, mat_c, map_mm, elev_m)],
                  envu[year == 0L, .(siteID, mat_nat = mat_c, map_nat = map_mm, elev_nat = elev_m)], by = "siteID")[order(mat_c)]
picks <- site_env[round(quantile(seq_len(.N), seq(0.1, 0.9, length.out = PARAMS$bare_locations)))]
picks[, `:=`(sim_site = .I, loc = sprintf("%.1f C, %.0f mm, %.0f m", mat_nat, map_nat, elev_nat))]
env_bare <- picks[, .(year = 0:PARAMS$bare_years), by = .(siteID = sim_site, mat_c, map_mm, elev_m)]
for (v in names(PARAMS$bare_natural)) env_bare[, (v) := zs(v, PARAMS$bare_natural[[v]])]
setcolorder(env_bare, c("siteID", "year", PARAMS$env_vars))
FINN.seed(PARAMS$seed)
sb <- simulateForest(model, env = env_bare, init_cohort = NULL, patches = PARAMS$bare_patches,
                     patch_size = PARAMS$patch_size_ha, device = "cpu")
SB <- as.data.table(sb$wide$site)[, .(siteID, year, species, ba, trees)]
for (cc in c("siteID", "year", "species")) SB[, (cc) := as.integer(as.character(get(cc)))]
SB <- merge(merge(SB, species, by = "species"), picks[, .(siteID = sim_site, loc)], by = "siteID")
SB[, `:=`(ba_ha = ba / PARAMS$patch_size_ha, stems_ha = trees / PARAMS$patch_size_ha, loc = factor(loc, picks$loc))]
EQ <- SB[, .(ba_last100 = mean(ba_ha[year > PARAMS$bare_years - 100]), ba_prev100 = mean(ba_ha[year > PARAMS$bare_years - 200 & year <= PARAMS$bare_years - 100])),
         by = .(loc, species_name)][
  , .(total_ba_ha = sum(ba_last100), drift_last_vs_previous_100 = sum(ba_last100) / sum(ba_prev100) - 1,
      n_species_above_1 = sum(ba_last100 > 1), dominant = species_name[which.max(ba_last100)], dominant_share = max(ba_last100) / sum(ba_last100)), by = loc]
EQ <- merge(EQ, SB[year == PARAMS$bare_years, .(stems_ha = sum(stems_ha)), by = loc], by = "loc")
fwrite(EQ, file.path(MODEL_DIR, "bare_ground_equilibrium.csv"))

## ---- 5. figures ------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 10))
## loss curves: stage 2 (baseline-scaled growth and mortality) and stage 3 (recruitment)
h2 <- fread(file.path(W_DIR, "loss_history_stage2.csv"))[, .(epoch, growth, mortality = mort)]
h3 <- fread(file.path(W_DIR, "loss_history_stage3.csv"))[, .(epoch, recruitment = reg)]
L <- rbind(melt(h2, id.vars = "epoch", variable.name = "loss")[, stage := "stage 2 (inverse calibration, 1000 epochs; stage 3 starts from epoch 700)"],
           melt(h3, id.vars = "epoch", variable.name = "loss")[, stage := "stage 3 (regeneration, 100 epochs)"])
g_loss <- ggplot(L, aes(epoch, value, colour = loss)) + geom_line(linewidth = 0.4) +
  geom_vline(data = data.table(stage = unique(L$stage)[1], x = 700), aes(xintercept = x), linetype = 2, colour = "grey40") +
  facet_wrap(~ stage, scales = "free", ncol = 1) + scale_colour_brewer(palette = "Dark2") +
  labs(x = "epoch", y = "loss (stage 2: scaled by the intercept-only baseline; stage 3: NB on stems per circle)", colour = NULL,
       title = "Training losses of the released model")
ggsave(file.path(FIG_DIR, "loss_curves.png"), g_loss, width = 8, height = 6, dpi = 120, bg = "white")

## observed vs predicted, species level, every response (density of site-period x species)
OP <- rbindlist(lapply(names(VARS), function(v) { x <- VARS[[v]]; r <- M[eval(x[[2]])]
  data.table(response = x[[1]], obs = r[[x[[3]]]], pred = r[[x[[4]]]]) }))[is.finite(obs) & is.finite(pred)]
OP[, response := factor(response, sapply(VARS, `[[`, 1))]
lab <- MET[level == "species", .(response = factor(label, sapply(VARS, `[[`, 1)), txt = sprintf("R2 %.2f, RMSE %.3g", r2, rmse))]
g_op <- ggplot(OP, aes(obs, pred)) + geom_bin2d(bins = 60) + geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.4) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt), hjust = -0.1, vjust = 1.5, size = 3) +
  scale_fill_viridis_c(trans = "log10", name = "site-periods\nx species") +
  facet_wrap(~ response, scales = "free", ncol = 3) +
  labs(x = "observed (end of decade; rates over the decade)", y = "predicted",
       title = sprintf("Holdout, species level: %d site-periods in 11 spatial blocks never fitted", uniqueN(M$siteID)))
ggsave(file.path(FIG_DIR, "holdout_obs_pred.png"), g_op, width = 10, height = 6.5, dpi = 120, bg = "white")

## stand totals per site-period
TT <- rbind(TOT[, .(response = "basal area (m2/ha)", obs = ba / PARAMS$patch_size_ha, pred = ba_m / PARAMS$patch_size_ha)],
            TOT[, .(response = "stems per ha", obs = trees / PARAMS$patch_size_ha, pred = trees_m / PARAMS$patch_size_ha)])
lab_t <- MET[level == "total" & response %in% c("ba", "trees"),
             .(response = fifelse(response == "ba", "basal area (m2/ha)", "stems per ha"), txt = sprintf("R2 %.2f (persistence %.2f)", r2,
               MET[level == "total, persistence reference"][match(response, response), r2]))]
g_tot <- ggplot(TT, aes(obs, pred)) + geom_point(alpha = 0.15, size = 0.6) + geom_abline(slope = 1, intercept = 0, colour = "red") +
  geom_text(data = lab_t, aes(x = -Inf, y = Inf, label = txt), hjust = -0.1, vjust = 1.5, size = 3.2) +
  facet_wrap(~ response, scales = "free") + labs(x = "observed, end of decade", y = "predicted", title = "Holdout, stand totals per site-period")
ggsave(file.path(FIG_DIR, "holdout_totals.png"), g_tot, width = 8, height = 4, dpi = 120, bg = "white")

## mortality ROC (tree-weighted)
roc <- MR[, .(pos = sum(pos), neg = sum(neg)), by = .(p = mort_period_m)][order(-p)][, .(fpr = cumsum(neg) / sum(neg), tpr = cumsum(pos) / sum(pos))]
g_roc <- ggplot(roc, aes(fpr, tpr)) + geom_path(colour = "#1a9e77") + geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  annotate("text", x = 0.6, y = 0.2, label = sprintf("AUC %.3f (tree-weighted)\ncalibration pred / obs %.3f", auc, calib), size = 3.3) +
  coord_equal() + labs(x = "false positive rate", y = "true positive rate", title = "Holdout: decade removal per tree at risk")
ggsave(file.path(FIG_DIR, "mortality_roc.png"), g_roc, width = 4.5, height = 4.5, dpi = 120, bg = "white")

## bare-ground succession: top species coloured, the rest grey
top <- SB[, .(mx = max(ba_ha)), by = species_name][order(-mx)][seq_len(PARAMS$n_top_species), species_name]
SB[, grp := fifelse(species_name %in% top, species_name, "other")]
g_bare <- ggplot() +
  geom_line(data = SB[grp == "other"], aes(year, ba_ha, group = species_name), colour = "grey80", linewidth = 0.3) +
  geom_line(data = SB[grp != "other"], aes(year, ba_ha, colour = species_name), linewidth = 0.5) +
  facet_wrap(~ loc, nrow = 1) + scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(x = "year", y = "basal area (m2/ha)",
       title = sprintf("Succession from bare ground (%d yr, %d patches of 0.1 ha, management tier 0, soil moisture at its long-term mean)", PARAMS$bare_years, PARAMS$bare_patches),
       subtitle = "five holdout sites at the 10-90 % quantiles of mean annual temperature; grey = the other species") +
  theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "bare_ground.png"), g_bare, width = 12, height = 4.2, dpi = 120, bg = "white")

## ---- 6. summary -------------------------------------------------------------------------
md_table <- function(dt, digits = 3) {
  f <- function(x) if (is.numeric(x)) ifelse(is.na(x), "", formatC(signif(x, digits), format = "fg", flag = "#")) else as.character(x)
  c(paste0("| ", paste(names(dt), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
    apply(as.data.frame(lapply(dt, f)), 1, function(r) paste0("| ", paste(sub("\\.$", "", r), collapse = " | "), " |")))
}
writeLines(c(
  "# Validation report, FINNmodelBWI (variant linear)", "",
  sprintf("Written %s by `calibration/05_validation.R`; FINN %s, torch %s; parameters `weights/%s/parameters.pt` sha256 `%s` rebuilt with `R/definition.R::build_model()`; data `prepared/%s/%s` (%d site-periods).",
          format(Sys.time(), "%Y-%m-%d %H:%M"), packageVersion("FINN"), packageVersion("torch"), PARAMS$release, sha, PARAMS$prep_version, PARAMS$holdout_set, uniqueN(obs$siteID)), "",
  "## Holdout metrics", "", md_table(MET[, .(response, level, n, rmse, r2, pearson_r, pred_over_obs, auc)]), "",
  "## Per species (holdout)", "", md_table(SPM[, .(species_name, n_trees_at_risk, ba_r2, ba_r2_persistence, stems_r2, growth_r, mort_auc, mort_pred_over_obs, stems_pred_over_obs, recruits_pred_over_obs)]), "",
  "## Bare-ground equilibrium (mean of the last 100 of 1000 years)", "", md_table(EQ), "",
  "## Figures", "", "| file | what |", "|---|---|",
  "| figures/loss_curves.png | stage-2 and stage-3 training losses of the released model |",
  "| figures/holdout_obs_pred.png | observed vs predicted per site-period x species, five responses |",
  "| figures/holdout_totals.png | observed vs predicted stand totals (basal area, stems) |",
  "| figures/mortality_roc.png | tree-weighted ROC of the decade removal probability |",
  "| figures/bare_ground.png | succession from bare ground at five sites |", "",
  "## Provenance", "", "| step | script | output |", "|---|---|---|",
  "| raw BWI, DWD and WorldClim layers -> two-interval-v8 | calibration/run_all.R steps 1-10 (bwi-dataprep, pinned commit) | prepared/two-interval-v8 (PREPARATION.md, PROVENANCE.md, expected_preparation.csv) |",
  "| staging + checksums | calibration/01_stage_data_and_weights.R | FINN-zoo-data/bwi/{prepared,weights} |",
  "| stage 1 prefit | calibration/02_fit_stage1.R | model_stage1.pt, stage1_coefficients.csv, dispersion_stage1.csv |",
  "| stage 2 inverse calibration | calibration/03_fit_stage2.R | epoch-700 checkpoint, species_parameters.csv, l2_deviation.csv |",
  "| stage 3 regeneration + cap | calibration/04_fit_stage3.R | model.pt |",
  "| parameter list (release asset) | calibration/06_extract_parameters.R | parameters.pt, parameters_summary.csv |",
  "| this report | calibration/05_validation.R | validation/metrics.csv, metrics_species.csv, bare_ground_equilibrium.csv, figures/ |"
), file.path(MODEL_DIR, "report.md"))
message("done")
