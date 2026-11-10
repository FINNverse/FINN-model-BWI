#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINNmodelBWI | calibration/02_fit_stage1.R
##
## ZOO COPY of the FINN-bwi study script
##   fits/2026-11-08-A-stage1-v8-drought2-mgmt-masked/fit.R  (branch inventory-simulator)
## with these changes only: FINN 0.2.0 (main 0848b60) in place of the release
## candidate 6705a50 (stage 1 reproduces bit for bit: identical
## stage1_coefficients.csv and dispersion_stage1.csv sha256 under both), data
## read from the zoo data folder, outputs written to DATA_DIR/fits/stage1.
## Runtime: minutes on a laptop (L-BFGS, full batch). A refit is NOT needed to
## use the model: the released weights already contain the result.
##
## WHAT   Stage 1 ALONE: the environmental response of growth, mortality and
##        regeneration, fitted per process against its own observations with
##        FINN's OWN process functions (FINN::competition, FINN::growth,
##        FINN::mortality, FINN::regeneration), every species sharing the same
##        allometry and light parameters (FINN's default centre values).
##          growth        each surviving tallied tree's mean annual relative
##                        increment against the mean of FINN::growth over the
##                        decade's years (MSE, as fit()'s growth loss)
##          mortality     each tallied tree's 0/1 removal over the decade
##                        against 1 - prod_t (1 - FINN::mortality) (binomial,
##                        as fit(loss_aggregation = c(mortality = "tree")))
##          regeneration  each site x species' recruits on the regeneration
##                        circle against the decade sum of FINN::regeneration,
##                        negative binomial with a per-species theta
##        Light comes from the OBSERVED stand at the decade start through
##        FINN::competition (the cohorts are the tallied trees with their
##        expansion weights; regeneration sees the light at 1 m), held over the
##        decade while the environment moves year by year: exactly the
##        per-tree path of FINN's fit() and predict_trees(), which section 4
##        checks numerically. Growth is fitted first; mortality then gets the
##        predicted growth, as in FINN's annual step.
## WHY    Stage 1 estimates the climate, management and drought responses
##        directly from the decadal observations; stage 2 pulls its env weights
##        towards these (L2 without the intercept). Env linear on climate,
##        management (tier 0-3, constant over the decade) and the single drought
##        index drought2 (mean of this and the previous year's soil-moisture
##        anomaly); no x, y; no sign constraint.
## SPECIES PARAMETERS (all species identical, FINN's init_species_parameters
##        centre values; nothing species-specific beyond the env response):
##          competition  height 0.5, strength 0.2
##          growth       light 0.5, size 0.05
##          mortality    light 0, size 0, growth 0
##          regeneration light 0.5
## DATA   DATA_DIR/prepared/two-interval-v8/train: tree_obs_dt.csv (tallied trees),
##        obs_dt.csv (regeneration and the annualised rates), env_dt.csv
## RUN    Rscript calibration/02_fit_stage1.R   (from the FINNmodelBWI repo root)
## OUTPUT DATA_DIR/fits/stage1: stage1_coefficients.csv, dispersion_stage1.csv,
##        model_stage1.pt (a FINN model holding the stage-1 parameters, all
##        frozen), fit_summary.md
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0 = main 0848b60
suppressPackageStartupMessages({
  library(FINN)
  library(data.table)
  library(torch)
})

DATA_DIR <- Filter(dir.exists, c(file.path("..", "..", "FINN-zoo-data", "bwi"), file.path("..", "FINN-zoo-data", "bwi")))[1]
FIT_DIR  <- file.path(DATA_DIR, "fits", "stage1")
dir.create(FIT_DIR, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  prep_version  = "two-interval-v12",      # v8 without env row year 0 (PREPARATION_v12.md): the per-tree window is env years 1..10 in both layouts
  train_set     = "train",
  env_vars      = c("mat_c", "map_mm", "elev_m", "management", "drought2"),   # single drought index drought2; the drought term must stay last
  mort_target   = "tree",                  # per-tree 0/1 over the decade (the final model's target)
  period        = 10L,                     # years per decade; the per-tree window is env years 1..10
  patch_size_ha = 0.1,                     # the expansion's patch area (PREPARATION.md)
  recruits_dbh  = 7,                       # BWI caliper threshold (cm)
  ## the common species parameters: FINN's init_species_parameters() centre values
  par_competition = c(0.5, 0.2),
  par_growth      = c(0.5, 0.05),
  par_mortality   = c(0, 0, 0),            # (light, size, growth): growth term 0, so growth_scale does not enter here
  par_regeneration = 0.5,
  n_quantiles   = 10L,                     # FINN's default height-quantile competition
  light_mask    = TRUE,                    # 0.2.0's masked competition: padding slots no longer set the height-quantile edges
  reg_floor     = 1e-3,                    # FINN::regeneration additive floor (recruits/ha/yr)
  ## species with too little information keep the intercept only
  growth_min_n  = 10L,                     # surviving trees with a growth value
  mort_min_nonzero = 20L,                  # site-periods with a removal
  reg_min_nonzero  = 20L,                  # site-periods with recruits
  ## the regeneration circle: obs * w_obs = stems on the circle
  reg_circle_r_m = 1.75, reg_incr_cm_yr = 0.30, reg_band_cm = 1.9,
  theta_floor   = 0.01,                    # FINN's p-clamp bites above 1e4 * theta
  lbfgs_iter    = 200L,                    # L-BFGS iterations per process (full batch, strong Wolfe)
  seed          = 123L
)
torch_manual_seed(PARAMS$seed)
P <- PARAMS$period
stopifnot(identical(formals(FINN::mortality)$growth_scale, 0.01))   # growth in percent (0.2.0 default)
stopifnot("light_mask" %in% names(formals(FINN::competition)), "light_mask" %in% names(formals(createProcess)))   # a FINN with masked competition (>= 0.2.0)

## ---- 1. data -----------------------------------------------------------------
train_dir <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs  <- fread(file.path(train_dir, "obs_dt.csv"))
env  <- fread(file.path(train_dir, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
tree <- fread(file.path(train_dir, "tree_obs_dt.csv"))
species_dt <- fread(file.path(train_dir, "species_dt.csv"))
n_species <- max(obs$species)
sites <- sort(unique(obs$siteID)); S <- length(sites)
message(sprintf("train: %d site-periods, %d species, %d tallied trees (%d with a fate, %d with growth)",
                S, n_species, nrow(tree), sum(!is.na(tree$died)), sum(!is.na(tree$growth))))

## the observed stand as padded [sites, 1, trees] tensors (empty sites: no trees)
tree[, site_i := match(siteID, sites)]
setorder(tree, site_i, -dbh)
tree[, slot := seq_len(.N), by = site_i]
Tm <- max(tree$slot)
pad <- function(col, fill) {
  a <- array(fill, c(S, 1L, Tm)); a[cbind(tree$site_i, 1L, tree$slot)] <- tree[[col]]; a
}
dbh_t     <- torch_tensor(pad("dbh", 0), dtype = torch_float())
trees_t   <- torch_tensor(pad("trees", 0), dtype = torch_float())
species_t <- torch_tensor(pad("species", 1L), dtype = torch_long())
growth_o  <- torch_tensor(pad("growth", NaN), dtype = torch_float())
died_o    <- torch_tensor(pad("died", NaN), dtype = torch_float())

## the environment of the per-tree window (env years 1..10) as [sites, years, 1 + vars]
E <- env[year >= 1L & year <= P][order(match(siteID, sites), year)]
stopifnot(nrow(E) == S * P, identical(unique(E$siteID), sites))
X <- array(as.matrix(E[, PARAMS$env_vars, with = FALSE]), c(P, S, length(PARAMS$env_vars)))
X <- aperm(X, c(2, 1, 3))
X_t <- torch_cat(list(torch_ones(S, P, 1), torch_tensor(X, dtype = torch_float())), dim = 3)
K <- X_t$shape[3]

## ---- 2. FINN's process functions and the common species parameters -------------
## The exported process functions read a few fields of the model they normally
## run inside (`self`); outside a model these are the defaults of finn().
in_model <- function(f) {
  environment(f) <- list2env(list(self = list(record_raws = FALSE, reg_floor = PARAMS$reg_floor)),
                             parent = environment(f))
  f
}
f_growth <- in_model(FINN::growth); f_mort <- in_model(FINN::mortality); f_reg <- in_model(FINN::regeneration)
common <- function(v) matrix(v, n_species, length(v), byrow = TRUE)
parComp <- torch_tensor(common(PARAMS$par_competition), dtype = torch_float())
parGrow <- torch_tensor(common(PARAMS$par_growth), dtype = torch_float())
parMort <- torch_tensor(common(PARAMS$par_mortality), dtype = torch_float())
parReg  <- torch_tensor(rep(PARAMS$par_regeneration, n_species), dtype = torch_float())

light_t <- FINN::competition(dbh = dbh_t, species = species_t, trees = trees_t, parComp = parComp,
                             h = NULL, patch_size_ha = PARAMS$patch_size_ha, n_quantiles = PARAMS$n_quantiles,
                             light_mask = PARAMS$light_mask)
light_reg <- FINN::competition(dbh = dbh_t, species = species_t, trees = trees_t, parComp = parComp,
                               h = 1, patch_size_ha = PARAMS$patch_size_ha, n_quantiles = PARAMS$n_quantiles,
                               light_mask = PARAMS$light_mask)

## linear env response: pred[s, t, tree] = [1, x_st] . W[species of the tree]
pred_trees <- function(W) {
  Pst <- torch_matmul(X_t, W$t())                                   # [S, P, species]
  idx <- species_t$expand(c(S, P, Tm))
  torch_gather(Pst, 3, idx)                                         # [S, P, trees]
}
growth_rates <- function(Wg) {                                       # annual rates [S, P, trees]
  pg <- pred_trees(Wg)
  torch_cat(lapply(seq_len(P), function(t)
    f_growth(dbh = dbh_t, species = species_t, parGrowth = parGrow,
             pred = pg[, t, ]$unsqueeze(2), light = light_t)), dim = 2)
}
mort_rates <- function(Wm, g) {                                      # annual probabilities [S, P, trees]
  pm <- pred_trees(Wm)
  torch_cat(lapply(seq_len(P), function(t)
    f_mort(dbh = dbh_t, species = species_t, trees = trees_t, parMort = parMort,
           pred = pm[, t, ]$unsqueeze(2), light = light_t, growth = g[, t, ]$unsqueeze(2))), dim = 2)
}
reg_sum <- function(Wr) {                                            # recruits/ha over the decade [S, species]
  Pst <- torch_matmul(X_t, Wr$t())                                  # [S, P, species]
  sp_all <- torch_arange(1, n_species, dtype = torch_long())$reshape(c(1, 1, -1))$expand(c(S, 1, n_species))
  Reduce(`+`, lapply(seq_len(P), function(t)
    f_reg(species = sp_all, parReg = parReg, pred = Pst[, t, ], light = light_reg)[, 1, ]))
}

## ---- 3. stage-1 fits, one process at a time -----------------------------------
## slopes of species with too little information are held at 0 (intercept only)
slope_mask <- function(keep) {
  m <- matrix(1, n_species, K); m[!keep, -1] <- 0; torch_tensor(m, dtype = torch_float())
}
fit_lbfgs <- function(W0, mask, loss_fn, extra = NULL) {
  W <- torch_tensor(W0, dtype = torch_float(), requires_grad = TRUE)
  params <- c(list(W), extra)
  opt <- optim_lbfgs(params, lr = 1, max_iter = PARAMS$lbfgs_iter, line_search_fn = "strong_wolfe",
                     tolerance_grad = 1e-9, tolerance_change = 1e-12, history_size = 20)
  closure <- function() { opt$zero_grad(); l <- loss_fn(W * mask); l$backward(); l }
  l0 <- as.numeric(loss_fn(W * mask)$item())
  opt$step(closure)
  list(W = as.matrix((W * mask)$detach()), loss0 = l0, loss1 = as.numeric(loss_fn(W * mask)$item()))
}
nanmask <- function(t) t$isnan()$logical_not()

## growth: MSE of the mean annual rate over the decade, survivors only
gm <- nanmask(growth_o)
g_scale <- as.numeric(growth_o[gm]$var()$item())
n_growth <- tree[!is.na(growth), .N, by = species][match(seq_len(n_species), species), N]
n_growth[is.na(n_growth)] <- 0L
Wg0 <- matrix(0, n_species, K); Wg0[, 1] <- log(mean(tree$growth, na.rm = TRUE) / 0.5)
fg <- fit_lbfgs(Wg0, slope_mask(n_growth >= PARAMS$growth_min_n), function(W) {
  gbar <- growth_rates(W)$mean(2)$unsqueeze(2)
  ((gbar[gm] - growth_o[gm])^2)$mean() / g_scale
})
Wg <- fg$W
G_pred <- growth_rates(torch_tensor(Wg, dtype = torch_float()))$detach()

## mortality: per tree over the decade
dm <- nanmask(died_o)
nz_mort <- obs[year == P & n_at_risk > 0 & is.finite(mort) & mort > 0, .N, by = species][match(seq_len(n_species), species), N]
nz_mort[is.na(nz_mort)] <- 0L
a0 <- qlogis(1 - (1 - mean(tree$died, na.rm = TRUE))^(1 / P))
Wm0 <- matrix(0, n_species, K); Wm0[, 1] <- a0
stopifnot(PARAMS$mort_target == "tree")
mort_loss <- function(W) {
  m <- mort_rates(W, G_pred)$clamp(1e-6, 1 - 1e-6)
  p <- (1 - (1 - m)$log()$sum(2)$exp())$unsqueeze(2)$clamp(1e-7, 1 - 1e-7)
  y <- died_o[dm]; pp <- p[dm]
  -(y * pp$log() + (1 - y) * (1 - pp)$log())$mean()
}
fm <- fit_lbfgs(Wm0, slope_mask(nz_mort >= PARAMS$mort_min_nonzero), mort_loss)
Wm <- fm$W

## regeneration: negative binomial on stems per circle, theta per species
w_obs <- (pi * PARAMS$reg_circle_r_m^2 / 1e4) / (PARAMS$reg_incr_cm_yr / PARAMS$reg_band_cm * P)
R <- dcast(obs[year == P, .(siteID, species, y = pmax(reg, 0) * w_obs)], siteID ~ species, value.var = "y")
stopifnot(identical(R$siteID, sites))
Y_t <- torch_tensor(as.matrix(R[, -1]), dtype = torch_float())
stopifnot(max(abs(as.matrix(R[, -1]) - round(as.matrix(R[, -1]))), na.rm = TRUE) < 1e-6)
ym <- nanmask(Y_t)
nz_reg <- colSums(as.matrix(R[, -1]) > 0, na.rm = TRUE)
Wr0 <- matrix(0, n_species, K)
Wr0[, 1] <- log(pmax(colMeans(as.matrix(R[, -1]), na.rm = TRUE) / w_obs / P, PARAMS$reg_floor))
log_theta <- torch_zeros(n_species, requires_grad = TRUE)
nb_nll <- function(y, mu, theta) -(torch_lgamma(y + theta) - torch_lgamma(theta) - torch_lgamma(y + 1) +
                                   theta * (theta / (theta + mu))$log() + y * (mu / (theta + mu))$log())
fr <- fit_lbfgs(Wr0, slope_mask(nz_reg >= PARAMS$reg_min_nonzero), function(W) {
  mu <- reg_sum(W) * w_obs + 0.001
  theta <- (log_theta$exp() + PARAMS$theta_floor)$unsqueeze(1)$expand_as(mu)
  nb_nll(Y_t[ym], mu[ym], theta[ym])$mean()
}, extra = list(log_theta))
Wr <- fr$W
theta_s1 <- as.numeric((log_theta$exp() + PARAMS$theta_floor)$detach())

coef_dt <- rbindlist(lapply(c("growth", "mortality", "regeneration"), function(p) {
  W <- switch(p, growth = Wg, mortality = Wm, regeneration = Wr)
  data.table(process = p, species = seq_len(n_species),
             intercept_only = switch(p, growth = n_growth < PARAMS$growth_min_n,
                                     mortality = nz_mort < PARAMS$mort_min_nonzero, regeneration = nz_reg < PARAMS$reg_min_nonzero),
             `colnames<-`(as.data.table(W), c("intercept", PARAMS$env_vars)))
}))
coef_dt <- merge(coef_dt, species_dt, by = "species")
fwrite(coef_dt[order(process, species)], file.path(FIT_DIR, "stage1_coefficients.csv"))
fwrite(data.table(species = seq_len(n_species), n_nonzero = nz_reg, theta_stage1 = theta_s1),
       file.path(FIT_DIR, "dispersion_stage1.csv"))
message("stage-1 mortality drought slopes: ",
        paste(sprintf("%s %.2f", species_dt$species_name[order(species_dt$species)], Wm[, K]), collapse = ", "))
message("stage-1 mortality management slopes: ",
        paste(sprintf("%s %.2f", species_dt$species_name[order(species_dt$species)],
                      Wm[, 1L + match("management", PARAMS$env_vars)]), collapse = ", "))

## ---- 4. the stage-1 model as a FINN model, checked against FINN's own tree path ----
env_formula <- as.formula(paste("~", paste(PARAMS$env_vars, collapse = " + ")))
model <- finn(
  N_species    = n_species,
  recruits_dbh = PARAMS$recruits_dbh,
  reg_floor    = PARAMS$reg_floor,
  recruit_obs_weight = w_obs,
  competition_process  = createProcess(~0, FINN::competition, initSpecies = common(PARAMS$par_competition),
                                       optimizeSpecies = FALSE, n_quantiles = PARAMS$n_quantiles,
                                       light_mask = PARAMS$light_mask),
  growth_process       = createProcess(env_formula, FINN::growth, initSpecies = common(PARAMS$par_growth),
                                       initEnv = list(Wg), optimizeSpecies = FALSE, optimizeEnv = FALSE),
  mortality_process    = createProcess(env_formula, FINN::mortality, initSpecies = common(PARAMS$par_mortality),
                                       initEnv = list(Wm), optimizeSpecies = FALSE, optimizeEnv = FALSE),
  regeneration_process = createProcess(env_formula, FINN::regeneration,
                                       initSpecies = matrix(PARAMS$par_regeneration, n_species, 1L),
                                       initEnv = list(Wr), dispersion_parameter = theta_s1,
                                       optimizeSpecies = FALSE, optimizeEnv = FALSE)
)
torch_save(model, file.path(FIT_DIR, "model_stage1.pt"))

## FINN's predict_trees() is the code path of fit()'s per-tree loss: it must
## reproduce this script's growth and interval mortality, on ALL sites in one batch.
chk_i <- seq_len(S)
chk <- tree[site_i %in% chk_i][, period_length := P]
pt <- model$predict_trees(chk[, .(siteID, year, patchID, species, dbh, trees, period_length)],
                          env = env[siteID %in% sites[chk_i]], patches = 1L, patch_size = PARAMS$patch_size_ha)
m_own <- (1 - (1 - mort_rates(torch_tensor(Wm, dtype = torch_float()), G_pred)$clamp(1e-6, 1 - 1e-6))$log()$sum(2)$exp())
own <- data.table(site_i = chk$site_i, slot = chk$slot)
own[, `:=`(g_own = as.array(G_pred$mean(2))[cbind(site_i, slot)], m_own = as.array(m_own)[cbind(site_i, slot)])]
check <- c(growth = max(abs(pt$growth_pred - own$g_own)), mortality = max(abs(pt$mort_pred - own$m_own)))
message(sprintf("vs FINN predict_trees on %d sites: max |diff| growth %.2g, interval mortality %.2g",
                length(chk_i), check[["growth"]], check[["mortality"]]))
stopifnot(all(check < 1e-4))

## ---- 5. what was fitted --------------------------------------------------------
slopes <- function(W, v) paste(sprintf("%s %.2f", species_dt$species_name[order(species_dt$species)],
                                       W[, 1L + match(v, PARAMS$env_vars)]), collapse = ", ")
writeLines(c(
  "# stage 1 (prefit), FINNmodelBWI", "",
  sprintf("Fitted %s by `calibration/02_fit_stage1.R` (FINN %s, torch %s).", format(Sys.time(), "%Y-%m-%d %H:%M"),
          packageVersion("FINN"), packageVersion("torch")), "",
  sprintf("Data: `prepared/%s/%s`, %d site-periods, %d tallied trees. Env `%s`, linear. Mortality target: %s.",
          PARAMS$prep_version, PARAMS$train_set, S, nrow(tree), paste(deparse(env_formula), collapse = ""), PARAMS$mort_target), "",
  sprintf("Common species parameters (FINN's default centres): competition (0.5, 0.2), growth (0.5, 0.05), mortality (0, 0, 0), regeneration 0.5. Competition light_mask = %s.", PARAMS$light_mask), "",
  "| process | loss start | loss end | intercept only |", "|---|---|---|---|",
  sprintf("| growth (MSE / var) | %.4f | %.4f | %s |", fg$loss0, fg$loss1, paste(species_dt$species_name[n_growth < PARAMS$growth_min_n], collapse = ", ")),
  sprintf("| mortality (%s) | %.4f | %.4f | %s |", PARAMS$mort_target, fm$loss0, fm$loss1, paste(species_dt$species_name[nz_mort < PARAMS$mort_min_nonzero], collapse = ", ")),
  sprintf("| regeneration (NB per circle) | %.4f | %.4f | %s |", fr$loss0, fr$loss1, paste(species_dt$species_name[nz_reg < PARAMS$reg_min_nonzero], collapse = ", ")), "",
  sprintf("Mortality drought2 slopes per pooled sd: %s.", slopes(Wm, "drought2")), "",
  sprintf("Growth drought2 slopes per pooled sd: %s.", slopes(Wg, "drought2")), "",
  sprintf("Mortality management slopes (per tier): %s.", slopes(Wm, "management")), "",
  sprintf("Theta (circle counts): median %.3g, range %.3g-%.3g.", median(theta_s1), min(theta_s1), max(theta_s1)), "",
  sprintf("Check against FINN's predict_trees() on %d sites: max |diff| growth %.2g, interval mortality %.2g.",
          length(chk_i), check[["growth"]], check[["mortality"]])
), file.path(FIT_DIR, "fit_summary.md"))
message("done")
