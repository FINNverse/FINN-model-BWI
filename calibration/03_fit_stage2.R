#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINNmodelBWI | calibration/03_fit_stage2.R
##
## ZOO COPY of the FINN-bwi study script
##   fits/2026-11-08-A-stage2-v8-drought2-ushape-lambda0.1-mortgrowth-masked/fit.R
## (branch inventory-simulator) with these changes only:
##   - FINN 0.2.0 (main 0848b60) in place of the release candidate 6705a50.
##     Under 0.2.0 a numeric `weights` vector multiplies the BASELINE-SCALED
##     losses (FINN#37); the study's weights were built for raw losses, so
##     fit(weight_scale = "raw") restores the objective of the released model
##     (the study's own 0848b60 refits do the same). The L2 optimizer builds its
##     tensors on the parameters' device (inert on CPU, required on a GPU).
##   - data read from the zoo data folder; the prefit (L2 anchor) read from the
##     released weights folder; outputs to DATA_DIR/fits/stage2.
## Runtime: 1000 epochs x ~22 s = ~6 h single-threaded on a laptop (Apple M3);
## ~4 h with 4 threads. A refit is NOT needed to use the model: the released
## weights contain the result (epoch-700 checkpoint, see 04_fit_stage3.R).
##
## WHAT   Stage 2, inverse calibration: every training site-period simulated
##        year by year from its observed stand, with the OBSERVED recruitment
##        fed in (env columns reg_obs_<s> through an identity env layer and a
##        custom regeneration returning exp(pred)), fitting
##          growth       species parameters + env weights
##          mortality    species parameters + env weights + U-shape (sen)
##          competition  species parameters
##        The saved model.pt carries the stage-1 regeneration for simulation;
##        stage 3 refits the regeneration from the epoch-700 checkpoint.
## MODEL  - mortality = FINN::mortality's linear predictor with growth in
##          percent (growth / 0.01) plus a U-shape term sen * (dbh/100)^2,
##          sen = 5 sigmoid(raw) in [0, 5] per species, starting at +0.05 on the
##          logit at 1 m dbh (custom parameters mort_sen_raw trainable,
##          mort_sen_max and mort_growth_scale fixed);
##        - L2 penalty lambda * sum(mask * (W - W_prefit)^2) on the env weights
##          of growth and mortality WITHOUT the intercept column, lambda 0.1;
##        - management = linear tier 0-3, constant over the decade;
##        - conditional per-tree mortality scoring (loss_aggregation tree, on
##          the observed stand), growth aggregated per site and species;
##        - Adam lr 0.01, batch 512, checkpoints every 100 epochs.
## HOW    (FINN interface only)
##        - L2 through fit(optimizer =): an Adam subclass whose step() adds
##          2 lambda mask (W - W1) to the gradients of nn_growth.0.weight and
##          nn_mortality.0.weight. FINN clips each parameter group's gradient
##          (clip_norm 2) BEFORE step(), so the penalty gradient is unclipped.
##        - Loss weights: growth and mortality scaled by their intercept-only
##          baselines (FINN's "auto"), read by a 1-epoch lr = 0 probe fit;
##          dbh, ba, trees, regeneration weighted 0; weight_scale = "raw".
##        - period_mode "per_site": FINN scores loss_aggregation "tree" there only.
## DATA   DATA_DIR/prepared/two-interval-v8/train; prefit DATA_DIR/weights/zoo-bwi-v1
## RUN    Rscript calibration/03_fit_stage2.R   (from the FINNmodelBWI repo root; OMP_NUM_THREADS=1 for a reproducible single-thread run)
## OUTPUT DATA_DIR/fits/stage2: model_init.pt, model_stage2.pt (training model,
##        observed recruitment), model.pt (stage-2 growth/mortality/competition
##        + stage-1 regeneration, for simulation), checkpoints/epoch_<n>model.pt,
##        loss_history.csv, l2_deviation.csv, species_parameters.csv,
##        mortality_size_curve.csv, fit_summary.md
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0 = main 0848b60
suppressPackageStartupMessages({
  library(FINN)
  library(data.table)
  library(torch)
})

DATA_DIR <- Filter(dir.exists, c(file.path("..", "..", "FINN-zoo-data", "bwi"), file.path("..", "FINN-zoo-data", "bwi")))[1]
FIT_DIR  <- file.path(DATA_DIR, "fits", "stage2")

PARAMS <- list(
  prep_version  = "two-interval-v8",       # the study's data version
  train_set     = "train",                 # 32,191 training site-periods
  holdout_set   = "holdout_blocked",       # not read here; recorded for the validation
  stage1_dir    = file.path(DATA_DIR, "weights", "zoo-bwi-v1"),   # the released prefit (model_stage1.pt, stage1_coefficients.csv, dispersion_stage1.csv); DATA_DIR/fits/stage1 from 02_fit_stage1.R is bit-identical
  env_vars      = c("mat_c", "map_mm", "elev_m", "management", "drought2"),   # = the prefit's: single drought index
  losses        = "mortgrowth",            # growth + mortality only ("all" would add dbh, ba, trees at the census)
  l2_lambda     = 0.1,                     # pull of the env slopes towards the prefit
  l2_intercept  = FALSE,                   # intercept column not penalised
  sen_max       = 5,                       # U-shape cap: at most +5 on the logit at 1 m
  sen_init      = 0.05,                    # U-shape start: +0.05 at 1 m, practically stage 1's linear size term
  growth_scale  = 0.01,                    # growth enters mortality as growth / 0.01 (growth in percent)
  period        = 10L,
  patch_size_ha = 0.1,                     # the expansion's patch area (PREPARATION.md)
  n_patches     = 25L,                     # passed to fit(); with an initial stand FINN uses its single plot
  n_quantiles   = 10L,                     # FINN's default height-quantile competition, as stage 1
  light_mask    = TRUE,                    # 0.2.0's masked competition
  recruits_dbh  = 7,                       # BWI caliper threshold
  dbh_binsize   = 2,                       # angle-count virtual trees -> cohorts
  reg_floor     = 1e-3,                    # FINN::regeneration additive floor
  reg_obs_eps   = 1e-8,                    # recruits/ha/yr standing in for an observed zero (log scale)
  reg_circle_r_m = 1.75, reg_incr_cm_yr = 0.30, reg_band_cm = 1.9,   # the regeneration circle
  epochs        = 1000L,                   # stage 3 starts from epoch 700 (the drought response is flat from 700 on)
  checkpoints   = 100L,                    # a model file and a parameter snapshot every 100 epochs
  batchsize     = 512L,
  lr            = 0.01,
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse",
                    mortality = "binomial", regeneration = "nbinom"),
  weight_scale  = "raw",                   # FINN 0.2.0 (#37): the numeric weights act on the raw losses, as in the study
  bound_tol     = 0.001,                   # within 0.1 % of its range from a bound = at the bound
  device        = "cpu",
  seed          = 123L
)
dir.create(file.path(FIT_DIR, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
stopifnot(identical(formals(FINN::mortality)$growth_scale, PARAMS$growth_scale))   # growth in percent (0.2.0 default)
stopifnot("light_mask" %in% names(formals(createProcess)))   # a FINN with masked competition (>= 0.2.0)
stopifnot("weight_scale" %in% names(formals(fit)))             # FINN >= 0.2.0 (main 0848b60)

## ---- 1. data -----------------------------------------------------------------
train_dir <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs  <- fread(file.path(train_dir, "obs_dt.csv"))
env  <- fread(file.path(train_dir, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init <- fread(file.path(train_dir, "init_cohort_dt.csv"))
tree_obs <- fread(file.path(train_dir, "tree_obs_dt.csv"))
species_dt <- fread(file.path(train_dir, "species_dt.csv"))
n_species <- max(obs$species)
init_cohort <- makeInitCohorts(init, dbh_binsize = PARAMS$dbh_binsize, Nspecies = n_species)
P <- PARAMS$period
stopifnot(all(obs$period_length == P))     # every observation is a decade

## the observed recruitment as env columns: log recruits per ha per year of each
## site-period and species, constant over the decade's years
reg_cols <- paste0("reg_obs_", seq_len(n_species))
REG <- dcast(obs[year == P, .(siteID, col = paste0("reg_obs_", species),
                              v = log(pmax(reg, 0) / P + PARAMS$reg_obs_eps))],
             siteID ~ col, value.var = "v")
setcolorder(REG, c("siteID", reg_cols))
env_obsreg <- merge(env, REG, by = "siteID")
setorder(env_obsreg, siteID, year)
message(sprintf("train: %d site-periods, %d species, %d tallied trees (%d at risk)",
                uniqueN(obs$siteID), n_species, nrow(tree_obs), sum(!is.na(tree_obs$died))))

## ---- 2. stage 1 (the prefit and the L2 anchor) ---------------------------------
## env weights from the coefficient table (a loaded model does not expose its
## frozen env layers), species parameters and dispersion from the model file
cf1 <- fread(file.path(PARAMS$stage1_dir, "stage1_coefficients.csv"))
W1 <- lapply(c(growth = "growth", mortality = "mortality", regeneration = "regeneration"), function(p)
  unname(as.matrix(cf1[process == p][order(species), c("intercept", PARAMS$env_vars), with = FALSE])))
s1 <- torch_load(file.path(PARAMS$stage1_dir, "model_stage1.pt"))
par1 <- list(competition = as.matrix(s1$par_competition$detach()), growth = as.matrix(s1$par_growth$detach()),
             mortality = as.matrix(s1$par_mortality$detach()), regeneration = as.matrix(s1$par_regeneration$detach()))
theta_s1 <- as.numeric(s1$par_theta_recruits$detach())
th_csv <- fread(file.path(PARAMS$stage1_dir, "dispersion_stage1.csv"))[order(species), theta_stage1]
stopifnot(max(abs(theta_s1 - th_csv) / th_csv) < 1e-4)
layer <- function(nn) as_array(nn$parameters[[1]]$detach()$cpu())
w_obs <- (pi * PARAMS$reg_circle_r_m^2 / 1e4) / (PARAMS$reg_incr_cm_yr / PARAMS$reg_band_cm * P)
env_formula <- as.formula(paste("~", paste(PARAMS$env_vars, collapse = " + ")))

## ---- 3. the training model -------------------------------------------------------
## custom regeneration: the observed annual recruits per ha, nothing else
regeneration_observed <- function(species, parReg, pred, light) {
  torch_exp(pred)[, NULL]$`repeat`(c(1, species$shape[2], 1))
}
obsreg_formula <- as.formula(paste("~ 0 +", paste(reg_cols, collapse = " + ")))
## U-shaped mortality with growth in percent: FINN::mortality's linear predictor
## (growth / mort_growth_scale) plus sen * (dbh/100)^2, sen = sen_max * sigmoid(raw)
## per species in [0, sen_max]
mortality_ushape <- function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL) {
  if (is.null(growth)) growth <- self$g
  d   <- dbh / 100
  sen <- (self$mort_sen_max * torch::torch_sigmoid(self$mort_sen_raw))[species]
  torch::torch_sigmoid(pred + growth / self$mort_growth_scale * parMort[, 3][species] + light * parMort[, 1][species] +
                       parMort[, 2][species] * d + sen * d * d)
}
sen_raw_init <- rep(qlogis(PARAMS$sen_init / PARAMS$sen_max), n_species)
build_training_model <- function() finn(
  N_species    = n_species,
  recruits_dbh = PARAMS$recruits_dbh,
  reg_floor    = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, initSpecies = par1$competition,
                                       optimizeSpecies = TRUE, n_quantiles = PARAMS$n_quantiles,
                                       light_mask = PARAMS$light_mask),
  growth_process       = createProcess(env_formula, FINN::growth, initSpecies = par1$growth,
                                       initEnv = list(W1$growth),
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(env_formula, mortality_ushape, initSpecies = par1$mortality,
                                       initEnv = list(W1$mortality),
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(mort_sen_raw = list(init = sen_raw_init, trainable = TRUE),
                                                                mort_sen_max = list(init = PARAMS$sen_max, trainable = FALSE),
                                                                mort_growth_scale = list(init = PARAMS$growth_scale, trainable = FALSE))),
  regeneration_process = createProcess(obsreg_formula, regeneration_observed,
                                       initEnv = list(diag(n_species)),
                                       optimizeSpecies = FALSE, optimizeEnv = FALSE,
                                       sample_regeneration = FALSE)
)
model <- build_training_model()
## the training model starts at stage 1 (its env layers are built inside fit() from initEnv)
stopifnot(max(abs(as.matrix(model$par_growth$detach()) - par1$growth)) < 1e-4,
          max(abs(as.matrix(model$par_mortality$detach()) - par1$mortality)) < 1e-4,
          max(abs(as.matrix(model$par_competition$detach()) - par1$competition)) < 1e-4)
torch_save(model, file.path(FIT_DIR, "model_init.pt"))
model$period_mode <- "per_site"            # FINN scores loss_aggregation = "tree" on its per-site path only

## ---- 4. loss weights ---------------------------------------------------------------
## Every loss is scaled by its intercept-only baseline (FINN's weights = "auto").
## For "mortgrowth" FINN needs a numeric vector, so a probe fit of one epoch at
## lr = 0 on a throwaway copy reads the baselines.
if (PARAMS$losses == "all") {
  weights <- "auto"
} else if (PARAMS$losses == "mortgrowth") {
  probe <- build_training_model(); probe$period_mode <- "per_site"
  fit(probe, env = env_obsreg, data = obs, tree_data = tree_obs, loss_aggregation = c(mortality = "tree"),
      init_cohort = init_cohort, epochs = 1L, batchsize = PARAMS$batchsize, lr = 0, weights = "auto",
      optimizer = torch::optim_adam,
      loss_family = PARAMS$loss_family, patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha,
      device = PARAMS$device, env_autoscale = FALSE, checkpoints = 1000L, plot_progress = FALSE)
  base <- probe$loss_baseline
  weights <- c(0, 0, 0, 1 / base[["growth"]], 1 / base[["mortality"]], 0)
  rm(probe); gc()
} else stop("unknown losses")

## ---- 5. stage 2 with the L2 pull towards stage 1 ---------------------------------
## Adam whose step() adds the gradient of lambda * sum(mask * (W - W1)^2) for the
## env weights of growth and mortality (matched by parameter name); mask is 0 in
## the intercept column (column 1 of the env weights, ~ 1 + env_vars)
l2_mask <- matrix(1, n_species, 1L + length(PARAMS$env_vars))
if (!PARAMS$l2_intercept) l2_mask[, 1] <- 0
optim_adam_l2 <- torch::optimizer("optim_adam_l2", inherit = torch::optim_adam,
  initialize = function(params, lr, anchors, lambda, mask) {
    stopifnot(all(names(anchors) %in% names(params)), all(dim(params[[names(anchors)[1]]]) == dim(anchors[[1]])),
              all(dim(mask) == dim(anchors[[1]])))
    self$l2 <- lapply(names(anchors), function(nm)
      list(p = params[[nm]],
           a = torch_tensor(anchors[[nm]], dtype = params[[nm]]$dtype, device = params[[nm]]$device),
           m = torch_tensor(mask, dtype = params[[nm]]$dtype, device = params[[nm]]$device)))
    self$lambda <- lambda
    super$initialize(params, lr = lr)
  },
  step = function(closure = NULL) {
    with_no_grad(for (x in self$l2) {
      g <- 2 * self$lambda * x$m * (x$p - x$a)
      if (is.null(x$p$grad) || !x$p$grad$numel()) x$p$grad <- g else x$p$grad$add_(g)
    })
    super$step(closure)
  })
anchors <- list(nn_growth.0.weight = W1$growth, nn_mortality.0.weight = W1$mortality)
optimizer_l2 <- function(params, lr, ...) optim_adam_l2(params, lr = lr, anchors = anchors, lambda = PARAMS$l2_lambda, mask = l2_mask)

FINN.seed(PARAMS$seed)
t0 <- Sys.time()
fit(model,
    env              = env_obsreg,
    data             = obs,
    tree_data        = tree_obs,
    loss_aggregation = c(mortality = "tree"),
    init_cohort      = init_cohort,
    epochs           = PARAMS$epochs,
    batchsize        = PARAMS$batchsize,
    lr               = PARAMS$lr,
    weights          = weights,
    weight_scale     = PARAMS$weight_scale,
    loss_family      = PARAMS$loss_family,
    patches          = PARAMS$n_patches,
    patch_size       = PARAMS$patch_size_ha,
    optimizer        = optimizer_l2,
    device           = PARAMS$device,
    env_autoscale    = FALSE,              # stage 1 and reg_obs are on the prepared scale
    folder           = file.path(FIT_DIR, "checkpoints"),
    checkpoints      = PARAMS$checkpoints,
    plot_progress    = FALSE)
min_fit <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
torch_save(model, file.path(FIT_DIR, "model_stage2.pt"))
hist <- as.data.table(do.call(rbind, model$history))[, 1:6]
setnames(hist, c("dbh", "ba", "trees", "growth", "mort", "reg"))
hist[, epoch := .I]
fwrite(hist, file.path(FIT_DIR, "loss_history.csv"))

## the penalty and the deviation from stage 1 per coefficient
W2 <- list(growth = layer(model$nn_growth), mortality = layer(model$nn_mortality))
penalty <- PARAMS$l2_lambda * sum(sapply(names(W2), function(p) sum(l2_mask * (W2[[p]] - W1[[p]])^2)))
DEV <- rbindlist(lapply(names(W2), function(p) {
  d <- data.table(process = p, species = seq_len(n_species),
                  `colnames<-`(as.data.table(W2[[p]] - W1[[p]]), c("intercept", PARAMS$env_vars)))
  merge(d, species_dt, by = "species")
}))
fwrite(DEV[order(process, species)], file.path(FIT_DIR, "l2_deviation.csv"))

## the fed recruitment must equal the observed period totals over fit()'s
## window: observation (year Y, period p) <- env rows Y-p+1..Y = sim years Y-p+2..Y+1
sim_chk <- simulateForest(model, env = env_obsreg, init_cohort = init_cohort,
                          patch_size = PARAMS$patch_size_ha, device = "cpu")
S <- as.data.table(sim_chk$wide$site)
for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
fed <- obs[, .(sim_year = (year - period_length + 2L):(year + 1L)), by = .(siteID, obs_year = year, species, reg)]
fed <- merge(fed, S[, .(siteID, sim_year = year, species, r = r_mean_ha)], by = c("siteID", "sim_year", "species"))[
  , .(reg = reg[1], reg_m = sum(r)), by = .(siteID, obs_year, species)]
fed_ratio <- sum(fed$reg_m) / sum(fed$reg)
stopifnot(abs(fed_ratio - 1) < 1e-3)

## ---- 6. the model for simulation: stage-2 processes + stage-1 regeneration -------
stage2 <- model
sen_raw2 <- as.numeric(stage2$mort_sen_raw$detach()$cpu())
model <- finn(
  N_species    = n_species,
  recruits_dbh = PARAMS$recruits_dbh,
  reg_floor    = PARAMS$reg_floor,
  recruit_obs_weight = w_obs,
  competition_process  = createProcess(~0, FINN::competition, initSpecies = as.matrix(stage2$par_competition$detach()),
                                       optimizeSpecies = FALSE, n_quantiles = PARAMS$n_quantiles,
                                       light_mask = PARAMS$light_mask),
  growth_process       = createProcess(env_formula, FINN::growth, initSpecies = as.matrix(stage2$par_growth$detach()),
                                       initEnv = list(W2$growth), optimizeSpecies = FALSE, optimizeEnv = FALSE),
  mortality_process    = createProcess(env_formula, mortality_ushape, initSpecies = as.matrix(stage2$par_mortality$detach()),
                                       initEnv = list(W2$mortality), optimizeSpecies = FALSE, optimizeEnv = FALSE,
                                       custom_parameters = list(mort_sen_raw = list(init = sen_raw2, trainable = FALSE),
                                                                mort_sen_max = list(init = PARAMS$sen_max, trainable = FALSE),
                                                                mort_growth_scale = list(init = PARAMS$growth_scale, trainable = FALSE))),
  regeneration_process = createProcess(env_formula, FINN::regeneration, initSpecies = par1$regeneration,
                                       initEnv = list(W1$regeneration), dispersion_parameter = theta_s1,
                                       optimizeSpecies = FALSE, optimizeEnv = FALSE)
)
stopifnot(max(abs(as.numeric(model$par_theta_recruits$detach()) - theta_s1)) < 1e-4 * max(theta_s1),
          max(abs(as.matrix(model$par_mortality$detach()) - as.matrix(stage2$par_mortality$detach()))) < 1e-4)
torch_save(model, file.path(FIT_DIR, "model.pt"))

## the fitted species parameters and which sit at a bound (FINN's default
## bounds: competition height [0.3, 0.7], strength [0, 2]; growth light
## [0.01, 0.99], size [0.01, 4]; mortality light, size, growth [-10, 10])
spar <- function(block, nm) `colnames<-`(as.data.table(as.matrix(stage2[[block]]$detach())), nm)
bnd  <- function(block) list(lo = as.numeric(stage2[[paste0(block, "_lower")]]$cpu()), hi = as.numeric(stage2[[paste0(block, "_upper")]]$cpu()))
SPP <- cbind(species_dt[order(species)], spar("par_competition", c("comp_height", "comp_strength")),
             spar("par_growth", c("growth_light", "growth_size")), spar("par_mortality", c("mort_light", "mort_size", "mort_growth")))
BL <- list(par_competition = c("comp_height", "comp_strength"), par_growth = c("growth_light", "growth_size"),
           par_mortality = c("mort_light", "mort_size", "mort_growth"))
AT <- rbindlist(lapply(names(BL), function(b) { bb <- bnd(b); rbindlist(lapply(seq_along(BL[[b]]), function(j) {
  v <- SPP[[BL[[b]][j]]]; tol <- PARAMS$bound_tol * (bb$hi[j] - bb$lo[j])
  data.table(parameter = BL[[b]][j], lower = bb$lo[j], upper = bb$hi[j], n_at_lower = sum(v <= bb$lo[j] + tol),
             n_at_upper = sum(v >= bb$hi[j] - tol),
             species_at_bound = paste(SPP$species_name[v <= bb$lo[j] + tol | v >= bb$hi[j] - tol], collapse = ", "))
})) }))
fwrite(SPP, file.path(FIT_DIR, "species_parameters.csv"))
## the fitted size curve per species: sen, the linear size slope, and the dbh of
## minimum mortality d* = -p_dbh / (2 sen) (Inf when sen ~ 0 or p_dbh >= 0)
SEN <- data.table(species = seq_len(n_species), p_dbh = SPP$mort_size, sen = PARAMS$sen_max * plogis(sen_raw2))
SEN[, dbh_min_mortality_cm := fifelse(sen > 1e-3 & p_dbh < 0, -p_dbh / (2 * sen) * 100, Inf)]
SEN <- merge(SEN, species_dt, by = "species")
fwrite(SEN, file.path(FIT_DIR, "mortality_size_curve.csv"))
stopifnot(max(abs(sen_raw2 - sen_raw_init)) > 0 || PARAMS$epochs < 1L)   # the U-shape parameter was fitted

## ---- 7. what was fitted --------------------------------------------------------
sum_l <- function(rows) { rows <- unique(pmin(pmax(rows, 1L), nrow(hist)))
  round(colMeans(as.matrix(hist[rows, .(dbh, ba, trees, growth, mort)])), 4) }
fmt <- function(v) paste(sprintf("%s %.4g", names(v), v), collapse = ", ")
dr <- function(p, v) paste(sprintf("%s %.2f -> %.2f", species_dt$species_name[order(species_dt$species)],
                                   W1[[p]][, 1L + match(v, PARAMS$env_vars)], W2[[p]][, 1L + match(v, PARAMS$env_vars)]), collapse = ", ")
writeLines(c(
  "# stage 2 (inverse calibration), FINNmodelBWI", "",
  sprintf("Fitted %s by `calibration/03_fit_stage2.R` (FINN %s, torch %s): %.1f min, %d epochs, on %s.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), packageVersion("FINN"), packageVersion("torch"), min_fit, PARAMS$epochs, PARAMS$device), "",
  sprintf("Data: `prepared/%s/%s` (%d site-periods, %d species, %d tallied trees). Holdout for validation: `%s`.",
          PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, nrow(tree_obs), PARAMS$holdout_set), "",
  sprintf("Prefit and L2 anchor: `%s`. Env `%s`, linear. Mortality = U-shape (sen in [0, %g], start %g at 1 m) with growth in percent (growth / %g). Competition light_mask = %s. Losses: %s (weights %s, weight_scale %s).",
          PARAMS$stage1_dir, paste(deparse(env_formula), collapse = ""), PARAMS$sen_max, PARAMS$sen_init, PARAMS$growth_scale, PARAMS$light_mask, PARAMS$losses,
          if (is.character(weights)) "auto" else paste(signif(weights, 3), collapse = ", "), PARAMS$weight_scale), "",
  sprintf("L2 penalty lambda %g on the env weights of growth and mortality (intercept %s): final penalty %.4g; max |W - W1| over the slopes growth %.3g, mortality %.3g; intercept growth %.3g, mortality %.3g (`l2_deviation.csv`).",
          PARAMS$l2_lambda, if (PARAMS$l2_intercept) "penalised" else "free", penalty,
          max(abs(W2$growth - W1$growth)[, -1]), max(abs(W2$mortality - W1$mortality)[, -1]),
          max(abs(W2$growth - W1$growth)[, 1]), max(abs(W2$mortality - W1$mortality)[, 1])), "",
  sprintf("Fitted: growth (species + env), mortality (species + env + U-shape), competition. Fixed: observed recruitment fed (ratio %.4f); the saved `model.pt` carries the stage-1 regeneration.", fed_ratio), "",
  sprintf("U-shape: sen %s; dbh of minimum mortality (cm) %s (`mortality_size_curve.csv`).",
          paste(sprintf("%.3f", range(SEN$sen)), collapse = "-"),
          paste(sprintf("%s %.0f", SEN$species_name, SEN$dbh_min_mortality_cm), collapse = ", ")), "",
  sprintf("Mortality drought2 slope (per pooled sd), stage 1 -> stage 2: %s.", dr("mortality", "drought2")), "",
  sprintf("Growth drought2 slope, stage 1 -> stage 2: %s.", dr("growth", "drought2")), "",
  sprintf("Mortality management slope, stage 1 -> stage 2: %s.", dr("mortality", "management")), "",
  sprintf("Mortality growth term (per 1 %% annual growth): %s (`species_parameters.csv`).",
          paste(sprintf("%s %.2f", SPP$species_name, SPP$mort_growth), collapse = ", ")), "",
  sprintf("Species parameters at a bound (within %g %% of the range):", 100 * PARAMS$bound_tol), "",
  "| parameter | lower | upper | at lower | at upper | species |", "|---|---|---|---|---|---|",
  sprintf("| %s | %g | %g | %d | %d | %s |", AT$parameter, AT$lower, AT$upper, AT$n_at_lower, AT$n_at_upper, AT$species_at_bound), "",
  sprintf("Baseline-scaled losses (mean of 10 epochs): start %s; end %s.", fmt(sum_l(1:10)), fmt(sum_l((nrow(hist) - 9):nrow(hist)))), "",
  sprintf("Fit: batch %d, lr %g, env_autoscale FALSE, checkpoints every %d epochs.", PARAMS$batchsize, PARAMS$lr, PARAMS$checkpoints)
), file.path(FIT_DIR, "fit_summary.md"))
message(sprintf("done: %.1f min", min_fit))
