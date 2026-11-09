#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINNmodelBWI | calibration/04_fit_stage3.R
##
## ZOO COPY of the FINN-bwi study script
##   fits/2026-11-08-A-stage3-v8-drought2-ushape-lambda0.1-mortgrowth-masked-regcap-ep700/fit.R
## (branch inventory-simulator) with these changes only: FINN 0.2.0 (main
## 0848b60) with fit(weight_scale = "raw") (FINN#37; the weights vector is
## recruitment-only, so the scale changes only the effective learning rate of
## the one active loss, but the study's 0848b60 refits use "raw" and so does
## this copy); the L2 optimizer builds its tensors on the parameters' device;
## data from the zoo data folder, stage 2 from DATA_DIR/fits/stage2, the
## prefit from the released weights folder; output to DATA_DIR/fits/stage3.
## Runtime: 100 epochs, ~30 min single-threaded on a laptop. The result of
## this script IS the released model.pt (stage-2 demography + stage-3
## regeneration).
##
## WHAT   Stage 3 (regeneration), from the EPOCH-700 checkpoint of stage 2:
##          - growth, mortality and competition FROZEN at stage 2;
##          - regeneration (species light parameter + env response + a fitted
##            Beverton-Holt cap) fitted against the recruitment loss alone,
##            negative binomial on stems per regeneration circle
##            (recruit_obs_weight), dispersion frozen at the prefit's theta;
##          - L2 pull (lambda 1) of the regeneration env SLOPES towards the
##            prefit, intercept free; 100 epochs, checkpoints every 25.
## WHY    The inverse calibration changed competition and therefore the light
##        at 1 m, so the regeneration must be refitted on the frozen demography.
##        The cap (FINN::regeneration_saturation, recruits = K m / (K + m),
##        K = exp(reg_logK) per species, start 500 recruits per ha per year)
##        stops the one-year recruitment bursts that a steep climate response
##        plus the NB draw produce at a few sites, which otherwise dominate the
##        holdout stem-count error. The cap parameter is not L2-pulled (it has
##        no prefit counterpart).
## DATA   DATA_DIR/prepared/two-interval-v8/train
## RUN    Rscript calibration/04_fit_stage3.R   (from the FINNmodelBWI repo root, after 03_fit_stage2.R)
## OUTPUT DATA_DIR/fits/stage3: source_stage2.pt, model_init.pt, model.pt (stage-2
##        demography + stage-3 regeneration), checkpoints/epoch_<n>model.pt,
##        loss_history.csv, regeneration_vs_stage1.csv, fit_summary.md
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0 = main 0848b60
suppressPackageStartupMessages({
  library(FINN)
  library(data.table)
  library(torch)
})

DATA_DIR <- Filter(dir.exists, c(file.path("..", "..", "FINN-zoo-data", "bwi"), file.path("..", "FINN-zoo-data", "bwi")))[1]
FIT_DIR  <- file.path(DATA_DIR, "fits", "stage3")

PARAMS <- list(
  prep_version  = "two-interval-v8",       # the study's data version
  train_set     = "train",
  holdout_set   = "holdout_blocked",
  stage2_dir    = file.path(DATA_DIR, "fits", "stage2"),          # written by 03_fit_stage2.R
  stage2_epoch  = 700L,                    # the checkpoint the study evaluates (drought response flat from 700 on)
  source_stage2 = file.path(FIT_DIR, "source_stage2.pt"),         # copied from stage2_dir on the first run
  stage1_dir    = file.path(DATA_DIR, "weights", "zoo-bwi-v1"),   # regeneration start and L2 anchor (= stage 2's prefit)
  env_vars      = c("mat_c", "map_mm", "elev_m", "management", "drought2"),   # = stage 2's: single drought index
  l2_lambda     = 1,                       # pull of the regeneration env slopes towards stage 1
  l2_free_cols  = 1L,                      # column 1 = intercept: no pull, the recruitment level is free
  period        = 10L,
  patch_size_ha = 0.1,
  n_patches     = 25L,
  n_quantiles   = 10L,
  light_mask    = TRUE,                    # 0.2.0's masked competition, as stage 2
  recruits_dbh  = 7,
  dbh_binsize   = 2,
  reg_floor     = 1e-3,
  reg_circle_r_m = 1.75, reg_incr_cm_yr = 0.30, reg_band_cm = 1.9,
  epochs        = 100L,
  checkpoints   = 25L,
  batchsize     = 512L,
  lr            = 0.01,
  weights       = c(0, 0, 0, 0, 0, 1),     # recruitment loss only (dbh, ba, trees, growth, mortality, regeneration)
  weight_scale  = "raw",                   # FINN 0.2.0 (#37): the vector acts on the raw losses, as in the study
  reg_cap_K_init = 500,                    # Beverton-Holt cap start, recruits per ha per year (observed site totals p99 492, p99.9 821)
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse",
                    mortality = "binomial", regeneration = "nbinom"),
  frozen_tol    = 1e-4,
  sen_max       = 5,                       # as stage 2: U-shape cap
  growth_scale  = 0.01,                    # as stage 2: growth enters mortality as growth / 0.01
  device        = "cpu",
  seed          = 123L
)
dir.create(file.path(FIT_DIR, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
stopifnot(identical(formals(FINN::mortality)$growth_scale, PARAMS$growth_scale))
stopifnot("light_mask" %in% names(formals(createProcess)), "weight_scale" %in% names(formals(fit)))   # FINN >= 0.2.0

## ---- 0. the stage-2 source: the epoch-700 checkpoint, copied once -------------------
if (!file.exists(PARAMS$source_stage2)) {
  ck <- list.files(file.path(PARAMS$stage2_dir, "checkpoints"), pattern = "^epoch_[0-9]+model\\.pt$", recursive = TRUE)
  ck_dt <- data.table(file = ck, n = as.integer(sub(".*epoch_([0-9]+)model\\.pt$", "\\1", ck)), offset = 0L)
  ck_dt[grepl("^resume_[0-9]+/", file), offset := as.integer(sub("^resume_([0-9]+)/.*", "\\1", file))]
  ck_dt[, epoch := offset + n]
  src <- ck_dt[epoch == PARAMS$stage2_epoch]
  if (nrow(src) != 1L) stop(sprintf("no epoch-%d checkpoint in %s/checkpoints yet", PARAMS$stage2_epoch, PARAMS$stage2_dir))
  stopifnot(file.copy(file.path(PARAMS$stage2_dir, "checkpoints", src$file), PARAMS$source_stage2))
  message("source: ", file.path(PARAMS$stage2_dir, "checkpoints", src$file))
}
src_md5 <- unname(tools::md5sum(PARAMS$source_stage2))

## ---- 1. data -----------------------------------------------------------------
train_dir <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs  <- fread(file.path(train_dir, "obs_dt.csv"))
env  <- fread(file.path(train_dir, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init <- fread(file.path(train_dir, "init_cohort_dt.csv"))
species_dt <- fread(file.path(train_dir, "species_dt.csv"))
n_species <- max(obs$species)
init_cohort <- makeInitCohorts(init, dbh_binsize = PARAMS$dbh_binsize, Nspecies = n_species)
P <- PARAMS$period
w_obs_p <- function(p) (pi * PARAMS$reg_circle_r_m^2 / 1e4) / (PARAMS$reg_incr_cm_yr / PARAMS$reg_band_cm * p)
w_obs <- w_obs_p(P)                         # FINN's recruit_obs_weight is one number: the decade's factor
stopifnot(max(abs(obs[reg > 0, reg * w_obs_p(period_length) - round(reg * w_obs_p(period_length))])) < 1e-6)   # whole stems per circle
stopifnot(all(obs$period_length == P))      # every observation is a decade, so w(10) is exact

## ---- 2. stage 2 (frozen) and stage 1 (regeneration start, L2 anchor) -----------
stage2 <- torch_load(PARAMS$source_stage2)
layer <- function(nn) as_array(nn$parameters[[1]]$detach()$cpu())
cf1 <- fread(file.path(PARAMS$stage1_dir, "stage1_coefficients.csv"))
W1_reg <- unname(as.matrix(cf1[process == "regeneration"][order(species), c("intercept", PARAMS$env_vars), with = FALSE]))
s1 <- torch_load(file.path(PARAMS$stage1_dir, "model_stage1.pt"))
par1_reg <- as.matrix(s1$par_regeneration$detach())
theta_s1 <- as.numeric(s1$par_theta_recruits$detach())
env_formula <- as.formula(paste("~", paste(PARAMS$env_vars, collapse = " + ")))
frozen <- list(competition = as.matrix(stage2$par_competition$detach()), growth = as.matrix(stage2$par_growth$detach()),
               mortality = as.matrix(stage2$par_mortality$detach()), growth_env = layer(stage2$nn_growth),
               mortality_env = layer(stage2$nn_mortality), sen_raw = as.numeric(stage2$mort_sen_raw$detach()$cpu()))
## the stage-2 mortality: U-shape with growth in percent (as 03_fit_stage2.R)
mortality_ushape <- function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL) {
  if (is.null(growth)) growth <- self$g
  d   <- dbh / 100
  sen <- (self$mort_sen_max * torch::torch_sigmoid(self$mort_sen_raw))[species]
  torch::torch_sigmoid(pred + growth / self$mort_growth_scale * parMort[, 3][species] + light * parMort[, 1][species] +
                       parMort[, 2][species] * d + sen * d * d)
}

## ---- 3. the stage-3 model ----------------------------------------------------------
model <- finn(
  N_species    = n_species,
  recruits_dbh = PARAMS$recruits_dbh,
  reg_floor    = PARAMS$reg_floor,
  recruit_obs_weight = w_obs,              # the loss scores stems on the circle
  competition_process  = createProcess(~0, FINN::competition, initSpecies = frozen$competition,
                                       optimizeSpecies = FALSE, n_quantiles = PARAMS$n_quantiles,
                                       light_mask = PARAMS$light_mask),
  growth_process       = createProcess(env_formula, FINN::growth, initSpecies = frozen$growth,
                                       initEnv = list(frozen$growth_env), optimizeSpecies = FALSE, optimizeEnv = FALSE),
  mortality_process    = createProcess(env_formula, mortality_ushape, initSpecies = frozen$mortality,
                                       initEnv = list(frozen$mortality_env), optimizeSpecies = FALSE, optimizeEnv = FALSE,
                                       custom_parameters = list(mort_sen_raw = list(init = frozen$sen_raw, trainable = FALSE),
                                                                mort_sen_max = list(init = PARAMS$sen_max, trainable = FALSE),
                                                                mort_growth_scale = list(init = PARAMS$growth_scale, trainable = FALSE))),
  regeneration_process = createProcess(env_formula, FINN::regeneration_saturation,
                                       custom_parameters = list(reg_logK = list(init = rep(log(PARAMS$reg_cap_K_init), n_species), trainable = TRUE)),
                                       initSpecies = par1_reg,
                                       initEnv = list(W1_reg), dispersion_parameter = theta_s1,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE)
)
## FINN always makes theta a trainable parameter; take it out of the gradient
invisible(model$par_theta_recruits_raw$requires_grad_(FALSE))
stopifnot(max(abs(as.numeric(model$par_theta_recruits$detach()) - theta_s1)) < 1e-4 * max(theta_s1))
torch_save(model, file.path(FIT_DIR, "model_init.pt"))

## Adam with the L2 pull of the regeneration env slopes towards stage 1; the mask
## (0 in the intercept column) leaves the intercept free
optim_adam_l2 <- torch::optimizer("optim_adam_l2", inherit = torch::optim_adam,
  initialize = function(params, lr, anchors, masks, lambda) {
    stopifnot(all(names(anchors) %in% names(params)), all(dim(params[[names(anchors)[1]]]) == dim(anchors[[1]])),
              identical(names(masks), names(anchors)))
    self$l2 <- lapply(names(anchors), function(nm)
      list(p = params[[nm]],
           a = torch_tensor(anchors[[nm]], dtype = params[[nm]]$dtype, device = params[[nm]]$device),
           m = torch_tensor(masks[[nm]], dtype = params[[nm]]$dtype, device = params[[nm]]$device)))
    self$lambda <- lambda
    super$initialize(params, lr = lr)
  },
  step = function(closure = NULL) {
    with_no_grad(for (x in self$l2) {
      g <- 2 * self$lambda * (x$p - x$a) * x$m
      if (is.null(x$p$grad) || !x$p$grad$numel()) x$p$grad <- g else x$p$grad$add_(g)
    })
    super$step(closure)
  })
mask_reg <- matrix(1, nrow(W1_reg), ncol(W1_reg))
mask_reg[, PARAMS$l2_free_cols] <- 0
anchors <- list(nn_regeneration.0.weight = W1_reg)
masks   <- list(nn_regeneration.0.weight = mask_reg)
optimizer_l2 <- function(params, lr, ...) optim_adam_l2(params, lr = lr, anchors = anchors, masks = masks, lambda = PARAMS$l2_lambda)

## ---- 4. fit ---------------------------------------------------------------------------
FINN.seed(PARAMS$seed)
t0 <- Sys.time()
fit(model,
    env           = env,
    data          = obs,
    init_cohort   = init_cohort,
    epochs        = PARAMS$epochs,
    batchsize     = PARAMS$batchsize,
    lr            = PARAMS$lr,
    weights       = PARAMS$weights,
    weight_scale  = PARAMS$weight_scale,
    loss_family   = PARAMS$loss_family,
    patches       = PARAMS$n_patches,
    patch_size    = PARAMS$patch_size_ha,
    optimizer     = optimizer_l2,
    device        = PARAMS$device,
    env_autoscale = FALSE,
    folder        = file.path(FIT_DIR, "checkpoints"),
    checkpoints   = PARAMS$checkpoints,
    plot_progress = FALSE)
min_fit <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
torch_save(model, file.path(FIT_DIR, "model.pt"))
hist <- as.data.table(do.call(rbind, model$history))[, 1:6]
setnames(hist, c("dbh", "ba", "trees", "growth", "mort", "reg"))
hist[, epoch := .I]
fwrite(hist, file.path(FIT_DIR, "loss_history.csv"))

## frozen blocks must still be stage 2's
dev <- c(growth_env    = max(abs(layer(model$nn_growth) - frozen$growth_env)),
         mortality_env = max(abs(layer(model$nn_mortality) - frozen$mortality_env)),
         growth_par    = max(abs(as.matrix(model$par_growth$detach()) - frozen$growth)),
         mortality_par = max(abs(as.matrix(model$par_mortality$detach()) - frozen$mortality)),
         competition   = max(abs(as.matrix(model$par_competition$detach()) - frozen$competition)),
         sen           = max(abs(as.numeric(model$mort_sen_raw$detach()$cpu()) - frozen$sen_raw)),
         theta_rel     = max(abs(as.numeric(model$par_theta_recruits$detach()) - theta_s1) / theta_s1))
message(sprintf("max |change| of the frozen blocks: %s", paste(sprintf("%s %.2g", names(dev), dev), collapse = ", ")))
stopifnot(all(dev < PARAMS$frozen_tol))

## the regeneration against stage 1
W3 <- layer(model$nn_regeneration)
REGDEV <- merge(data.table(species = seq_len(n_species), light_stage1 = par1_reg[, 1],
                           light_stage3 = as.numeric(as.matrix(model$par_regeneration$detach())[, 1]),
                           `colnames<-`(as.data.table(W1_reg), paste0("s1_", c("intercept", PARAMS$env_vars))),
                           `colnames<-`(as.data.table(W3), paste0("s3_", c("intercept", PARAMS$env_vars)))),
                species_dt, by = "species")
REGDEV[, K_fit := exp(as.numeric(model$reg_logK$detach()$cpu()))[species]]   # the fitted cap, recruits per ha per year
fwrite(REGDEV, file.path(FIT_DIR, "regeneration_vs_stage1.csv"))

## predicted / observed recruitment on train over fit()'s windows: observation
## (year Y, period p) <- sim years Y-p+2..Y+1; start (stage-1 regeneration) and end
rec_ratio <- function(m) {
  FINN.seed(PARAMS$seed)
  S <- as.data.table(simulateForest(m, env = env, init_cohort = init_cohort, patch_size = PARAMS$patch_size_ha, device = "cpu")$wide$site)
  for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
  W <- unique(obs[, .(siteID, obs_year = year, period_length)])[
    , .(year = (obs_year - period_length + 2L):(obs_year + 1L)), by = .(siteID, obs_year)]
  sum(merge(W, S[, .(siteID, year, r_mean_ha)], by = c("siteID", "year"))$r_mean_ha) / sum(pmax(obs$reg, 0))
}
ratio_start <- rec_ratio(torch_load(file.path(FIT_DIR, "model_init.pt")))
ratio_end   <- rec_ratio(model)

## ---- 5. what was fitted --------------------------------------------------------
writeLines(c(
  "# stage 3 (regeneration with a fitted cap), FINNmodelBWI", "",
  sprintf("Fitted %s by `calibration/04_fit_stage3.R` (FINN %s, torch %s): %.1f min, %d epochs, on %s.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), packageVersion("FINN"), packageVersion("torch"), min_fit, PARAMS$epochs, PARAMS$device), "",
  sprintf("Stage-2 source (frozen): `%s` (md5 %s) = epoch-%d checkpoint of `%s`. Regeneration start and L2 anchor: `%s`. Competition light_mask = %s.",
          PARAMS$source_stage2, substr(src_md5, 1, 8), PARAMS$stage2_epoch, PARAMS$stage2_dir, PARAMS$stage1_dir, PARAMS$light_mask), "",
  sprintf("Recruitment loss only (NB on stems per circle, w %.4g, theta frozen at stage 1; weight_scale %s); L2 lambda %g on the regeneration env slopes, intercept free: slopes sum((W - W1)^2) %.4g, max |change| %.3g; intercept change %s.",
          w_obs, PARAMS$weight_scale, PARAMS$l2_lambda, sum(((W3 - W1_reg) * mask_reg)^2), max(abs((W3 - W1_reg) * mask_reg)),
          paste(sprintf("%.2f", range(W3[, 1] - W1_reg[, 1])), collapse = " to ")), "",
  sprintf("Max change of the frozen blocks: %s.", paste(sprintf("%s %.2g", names(dev), dev), collapse = ", ")), "",
  sprintf("Predicted / observed recruitment on train: %.2f with the stage-1 regeneration, %.2f after stage 3.", ratio_start, ratio_end), "",
  sprintf("Recruitment cap K (Beverton-Holt, recruits per ha per year), start %g, fitted: %s.", PARAMS$reg_cap_K_init,
          paste(sprintf("%s %.0f", REGDEV$species_name, REGDEV$K_fit), collapse = ", ")), "",
  sprintf("Regeneration light parameter, stage 1 -> stage 3: %s.",
          paste(sprintf("%s %.2f -> %.2f", REGDEV$species_name, REGDEV$light_stage1, REGDEV$light_stage3), collapse = ", ")), "",
  sprintf("Recruitment loss: start %.4g, end %.4g (mean of 10 epochs).", mean(hist$reg[1:min(10, nrow(hist))]), mean(tail(hist$reg, 10)))
), file.path(FIT_DIR, "fit_summary.md"))
message(sprintf("done: %.1f min", min_fit))
