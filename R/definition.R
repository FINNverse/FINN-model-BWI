#' Build the model structure and set its fitted parameters
#'
#' The model definition of this package: the same `createProcess()` and
#' `finn()` calls as `calibration/04_fit_stage3.R` (which produced the released
#' parameters), with every parameter frozen at its fitted value. The model is
#' rebuilt with the installed FINN rather than restored from a saved model
#' object, so it always runs with the installed package's methods.
#'
#' Variant `"linear"`: FINN's process functions with species-specific LINEAR
#' environmental responses of growth, mortality and regeneration on the link
#' scale, over five covariates in this order: `mat_c`, `map_mm`, `elev_m`,
#' `management`, `drought2` (see `covariates.csv` and [scale_env()] for their
#' units and scaling). Competition is FINN's masked height-quantile light
#' competition (`light_mask = TRUE`, 10 quantiles). Mortality is FINN's logit
#' with growth in percent (growth / 0.01) plus a U-shaped size term
#' `sen * (dbh / 100)^2`, `sen = 5 * sigmoid(mort_sen_raw)` per species.
#' Regeneration is `FINN::regeneration_saturation`: recruits enter at 7 cm dbh,
#' negative-binomially distributed, with a Beverton-Holt cap `K = exp(reg_logK)`
#' per species (recruits per ha per year).
#'
#' @param variant Name of the model variant (see [model_info()]); only
#'   `"linear"` exists.
#' @param parameters Named list of parameter tensors, as saved by
#'   `calibration/06_extract_parameters.R`: `competition` (species x 2: height
#'   allometry, competition strength), `growth` (species x 2: light, size),
#'   `mortality` (species x 3: light, size, growth), `regeneration_light`
#'   (species x 1), `theta` (recruit dispersion per species, as fitted on the
#'   regeneration circles), `env_growth`, `env_mortality`, `env_regeneration`
#'   (species x 6: intercept + the five covariates), `mort_sen_raw` (species),
#'   `mort_sen_max` (1), `mort_growth_scale` (1), `reg_logK` (species).
#' @param theta_sim_factor Multiplier applied to the fitted recruit dispersion
#'   `theta` for simulation. The fitted `theta` describes counts on the 9.6 m2
#'   regeneration circles; FINN draws one recruit count per 0.1 ha patch and
#'   year, which averages part of that overdispersion out. The calibration's
#'   validation uses 3.5 (estimated from the between-circle covariance of the
#'   observed counts); pass 1 to keep the fitted value.
#' @return A fitted `FINN` model (all parameters frozen) for
#'   [FINN::simulateForest()] or [simulate_model()].
#' @keywords internal
build_model = function(variant = "linear", parameters, theta_sim_factor = 3.5) {
  if (!identical(variant, "linear"))
    stop("Unknown variant '", variant, "'; this package has the variant 'linear'.", call. = FALSE)
  need = c("competition", "growth", "mortality", "regeneration_light", "theta", "env_growth", "env_mortality",
           "env_regeneration", "mort_sen_raw", "mort_sen_max", "mort_growth_scale", "reg_logK")
  miss = setdiff(need, names(parameters))
  if (length(miss)) stop("parameters lack: ", paste(miss, collapse = ", "), call. = FALSE)
  m = function(x) as.matrix(torch::as_array(x$detach()$cpu()))
  v = function(x) as.numeric(torch::as_array(x$detach()$cpu()))
  n_species = nrow(m(parameters$competition))
  ## the env formula is written out in every createProcess() call below (FINN
  ## evaluates the formula argument itself, so it cannot be held in a variable)
  stopifnot(ncol(m(parameters$env_growth)) == 6L, nrow(m(parameters$env_growth)) == n_species)

  ## U-shaped mortality with growth in percent: FINN::mortality's linear predictor
  ## (growth / mort_growth_scale) plus sen * (dbh / 100)^2, sen in [0, mort_sen_max]
  mortality_ushape = function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL) {
    if (is.null(growth)) growth = self$g
    d   = dbh / 100
    sen = (self$mort_sen_max * torch::torch_sigmoid(self$mort_sen_raw))[species]
    torch::torch_sigmoid(pred + growth / self$mort_growth_scale * parMort[, 3][species] + light * parMort[, 1][species] +
                         parMort[, 2][species] * d + sen * d * d)
  }

  model = FINN::finn(
    N_species    = n_species,
    recruits_dbh = 7,                      # BWI caliper threshold (cm): recruits enter the stand at 7 cm
    reg_floor    = 1e-3,                   # additive floor of the recruitment mean (recruits per ha per year)
    competition_process  = FINN::createProcess(~0, FINN::competition, initSpecies = m(parameters$competition),
                                               optimizeSpecies = FALSE, n_quantiles = 10L, light_mask = TRUE),
    growth_process       = FINN::createProcess(~ mat_c + map_mm + elev_m + management + drought2, FINN::growth,
                                               initSpecies = m(parameters$growth),
                                               initEnv = list(m(parameters$env_growth)),
                                               optimizeSpecies = FALSE, optimizeEnv = FALSE),
    mortality_process    = FINN::createProcess(~ mat_c + map_mm + elev_m + management + drought2, mortality_ushape,
                                               initSpecies = m(parameters$mortality),
                                               initEnv = list(m(parameters$env_mortality)),
                                               optimizeSpecies = FALSE, optimizeEnv = FALSE,
                                               custom_parameters = list(
                                                 mort_sen_raw      = list(init = v(parameters$mort_sen_raw), trainable = FALSE),
                                                 mort_sen_max      = list(init = v(parameters$mort_sen_max), trainable = FALSE),
                                                 mort_growth_scale = list(init = v(parameters$mort_growth_scale), trainable = FALSE))),
    regeneration_process = FINN::createProcess(~ mat_c + map_mm + elev_m + management + drought2, FINN::regeneration_saturation,
                                               initSpecies = m(parameters$regeneration_light),
                                               initEnv = list(m(parameters$env_regeneration)),
                                               dispersion_parameter = theta_sim_factor * v(parameters$theta),
                                               custom_parameters = list(reg_logK = list(init = v(parameters$reg_logK), trainable = FALSE)),
                                               optimizeSpecies = FALSE, optimizeEnv = FALSE)
  )
  model$eval()
  model
}
