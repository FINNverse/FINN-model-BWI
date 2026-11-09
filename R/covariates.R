#' Covariates of the model and how to derive them for new sites
#'
#' Reads `inst/model/covariates.csv`: one row per environmental covariate the
#' model expects in `env`, with its unit, source and derivation.
#'
#' @return A data.frame with columns `name`, `unit`, `source`, `derivation`.
#' @export
covariate_table = function() {
  utils::read.csv(system.file("model", "covariates.csv", package = utils::packageName()), stringsAsFactors = FALSE)
}

#' Put an environment table in natural units on the model's scale
#'
#' The model was fitted on z-scaled covariates (one mean and sd per covariate,
#' pooled over all calibration site-periods; `inst/model/env_scales.csv`).
#' `management` was left on its natural 0-3 tier scale and `drought2` is already
#' a standardised anomaly, but both are also listed in the scaling table (mean
#' 0 / sd 1 and the pooled mean / sd of the index), so every column goes
#' through the same formula `(x - mean) / sd`.
#'
#' @param env A data.frame with `siteID`, `year` and the covariates in natural
#'   units: `mat_c` (mean annual temperature, deg C), `map_mm` (annual
#'   precipitation, mm), `elev_m` (elevation, m), `management` (regular-harvest
#'   tier 0, 1, 2 or 3 per site, see [covariate_table()]) and `drought2`
#'   (mean of the current and the previous year's growing-season soil-moisture
#'   z-anomaly; 0 = long-term mean, negative = drier). Row `year = 0` is the
#'   initial state; row `k` drives simulated year `k + 1`.
#' @return A `data.table` (what [FINN::simulateForest()] expects) with the
#'   covariate columns scaled, ordered as the model expects.
#' @export
scale_env = function(env) {
  sc = utils::read.csv(system.file("model", "env_scales.csv", package = utils::packageName()), stringsAsFactors = FALSE)
  vars = c("mat_c", "map_mm", "elev_m", "management", "drought2")
  miss = setdiff(c("siteID", "year", vars), names(env))
  if (length(miss)) stop("env lacks the columns: ", paste(miss, collapse = ", "), call. = FALSE)
  if (!all(env$management %in% 0:3)) stop("management must be a tier 0, 1, 2 or 3", call. = FALSE)
  out = as.data.frame(env)[, c("siteID", "year", vars)]
  for (v in vars) out[[v]] = (out[[v]] - sc$mean[sc$variable == v]) / sc$sd[sc$variable == v]
  out = out[order(out$siteID, out$year), ]
  rownames(out) = NULL
  data.table::as.data.table(out)
}
