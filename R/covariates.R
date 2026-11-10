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

#' Derive the model's drought index from monthly soil moisture
#'
#' The calibration's `drought2` is built from the DWD Climate Data Center's
#' monthly 1 km grids of modelled soil moisture (percent of usable field
#' capacity under grass over sandy loam,
#' `https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/soil_moist/`,
#' CC BY 4.0): per site, the growing-season mean (April to September) of each
#' year is z-standardised with the site's own mean and standard deviation over
#' the baseline years 2000-2022, and `drought2` of year `t` is the mean of the
#' anomalies of `t` and `t - 1`. This function does exactly that for monthly
#' values a user has extracted at new sites (or for any other soil-moisture or
#' water-balance series, in which case the index is defined analogously but
#' is not the calibration's).
#'
#' @param monthly A data.frame with `siteID`, `year`, `month` (1-12) and
#'   `soil_moist` (any unit; only the anomaly is used).
#' @param baseline Years over which each site's mean and sd are taken; the
#'   calibration used 2000:2022. The years must be present in `monthly`.
#' @param months Growing-season months, default `4:9`.
#' @return A data.frame with `siteID`, `year`, `z` (that year's anomaly) and
#'   `drought2`; the first year of each site has no previous year and gets
#'   `drought2 = z` (the calibration's first year, 2002 or 2012, had a
#'   previous year available).
#' @export
drought_index = function(monthly, baseline = 2000:2022, months = 4:9) {
  m = as.data.frame(monthly)
  need = c("siteID", "year", "month", "soil_moist")
  if (!all(need %in% names(m))) stop("monthly needs the columns ", paste(need, collapse = ", "), call. = FALSE)
  m = m[m$month %in% months & is.finite(m$soil_moist), ]
  gs = stats::aggregate(soil_moist ~ siteID + year, m, mean)
  gs = gs[order(gs$siteID, gs$year), ]
  base = gs[gs$year %in% baseline, ]
  if (!nrow(base)) stop("no baseline years in monthly", call. = FALSE)
  mu = stats::aggregate(soil_moist ~ siteID, base, mean); names(mu)[2] = "mu"
  sd = stats::aggregate(soil_moist ~ siteID, base, stats::sd); names(sd)[2] = "sd"
  gs = merge(merge(gs, mu, by = "siteID"), sd, by = "siteID")
  gs$z = (gs$soil_moist - gs$mu) / gs$sd
  gs = gs[order(gs$siteID, gs$year), ]
  prev = ave(gs$z, gs$siteID, FUN = function(z) c(NA, z[-length(z)]))
  gs$drought2 = ifelse(is.na(prev), gs$z, (gs$z + prev) / 2)
  gs[, c("siteID", "year", "z", "drought2")]
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
#'   z-anomaly; 0 = long-term mean, negative = drier). One row per annual
#'   transition: row `year = k` drives simulated year `k` (years 1..n); the
#'   initial state is the init cohort, there is no year-0 row.
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
