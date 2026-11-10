#' Simulate the same sites under the full management range (tiers 0-3)
#'
#' The management tier of this model is derived from the same decade's
#' removals in the calibration data, so it is a post-treatment covariate: a
#' forced tier is a sensitivity to the harvest confounder, not a harvest
#' scenario, and results should be reported as the range over tiers 0-3
#' rather than for one tier. This helper runs [FINN::simulateForest()] once
#' per tier with the same seed and returns the stand totals per site, year
#' and tier (and, optionally, the per-species table), so a range can be drawn
#' directly.
#'
#' @param model A fitted model from [load_model()].
#' @param env Environment table in natural units (see [scale_env()]); its
#'   `management` column is overwritten by each tier in turn.
#' @param init_cohort Initial cohorts from [FINN::makeInitCohorts()], or
#'   `NULL` for bare ground.
#' @param tiers Integer vector of tiers to run, default `0:3`.
#' @param seed Seed passed to [FINN::FINN.seed()] before every run, so the
#'   tiers differ only through the covariate (the random draws still differ
#'   between tiers because the stands diverge).
#' @param patch_size,patches Passed to [FINN::simulateForest()]; the model's
#'   fitted patch size is 0.1 ha.
#' @param species Return the per-species table instead of the stand totals.
#' @param ... Further arguments for [FINN::simulateForest()].
#' @return A `data.table` with `tier`, `siteID`, `year` and, for the totals,
#'   `ba_ha` (m2 per ha), `stems_ha`, `removal` (stem-weighted annual removal
#'   probability), `recruits_ha` (recruits per ha in that year); per species
#'   additionally `species`, `name`, `dbh`, `growth`, `mort`.
#' @export
simulate_tiers = function(model, env, init_cohort = NULL, tiers = 0:3, seed = 1L, patch_size = 0.1, patches = 20L,
                          species = FALSE, ...) {
  stopifnot(all(tiers %in% 0:3))
  sp = species_table()
  out = lapply(tiers, function(t) {
    e = as.data.frame(env); e$management = t
    FINN::FINN.seed(seed)
    s = FINN::simulateForest(model, env = scale_env(e), init_cohort = init_cohort, patch_size = patch_size,
                             patches = patches, device = "cpu", ...)$wide$site
    s = data.table::as.data.table(lapply(s, function(x) if (is.factor(x)) as.integer(as.character(x)) else x))
    s[, tier := t]
    s
  })
  S = data.table::rbindlist(out)
  if (species) {
    S = merge(S, data.table::as.data.table(sp), by = "species")
    return(S[order(tier, siteID, year, species), c("tier", "siteID", "year", "species", "name", "dbh", "ba", "trees", "growth", "mort", "r_mean_ha"), with = FALSE])
  }
  S[, .(ba_ha = sum(ba) / patch_size, stems_ha = sum(trees) / patch_size,
        removal = if (sum(trees) > 0) sum(mort * trees) / sum(trees) else NA_real_,
        recruits_ha = sum(r_mean_ha, na.rm = TRUE)), by = .(tier, siteID, year)][order(tier, siteID, year)]
}
