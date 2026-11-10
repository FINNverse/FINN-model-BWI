# FINNmodelBWI 0.1.0

First fitted release (tag `v0.1.0`), variant `linear`.

* Data: German National Forest Inventory (BWI) 2002 -> 2012 -> 2022, prepared
  set `two-interval-v8` (40,787 corner-decades; 32,191 fitted, 8,596 in a
  spatially blocked holdout of 11 of 53 100-km blocks), 15 species. Holdout
  metrics computed on the corrected env layout `two-interval-v12` (see below).
* Model: FINN process functions with species-specific linear environmental
  responses (mean annual temperature, annual precipitation, elevation,
  management tier 0-3, two-year soil-moisture drought index), masked
  height-quantile competition, U-shaped size mortality with growth in percent,
  recruitment with a fitted Beverton-Holt cap; 0.1 ha patches, annual step.
* Fit: three stages (prefit per process; inverse calibration of growth,
  mortality and competition against per-tree removals and survivor growth with
  an L2 pull to the prefit, Adam, 1000 epochs, epoch 700 used; regeneration
  refit on the frozen demography, 100 epochs). Fitted under FINN 0.2.0 release
  candidate 6705a50 (fast-setup + F1); the parameter list is extracted,
  rebuilt and validated under FINN 0.2.0 main 0848b60 (torch 0.15.1), where
  `fit(weight_scale = "raw")` reproduces the fit objective.
* Holdout (spatially blocked, both decades, `calibration/validation/report.md`):
  species-level basal area R2 0.86 (persistence of the start inventory 0.83),
  stem density R2 0.68 (persistence 0.44), species-mean dbh R2 0.38, survivor
  growth R2 0.47, recruitment R2 0.02; stand totals: basal area R2 0.44
  (persistence 0.28), stems 0.48 (0.11); mortality (total removal, tree-weighted)
  AUC 0.76, predicted/observed removals 0.93; stems predicted/observed 1.23.
* Release assets (GitHub release `v0.1.0`): `parameters.pt`, a named list of
  parameter tensors rebuilt by `build_model()` (what `load_model()` uses), and
  `model.pt`, the whole stage-3 model object for `torch::torch_load()` under the
  same FINN version; sha256 of both in `inst/model/weights.json`. Pinned for
  simulation: one torch thread and `batchsize = 1000`; a different batchsize
  changes the random draws but not the model.
* Env layout: the v0.1.0 parameters come from the study's fit on
  `two-interval-v8`, whose env grid carried an extra year-0 row, so that fit
  simulated 11 annual transitions per 10-year interval. The fitted rates are
  unaffected (their windows were 10 years); the holdout state metrics above are
  computed on the corrected layout `two-interval-v12` (10 env rows, census at
  step 10), which all calibration scripts now read. One env row = one annual
  transition; there is no year-0 row.
* Helpers: `scale_env()`, `covariate_table()`, `simulate_tiers()` (the same
  sites under tiers 0-3), `drought_index()` (the calibration's `drought2` from
  monthly soil moisture).

Known limitations:

* The management tier is derived from the same decade's removals in the
  calibration, i.e. a post-treatment covariate: report results as the range
  over tiers 0-3, not for one tier. Mortality is total removal (natural death
  plus harvest); the tier raises the removal rate, it does not cut trees.
* FINN's simulated mortality has a floor of 0.01 per model step
  (FINNverse/FINN#51), here 1 % per year. The fit did not see the floor;
  it binds for species whose fitted removal rate lies below it (beech at tier
  0: 0.55 % per year fitted, 1 % simulated; 51 % of all stem-years at tier 0,
  29 % at the recorded tiers, where the mean rate moves from 3.05 to 3.16 %
  per year). Long tier-0 and bare-ground runs (plateau 40-45 m2/ha, beech
  dominated) are therefore denser and younger than the fitted rates imply.
* Stem numbers run high (predicted/observed 1.23 on the holdout): the
  recruitment level is about twice the ingrowth the angle-count sample implies
  for the first decade (a documented choice of the data preparation);
  recruitment itself is not predicted (R2 0.02).
* No density-dependent mortality beyond light competition, no bark-beetle
  process; the drought response is logit-linear and extrapolates beyond
  `drought2` of about -2.3.

Licences: code GPL (>= 3); BWI base data Thuenen-Institut CC BY-SA 4.0, so
the fitted parameters (both assets) and the four BWI-derived example stands in
`inst/extdata` are released share-alike under CC BY-SA 4.0 with the
attribution in `inst/model/LICENSE-parameters.md`. WorldClim 2.1 values at the
site locations (CC BY-SA 4.0, Fick & Hijmans 2017); DWD soil-moisture grids
CC BY 4.0. The quickstart runs synthetic stand types and the real example
stands.

# FINNmodelBWI (development)

* Repository created from the FINN model template.
