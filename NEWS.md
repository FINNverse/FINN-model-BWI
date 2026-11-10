# FINNmodelBWI 0.1.0 (in preparation, not released)

First fitted model, variant `linear`:

* Data: German National Forest Inventory (BWI) 2002 -> 2012 -> 2022, prepared
  set `two-interval-v8` (40,787 corner-decades; 32,191 fitted, 8,596 in a
  spatially blocked holdout of 11 of 53 100-km blocks), 15 species.
* Model: FINN process functions with species-specific linear environmental
  responses (mean annual temperature, annual precipitation, elevation,
  management tier 0-3, two-year soil-moisture drought index), masked
  height-quantile competition, U-shaped size mortality with growth in percent,
  recruitment with a fitted Beverton-Holt cap; 0.1 ha patches, annual step.
* Fit: three stages (prefit per process; inverse calibration of growth,
  mortality and competition against per-tree removals and survivor growth with
  an L2 pull to the prefit, Adam, 1000 epochs, epoch 700 used; regeneration
  refit on the frozen demography, 100 epochs). FINN 0.2.0 release candidate
  6705a50 (fast-setup + F1); reproduced under FINN 0.2.0 main 0848b60 with
  `fit(weight_scale = "raw")`.
* Holdout (spatially blocked, both decades): basal area R2 0.86 (species level;
  persistence reference 0.83), stem density R2 0.68, survivor growth R2 0.47,
  mortality AUC 0.76 (tree-weighted), predicted/observed removals 0.93. See
  `calibration/validation/report.md`.
* Env layout: the v0.1.0 parameters come from the study's fit on
  `two-interval-v8`, whose env grid carried an extra year-0 row, so that fit
  simulated 11 annual transitions per 10-year interval. The fitted rates are
  unaffected (their windows were 10 years); the holdout state metrics above are
  computed on the corrected layout `two-interval-v12` (10 env rows, census at
  step 10), which all calibration scripts now read. One env row = one annual
  transition; there is no year-0 row.
* Release asset `parameters.pt` (sha256 in `inst/model/weights.json`), a named
  list of parameter tensors rebuilt by `build_model()`.
* Licences: code GPL (>= 3); BWI base data Thuenen-Institut CC BY-SA 4.0, so
  the fitted parameters and the four BWI-derived example stands in
  `inst/extdata` are released share-alike under CC BY-SA 4.0 with the
  attribution in `inst/model/LICENSE-parameters.md`. The quickstart runs
  synthetic stand types and the real example stands.

# FINNmodelBWI (development)

* Repository created from the FINN model template.
