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
* Holdout (spatially blocked, both decades): basal area R2 0.85 (species level;
  persistence reference 0.82), stem density R2 0.64, survivor growth R2 0.46,
  mortality AUC 0.76 (tree-weighted), predicted/observed removals 0.93. See
  `calibration/validation/report.md`.
* Release asset `parameters.pt` (sha256 in `inst/model/weights.json`), a named
  list of parameter tensors rebuilt by `build_model()`.
* No plot-level BWI rows are shipped (species table and summary counts only);
  the quickstart uses a synthetic stand, `calibration/09_example_sites.R`
  builds real stands locally after the data download. Data citation as the
  Thuenen Institute requests (README, PROVENANCE.md).

# FINNmodelBWI (development)

* Repository created from the FINN model template.
