# FINN-model-BWI

A fitted [FINN](https://github.com/FINNverse/FINN) forest model calibrated to
the German National Forest Inventory (Bundeswaldinventur, BWI: the 2002, 2012
and 2022 surveys, two consecutive decades on the same sample corners), packaged
as a versioned simulation tool (R package `FINNmodelBWI`). The model card with
data, model, fit and holdout validation is in the
[FINN-zoo](https://finnverse.github.io/FINN-zoo/models.html).

**Status: v0.1.0 in preparation, no public release yet.** Until the release
asset exists, pass the parameter file by hand:
`load_model(path = "<...>/parameters.pt")`.

```r
remotes::install_github("FINNverse/FINN-model-BWI")
library(FINNmodelBWI)
m = load_model()          # downloads the released parameters once, checks sha256
species_table()           # 15 species, index -> Latin name
covariate_table()         # the five covariates, units, sources, how to derive them
vignette("quickstart", package = "FINNmodelBWI")
```

What the model needs per site and year (`env`, natural units; `scale_env()`
puts it on the model's scale): mean annual temperature, annual precipitation,
elevation, a **management tier** 0-3 (regular harvest intensity; 0 = none) and
a **drought index** (`drought2`, the mean of the current and the previous
year's growing-season soil-moisture anomaly; 0 = average, -1 = a dry year such
as 2018). Simulations start from a user-given stand (cohorts of species, dbh,
stems on a 0.1 ha patch) or from bare ground.

This repository is a tool, not the code of a manuscript: it is versioned for
users (see `NEWS.md`) and may differ from the models in our papers. The
calibration pipeline that produced each release is in `calibration/`
(`run_all.R` runs it end to end from the public sources; `PROVENANCE.md` lists
every raw layer with licence and checksum).

License: GPL (>= 3) for code and fitted parameters; the calibration data keep
the licence of their source (see `calibration/PROVENANCE.md` and the model card).
