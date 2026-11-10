# FINN-model-BWI

A fitted [FINN](https://github.com/FINNverse/FINN) forest model calibrated to
the German National Forest Inventory (Bundeswaldinventur, BWI: the 2002, 2012
and 2022 surveys, two consecutive decades on the same sample corners), packaged
as a versioned simulation tool (R package `FINNmodelBWI`). The model card with
data, model, fit and holdout validation is in the
[FINN-zoo](https://finnverse.github.io/FINN-zoo/models.html).

**Current release: v0.1.0** (see `NEWS.md`). The fitted parameters are release
assets of this repository (`parameters.pt`, the named parameter list
`load_model()` rebuilds from, and `model.pt`, the whole fitted model object);
`load_model()` downloads `parameters.pt` once and checks its sha256. A local
copy can be passed by hand: `load_model(path = "<...>/parameters.pt")`.

```r
remotes::install_github("FINNverse/FINN@0848b60")           # FINN 0.2.0, the version the parameters were validated with
remotes::install_github("FINNverse/FINN-model-BWI@v0.1.0")
library(FINNmodelBWI)
torch::torch_set_num_threads(1L)   # pinned: the validation ran single-threaded
m = load_model()          # downloads the released parameters once, checks sha256
species_table()           # 15 species, index -> Latin name
covariate_table()         # the five covariates, units, sources, how to derive them
vignette("quickstart", package = "FINNmodelBWI")
```

Simulations are run with one torch thread and an explicit `batchsize` (1000
sites per `simulateForest()` batch, the value pinned in `model_info()`;
`simulate_model()` and `simulate_tiers()` use it by default). A different
batchsize changes the random draws but not the model.

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

## License

- **Code** (R package, calibration scripts): GPL (>= 3).
- **BWI base data: Thünen-Institut, licensed CC BY-SA 4.0.** Cite as Thünen
  requests: "Thünen-Institut, Dritte Bundeswaldinventur – Basisdaten (Stand
  20.03.2015)", https://bwi.info/Download/de/BWI-Basisdaten/ACCESS2003/ (BWI
  2002 and 2012), and "Thünen-Institut, Bundeswaldinventur 2022 und
  Kohlenstoffinventur 2017 – Inventurarchiv", https://ci.thuenen.de/rest/v1/
  (BWI 2022). The raw data are not in this repository; `calibration/`
  downloads them from the source.
- **Fitted parameters** (the release assets `parameters.pt` and `model.pt`)
  **and the BWI-derived example stands** in `inst/extdata`:
  share-alike, released under **CC BY-SA 4.0** with the attribution above
  (`inst/model/LICENSE-parameters.md`).
- Climate normals: WorldClim 2.1, Fick, S.E. & Hijmans, R.J. (2017) Int. J.
  Climatol. 37, 4302–4315, doi:10.1002/joc.5086 (values at the site
  locations only). Drought grids: Deutscher Wetterdienst, Climate Data Center,
  CC BY 4.0. Details in `calibration/PROVENANCE.md`.
