First release of FINNmodelBWI: a fitted FINN forest model for German forests, calibrated to the National Forest Inventory (Bundeswaldinventur: 2002, 2012 and 2022 surveys) and packaged as an R simulation tool.

Variant `linear`: 15 species, annual step on 0.1 ha patches, FINN's process functions with species-specific linear responses to mean annual temperature, annual precipitation, elevation, a management tier 0-3 and a two-year soil-moisture drought index; U-shaped size mortality, masked light competition, recruitment with a fitted Beverton-Holt cap. Staged fit (prefit, inverse calibration against per-tree removals and survivor growth, regeneration refit) under FINN 0.2.0 release candidate 6705a50; the parameters are extracted, rebuilt and validated under FINN 0.2.0 main 0848b60 (torch 0.15.1), where `fit(weight_scale = "raw")` reproduces the fit objective.

Install and use:

```r
remotes::install_github("FINNverse/FINN@0848b60")
remotes::install_github("FINNverse/FINN-model-BWI@v0.1.0")
library(FINNmodelBWI)
torch::torch_set_num_threads(1L)
m <- load_model()   # downloads parameters.pt once, checks its sha256, rebuilds the model with the installed FINN
```

Simulations are pinned to one torch thread and `batchsize = 1000`; a different batchsize changes the random draws but not the model. `simulate_tiers()` runs the same sites under tiers 0-3, `drought_index()` derives the drought covariate from monthly soil moisture, `scale_env()` puts natural units on the model's scale.

Assets:

| file | what | sha256 |
|---|---|---|
| parameters.pt | named list of parameter tensors; `load_model()` rebuilds from it | dadd5d300cee0976e6dc7874b75070c50a12f2da3fb1c48e95e38f896a156da7 |
| model.pt | the whole stage-3 model object (`torch_save`) for `torch::torch_load()` | 5c046997df7f882b0cd7d12e468f1eb9a2c89f61542d557a5f089611a8f2b857 |

Holdout skill (8,596 corner-decades in 11 spatially blocked 100-km blocks the fit never saw; the persistence of the start inventory in brackets): species-level basal area R2 0.86 (0.83), stems R2 0.68 (0.44), stand basal area R2 0.44 (0.28), survivor growth R2 0.47, mortality (total removal) AUC 0.76, predicted/observed removals 0.93, stems predicted/observed 1.23.

Limitations: the management tier is a post-treatment covariate, so report results as the range over tiers 0-3; mortality is total removal (natural death plus harvest); FINN's simulated mortality has a floor of 0.01 per annual step (FINNverse/FINN#51), which the fit did not see and which binds for low-mortality species such as beech at tier 0; stem numbers run high; no bark-beetle process. Details in NEWS.md and `calibration/validation/report.md`.

Licences: code GPL (>= 3). BWI base data: Thuenen-Institut, CC BY-SA 4.0, so the fitted parameters (both assets) and the four BWI-derived example stands in the package are released share-alike under CC BY-SA 4.0 with the attribution in `inst/model/LICENSE-parameters.md`. Climate normals: WorldClim 2.1 (Fick & Hijmans 2017, doi:10.1002/joc.5086), free for academic and other non-commercial use; the rasters are not redistributed, only values extracted at the plot locations are included. DWD soil-moisture grids: CC BY 4.0. The raw data are not redistributed; `calibration/` rebuilds the prepared data from the public sources.
