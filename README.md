# FINN-model-BWI

A fitted [FINN](https://github.com/FINNverse/FINN) forest model calibrated to the German National Forest Inventory (Bundeswaldinventur, BWI),
packaged as a versioned simulation tool (R package `FINNmodelBWI`). The model card is in the
[FINN-zoo](https://finnverse.github.io/FINN-zoo/models.html).

**Status: in preparation, no release yet.**

```r
remotes::install_github("FINNverse/FINN-model-BWI")
library(FINNmodelBWI)
m = load_model()          # downloads the released parameters once, checks sha256
species_table()
model_info()
```

This repository is a tool, not the code of a manuscript: it is versioned for
users (see NEWS.md) and may differ from the models in our papers. The
calibration pipeline that produced each release is in `calibration/`.

License: GPL (>= 3) for code and fitted parameters; the calibration data keep
the license of their source (see the model card).
