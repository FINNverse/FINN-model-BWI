# Licence of the fitted parameters and the BWI-derived data in this package

The code of the package and of `calibration/` is licensed GPL (>= 3) (file
`LICENSE`). This file covers everything derived from the German National
Forest Inventory (Bundeswaldinventur, BWI) data:

- the fitted model parameters: the release asset `parameters.pt` (named list
  of parameter tensors, sha256 in `weights.json`), the stage-3 model object
  `model.pt` of the calibration, and their human-readable tables
  (`calibration/parameters_summary.csv`, `calibration/validation/`);
- the example stands in `inst/extdata` (`example_env.csv`,
  `example_init.csv`: four inventory corners of the 2012 survey, environment
  in natural units, initial cohorts; positions as 1 km INSPIRE grid cells
  rounded to 0.05 degrees, no corner or tract identifiers, no personal names);
- the species table and the summary counts in `calibration/expected_preparation.csv`.

**BWI base data: Thuenen-Institut, licensed CC BY-SA 4.0**
(Creative Commons Attribution-ShareAlike 4.0 International,
https://creativecommons.org/licenses/by-sa/4.0/). Under its share-alike term
the items above are released under **CC BY-SA 4.0** as well, with this
attribution:

- Thuenen-Institut, Dritte Bundeswaldinventur - Basisdaten (Stand 20.03.2015),
  https://bwi.info/Download/de/BWI-Basisdaten/ACCESS2003/ (BWI 2002 and 2012,
  the 2002-2012 change tables, metadata, INSPIRE grid cells).
- Thuenen-Institut, Bundeswaldinventur 2022 und Kohlenstoffinventur 2017 -
  Inventurarchiv, https://ci.thuenen.de/rest/v1/ (BWI 2022 and the 2017
  carbon inventory).

The parameters also depend on climate normals from WorldClim 2.1 (Fick, S.E. &
Hijmans, R.J. (2017) WorldClim 2: new 1-km spatial resolution climate surfaces
for global land areas, Int. J. Climatol. 37, 4302-4315, doi:10.1002/joc.5086),
free for academic and other non-commercial use; the rasters are not
redistributed, only values extracted at the plot locations are included. They
also depend on Deutscher Wetterdienst Climate Data Center soil-moisture grids
(CC BY 4.0). Model: FINNverse/FINN-model-BWI,
Kaeber & Pichler (2026).
