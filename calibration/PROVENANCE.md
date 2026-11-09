# Provenance of the raw layers behind `two-interval-v8` and the v0.1.0 parameters

Every raw file the calibration reads, with source, licence, date pulled and
SHA-256 (full list, one line per file: `provenance_sha256.txt`; 1,124 files).
The preparation code is `github.com/ykaeber/bwi-dataprep` at commit
`3dc3c059f190acda92421b06fbe595cd2d1f02e6` (branch `port-two-interval-v8`,
GPL-3), run in the order of `run_all.R`; its `scripts/91` + `scripts/71`
reproduce the study's 43 FINN tables byte for byte (`two-interval-v8.sha256`).

## 1. BWI 2002 and 2012 base data: ACCESS2003 release (Thuenen Institute)

- Source: https://bwi.info/Download/de/BWI-Basisdaten/ACCESS2003/ (data status
  20.03.2015, INSPIRE shapefile 02.03.2016). Pulled 2026-07-28.
- Licence: **BWI base data: Thuenen-Institut, licensed CC BY-SA 4.0** (stated
  by the Thuenen BWI support in the API-access correspondence of 30 July
  2026, on which access was granted; the download page itself prints no
  licence name). Anonymised: tract positions only as 1 km INSPIRE grid cells
  "to protect the anonymity of the sample"
  (`_Hinweis_BWI-DatenUnterACCESS2003.pdf`). Citation as Thuenen requests
  (`Zitiervorschlag-BWI-Basisdaten-ACCESS2003-(2015).txt`):
  **"Thuenen-Institut, Dritte Bundeswaldinventur - Basisdaten (Stand
  20.03.2015)", https://bwi.info/Download/de/BWI-Basisdaten/ACCESS2003/**.
  The raw tables are not in this repository (`run_all.R` downloads them);
  everything derived from them (fitted parameters, example stands, summary
  tables) is released share-alike under CC BY-SA 4.0 with that attribution
  (`inst/model/LICENSE-parameters.md`). (Thuenen's field-database server
  software, github.com/Thuenen-Forest-Ecosystems/TFM-Server, is AGPL-3.0; that
  licence covers the software, not the data.)
- Files (SHA-256):

| file | bytes | sha256 |
|---|---|---|
| bwi20150320_alle_daten2002.zip | 111,970,463 | `0f2996171c3ed85ca81379d65ccae3c61da0961aeb33af8ba311de1bbfb07ac7` |
| bwi20150320_alle_daten2012.zip | 79,999,650 | `75569cff6e5aeaa9d3d1446b064e77c2fb3e37df21597367559a872329ae4551` |
| bwi20150320_alle_daten2002-2012.zip | 79,722,836 | `2a48b5b6c5a4ffaffdca0d4edf63dcf41a99c625978ed5e194aadea26d687b25` |
| bwi20150320_alle_metadaten.zip | 4,342,553 | `28ebcc22b2784c00d36c3c526c33db878901a43b2709038b91add72cc5b73e7f` |
| ShapeFile_Tnr_INSPIRE_Poly1000.zip | 1,583,409 | `263ded9d29628378a43f6680f2472323afabc1fc975c3771f00c11465bfcb47f` |
| Zitiervorschlag-BWI-Basisdaten-ACCESS2003-(2015).txt | 396 | `3496920bec8469fc6a88094d0da16e8317ade6c1031bd06524ed10137873e6d9` |
| _Hinweis_BWI-DatenUnterACCESS2003.pdf | 56,860 | `c47237dd5bae8954cbd1a37b83b503d578a3dac41c8ac3cae8be02d51f0d48cf` |
| Beschreibung_z_Tabellen.xlsx | 15,061 | `b0384b5c33728e9941282145f17629c85c8a59c705a2192f15a9a99ea72814da` |

The Access tables are exported to csv with mdbtools (`scripts/00`); the
field-crew name column (`TruppFuehrer`) and every other person-related column
are never used or redistributed.

## 2. BWI 2022 and the 2017 carbon inventory: Thuenen inventory-archive API

- Source: https://ci.thuenen.de/rest/v1/ (PostgREST, schema
  `inventory_archive`; documentation
  https://thuenen-forest-ecosystems.github.io/TFM-Documentation/), public
  anonymous key documented by Thuenen. Pulled 2026-07-30 with
  `bwi-dataprep/scripts/fetch_bwi4.sh` and `fetch_bwi4_tree.sh` (all
  intervals bwi2012, ci2017, bwi2022, ci2027; the panel uses bwi2012, ci2017
  and bwi2022).
- Licence and citation: **Thuenen-Institut, CC BY-SA 4.0** (as above); the API
  is meant for filtered queries rather than bulk extraction, Thuenen offers a
  self-hosted SQL dump for bulk use. Not redistributed here. Cite:
  **"Thuenen-Institut, Bundeswaldinventur 2022 und Kohlenstoffinventur 2017 -
  Inventurarchiv", https://ci.thuenen.de/rest/v1/** (documentation
  https://thuenen-forest-ecosystems.github.io/TFM-Documentation/). A live API:
  the data "may change subsequently due to corrections for quality assurance",
  so a rerun compares the sha256 below and reports differences (`run_all.R`
  step 1b).
- Files (SHA-256), `data/raw/bwi4/`:

| file | rows | sha256 |
|---|---|---|
| plot.csv | 783,532 | `edb81896e2c6ad567bbe85014b1a4602246132992edfc2b0f323736746cf1190` |
| tree.csv | 1,722,072 | `00466cc51bf357f325da26cc34d2da2ec7b05a254b23479995442403587646fe` |
| regeneration.csv | | `fcd98425e0c222c366dbacb0006688006ba98a58d591bf93052c45d838cc3a97` |
| deadwood.csv | | `01a7b41224e57ffc75ce6a91cab18ff50dc7eaf6c571a5ace4fc05319c0f697c` |
| cluster.csv | | `1378e49fe6c19443cd56948eca58b293118daca62867a8014d7654f7c3bc9f15` |
| edges.csv | | `01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b` |
| lookup_tree_species.csv | | `a6b7432013da286902065f1e454434f9aeca1fce5552c9be217c283645c1d378` |
| lookup_tree_status.csv | | `c3d916f63b5411c9505d60184f03b1ff03530ee4d791d951ac0cd6525a584660` |
| lookup_federal_state.csv | | `4ceb7940c9a633dbc8816b3f50649ef77af5a1d2c3a4aa956db0166d845390d1` |
| lookup_interval.csv (empty download) | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |

## 3. Drought index: DWD Climate Data Center, monthly 1 km grids

- Source: https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/
  (`soil_moist`, `evapo_p`, `evapo_r`, `drought_index`), 2000-01 to 2022-12,
  276 files each, 1,104 files (`grids_germany_monthly_<type>_<yyyymm>.asc.gz`).
  Pulled 2026-07-29/30 with `bwi-dataprep/scripts/fetch_dwd_drought.sh`.
- Licence: CC BY 4.0 ("Terms of use for data on the CDC-OpenData area: The
  Creative Commons BY 4.0 Licence applies", DWD, status May 2024;
  https://opendata.dwd.de/climate_environment/CDC/Terms_of_use.pdf). Attribution:
  Deutscher Wetterdienst.
- Only `soil_moist` enters the model (growing-season mean April-September,
  per-site z-anomaly 2000-2022; `drought2` = mean of the current and the
  previous year); the other three grid types are extracted for comparison.
- SHA-256 of all 1,104 files: `provenance_sha256.txt`.

## 4. Climate normals and elevation fallback: WorldClim 2.1

- Source: WorldClim 2.1, https://www.worldclim.org, 30 arc-second bioclimatic
  variables BIO1 (annual mean temperature) and BIO12 (annual precipitation),
  1970-2000; `wc2.1_30s_elev.tif` as the elevation fallback where the
  inventory has no tract elevation. Rasters dated 2019-02-28; the copy used
  was downloaded 2025-09-09.
- Citation: Fick, S.E. & Hijmans, R.J. (2017) WorldClim 2: new 1-km spatial
  resolution climate surfaces for global land areas. International Journal of
  Climatology 37, 4302-4315, doi:10.1002/joc.5086.
- Terms: WorldClim states "freely available for academic use and other
  non-commercial use. Redistribution or commercial use is not allowed without
  prior permission." The rasters are therefore not shipped; `run_all.R`
  downloads them. The values extracted at the site locations (two numbers per
  site) are part of the prepared data.
- SHA-256: `wc2.1_30s_bio_1.tif`
  `382c991cea4a20b42e345fcd5a897788feee05bc0d613dec032ff0cad1f4c375`,
  `wc2.1_30s_bio_12.tif`
  `6df0ef2f57f900733a569fcc4c8debd27b490978c5d82c2970485f0e66abbf28`.

## 5. Coordinates

Tract positions are the public 1 km INSPIRE (ETRS89-LAEA) grid cells of the
ACCESS2003 release (`ShapeFile_Tnr_INSPIRE_Poly1000.zip`); the cell centroid is
used for the WorldClim and DWD extractions (positional error <= ~700 m). No
exact coordinates exist in any file of this pipeline. Elevation comes from the
inventory's own per-tract value (`b0_tab.HoeheNN`), not from a DEM at the
anonymised position.

## 6. Prepared data and weights

- `two-interval-v8`: `PREPARATION.md` (copy of the version folder's) and the
  46-file manifest `two-interval-v8.sha256`; the committed summary
  `expected_preparation.csv` (written by `07_expected_preparation.R`).
- Release asset `parameters.pt` (v0.1.0): sha256
  `dadd5d300cee0976e6dc7874b75070c50a12f2da3fb1c48e95e38f896a156da7`, extracted
  from the stage-3 model object (sha256
  `5c046997df7f882b0cd7d12e468f1eb9a2c89f61542d557a5f089611a8f2b857`) by
  `06_extract_parameters.R`; `parameters_extraction.txt` records the check.
  Licence of the parameters and of every BWI-derived file: CC BY-SA 4.0
  (`inst/model/LICENSE-parameters.md`).
