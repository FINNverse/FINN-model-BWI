# Validation report, FINNmodelBWI (variant linear)

Written 2026-11-10 06:00 by `calibration/05_validation.R`; FINN 0.2.0, torch 0.15.1; parameters `weights/v0.1.0/parameters.pt` sha256 `dadd5d300cee0976e6dc7874b75070c50a12f2da3fb1c48e95e38f896a156da7` rebuilt with `R/definition.R::build_model()`; data `prepared/two-interval-v12/holdout_blocked` (8596 site-periods).

Layout note: the v0.1.0 parameters come from the study's fit on `two-interval-v8`, whose env grid carried an extra year-0 row, so that fit simulated 11 annual transitions per 10-year interval; the fitted rates are unaffected (the rate windows were 10 years), and the state metrics below are computed on the corrected layout `two-interval-v12` (10 transitions, census at step 10).

## Holdout metrics

| response | level | n | rmse | r2 | pearson_r | pred_over_obs | auc |
|---|---|---|---|---|---|---|---|
| ba | species | 129000 | 0.2730 | 0.8580 | 0.9320 | 1.070 |  |
| trees | species | 129000 | 11.90 | 0.6800 | 0.8400 | 1.230 |  |
| dbh | species | 14400 | 11.80 | 0.3800 | 0.7130 | 0.8660 |  |
| growth | species | 13600 | 0.007370 | 0.4740 | 0.6900 | 1.030 |  |
| reg | species | 129000 | 262.0 | 0.01640 | 0.1310 | 0.9500 |  |
| ba | total | 8600 | 1.050 | 0.4420 | 0.7070 | 1.070 |  |
| trees | total | 8600 | 47.70 | 0.4780 | 0.7500 | 1.230 |  |
| reg | total | 8600 | 1020 | 0.01850 | 0.1500 | 0.9500 |  |
| ba | species, persistence reference | 129000 | 0.3030 | 0.8250 | 0.9140 | 1.020 |  |
| trees | species, persistence reference | 129000 | 15.90 | 0.4350 | 0.8150 | 1.290 |  |
| ba | total, persistence reference | 8600 | 1.190 | 0.2820 | 0.6300 | 1.020 |  |
| trees | total, persistence reference | 8600 | 62.40 | 0.1070 | 0.6990 | 1.290 |  |
| mortality | tree-weighted | 65600 |  |  |  | 0.9320 | 0.7570 |

## Per species (holdout)

| species_name | n_trees_at_risk | ba_r2 | ba_r2_persistence | stems_r2 | growth_r | mort_auc | mort_pred_over_obs | stems_pred_over_obs | recruits_pred_over_obs |
|---|---|---|---|---|---|---|---|---|---|
| Pinus sylvestris | 24000 | 0.8870 | 0.8320 | 0.7210 | 0.7620 | 0.7700 | 0.9580 | 1.130 | 0.5130 |
| Picea abies | 12800 | 0.6500 | 0.5900 | 0.5560 | 0.7140 | 0.6320 | 0.8770 | 1.230 | 1.300 |
| Fagus sylvatica | 10200 | 0.8770 | 0.8430 | 0.5360 | 0.5510 | 0.7910 | 1.040 | 1.500 | 1.290 |
| Quercus petraea | 5320 | 0.9020 | 0.9030 | 0.7330 | 0.5900 | 0.7580 | 0.7660 | 1.100 | 0.4860 |
| Pseudotsuga menziesii | 3080 | 0.8800 | 0.8710 | 0.7100 | 0.6830 | 0.7780 | 0.8840 | 1.090 | 1.130 |
| Quercus robur | 2760 | 0.9000 | 0.9010 | 0.6290 | 0.6530 | 0.7550 | 0.9850 | 1.250 | 1.090 |
| Betula pendula | 2000 | 0.7780 | 0.7440 | 0.6010 | 0.7300 | 0.7220 | 0.9380 | 1.200 | 0.5940 |
| Larix decidua | 1420 | 0.8000 | 0.7950 | 0.6250 | 0.6930 | 0.7540 | 0.8850 | 1.000 | 0.01050 |
| Carpinus betulus | 1060 | 0.8490 | 0.8500 | 0.7020 | 0.4750 | 0.7420 | 1.060 | 1.360 | 0.8300 |
| Alnus glutinosa | 805.0 | 0.9230 | 0.9150 | 0.8740 | 0.6270 | 0.7520 | 0.9490 | 1.110 | 3.110 |
| Larix kaempferi | 739.0 | 0.8830 | 0.8200 | 0.6860 | 0.6550 | 0.7780 | 1.030 | 0.9860 | 0.01750 |
| Fraxinus excelsior | 557.0 | 0.8050 | 0.7960 | 0.5430 | 0.6380 | 0.7130 | 1.290 | 1.180 | 0.9560 |
| Acer pseudoplatanus | 365.0 | 0.7100 | 0.7920 | 0.2920 | 0.7220 | 0.7400 | 0.5430 | 2.040 | 1.210 |
| Prunus avium | 320.0 | 0.8310 | 0.8490 | 0.6840 | 0.6690 | 0.7020 | 0.7010 | 1.150 | 0.4550 |
| Abies alba | 141.0 | 0.8970 | 0.9290 | -0.1520 | 0.6530 | 0.5930 | 1.010 | 4.600 | 3.370 |

## Bare-ground equilibrium (mean of the last 100 of 1000 years)

| loc | total_ba_ha | drift_last_vs_previous_100 | n_species_above_1 | dominant | dominant_share | stems_ha |
|---|---|---|---|---|---|---|
| 7.5 C, 998 mm, 504 m | 40.10 | -0.05180 | 2.000 | Fagus sylvatica | 0.9180 | 552.0 |
| 8.4 C, 693 mm, 281 m | 43.60 | 0.03500 | 3.000 | Fagus sylvatica | 0.8630 | 495.0 |
| 8.7 C, 578 mm, 84 m | 42.70 | 0.08160 | 3.000 | Fagus sylvatica | 0.9000 | 487.0 |
| 8.9 C, 615 mm, 244 m | 42.00 | 0.004930 | 3.000 | Fagus sylvatica | 0.7810 | 593.0 |
| 9.3 C, 760 mm, 292 m | 44.60 | 0.09250 | 3.000 | Fagus sylvatica | 0.8140 | 597.0 |

## Figures

| file | what |
|---|---|
| figures/loss_curves.png | stage-2 and stage-3 training losses of the released model |
| figures/holdout_obs_pred.png | observed vs predicted per site-period x species, five responses |
| figures/holdout_totals.png | observed vs predicted stand totals (basal area, stems) |
| figures/mortality_roc.png | tree-weighted ROC of the decade removal probability |
| figures/bare_ground.png | succession from bare ground at five sites |

## Provenance

| step | script | output |
|---|---|---|
| raw BWI, DWD and WorldClim layers -> two-interval-v8 | calibration/run_all.R steps 1-10 (bwi-dataprep, pinned commit) | prepared/two-interval-v8 (PREPARATION.md, PROVENANCE.md, expected_preparation.csv) |
| staging + checksums | calibration/01_stage_data_and_weights.R | FINN-zoo-data/bwi/{prepared,weights} |
| stage 1 prefit | calibration/02_fit_stage1.R | model_stage1.pt, stage1_coefficients.csv, dispersion_stage1.csv |
| stage 2 inverse calibration | calibration/03_fit_stage2.R | epoch-700 checkpoint, species_parameters.csv, l2_deviation.csv |
| stage 3 regeneration + cap | calibration/04_fit_stage3.R | model.pt |
| parameter list (release asset) | calibration/06_extract_parameters.R | parameters.pt, parameters_summary.csv |
| this report | calibration/05_validation.R | validation/metrics.csv, metrics_species.csv, bare_ground_equilibrium.csv, figures/ |
