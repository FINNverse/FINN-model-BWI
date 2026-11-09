# Validation report, FINNmodelBWI (variant linear)

Written 2026-11-09 22:20 by `calibration/05_validation.R`; FINN 0.2.0, torch 0.15.1; parameters `weights/zoo-bwi-v1/parameters.pt` sha256 `dadd5d300cee0976e6dc7874b75070c50a12f2da3fb1c48e95e38f896a156da7` rebuilt with `R/definition.R::build_model()`; data `prepared/two-interval-v8/holdout_blocked` (8596 site-periods).

## Holdout metrics

| response | level | n | rmse | r2 | pearson_r | pred_over_obs | auc |
|---|---|---|---|---|---|---|---|
| ba | species | 129000 | 0.2760 | 0.8550 | 0.9320 | 1.080 |  |
| trees | species | 129000 | 12.00 | 0.6770 | 0.8370 | 1.240 |  |
| dbh | species | 14400 | 12.20 | 0.3450 | 0.6930 | 0.8660 |  |
| growth | species | 13600 | 0.007460 | 0.4610 | 0.6790 | 1.020 |  |
| reg | species | 129000 | 262.0 | 0.01600 | 0.1300 | 0.9410 |  |
| ba | total | 8600 | 1.060 | 0.4310 | 0.7050 | 1.080 |  |
| trees | total | 8600 | 48.10 | 0.4700 | 0.7440 | 1.240 |  |
| reg | total | 8600 | 1020 | 0.01790 | 0.1470 | 0.9410 |  |
| ba | species, persistence reference | 129000 | 0.3030 | 0.8250 | 0.9140 | 1.020 |  |
| trees | species, persistence reference | 129000 | 15.90 | 0.4350 | 0.8150 | 1.290 |  |
| ba | total, persistence reference | 8600 | 1.190 | 0.2820 | 0.6300 | 1.020 |  |
| trees | total, persistence reference | 8600 | 62.40 | 0.1070 | 0.6990 | 1.290 |  |
| mortality | tree-weighted | 65600 |  |  |  | 0.9310 | 0.7580 |

## Per species (holdout)

| species_name | n_trees_at_risk | ba_r2 | ba_r2_persistence | stems_r2 | growth_r | mort_auc | mort_pred_over_obs | stems_pred_over_obs | recruits_pred_over_obs |
|---|---|---|---|---|---|---|---|---|---|
| Pinus sylvestris | 24000 | 0.8850 | 0.8320 | 0.7140 | 0.7540 | 0.7700 | 0.9550 | 1.120 | 0.5060 |
| Picea abies | 12800 | 0.6420 | 0.5900 | 0.5630 | 0.7080 | 0.6320 | 0.8810 | 1.240 | 1.290 |
| Fagus sylvatica | 10200 | 0.8750 | 0.8430 | 0.5250 | 0.5170 | 0.7910 | 1.040 | 1.550 | 1.280 |
| Quercus petraea | 5320 | 0.9000 | 0.9030 | 0.7360 | 0.5860 | 0.7600 | 0.7650 | 1.100 | 0.4810 |
| Pseudotsuga menziesii | 3080 | 0.8740 | 0.8710 | 0.7190 | 0.6830 | 0.7780 | 0.8790 | 1.050 | 1.130 |
| Quercus robur | 2760 | 0.8990 | 0.9010 | 0.6380 | 0.6460 | 0.7500 | 0.9860 | 1.260 | 1.080 |
| Betula pendula | 2000 | 0.7750 | 0.7440 | 0.6080 | 0.7220 | 0.7230 | 0.9390 | 1.220 | 0.5850 |
| Larix decidua | 1420 | 0.8130 | 0.7950 | 0.6310 | 0.6910 | 0.7540 | 0.8800 | 0.9720 | 0.01050 |
| Carpinus betulus | 1060 | 0.8470 | 0.8500 | 0.7110 | 0.4750 | 0.7390 | 1.060 | 1.340 | 0.8280 |
| Alnus glutinosa | 805.0 | 0.9320 | 0.9150 | 0.8750 | 0.6320 | 0.7620 | 0.9370 | 1.090 | 3.080 |
| Larix kaempferi | 739.0 | 0.8680 | 0.8200 | 0.6860 | 0.6620 | 0.7790 | 1.020 | 0.9710 | 0.01750 |
| Fraxinus excelsior | 557.0 | 0.8200 | 0.7960 | 0.5360 | 0.6400 | 0.7030 | 1.280 | 1.270 | 0.9420 |
| Acer pseudoplatanus | 365.0 | 0.6930 | 0.7920 | 0.2850 | 0.7180 | 0.7420 | 0.5460 | 2.170 | 1.200 |
| Prunus avium | 320.0 | 0.8080 | 0.8490 | 0.6790 | 0.6520 | 0.7090 | 0.7060 | 1.180 | 0.4470 |
| Abies alba | 141.0 | 0.8850 | 0.9290 | -0.08820 | 0.6340 | 0.5970 | 1.020 | 4.530 | 3.370 |

## Bare-ground equilibrium (mean of the last 100 of 1000 years)

| loc | total_ba_ha | drift_last_vs_previous_100 | n_species_above_1 | dominant | dominant_share | stems_ha |
|---|---|---|---|---|---|---|
| 7.5 C, 998 mm, 504 m | 40.10 | -0.05150 | 2.000 | Fagus sylvatica | 0.9180 | 552.0 |
| 8.4 C, 693 mm, 281 m | 43.50 | 0.03410 | 3.000 | Fagus sylvatica | 0.8630 | 495.0 |
| 8.7 C, 578 mm, 84 m | 42.70 | 0.08160 | 3.000 | Fagus sylvatica | 0.9000 | 487.0 |
| 8.9 C, 615 mm, 244 m | 41.90 | 0.003970 | 3.000 | Fagus sylvatica | 0.7810 | 593.0 |
| 9.3 C, 760 mm, 292 m | 44.60 | 0.09240 | 3.000 | Fagus sylvatica | 0.8140 | 597.0 |

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
