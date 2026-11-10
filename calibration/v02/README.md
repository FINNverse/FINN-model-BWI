# v0.2 round (zoo brief of 2026-11-10)

Audit: `AUDIT_v02.md` (issues, how the pipeline handles them, what is open),
evidence from `audit_residuals.R` in `audit/`.

## Where v0.2's parameters come from

Not from a separate zoo fit. The study (FINN-bwi) refits its three candidate
chains on `two-interval-v12` (the corrected env layout) under FINN 0.2.0 main
0848b60 with `fit(weight_scale = "raw")`; the laptop and Max's node run them.
Duplicate compute is not wanted. The user chooses the candidate; the chain with
v0.1.0's structure (linear env responses, single drought index `drought2`,
fitted recruitment cap, stage-2 epoch 700) is

- stage 1 `fits/2026-11-10-A-stage1-v12-drought2-mgmt-masked-finn0848b60`
- stage 2 `fits/2026-11-10-A-stage2-v12-drought2-ushape-lambda0.1-mortgrowth-masked-finn0848b60-raw`
- stage 3 `fits/2026-11-10-A-stage3-v12-drought2-ushape-lambda0.1-mortgrowth-masked-regcap-ep700-finn0848b60-raw`

and the two others are the two-term drought chain (`...-v12-ushape-lambda0.1-...`,
env columns `drought` + `drought_lag1`; needs a `build_model()` variant with six
env columns) and the DNN-climate + linear-management chain
(`...-v12-drought2-ushape-dnnclim-linmgmt-...`; needs a network env response in
`build_model()`, i.e. a new variant and a larger parameter list).

## Steps once the chosen chain has its stage-3 `model.pt`

1. `Rscript calibration/v02/01_stage_v12_weights.R` (edit `stage*_dir` in its
   PARAMS if the user picks another chain) -> `FINN-zoo-data/bwi/weights/v0.2.0/`.
2. Set `"tag": "v0.2.0"` in `inst/model/weights.json` (and `finn_commit_fitted`
   to `0848b60`, `data_version` to `two-interval-v12`); `05_validation.R` and
   `06_extract_parameters.R` read the tag from there.
3. `OMP_NUM_THREADS=1 Rscript calibration/06_extract_parameters.R` ->
   `parameters.pt` + sha256 (goes into `weights.json`) + reproduction check on
   the pilot set.
4. `OMP_NUM_THREADS=1 Rscript calibration/05_validation.R` -> `validation/`
   (metrics, figures); copy `metrics.csv`, `metrics_species.csv` and the figures
   into the zoo folder `BWI-DE_staged-linear/` and update the card numbers.
5. `Rscript calibration/v02/audit_residuals.R` for the v0.2 residual tables.
6. Bump `DESCRIPTION` to 0.2.0, write the NEWS entry (data version, FINN commit,
   holdout metrics, every change and its reason), commit; the coordinator tags
   and uploads the release asset.

No FINN package change is needed for the default candidate.

## Zoo-specific improvements in this round (no refit)

- `simulate_tiers()` (R/tiers.R): runs the same sites under tiers 0-3 with the
  same seed and returns stand totals (or the species table) per tier, so the
  management range the card asks for is one call.
- `drought_index()` (R/covariates.R): the calibration's `drought2` from monthly
  soil moisture a user extracts at new sites (growing-season mean, per-site
  z over 2000-2022, two-year mean), with the DWD source and licence in its help.
- `covariates.csv` and the vignette say how each covariate is derived; the
  quickstart has a tier-range example.
