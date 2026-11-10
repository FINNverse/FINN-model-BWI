# two-interval-v12

Written by `scripts/148_prepare_two_interval_v12_no_year0.R` at commit `d0c5488` on 2026-11-09.

v12 = v8 without env row `year 0`. Nothing else changes.

**Why.** Every earlier two-interval version laid the env grid out as year 0..10 (calendar t0..t0+10)
with the census at obs year 10. Year 0 was meant as the initial state, but FINN has no initial-state
env row: every env row drives one annual transition, and the initial state is `init_cohort_dt`.
`fit()` therefore scored dbh / ba / trees after 11 transitions for a 10-year interval; the rate
windows were already the correct 10 years. In v12 the census is at step 10.

**Calendar.** env year k = calendar t0 + k, k = 1..10. Simulated year k = calendar t0 + k.
Scripts reading v12 take the scored state at `sim_year = obs_year` (v8: `obs_year + 1`).

**Scaling.** `env_scales.csv` is v8's; every kept env value is v8's value unchanged.

**Files.** `env_dt.csv` and `env_unscaled_dt.csv` lose the year-0 rows; every other file is a byte
copy of v8. Checks: `v8_v12_env_check.csv`. Manifest: `../two-interval-v12.sha256`.

## PARAMS

- `version` = `two-interval-v12`
- `base_dir` = `outputs/prepared/two-interval-v8`
- `drop_year` = `0`
- `P` = `10`
- `env_files` = `env_dt.csv, env_unscaled_dt.csv`
