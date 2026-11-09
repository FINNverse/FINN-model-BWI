# Prepared data: two-interval-v8

Written by `scripts/131_prepare_two_interval_v8.R` (a copy of `scripts/92_prepare_two_interval_v5.R`; decade-2 panel: `scripts/91_import_new_schema_chained.R`) at commit `3fcc67a` on 2026-11-08 18:19.

## What changed vs two-interval-v7 (and v5)

v8 = v7 with the lateral-ingrowth double count removed (FINN-bwi issue #1; assessment
`outputs/assessment/2026-11-08-ingrowth-issue1/README.md`, option (a) of section 6; numbers
`scripts/129_ingrowth_issue1_check.R`). v5 and v7 are not modified.

1. **Initial stand = the start tally only.** Decade 1's ACCESS2003 panel gives every tree that
   grew into the 2012 angle-count sample laterally (`ingrowth_lateral`, back-modelled 2002 dbh
   >= 7 cm) a back-modelled 2002 state with N_ha > 0. v5 expanded these trees into the 2002
   stand on top of the 2002 tally that already represents them. v8 drops them from the start
   allocation (`init_cohort_dt.csv`) and from `tree_obs_dt.csv` (v5's shading-only rows,
   `died = NA`). Decade 2's start stand was already the plain 2012 tally and is unchanged.
2. **End state = fresh expansion of the living end tally.** v5 built `ba`, `trees`, `dbh` at
   year 10 from the start-tally expansion carried forward plus every new tree at its end
   expansion, so laterally ingrown trees were counted twice. v8 expands every living tree of
   the end tally at its own end-survey n_ha (4 / g(t1); decade 1 Thuenen's `N_ha` of BWI3 with
   the edge correction), largest-remainder allocation per patch x species, in both decades.
   No recruit label and no dbh ceiling enter. The t0 and t1 virtual trees are no longer chained
   one-to-one; FINN's chained growth / mortality from `makeObsData()` was already replaced by
   the tallied-tree values, so nothing reads the chain.
3. **`drought2` (v7)** = (drought + drought_lag1) / 2 on the natural scale, pooled-scaled over
   the national set, last column of the env tables and last row of `env_scales.csv`, exactly
   as `scripts/124`; recomputed in this script (not copied from v7), checked against v7 below.
4. **Consequence: 617 decade-1 site-periods drop out** (v5 / v7: 41404 site-periods, v8: 40787). Their start stand expands to no tree at all once the lateral ingrowth is removed: 616 corners whose 2002 angle count holds no tree of the 15 species (v5's start stand there was back-modelled lateral ingrowth only, mean 16.1 m2/ha; end tally mean 20.7 m2/ha) and 1 whose few tallied trees round to 0 stems on the 0.1 ha patch. `makeObsData()` needs both censuses, so they leave like any site without a start stand. Nothing is added; every other site-period keeps v5's siteName, corner, years and set, but **siteIDs are renumbered 1..n** over the kept site-periods in v5's order (interval-2 siteIDs are now offset by 20092 instead of v5's 20709).

Kept exactly as v5 on the kept site-periods (checked below): corner filters and the dbh screen
(run on the full panel, lateral ingrowth included), fates, freshly dead / CI2017 / calamity
rules, management tiers, the regeneration target `reg` (circles), `growth`, `growth_n`, `mort`,
`n_at_risk`, `n_died` and the at-risk rows of `tree_obs_dt` (tallied trees), the split
(v5's `split_blocks.csv`, copied; every kept site-period stays in its v5 set), the natural-scale
env and the pooled-scaling rule (recomputed over v8's national set, so the scaled env values
shift slightly, table below). Two implementation details make
that hold: (i) the tallied-tree sample (`n_trees > 0` in the start allocation) is selected with
v5's start allocation (lateral ingrowth included) -- with the v8 allocation alone, 35 (decade 1) / 0 (decade 2) at-risk tallied trees would move in or out of the sample through the rounding of near-zero expansions; (ii) end-survey rows of dead and harvested trees keep their start expansion (they are not part of the end stand; only `makeObsData()`'s overwritten mort reads them).

## Checks

### Initial stand and end state, v8 vs v5, per set and decade

Per hectare, from the prepared tables: `init` = living trees of `init_cohort_dt` (expanded),
`tally` = `tree_obs_dt` (tallied trees at their stems), `end` = `obs_dt` at year 10, summed over
species; `design` = the end tally expanded per tallied tree at its n_ha (no rounding);
`b2_baeume` = BWI2's own 2002 tally (Av 0, dbh >= 7 cm, the 15 species) on the same corners.
Ratios are ratios of sums over the site-periods; `*_ba_*` columns are means in m2/ha.

| set | interval | site_periods | init_ba_v5 | init_ba_v8 | init BA v8/v5 | init stems v8/v5 | init BA v8/b2_baeume | tally BA v8/b2_baeume | tally BA v5/b2_baeume | end_ba_v5 | end_ba_v8 | end BA v5/v8 | end stems v5/v8 | end BA v8/design | end stems v8/design | median site abs(end BA v8/design - 1) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| national | 1 | 20092 | 40.46 | 32.58 | 0.8052 | 0.7380 | 0.9984 | 0.9982 | 1.240 | 44.80 | 33.58 | 1.334 | 1.404 | 0.9999 | 1.000 | 0.01153 |
| national | 2 | 20695 | 32.58 | 32.58 | 1.000 | 1.000 |  |  |  | 35.91 | 29.89 | 1.201 | 1.212 | 1.001 | 1.000 | 0.01488 |
| train | 1 | 15917 | 41.17 | 33.03 | 0.8022 | 0.7303 | 0.9985 | 0.9984 | 1.245 | 45.50 | 34.11 | 1.334 | 1.405 | 0.9999 | 0.9999 | 0.01202 |
| train | 2 | 16274 | 32.97 | 32.97 | 1.000 | 1.000 |  |  |  | 36.16 | 30.07 | 1.203 | 1.215 | 1.001 | 1.000 | 0.01511 |
| holdout_blocked | 1 | 4175 | 37.79 | 30.88 | 0.8173 | 0.7643 | 0.9979 | 0.9977 | 1.221 | 42.14 | 31.52 | 1.337 | 1.399 | 0.9998 | 1.000 | 0.009821 |
| holdout_blocked | 2 | 4421 | 31.13 | 31.13 | 1.000 | 1.000 |  |  |  | 34.98 | 29.25 | 1.196 | 1.204 | 0.9994 | 0.9997 | 0.01412 |
| train_quarter | 1 | 4029 | 40.94 | 32.84 | 0.8022 | 0.7281 | 0.9987 | 0.9984 | 1.245 | 45.12 | 33.81 | 1.334 | 1.411 | 0.9998 | 0.9998 | 0.01257 |
| train_quarter | 2 | 4011 | 33.14 | 33.14 | 1.000 | 1.000 |  |  |  | 36.75 | 30.67 | 1.198 | 1.208 | 1.001 | 1.000 | 0.01508 |
| pilot | 1 | 988 | 40.53 | 32.55 | 0.8032 | 0.7451 | 0.9985 | 0.9973 | 1.243 | 44.81 | 33.46 | 1.339 | 1.416 | 0.9993 | 0.9998 | 0.01240 |
| pilot | 2 | 987 | 32.60 | 32.60 | 1.000 | 1.000 |  |  |  | 36.33 | 30.52 | 1.190 | 1.192 | 1.002 | 1.000 | 0.01434 |

National: decade-1 init BA -19.5 % and stems -26.2 % against v5 (README expected about -21 % / -29 %); decade-1 init / BWI2's 2002 tally 0.998 (tallied 0.998; v5 tallied 1.240). End-state BA target: v5 / v8 = 1.334 (decade 1) and 1.201 (decade 2) (README: 1.34 / 1.20); v8 / design estimate 0.9999 / 1.0006 -- identical by construction up to the integer rounding of the expansion.

Open decision, not applied: the decade-1 corner filter of `scripts/02` (start BA >= 5 m2/ha, the 15 species >= 99 % of it) was evaluated on v5's lateral-inflated start stand. Evaluated on the start tally, 483 of v8's 20092 decade-1 site-periods hold < 5 m2/ha of the 15 species (decade 2's filter already sees its plain 2012 tally). v8 keeps the corner selection of v5.

### What must be unchanged, per set

On the kept site-periods; v5 / v7 rows are mapped onto v8's siteIDs by siteName (row order kept).

| check | national | train | holdout_blocked | train_quarter | pilot |
|---|---|---|---|---|---|
| site-periods v7 (= v5), decade 1 / decade 2 | 20709 / 20695 | 16412 / 16274 | 4297 / 4421 | 4154 / 4011 | 1011 / 987 |
| site-periods v8, decade 1 / decade 2 | 20092 / 20695 | 15917 / 16274 | 4175 / 4421 | 4029 / 4011 | 988 / 987 |
| v8 siteID_dt rows = v5 rows of the kept site-periods (all columns but siteID, siteID_national) | TRUE | TRUE | TRUE | TRUE | TRUE |
| species_dt.csv = v5 (md5) | TRUE | TRUE | TRUE | TRUE | TRUE |
| obs_dt rows: v5 on the kept site-periods = v8 | TRUE | TRUE | TRUE | TRUE | TRUE |
| obs_dt year, species identical to v5 | TRUE | TRUE | TRUE | TRUE | TRUE |
| obs_dt growth, growth_n identical to v5 | TRUE | TRUE | TRUE | TRUE | TRUE |
| obs_dt mort, n_at_risk, n_died identical to v5 | TRUE | TRUE | TRUE | TRUE | TRUE |
| obs_dt reg, period_length identical to v5 | TRUE | TRUE | TRUE | TRUE | TRUE |
| obs_dt ba, dbh, trees identical to v5 (FALSE expected) | FALSE | FALSE | FALSE | FALSE | FALSE |
| tree_obs_dt rows v5: at risk (died not NA) / shading only (died NA) | 325682 / 40137 | 260003 / 32749 | 65679 / 7388 | 64934 / 8249 | 15746 / 1974 |
| tree_obs_dt rows: v5 at risk on the kept site-periods / v8 | 325681 / 325681 | 260002 / 260002 | 65679 / 65679 | 64934 / 64934 | 15746 / 15746 |
| tree_obs_dt v8 = v5 at-risk rows of the kept site-periods (all columns, identical) | TRUE | TRUE | TRUE | TRUE | TRUE |
| env_unscaled_dt = v7 on the kept site-periods (identical; max abs diff) | TRUE; 0 | TRUE; 0 | TRUE; 0 | TRUE; 0 | TRUE; 0 |
| env_dt vs v7 on the kept site-periods (max abs diff; pooled rescaling) | 0.00791 | 0.00791 | 0.00787 | 0.00791 | 0.00791 |
| init_cohort_dt decade 2 = v5 (all columns, siteID mapped; national only) | TRUE | NA | NA | NA | NA |

The pooled scaling, v7 vs v8 (natural-scale env identical on the kept site-periods; mean / sd
recomputed over v8's national set, hence the small shift of every scaled column):

| variable | mean_v7 | sd_v7 | mean_v8 | sd_v8 |
|---|---|---|---|---|
| mat_c | 8.16376 | 0.955826 | 8.16498 | 0.956352 |
| map_mm | 887.161 | 267.584 | 886.859 | 267.601 |
| elev_m | 391.886 | 268.044 | 391.512 | 268.104 |
| x | 9.92530 | 2.06785 | 9.92610 | 2.06939 |
| y | 50.2273 | 1.87950 | 50.2289 | 1.87983 |
| drought | -0.0458773 | 0.958011 | -0.0508716 | 0.957729 |
| drought_lag1 | 0.0443740 | 0.952997 | 0.0391739 | 0.952817 |
| management | 0 | 1.00000 | 0 | 1.00000 |
| drought2 | -0.000751635 | 0.680639 | -0.00584884 | 0.680540 |

### Basal-area conservation of the expansions (national)

`wzp_qc_i{1,2}.csv` hold the start expansion (as in v5); the end expansion is new in v8.

| interval | stand | patches | BA expanded / design | median abs rel err | p95 abs rel err |
|---|---|---|---|---|---|
| 1 | start (2002 tally) | 38787 | 1.000 | 0.008790 | 0.06004 |
| 1 | end (2012 tally, fresh) | 39938 | 0.9999 | 0.01135 | 0.07164 |
| 2 | start (2012 tally) | 37939 | 1.000 | 0.01147 | 0.07200 |
| 2 | end (2022 tally, fresh) | 36464 | 1.000 | 0.01450 | 0.09240 |

### Per-tree removal vs obs_dt (tallied trees; unchanged from v5)

| interval | at-risk tallied trees | removal share, trees | removal share, obs_dt (per stem) |
|---|---|---|---|
| 1 | 157147 | 0.229 | 0.263 |
| 2 | 168534 | 0.264 | 0.290 |

## Which results this moves

Every fit on v5 / v6 / v7. The initial stand changes in decade 1 (it feeds the light of stage 1,
the simulation of stage 2 and the frozen demography of stage 3), and every end-state target
(`ba`, `trees`, `dbh` at year 10) changes in both decades, so the whole chain has to be refitted
and every state metric of the reports changes (the decade-2 holdout BA ratio of ~0.905 was the
double count). Mortality, growth, the regeneration target and the drought columns are the
same data on the kept site-periods (the dropped decade-1 site-periods had no at-risk tallied
tree of the 15 species, or very few); the scaled env columns shift slightly with the pooled
mean / sd; fitted values move mainly through the decade-1 light. v6 (CI2017 halves) has
the same construction and would need the same fix.

## What this is

FINN inputs for **two decades on the same corners**: interval 1 = BWI2 -> BWI3
(2002 -> 2012, ACCESS2003 release) and interval 2 = BWI3 -> BWI4 (2012 -> 2022,
Thuenen inventory archive). Annual timestep, observation at year 10 of each
decade, per-year drought anomaly in `env_dt$drought` (+ `drought_lag1`, `drought2`).
Each decade is its own FINN site-period; `siteID_dt.csv` maps `siteID` to `interval`
and the physical `corner` (Tnr_Enr). Interval-2 siteIDs are offset by 20092.

## Contents (`national/`; every set folder has the same FINN tables)

| file | content |
|---|---|
| obs_dt.csv | siteID, year, species, ba, dbh, trees (end state: fresh expansion of the living end tally), growth, growth_n, mort, n_at_risk, n_died (tallied trees), reg (circles), period_length |
| env_dt.csv | siteID, year 0..10, mat_c, map_mm, elev_m, x, y, management, drought, drought_lag1, drought2 -- z-scaled ONCE over both decades; management left on 0..3 |
| env_unscaled_dt.csv / env_scales.csv | natural scale and the pooled mean/sd |
| init_cohort_dt.csv | start-of-decade stand = the start tally (BWI2 for i1, BWI3 for i2), for makeInitCohorts() |
| species_dt.csv | the 15 species, one code map for both decades |
| siteID_dt.csv | siteID -> interval, corner, qualified siteName (`Tnr_Enr@i1/@i2`) |
| tree_obs_dt.csv | tallied trees of the start tally: siteID, year, patchID, species, dbh (start, cm), trees (stems on the 0.1 ha patch), growth (/yr, survivors), died (removal 0/1 over the decade; no NA rows in v8) |
| reg_target.csv | regeneration target per survey x corner x species (diagnostic) |
| wzp_qc_i{1,2}.csv | basal-area conservation of the start expansion |

## Per decade

| interval | sites | obs rows | mean ba (per 0.1 ha patch) | mean mort (/yr) | median growth (/yr) | mean reg (/ha/period) | drought coverage |
|---|---|---|---|---|---|---|---|
| 1 | 20092 | 301380 | 0.22 | 0.1088 | 0.0113 | 15.24 | 100.0 % |
| 2 | 20695 | 310425 | 0.20 | 0.1433 | 0.0094 | 14.13 | 100.0 % |

Management tiers (sites): interval 1 M0 9100, M1 38, M2 4451, M3 6503; interval 2 M0 12068, M1 8, M2 2882, M3 5737.

Lateral ingrowth that v5 carried in the decade-1 start stand: 20.6 % of its start basal area (all species of the panel, before the restriction to the 15).

Interval-2 corner funnel (as v5):

| step | corners |
|---|---|
| corners in the BWI3->BWI4 panel | 66657 |
| no untraceable trees | 65669 |
| basal area >= 5 m2/ha | 63650 |
| 15 species >= 99 % of BA | 55161 |
| realised period 8-12 yr | 54093 |
| environment available | 37939 |

Calamity salvage booked as natural death in interval 2 (harvest_reason 1,2,3,4): 36663 of 105670 harvested stems (34.7 %).
Regeneration from the archive circles, classes 5+6, r = 1.75 m, stock -> flux at 0.30 cm/yr over a 1.9 cm band (factor 1.579 per stem/ha). Off-V share of individuals: bwi2012 15.3 %, bwi2022 16.0 %.

Decisions, the freshly dead / CI2017 / calamity details and their counts are documented in
`two-interval-v5/PREPARATION.md`; they are unchanged here.

## Sets (spatially blocked, v5's split)

Every tract keeps its v5 set, block and train_quarter / pilot membership (`split_blocks.csv` is a byte copy of v5's): whole tracts, both decades together, 11 of 53 blocks held out. Every set folder has the national file layout with siteIDs renumbered 1..n (`siteID_dt$siteID_national` maps back) and the national pooled scaling.

| set | site-periods | interval 1 | interval 2 | tracts | blocks |
|---|---|---|---|---|---|
| train | 32191 | 15917 | 16274 | 9698 | 42 |
| holdout_blocked | 8596 | 4175 | 4421 | 2690 | 11 |
| train_quarter | 8040 | 4029 | 4011 | 2451 | 42 |
| pilot | 1975 | 988 | 987 | 594 | 41 |

Distance from each holdout tract to the nearest training tract: median 22.1 km, minimum 2.0 km, 10th percentile 6.0 km.

## Caveats

- **The end state is a design estimate, not FINN's accounting of the start stand.** A survivor
  stands for 4 / g(t0) stems at the start and 4 / g(t1) at the end (it grew, so fewer); the
  expected stand still balances (start stems - deaths + recruits = end stems in expectation,
  assessment section 4), but no longer per site. The state targets are therefore noisier per
  site-period than v5's carried construction and unbiased instead of inflated.
- **The regeneration target is unchanged** (circle stock -> flux at 0.30 cm/yr), and it is about
  2x the decade-1 WZP true-recruit flux (assessment section 3). With the inflated end state gone,
  an over-prediction of stems will become visible; that is the increment assumption, not ingrowth.
- **Regeneration geometry for 2022** is assumed unchanged (r = 1.75 m for classes 5+6), as v5.
- **Tier M1 is near-empty in both decades** (BAF 4: one tallied tree carries 4 m2/ha), as v5.
- Environment exists only for interval-1 corners, as v5.

## PARAMS

```
List of 22
 $ version          : chr "two-interval-v8"
 $ panel2_dir       : chr "outputs/new_schema_chained"
 $ v5_dir           : chr "outputs/prepared/two-interval-v5"
 $ v7_dir           : chr "outputs/prepared/two-interval-v7"
 $ patch_size_ha    : num 0.1
 $ period           : int 10
 $ min_ba_ha_t0     : num 5
 $ species_X        : num 0.99
 $ period_range     : int [1:2] 8 12
 $ dbh_growth_thresh: num [1:2] -10 60
 $ env_vars         : chr [1:5] "mat_c" "map_mm" "elev_m" "x" ...
 $ drought_var      : chr "z_soil_moist"
 $ drought2_parents : chr [1:2] "drought" "drought_lag1"
 $ tier_cut         : Named num [1:2] 0.05 0.2
 $ calamity_reasons : int [1:4] 1 2 3 4
 $ reg_classes      : int [1:2] 5 6
 $ reg_circle_r_m   : num 1.75
 $ reg_incr_cm_yr   : num 0.3
 $ reg_band_cm      : num 1.9
 $ sets             : chr [1:4] "train" "holdout_blocked" "train_quarter" "pilot"
 $ bwi4_dir         : chr "data/raw/bwi4"
 $ seed             : int 42
```
