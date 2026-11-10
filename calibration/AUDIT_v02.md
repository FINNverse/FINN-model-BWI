# AUDIT v0.2 — BWI data issues and how the zoo pipeline handles them

Written 2026-11-10 for the zoo v0.2 round (`FINN-zoo-data/V02_BRIEF.md`). Sources:
the FINN-bwi study (branch `inventory-simulator`: `docs/bwi_data_notes.md`,
`docs/STATUS_REVIEW_2026-11-08.md`, `docs/HANDOVER_2026-11-08.md`, the v5/v8/v12
`PREPARATION.md`, `scripts/148`) and the study's memory notes, read-only. Evidence
for the released v0.1.0 parameters: `validation/` (holdout, corrected v12 layout)
and `v02/audit/` (`v02/audit_residuals.R`: residuals by decade x tier, start-BA
class and species, change metrics against persistence, loss trends; sections 12
and 13; `v02/audit_mortality_floor.R`: where FINN's 1 % mortality floor binds). Status flags: **handled** (in
the pipeline and the released parameters), **handled by v0.2** (needs the v12
refit, which the study runs), **open** (no fix yet; stated as a limitation).

## 1. Env rows vs observation windows — handled by v0.2 (metrics already corrected)

FINN runs one annual transition per env row; the initial state is the init
cohort. The v8 grid (years 0..10) made `fit()` score the census after 11
transitions (`year_sequence` 11 for a 10-year interval); the rate windows
(growth, per-tree removal, recruit sum) were 10 years, so the stage-2 objective
(growth + per-tree removal only; dbh/ba/trees weighted 0) and the stage-3
recruit loss were unaffected. Found by the study on 2026-11-09 (`scripts/148`,
a certain-death pulse in env row 10 killed the stand at step 11).

- Pipeline: every script reads `two-interval-v12` (env years 1..10;
  `run_all.R` step 2m = local copy of `scripts/148`, 2n manifest check); the
  scoring windows are sim years 1..10 and the state at step 10; the quickstart,
  vignette, `covariates.csv` and the example stands use one env row per
  transition; `05_validation.R` asserts `nrow(env) == n_sites * 10` and
  `setequal(year, 1:10)`.
- Evidence: the same v0.1.0 parameters on v12 vs the v8 layout: BA R² 0.858 vs
  0.855, stems 0.680 vs 0.677, dbh 0.380 vs 0.345, growth 0.474 vs 0.461, AUC
  0.757 vs 0.758 — the extra step mostly cost dbh and growth.
- v0.2: parameters from the study's v12 refit (stage 2 was at epoch ~520/1000
  at 05:00 UTC on 2026-11-10). Check at extraction: `fit_summary_stage2.md`
  names `two-interval-v12`.
- Proposal for FINN (no change made): `fit()` could warn when the number of
  env rows before the first observation of a site differs from its
  `period_length`.

## 2. Ingrowth vs recruitment — handled (v8 fix carried into v12)

Angle-count inclusion grows with d², so 91 % of trees first tallied in 2012 were
established trees (median 2002 dbh 23.6 cm) that grew into the sample
("ingrowth trap", `bwi_data_notes.md` §5). v3-v7 counted the back-modelled
start state of laterally ingrown trees on top of the start tally (initial BA
1.26x BWI2's own tally) and the end state carried the start expansion forward
plus every new tree (end BA 1.34 / 1.20x the design estimate).

- Pipeline (v8 = v12 tables): initial stand = the start tally only; end state =
  a fresh expansion of the living end tally (4 / g(t1)); 617 decade-1
  site-periods whose start stand was only back-modelled ingrowth drop out
  (40,787 site-periods). The recruitment target is NOT the new-tree list but the
  regeneration circles (see 6), so recruits are never confused with ingrowth.
- Open: in 2012-2022 the archive has no back-modelled diameters, so "new tree"
  there is every new sample tree (not used as a target; FINN-bwi issue #1
  closed by v8's construction). The end-state target is a design estimate per
  site-period (noisier than the carried construction, unbiased).

## 3. Sample sizes vs expansion — handled

Densities are expanded (n_ha = 4 / g, largest-remainder allocation on a 0.1 ha
patch; basal-area conservation median |error| 1.1 %); `n_at_risk`, `growth_n`
count tallied trees; the per-tree mortality loss (`tree_obs_dt`, `fit(tree_data
=, loss_aggregation = c(mortality = "tree"))`) scores each tallied tree once
while its expansion weight sets the light. Evidence: 260,002 tallied trees in
train vs ~3.3 M expanded stems; AUC 0.76 is the realistic ceiling for decadal
NFI removal (`warmstart-freeze-calibration`).

## 4. Harvest and the management tier — handled as a covariate; open as a scenario tool

Removals (1.56 % yr⁻¹ of BA) exceed natural mortality (0.26) six-fold. The
target is total removal with management as an unspecific confounder (user's
rule; censoring regular harvest would bias natural mortality low because
pre-emptive felling takes dying trees). The tier (0 none, 1 <= 5 %, 2 5-20 %, 3
> 20 % of the start BA removed by regular harvest over the decade; 2012-2022
calamity salvage counted as death) is computed from the SAME decade's removals.

- Pipeline: `management` is a linear per-species covariate in growth, mortality
  and regeneration (`prep/02`: all tiers; `prep/04`: `USE_MGMT <- TRUE`, literal).
  The card reports every removal / drought result as the range over tiers 0-3;
  `simulate_tiers()` makes that one call for users.
- Evidence: mortality management slope 0.37 (spruce) to 0.85 (beech) logit per
  tier; forcing all corners to tier 0-3 moves the national spruce removal ratio
  drought decade / previous decade from 1.42 (recorded) to 1.64-1.73.
- Open: a forced tier is a sensitivity to a post-treatment covariate, not a
  harvest process (no trees are cut; the removal rate shifts). Tier M1 is
  near-empty (38 + 8 corners, BAF 4). The 2007 Kyrill salvage sits inside the
  2002-2012 tier (no harvest reason in BWI 2012). The earlier "no-management"
  recipe made big trees immortal (size parameter at the bound) because the tier
  had absorbed harvest: never drop the covariate without a harvest process.
- v0.2 candidate: unchanged (linear tier). A DNN-climate + linear-management
  chain exists in the study; the linear one is the default for the tool because
  forced tiers stay a pure shift.

## 5. Fates: freshly dead, CI2017 chaining, calamity salvage — handled

- Freshly dead trees (dead < 12 months, recorded like live sample trees with a
  flag: ACCESS `tot`, archive `damage_dead`) are deaths; a freshly dead new tree
  is not ingrowth; dead at the start survey = not at risk (v5 fix, in
  `scripts/91` + `71`).
- 25 % of 2012-2022 corners were remeasured in the 2017 carbon inventory, so
  trees lost by 2017 had no 2022 row: fates chained 2012 -> 2017 -> 2022
  (4,831 removals, 647 deaths restored). Decade-2 removal 2.64 -> 3.02 % yr⁻¹;
  observed spruce spike 1.22 -> 1.37.
- Harvest reasons 1-4 (storm, insects/drought, fire, other calamity) in 2022
  count as natural death (36,663 of 105,670 harvested stems); BWI 2012 has no
  harvest reason.
- Evidence of correctness in the pipeline: v8 tables reproduce the study byte
  for byte (`two-interval-v8.sha256`, 35 of 46 files verified for the staged
  sets; the remaining 11 are train_quarter and QC files not staged).

## 6. Recruitment target and dispersion — handled; level open

- Target: stems of 5.0-6.9 cm on the 1.75 m regeneration circles (archive
  classes 5+6), stock -> flux at 0.30 cm yr⁻¹ over the 1.9 cm band, recruits
  per ha per PERIOD (`makeObsData` units; factor 1.579 per stem ha⁻¹). The
  tally-threshold "new trees" are not the target (ingrowth trap). The off-V
  share of circle individuals is 15-16 %.
- Dispersion: a free NB theta collapsed (1 -> 0.003) and recruitment exploded in
  recruitment-only stages; fix in the pipeline = NB in circle-count space
  (`recruit_obs_weight = 6.09e-4`) with theta frozen at the stage-1 ML value
  floored at 0.01 (8 of 15 species floored); for simulation theta x 3.5
  (`theta_sim_factor`, from the between-circle covariance; `scripts/79`).
- Cap: `regeneration_saturation()` with K fitted per species (start 500
  recruits ha⁻¹ yr⁻¹) stops one-year bursts from the steep stage-1 climate
  slopes (pine +3.4 per sd previous-year soil moisture); without it ~5
  site-periods set the holdout stems R² (0.25 without, 0.64 with cap on v5).
- Evidence (holdout, v12): recruits pred/obs 0.95 overall; by species 0.5
  (pine, sessile oak, wild cherry), 1.3 (beech, spruce), 3 (alder, fir),
  0.01-0.02 (both larches); recruitment R² 0.016; stems pred/obs 1.23.
- Open: the level rests on the assumed 0.30 cm yr⁻¹ increment and is about 2x
  the decade-1 recruit flux implied by the angle-count sample (recruitment
  target variants v9-v11 live on the study's branch
  `recruitment-target-variants`, issue #2, not adopted); larch recruitment is
  missed (few circles); recruit composition is weak. v0.2 keeps the target
  (the refit uses v12's tables = v8's); a changed target would be a new data
  version and a minor bump.

## 7. Loss competition and the staged recipe — handled

A joint fit of all losses pulled the mortality climate slopes to the wrong sign
(state losses vs a weakly identified process; spruce drought slope -2.6 -> +1.2),
drove the mortality growth term to the +-10 bound and made the U-shape
collapse under "all losses". The pipeline is the staged recipe: stage 1 per
process on its own loss with common species parameters (L-BFGS); stage 2
growth + per-tree removal only (weights 0 on dbh/ba/trees/recruits), observed
recruitment fed in, L2 pull λ 0.1 of the env slopes to the prefit (intercepts
free), epoch 700 of 1000; stage 3 recruitment alone on the frozen demography
with theta frozen and λ 1. Evidence: spruce mortality drought2 slope stage 1
-0.71 -> stage 2 -0.69; management slopes within 0.02 of the prefit; sen
(U-shape) 1.5-5.0 with dbh of minimum mortality 19-150 cm.

- `weight_scale`: the stage-2 weights are 1/baseline from a probe fit; under
  FINN >= 0.2.0 `fit(weight_scale = "raw")` keeps them acting on raw losses
  (named in `03_fit_stage2.R`, `04_fit_stage3.R`, `weights.json`, the card).
- Open: `growth_light` at its upper bound 0.99 for alder, birch, larch, spruce,
  pine; `comp_height` at 0.70 for sessile oak (FINN #30). The mortality loss
  itself never improves on the stage-1 GLM, so the drought response is
  essentially the prefit's (defensible: it is the per-tree ML estimate).

## 8. Convergence — handled (epoch-700 convention); flagged

Loss trends of the released stage 2 (`loss_history_stage2.csv`, baseline-scaled,
100-epoch windows): growth 0.628 (epoch 1) -> 0.487 (700) -> 0.496 (1000),
mortality 0.887 -> 0.804 -> 0.856; last-100 vs previous-100 at epoch 1000:
growth -0.5 %, mortality +0.02 %; at epoch 700: growth +2.3 %, mortality -0.5 %
(mini-batch noise; the study's `scripts/122` found the drought response flat
from epoch 700 to 2000 and the holdout BA R² within 0.005). Stage 3: recruit
loss 0.0388 -> 0.0360, last-25 vs previous-25 +3.9 % (noise around a flat
level). No `param_history` is shipped (it sits in the study's stage-2
model objects); the stage-2 slope trajectories are in the card figure
`stage2_mortality_slopes.png`.

- Open / v0.2 option without refit cost: a cosine lr decay settled stage 2
  reproducibly in the study but at a worse holdout calibration (1.28), so it is
  not adopted. The v12 refit uses the same 1000-epoch / epoch-700 convention.

## 9. Evaluation — handled

Spatially blocked holdout (11 of 53 100-km blocks, whole tracts, both decades;
median 22 km to the nearest training tract), persistence as the null (BA R²
0.858 vs 0.825; stems 0.680 vs 0.435; stand totals 0.442 vs 0.282 and 0.478 vs
0.107), tree-weighted removal AUC with ties at half (0.757) and calibration
(0.932), bare-ground runs converge (40-45 m² ha⁻¹, drift -5 to +9 % per
century, beech-dominated). Change metrics (ΔBA, Δstems vs persistence) are in
section 12. Mortality target: `obs_dt$mort` is the
annualised TOTAL removal rate; `n_died` counts natural deaths only — removals
must be derived from `mort` (the first zoo validation used `n_died` and got AUC
0.726 / calibration 2.9; corrected 2026-11-09).

## 10. FINN interface gotchas — handled / for Max's side

Handled in the pipeline: `env_autoscale = FALSE` (the prefit coefficients are
on the prepared scale); `period_mode <- "per_site"` so `loss_aggregation =
"tree"` is scored; `optimizer = torch::optim_adam` subclass for the L2 pull,
tensors on the parameter device; `light_mask = TRUE` written out; `torch_load`
keeps old methods, so the release is a parameter list rebuilt by
`build_model()`; `makeInitCohorts()` needs a data.table (vignette uses
`fread`); `scale_env()` returns a data.table because FINN's env extraction
indexes columns the data.table way.

Open, for Max's side (no FINN change made here):
1. A fresh `finn()` model builds its env layers at the first
   `simulateForest()`/`fit()` and consumes random numbers there, so
   `FINN.seed()` reproduces a `torch_load()`ed model's draws but not a rebuilt
   one's (deterministic rates identical; decade totals within 0.3 %).
2. `makeInitCohorts()` fails on a data.frame ("unused argument by =").
3. `fit()` could check env rows against `period_length` (issue 1).
4. `recruit_obs_weight` is a model field (1 in a rebuilt model); confirm
   `simulateForest()` never reads it.
5. `regeneration()`'s hard +0.2 floor per species (`reg_floor` 1e-3 here) and
   `regen_mask`: not needed for BWI (no lumped "other" class; off-V trees are
   dropped at preparation), noted for the other zoo models.

## 11. Coordinates, species, covariates — handled

Positions = 1 km INSPIRE cells (positional error <= 700 m); elevation from the
inventory (`b0_tab.HoeheNN`), not a DEM at the anonymised position; WorldClim
1970-2000 normals at 30"; drought = DWD soil moisture (growing-season mean,
per-site z over 2000-2022, two-year mean `drought2`; current vs previous year
not separable in decadal data, r 0.95, hence one index). 15 species = 99 % of
the basal area on the kept corners; no "other" class. No person-related
columns anywhere (the field-crew column is never read).

## 12. Residuals of v0.1.0 on the holdout (v12 layout) — evidence

From `v02/audit_residuals.R` (`v02/audit/*.csv`, `fig_residuals.png`); 8,596
site-periods, ratios are sums predicted / observed.

**By decade and recorded tier** (`residuals_by_decade_tier.csv`; M1 has 4 + 1
corners and is ignored):

| decade | tier | n | removals | growth | recruits | BA end | stems end |
|---|---|---|---|---|---|---|---|
| 2002-12 | 0 | 2,100 | 1.27 | 0.98 | 0.66 | 0.99 | 1.15 |
| 2002-12 | 2 | 832 | 1.19 | 0.99 | 1.13 | 1.00 | 1.18 |
| 2002-12 | 3 | 1,239 | **0.79** | 1.00 | 1.48 | **1.31** | **1.57** |
| 2012-22 | 0 | 2,775 | **0.71** | 1.06 | 0.94 | 1.09 | 1.26 |
| 2012-22 | 2 | 518 | **1.57** | 1.10 | 1.30 | 0.96 | 1.15 |
| 2012-22 | 3 | 1,127 | 0.99 | 1.10 | 1.24 | 1.10 | 1.21 |

The linear tier term averages over two decades with different harvest
regimes: it over-predicts removals at tier 0 and under-predicts at tier 3 in
2002-2012 (Kyrill salvage sits in the 2002-2012 tiers), and under-predicts the
tier-0 removals of the drought decade by 29 % (calamity salvage is booked as
death on tier-0 corners) while over-predicting tier 2 by 57 %. Growth is within
2 % in decade 1 and 6-10 % too high in decade 2 for every tier (the drought
decade's growth loss is under-represented by one two-year index). End-of-decade
stems are 15-26 % too many at every tier (recruitment level, section 6) and 57 %
too many at tier 3 in decade 1 (too few removals there).

**By start basal area** (`residuals_by_ba_class.csv`, m² ha⁻¹ at the start):

| class | n | removals | growth | recruits | BA end | stems end | BA R² (persistence) |
|---|---|---|---|---|---|---|---|
| < 10 | 417 | 0.76 | 0.85 | 0.40 | 0.69 | 0.76 | 0.57 (0.57) |
| 10-20 | 1,485 | 0.85 | 1.05 | 0.76 | 0.98 | 1.19 | 0.81 |
| 20-30 | 2,391 | 0.89 | 1.08 | 1.24 | 1.06 | 1.30 | 0.86 |
| 30-40 | 2,149 | 0.95 | 1.08 | 1.57 | 1.09 | 1.31 | 0.88 |
| >= 40 | 2,154 | 0.96 | 0.94 | 1.02 | 1.11 | 1.24 | 0.86 |

Sparse stands (< 10 m² ha⁻¹, 5 % of the holdout; young or just-harvested
corners) get too little growth, recruitment and removal and end 31 % too
light; dense stands get 9-11 % too much basal area and 24-31 % too many stems.
The recruitment ratio runs from 0.4 to 1.6 across the classes: the light
response of regeneration (stage 3) is steeper than the data (more recruits
in dense stands, fewer in open ones than observed).

**By species and decade** (`residuals_by_species_decade.csv`, removals pred /
obs decade 1 | decade 2): spruce 0.81 | 0.92, sessile oak 0.75 | 0.78,
sycamore 0.60 | 0.49, wild cherry 0.88 | 0.61, hornbeam 1.42 | 0.86, ash
1.49 | 1.21, beech 1.05 | 1.03, pine 0.99 | 0.93. Growth decade 1 is 3-17 %
under for 12 species and decade 2 5-16 % over for 10 species (beech 0.97 |
1.16, pine 1.04 | 1.12, Douglas fir 0.91 | 1.14): the decade contrast in
growth is not carried by the covariates. Recruits: larches 0.005-0.03, pine
0.50 | 0.60, sessile oak 0.61 | 0.37, birch 0.75 | 0.38 (under); beech 1.60 |
1.12, spruce 1.10 | 1.49, pedunculate oak 0.90 | 1.86, fir and alder 2-3 x
(over). End stems: fir 4.3-4.7 x, sycamore 1.9-2.2 x, beech 1.46-1.56 x.

**Change metrics** (`change_metrics.csv`; per site-period totals, per ha per
decade; persistence predicts zero change):

| group | n | ΔBA obs | ΔBA pred | ΔBA R² | Δstems obs | Δstems pred | Δstems R² |
|---|---|---|---|---|---|---|---|
| all | 8,596 | -0.66 | +1.41 | 0.22 | -181 | -35 | 0.36 |
| decade 1, tier 0 | 2,100 | +6.45 | +6.00 | 0.01 | -13 | +111 | 0.02 |
| decade 1, tier 3 | 1,239 | -10.17 | -2.92 | -0.36 | -506 | -252 | 0.56 |
| decade 2, tier 0 | 2,775 | +0.56 | +3.40 | 0.04 | -104 | +51 | 0.18 |
| decade 2, tier 2 | 518 | -0.27 | -1.75 | -0.54 | -187 | -104 | 0.62 |
| decade 2, tier 3 | 1,127 | -8.62 | -6.31 | 0.12 | -300 | -206 | 0.55 |

The model explains 22 % of the between-site variance of the decade's
basal-area change and 36 % of the stem change (persistence: 0 by
construction); the mean BA change is +1.4 m² ha⁻¹ per decade predicted against
-0.7 observed, the mean stem change -35 ha⁻¹ against -181. At tier 0 the model
has stems increasing (+111, +51) where the data have them falling (-13,
-104): too much recruitment and too little removal in unmanaged stands;
at tier 3 it removes only 30-70 % of the observed basal-area loss. The
within-tier ΔBA R² is near zero or negative for tiers 0-2: the tier carries
most of the explained change, the species-level responses little.

**Loss trends** (`loss_trend.csv`, from the shipped histories): stage-2
growth loss 0.628 (epoch 1) -> 0.487 (700) -> 0.496 (1000), mortality 0.887 ->
0.804 -> 0.856 (baseline-scaled); last-100 vs previous-100 epochs: growth
-0.5 %, mortality +0.02 % at 1000; +2.3 % and -0.5 % at 700. Stage 3: 0.0388
-> 0.0360, last-25 vs previous-25 +3.9 %. Flat within mini-batch noise; the
study's convergence check (`scripts/122`, drought response flat from epoch
700 to 2000) stands.

What this suggests for v0.2 and later (no FINN change): the tier x decade
pattern is the known limit of a linear, decade-constant tier (a cause-resolved
removal target or a tier that counts salvage would be a data change, major
bump); the stems surplus and the recruitment-by-density gradient point at the
regeneration light response and the recruitment level (section 6), both
stage-3 matters; the decade-2 growth surplus suggests a growth drought response
that the two-year index under-states (the study's two-term chain is the
alternative candidate).

## 13. FINN's simulated-mortality floor (0.01 per step, FINNverse/FINN#51)

`simulateForest()` and the state losses of `fit()` draw deaths with the
predicted annual probability clamped to [0.01, 0.99] (`binomial_from_gamma`
guard); the per-tree loss is not clamped, so the released parameters were
fitted on the unclamped probabilities but every simulation applies the floor.
`v02/audit_mortality_floor.R` (`v02/audit/mortality_floor.csv`, pilot set,
stem-weighted share of species-site-years whose predicted annual removal lies
below 0.01):

| run | all species | beech | spruce | pine | sessile oak |
|---|---|---|---|---|---|
| pilot, management as recorded | 29 % (mean rate 3.05 % -> 3.16 % floored) | 54 % | 13 % | 23 % | 74 % |
| pilot, every site tier 0 | 51 % (1.50 % -> 1.68 %) | 100 % (0.55 % -> 1.00 %) | 26 % | 35 % | 88 % |
| bare ground, tier 0, years 101-300 | 83 % (1.08 % -> 1.39 %) | 100 % (0.62 % -> 1.00 %) | 0 % | 0 % | 0 % |

At the recorded tiers the floor hardly matters (mean removal 3.05 vs 3.16 %
yr⁻¹). In tier-0 simulations it binds for every beech stem-year and most oak
stem-years: the fitted unmanaged beech removal is 0.55 % yr⁻¹ (U-shaped
mortality with the minimum at 45 cm dbh), and the simulation turns it into
1.0 %. Consequences for users: tier-0 and tier-range results for beech-, oak-
and fir-dominated stands carry removal rates about twice the fitted ones
(the range over tiers 0-3 narrows from below); bare-ground equilibria are
beech-dominated under a 1 % yr⁻¹ floor and would be denser and older without
it (the 40-45 m² ha⁻¹ plateau is partly the floor, partly the absence of
density-dependent mortality, FINN #28). The per-tree AUC and calibration of
the validation are unaffected (they use the predicted probability), the
simulated stems and basal area of unmanaged stands are. When FINN#51's
`createProcess(mortality_floor =)` is available, `build_model()` should pass
a small floor with `mortality_draw = "poisson"` and the bare-ground and tier-0
figures be redone; until then the card states the floor as a limitation.

## 14. What v0.2 changes and what it does not

Changes: parameters from the v12 refit (10 transitions per decade), the
layout-corrected recipe, `simulate_tiers()` and `drought_index()` for users,
this audit. Unchanged: data tables (v12 = v8 tables), species set, covariates,
model structure, staged recipe, theta freeze, cap. Therefore a minor bump
(v0.2.0): same inputs, refit on corrected data.
