# Release checklist: FINNmodelBWI v0.1.0

Prepared 2026-11-10. Nothing here has been executed; tag, push and release only
after the go. All commands from the repository root.

## Commit to tag

`da59470` on `main` ("FINNmodelBWI 0.1.0: FINN pinned at 0848b60, both release
assets ..."). It is the v0.1 package (`R/`, `inst/model/`, `DESCRIPTION` 0.1.0,
vignette); `calibration/` (incl. the v0.2 audit and the v12 staging under
`calibration/v02/`), `hpc/` and this file are excluded from the build by
`.Rbuildignore`. `R/definition.R` and `inst/model/` were not changed for v0.2
(`simulate_tiers()`, `drought_index()` and the v12 env-layout wording were
added for users of v0.1.0 and are part of it).

## Assets (outside git), sha256 verified 2026-11-10 with `shasum -a 256`

| asset | local path (relative to this repo) | sha256 | size |
|---|---|---|---|
| parameters.pt | ../../FINN-zoo-data/bwi/weights/v0.1.0/parameters.pt | dadd5d300cee0976e6dc7874b75070c50a12f2da3fb1c48e95e38f896a156da7 | 3,978 B |
| model.pt | ../../FINN-zoo-data/bwi/weights/v0.1.0/model.pt | 5c046997df7f882b0cd7d12e468f1eb9a2c89f61542d557a5f089611a8f2b857 | 54.5 MB |

`model.pt` is the stage-3 model object of the study fit under FINN 6705a50
(the fast-setup release candidate of 0.2.0), from which `parameters.pt` was
extracted (`calibration/06_extract_parameters.R`); `weights.json` records both
commits (`finn_commit_fitted` 6705a50, `finn_commit` 0848b60). The same sha256
values are in `inst/model/weights.json`; `load_model()` refuses a file whose
sha256 differs.

## Checks done (2026-11-10, one thread, FINN 0.2.0 = 0848b60, torch 0.15.1, R 4.5.0)

- `R CMD build` with the vignette (chunks not evaluated at build time) and
  `R CMD INSTALL` into a temporary library: OK.
- `load_model(path = <parameters.pt>)`: sha256 match, model rebuilt, 3-year
  bare-ground simulation with the pinned batchsize runs. `model.pt` sha256 match.
- Vignette renders without evaluation (as built) and, with
  `FINNMODELBWI_PARAMETERS=<parameters.pt>`, with every chunk evaluated
  (synthetic stands, the four real example stands, tiers 0-3, drought
  scenario; 7 s, no errors).

## Steps

```sh
# 0. sanity: the sha256 of the files to upload equal weights.json
shasum -a 256 ../../FINN-zoo-data/bwi/weights/v0.1.0/parameters.pt ../../FINN-zoo-data/bwi/weights/v0.1.0/model.pt

# 1. tag the release commit
git tag -a v0.1.0 da59470 -m "FINNmodelBWI 0.1.0: first fitted release (variant linear), parameters validated under FINN 0.2.0 commit 0848b60"

# 2. push main and the tag
git push origin main
git push origin v0.1.0

# 3. create the GitHub release with both assets
gh release create v0.1.0 --repo FINNverse/FINN-model-BWI \
  --title "FINNmodelBWI 0.1.0" \
  --notes-file RELEASE_NOTES_v0.1.0.md \
  ../../FINN-zoo-data/bwi/weights/v0.1.0/parameters.pt \
  ../../FINN-zoo-data/bwi/weights/v0.1.0/model.pt

# 4. verify the uploaded assets byte for byte
tmp=$(mktemp -d) && gh release download v0.1.0 --repo FINNverse/FINN-model-BWI -D "$tmp" && shasum -a 256 "$tmp"/*.pt

# 5. (user's go) make the repository public
gh repo edit FINNverse/FINN-model-BWI --visibility public --accept-visibility-change-consequences
```

## After the release

```r
remotes::install_github("FINNverse/FINN-model-BWI@v0.1.0")
library(FINNmodelBWI); torch::torch_set_num_threads(1L)
m <- load_model()             # downloads from the release URL in weights.json and checks the sha256
```

Then the FINN-zoo card (`BWI-DE_staged-linear/model-card.qmd`) cites tag
`v0.1.0` and these sha256 values.

## Release notes

The text of `RELEASE_NOTES_v0.1.0.md` (used by `--notes-file` above).
