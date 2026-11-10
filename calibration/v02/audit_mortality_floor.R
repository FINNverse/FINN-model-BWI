#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINNmodelBWI | calibration/v02/audit_mortality_floor.R
##
## WHAT   Where FINN's simulated-mortality floor binds for this model. The
##        forward pass draws deaths with the predicted annual probability
##        clamped to [0.01, 0.99] (FINNverse/FINN#51; the per-tree loss is not
##        clamped), so no cohort loses fewer than 1 % of its stems per year in
##        simulateForest(). This script simulates (a) the pilot set at the
##        recorded tier and forced to tier 0 and (b) a bare-ground run at one
##        average site, and reports the share of species-site-years (weighted
##        by stems) whose predicted annual removal probability lies below 0.01,
##        i.e. where the floor replaces the model's rate.
## WHY    AUDIT_v02.md, section 13 (mortality floor).
## INPUT  DATA_DIR/weights/<tag>/parameters.pt; DATA_DIR/prepared/two-interval-v12/pilot/*
## OUTPUT calibration/v02/audit/mortality_floor.csv
## RUN    OMP_NUM_THREADS=1 Rscript calibration/v02/audit_mortality_floor.R   (from the FINNmodelBWI repo root, ~20 s, one thread)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0 = main 0848b60
suppressPackageStartupMessages({ library(FINN); library(torch); library(data.table) })
torch_set_num_threads(1L)
source("R/definition.R")
DATA_DIR <- Filter(dir.exists, c(file.path("..", "..", "FINN-zoo-data", "bwi"), file.path("..", "FINN-zoo-data", "bwi")))[1]
OUT <- file.path("calibration", "v02", "audit"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PARAMS <- list(
  release      = jsonlite::read_json(file.path("inst", "model", "weights.json"))$variants$linear$tag,
  prep_version = "two-interval-v12", set = "pilot",
  env_vars     = c("mat_c", "map_mm", "elev_m", "management", "drought2"),
  floor        = 0.01,                   # FINN's clamp of the drawn mortality probability (FINN#51)
  patch_size_ha = 0.1, dbh_binsize = 2, theta_sim_k = 3.5,
  bare_years   = 300L, bare_patches = 15L,
  bare_site    = c(mat_c = 8.2, map_mm = 890, elev_m = 390),   # about the calibration means (natural units)
  seed         = 1L
)
model <- build_model("linear", torch_load(file.path(DATA_DIR, "weights", PARAMS$release, "parameters.pt")), theta_sim_factor = PARAMS$theta_sim_k)
d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$set)
env <- fread(file.path(d, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init <- fread(file.path(d, "init_cohort_dt.csv")); species <- fread(file.path(d, "species_dt.csv")); scales <- fread(file.path(d, "env_scales.csv"))
zs <- function(v, x) (x - scales[variable == v, mean]) / scales[variable == v, sd]
ic <- makeInitCohorts(init, dbh_binsize = PARAMS$dbh_binsize, Nspecies = nrow(species))
share <- function(S, label) {   # stem-weighted share of species-site-years with predicted mort < floor, by species and overall
  S <- S[trees > 0 & is.finite(mort)]
  rbind(S[, .(run = label, species_name = "all species", stem_years = sum(trees), share_below_floor = sum(trees * (mort < PARAMS$floor)) / sum(trees),
              mean_mort = sum(trees * mort) / sum(trees), mean_mort_floored = sum(trees * pmax(mort, PARAMS$floor)) / sum(trees))],
        merge(S, species, by = "species")[, .(run = label, stem_years = sum(trees), share_below_floor = sum(trees * (mort < PARAMS$floor)) / sum(trees),
              mean_mort = sum(trees * mort) / sum(trees), mean_mort_floored = sum(trees * pmax(mort, PARAMS$floor)) / sum(trees)), by = species_name][order(-stem_years)])
}
sim <- function(e, init_cohort, patches = NULL) {
  FINN.seed(PARAMS$seed)
  args <- list(model, env = e, init_cohort = init_cohort, patch_size = PARAMS$patch_size_ha, device = "cpu")
  if (!is.null(patches)) args$patches <- patches
  S <- as.data.table(do.call(simulateForest, args)$wide$site)
  for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
  S
}
R1 <- share(sim(env, ic), "pilot, management as recorded")
R2 <- share(sim(copy(env)[, management := zs("management", 0)], ic), "pilot, every site tier 0")
eb <- data.table(siteID = 1L, year = seq_len(PARAMS$bare_years))
for (v in names(PARAMS$bare_site)) eb[, (v) := zs(v, PARAMS$bare_site[[v]])]
eb[, `:=`(management = zs("management", 0), drought2 = zs("drought2", 0))]
SB <- sim(eb[, c("siteID", "year", PARAMS$env_vars), with = FALSE], NULL, patches = PARAMS$bare_patches)
R3 <- share(SB[year > 100L], "bare ground, tier 0, years 101-300")
RES <- rbind(R1, R2, R3)
fwrite(RES, file.path(OUT, "mortality_floor.csv"))
print(RES[species_name %in% c("all species", "Fagus sylvatica", "Picea abies", "Pinus sylvestris", "Quercus petraea")], digits = 3)
message("done")
