#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINNmodelBWI | calibration/08_quickstart_check.R
##
## WHAT   The model card's Quickstart as a plain script, without installing the
##        package: rebuild the model from parameters.pt through
##        R/definition.R::build_model(), scale an environment table given in
##        NATURAL units with the shipped scaling constants, set the management
##        tier and the drought index, and run simulateForest() (a) from a given
##        stand and (b) from bare ground. Prints basal area per species.
## WHY    Checks that the released parameters load and run on the stated FINN
##        version and shows what a user has to supply. Seconds on one thread.
## INPUT  DATA_DIR/weights/v0.1.0/parameters.pt; inst/model/{env_scales.csv, species.csv}
## OUTPUT console only
## RUN    Rscript calibration/08_quickstart_check.R   (from the FINNmodelBWI repo root)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0 = main 0848b60
suppressPackageStartupMessages({ library(FINN); library(torch); library(data.table) })
torch_set_num_threads(1L)
source("R/definition.R")                 # build_model()

DATA_DIR <- Filter(dir.exists, c(file.path("..", "..", "FINN-zoo-data", "bwi"), file.path("..", "FINN-zoo-data", "bwi")))[1]
PARAMS <- list(
  release       = "v0.1.0",                # weights folder = release tag
  env_vars      = c("mat_c", "map_mm", "elev_m", "management", "drought2"),   # the model's env columns, in the formula order
  theta_sim_k   = 3.5,                     # recruit-draw dispersion factor (build_model default)
  patch_size_ha = 0.1,                     # the fitted patch area; keep it
  years         = 30L,
  site          = c(mat_c = 8.0, map_mm = 800, elev_m = 400),   # a lowland site in natural units (deg C, mm, m)
  tiers         = c(0, 3),                 # management tier: 0 = no regular harvest ... 3 = > 20 % of the start basal area per decade
  drought_years = 18:20,                   # a three-year drought (simulated years), drought2 = -1.5 sd of soil moisture
  seed          = 1L
)

## ---- 1. model from the released parameters -------------------------------------
model <- build_model("linear", torch_load(file.path(DATA_DIR, "weights", PARAMS$release, "parameters.pt")), theta_sim_factor = PARAMS$theta_sim_k)
species <- fread(file.path("inst", "model", "species.csv"))
scales  <- fread(file.path("inst", "model", "env_scales.csv"))
zs <- function(v, x) (x - scales[variable == v, mean]) / scales[variable == v, sd]   # natural -> model scale (what scale_env() does)

## ---- 2. environment: one site, 30 years, two management tiers ----------------------
## one env row per annual transition: row k drives simulated year k; the initial
## state is the init cohort (there is no year-0 env row).
env <- CJ(siteID = seq_along(PARAMS$tiers), year = seq_len(PARAMS$years))
for (v in names(PARAMS$site)) env[, (v) := zs(v, PARAMS$site[[v]])]
env[, management := zs("management", PARAMS$tiers[siteID])]
env[, drought2 := zs("drought2", fifelse(year %in% PARAMS$drought_years, -1.5, 0))]
setcolorder(env, c("siteID", "year", PARAMS$env_vars))

## ---- 3a. from a given stand: a spruce-beech stand, cohorts in cm dbh with stems per 0.1 ha
sp <- function(nm) species[name == nm, species]
init <- rbindlist(lapply(seq_along(PARAMS$tiers), function(s) data.table(
  siteID = s, patchID = 1L, species = c(sp("Picea abies"), sp("Picea abies"), sp("Fagus sylvatica")),
  dbh = c(35, 20, 25), trees = c(20, 15, 10))))
init[, treeName := paste(siteID, species, dbh, sep = "_")]
FINN.seed(PARAMS$seed)
sim <- simulateForest(model, env = env, init_cohort = makeInitCohorts(init, Nspecies = nrow(species)), patch_size = PARAMS$patch_size_ha, device = "cpu")
S <- as.data.table(sim$wide$site)
for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
S <- merge(S, species, by = "species")[, tier := PARAMS$tiers[siteID]]
cat("\nBasal area (m2/ha) of the given stand, years 1, 10, 20, 30, by management tier:\n")
print(dcast(S[year %in% c(1, 10, 20, 30), .(ba_ha = round(sum(ba) / PARAMS$patch_size_ha, 1)), by = .(tier, year)], tier ~ year, value.var = "ba_ha"))
cat("\nMean annual removal probability of Picea abies in the drought years vs before, by tier:\n")
print(S[name == "Picea abies" & trees > 0, .(before = mean(mort[year %in% 11:17]), drought = mean(mort[year %in% PARAMS$drought_years])), by = tier])

## ---- 3b. from bare ground: 20 patches, no initial stand ------------------------------
FINN.seed(PARAMS$seed)
bare <- simulateForest(model, env = env[siteID == 1L], init_cohort = NULL, patches = 20L, patch_size = PARAMS$patch_size_ha, device = "cpu")
B <- as.data.table(bare$wide$site)
for (cc in c("year", "species")) B[, (cc) := as.integer(as.character(get(cc)))]
B <- merge(B, species, by = "species")
cat(sprintf("\nBare ground after %d years (tier 0): total %.1f m2/ha; top species:\n", PARAMS$years, B[year == PARAMS$years, sum(ba)] / PARAMS$patch_size_ha))
print(B[year == PARAMS$years][order(-ba)][1:4, .(name, ba_ha = round(ba / PARAMS$patch_size_ha, 2), stems_ha = round(trees / PARAMS$patch_size_ha))])
cat(sprintf("\nFINN %s, torch %s\n", packageVersion("FINN"), packageVersion("torch")))
