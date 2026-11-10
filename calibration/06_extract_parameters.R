#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINNmodelBWI | calibration/06_extract_parameters.R
##
## WHAT   Turns the fitted stage-3 model object (model.pt, a whole FINN model
##        saved with torch_save) into the release asset of this package:
##        parameters.pt, a NAMED LIST of parameter tensors that
##        R/definition.R::build_model() rebuilds the model from with whatever
##        FINN version is installed. Writes its sha256 next to it and checks
##        that the rebuilt model reproduces the saved one (same seed, same
##        simulation on the pilot set: identical growth / mortality / state,
##        recruitment identical in distribution because the same theta and seed
##        are drawn).
## WHY    A torch_load()ed FINN model carries the methods of the FINN version
##        it was saved under; a parameter list does not (model repo convention).
## INPUT  SOURCE_MODEL: DATA_DIR/fits/stage3/model.pt if the calibration chain
##        was run here, else DATA_DIR/weights/<release>/model.pt (staged copy)
##        DATA_DIR/prepared/two-interval-v8/pilot (the reproduction check)
## OUTPUT DATA_DIR/weights/<release>/parameters.pt (+ parameters.pt.sha256)
##        calibration/parameters_summary.csv  every parameter per species, human-readable
## RUN    Rscript calibration/06_extract_parameters.R   (from the FINNmodelBWI repo root, < 1 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0 = main 0848b60
suppressPackageStartupMessages({ library(FINN); library(torch); library(data.table) })
torch_set_num_threads(1L)
source("R/definition.R")                 # build_model(): the model definition this package ships

DATA_DIR <- Filter(dir.exists, c(file.path("..", "..", "FINN-zoo-data", "bwi"), file.path("..", "FINN-zoo-data", "bwi")))[1]
PARAMS <- list(
  release       = "v0.1.0",              # weights folder = the model repo's semver release tag (weights.json)
  prep_version  = "two-interval-v12",    # corrected layout (10 env rows per decade); the comparison needs only identical inputs for both models
  check_set     = "pilot",               # 1,975 site-periods: enough to compare the rebuilt and the saved model
  patch_size_ha = 0.1,
  dbh_binsize   = 2,
  theta_sim_k   = 3.5,                   # the simulation dispersion factor build_model() applies by default
  tol           = 1e-5,                  # max |difference| of the first year's growth and mortality rates (deterministic)
  tol_total     = 0.05,                  # relative difference of the end-of-decade basal-area and stem totals (random draws differ)
  seed          = 1L
)
W_DIR <- file.path(DATA_DIR, "weights", PARAMS$release)
src <- c(file.path(DATA_DIR, "fits", "stage3", "model.pt"), file.path(W_DIR, "model.pt"))
src <- src[file.exists(src)][1]
if (is.na(src)) stop("no stage-3 model.pt found (run 04_fit_stage3.R or 01_stage_data_and_weights.R)")
message("source model: ", src)
saved <- torch_load(src)
saved$eval()

## ---- 1. the named parameter list ---------------------------------------------------
cpu <- function(x) x$detach()$cpu()$clone()
layer <- function(nn) cpu(nn$parameters[[1]])
parameters <- list(
  competition        = cpu(saved$par_competition),          # species x 2: height allometry, competition strength
  growth             = cpu(saved$par_growth),               # species x 2: light response, size decline
  mortality          = cpu(saved$par_mortality),            # species x 3: light, size (linear), growth (per 1 % growth)
  regeneration_light = cpu(saved$par_regeneration),         # species x 1: light response of recruitment at 1 m
  theta              = cpu(saved$par_theta_recruits),       # recruit dispersion per species, as fitted on the circles
  env_growth         = layer(saved$nn_growth),              # species x 6: intercept + mat_c, map_mm, elev_m, management, drought2
  env_mortality      = layer(saved$nn_mortality),
  env_regeneration   = layer(saved$nn_regeneration),
  mort_sen_raw       = cpu(saved$mort_sen_raw),             # U-shape: sen = mort_sen_max * sigmoid(raw)
  mort_sen_max       = cpu(saved$mort_sen_max),
  mort_growth_scale  = cpu(saved$mort_growth_scale),        # growth enters the mortality logit as growth / 0.01
  reg_logK           = cpu(saved$reg_logK)                  # log of the Beverton-Holt recruitment cap, recruits per ha per year
)
stopifnot(all(sapply(parameters, inherits, "torch_tensor")))
dir.create(W_DIR, recursive = TRUE, showWarnings = FALSE)
## torch_save() is not byte-reproducible, so an existing file with the same
## values is kept (its sha256 is what weights.json and the card cite)
p_file <- file.path(W_DIR, "parameters.pt")
same_values <- function(a, b) identical(names(a), names(b)) &&
  all(mapply(function(x, y) all(dim(x) == dim(y)) && max(abs(as_array(x) - as_array(y))) == 0, a, b))
if (file.exists(p_file) && same_values(torch_load(p_file), parameters)) {
  message("parameters.pt exists with identical values: kept")
} else torch_save(parameters, p_file)
sha <- digest::digest(p_file, algo = "sha256", file = TRUE)
writeLines(sprintf("%s  parameters.pt", sha), file.path(W_DIR, "parameters.pt.sha256"))
message(sprintf("parameters.pt: %d bytes, sha256 %s", file.size(file.path(W_DIR, "parameters.pt")), sha))

## ---- 2. the rebuilt model reproduces the saved one -----------------------------------
## The deterministic rates of the first simulated year (growth and mortality
## probability from the observed initial stand) must agree exactly. Later years
## and the recruit counts differ by the random draws only: a freshly built model
## creates its env layers at the first simulation, which consumes random numbers,
## so the same seed does not give the same draw sequence as the saved object.
## The stand totals after the decade must agree within the draw-to-draw noise.
rebuilt <- build_model("linear", torch_load(p_file), theta_sim_factor = PARAMS$theta_sim_k)
saved$par_theta_recruits <- PARAMS$theta_sim_k * as.numeric(saved$par_theta_recruits$detach())   # the validation's simulation setting
d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$check_set)
env  <- fread(file.path(d, "env_dt.csv"))[, .(siteID, year, mat_c, map_mm, elev_m, management, drought2)]
init <- fread(file.path(d, "init_cohort_dt.csv"))
species <- fread(file.path(d, "species_dt.csv"))
ic <- makeInitCohorts(init, dbh_binsize = PARAMS$dbh_binsize, Nspecies = nrow(species))
run <- function(m) { FINN.seed(PARAMS$seed)
  S <- as.data.table(simulateForest(m, env = env, init_cohort = ic, patch_size = PARAMS$patch_size_ha, device = "cpu")$wide$site)
  for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
  S[order(siteID, year, species)] }
A <- run(saved); B <- run(rebuilt)
stopifnot(nrow(A) == nrow(B))
diff <- c(growth_year1 = max(abs(A[year == 1L, growth] - B[year == 1L, growth]), na.rm = TRUE),
          mort_year1   = max(abs(A[year == 1L, mort] - B[year == 1L, mort]), na.rm = TRUE),
          ba_total_end_rel    = abs(B[year == max(year), sum(ba)] / A[year == max(year), sum(ba)] - 1),
          trees_total_end_rel = abs(B[year == max(year), sum(trees)] / A[year == max(year), sum(trees)] - 1))
message("saved vs rebuilt on ", PARAMS$check_set, ": ", paste(sprintf("%s %.3g", names(diff), diff), collapse = ", "))
stopifnot(diff[["growth_year1"]] < PARAMS$tol, diff[["mort_year1"]] < PARAMS$tol,   # deterministic rates: exact
          diff[["ba_total_end_rel"]] < PARAMS$tol_total, diff[["trees_total_end_rel"]] < PARAMS$tol_total)

## ---- 3. human-readable parameter table -------------------------------------------------
mt <- function(x) as.matrix(as_array(x))
env_names <- c("intercept", "mat_c", "map_mm", "elev_m", "management", "drought2")
TAB <- cbind(species[order(species)],
             `colnames<-`(as.data.table(mt(parameters$competition)), c("comp_height", "comp_strength")),
             `colnames<-`(as.data.table(mt(parameters$growth)), c("growth_light", "growth_size")),
             `colnames<-`(as.data.table(mt(parameters$mortality)), c("mort_light", "mort_size", "mort_growth")),
             data.table(reg_light = as.numeric(mt(parameters$regeneration_light)),
                        theta_circle = as.numeric(as_array(parameters$theta)),
                        mort_sen = as.numeric(as_array(parameters$mort_sen_max)) * plogis(as.numeric(as_array(parameters$mort_sen_raw))),
                        reg_cap_K_per_ha_yr = exp(as.numeric(as_array(parameters$reg_logK)))),
             `colnames<-`(as.data.table(mt(parameters$env_growth)), paste0("growth_", env_names)),
             `colnames<-`(as.data.table(mt(parameters$env_mortality)), paste0("mort_", env_names)),
             `colnames<-`(as.data.table(mt(parameters$env_regeneration)), paste0("reg_", env_names)))
fwrite(TAB, file.path("calibration", "parameters_summary.csv"))
writeLines(c(sprintf("source_model: %s", basename(src)), sprintf("source_model_sha256: %s", digest::digest(src, algo = "sha256", file = TRUE)),
             sprintf("parameters_pt_sha256: %s", sha), sprintf("finn_version: %s", packageVersion("FINN")),
             sprintf("torch_version: %s", packageVersion("torch")), sprintf("extracted: %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
             sprintf("reproduction_check_%s: %s", names(diff), signif(diff, 3))),
           file.path("calibration", "parameters_extraction.txt"))
message("done")
