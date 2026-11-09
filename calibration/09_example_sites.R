#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINNmodelBWI | calibration/09_example_sites.R
##
## WHAT   Writes the example sites the package ships in inst/extdata: a few
##        holdout site-periods (never fitted) with their environment in NATURAL
##        units (year 0..10, the decade's observed drought index) and their
##        initial stand (cohorts: species, dbh, stems on the 0.1 ha patch), for
##        the quickstart vignette and load_model() tests.
## WHY    Users need a worked example of the input tables. Only what the BWI
##        licence allows is redistributed: no corner or tract identifiers, no
##        coordinates beyond the 1 km INSPIRE cell centroid (here rounded to
##        0.05 degrees), no person-related columns.
## INPUT  DATA_DIR/prepared/two-interval-v8/holdout_blocked/{env_unscaled_dt, init_cohort_dt, siteID_dt, species_dt}.csv
## OUTPUT inst/extdata/example_env.csv    siteID, year, mat_c, map_mm, elev_m, management, drought2, lon, lat
##        inst/extdata/example_init.csv   siteID, patchID, species, species_name, dbh, trees
## RUN    Rscript calibration/09_example_sites.R   (from the FINNmodelBWI repo root, seconds)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))
DATA_DIR <- Filter(dir.exists, c(file.path("..", "..", "FINN-zoo-data", "bwi"), file.path("..", "FINN-zoo-data", "bwi")))[1]
PARAMS <- list(
  prep_version = "two-interval-v8",
  set          = "holdout_blocked",       # never fitted
  interval     = 2L,                      # the 2012-2022 decade: its drought index carries 2018-2020
  n_sites      = 4L,                      # one per dominant species below
  dominants    = c("Picea abies", "Fagus sylvatica", "Pinus sylvestris", "Quercus petraea"),
  min_ba_ha    = 20, max_ba_ha = 45,      # ordinary closed stands
  min_tier     = 0L,                      # any management tier
  coord_round  = 0.05,                    # degrees; coarser than the 1 km INSPIRE cell the BWI publishes
  seed         = 7L
)
d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$set)
envu <- fread(file.path(d, "env_unscaled_dt.csv")); init <- fread(file.path(d, "init_cohort_dt.csv"))
sid  <- unique(fread(file.path(d, "siteID_dt.csv"))[, .(siteID, interval)]); sp <- fread(file.path(d, "species_dt.csv"))
init <- init[living == TRUE, .(siteID, patchID, species, dbh, trees)]   # init_cohort_dt carries many bookkeeping columns; keep what the example needs
stand <- merge(init, sp, by = "species")[, .(ba_ha = sum(trees * pi * (dbh / 200)^2) / 0.1, stems = sum(trees)), by = .(siteID, species_name)]
tot <- stand[, .(ba_ha = sum(ba_ha), dominant = species_name[which.max(ba_ha)], share = max(ba_ha) / sum(ba_ha)), by = siteID]
tot <- merge(tot, sid, by = "siteID")[interval == PARAMS$interval & ba_ha >= PARAMS$min_ba_ha & ba_ha <= PARAMS$max_ba_ha & share >= 0.6]
tot <- merge(tot, envu[year == 0L, .(siteID, management)], by = "siteID")[management >= PARAMS$min_tier]
set.seed(PARAMS$seed)
pick <- rbindlist(lapply(PARAMS$dominants, function(s) tot[dominant == s][sample(.N, 1L)]))
pick[, example := .I]
ENV <- merge(envu[siteID %in% pick$siteID], pick[, .(siteID, example)], by = "siteID")[
  , .(siteID = example, year, mat_c = round(mat_c, 2), map_mm = round(map_mm), elev_m = round(elev_m), management = as.integer(management),
      drought2 = round(drought2, 3), lon = round(x / PARAMS$coord_round) * PARAMS$coord_round, lat = round(y / PARAMS$coord_round) * PARAMS$coord_round)][order(siteID, year)]
INIT <- merge(merge(init[siteID %in% pick$siteID], pick[, .(siteID, example)], by = "siteID"), sp, by = "species")[
  , .(siteID = example, patchID = 1L, species, species_name, dbh = round(dbh, 1), trees = as.integer(trees))][order(siteID, species, -dbh)]
dir.create(file.path("inst", "extdata"), showWarnings = FALSE, recursive = TRUE)
fwrite(ENV, file.path("inst", "extdata", "example_env.csv"))
fwrite(INIT, file.path("inst", "extdata", "example_init.csv"))
print(merge(pick[, .(example, dominant, ba_ha = round(ba_ha, 1), management)], ENV[year == 0L, .(siteID, mat_c, map_mm, elev_m, lon, lat)], by.x = "example", by.y = "siteID"))
message(sprintf("wrote %d sites: %d env rows, %d cohorts", nrow(pick), nrow(ENV), nrow(INIT)))
