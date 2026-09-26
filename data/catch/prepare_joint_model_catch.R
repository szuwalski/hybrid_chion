# Prepare snow and Tanner crab fishery data for the joint model.
#
# This script uses only base R. It writes:
#   1. annual_catch.csv       - annual numbers and biomass (metric tonnes)
#   2. size_composition.csv   - observed counts and proportions in model bins
#   3. joint_model_catch.rds  - the same two tables plus model-ready matrices
#
# In the source files, the folder identifies the species that was caught and
# target_stock identifies the fishery:
#   BSSC      = snow crab fishery
#   EBT/WBT   = Tanner crab fishery (east and west are combined here)
# Catch of a species in the other species' fishery is labelled "incidental".

input_dir <- "data/catch"
output_dir <- file.path(input_dir, "joint_model_ready")

# These are the male size-bin midpoints used by the current joint model.
model_sizes <- list(
  snow = seq(27.5, 132.5, by = 5),
  tanner = seq(27.5, 182.5, by = 5)
)

read_catch_file <- function(...) {
  read.csv(file.path(input_dir, ...), stringsAsFactors = FALSE)
}

fishery_name <- function(target_stock) {
  ifelse(target_stock == "BSSC", "snow",
         ifelse(target_stock %in% c("EBT", "WBT"), "tanner", NA_character_))
}

catch_type <- function(species, fishery) {
  ifelse(species == fishery, "directed", "incidental")
}

sum_or_na <- function(x) {
  if (length(x) == 0L || all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
}

aggregate_values <- function(x, value_names, group_names) {
  groups <- interaction(x[group_names], drop = TRUE, lex.order = TRUE)
  pieces <- split(x, groups)
  out <- lapply(pieces, function(z) {
    ans <- z[1, group_names, drop = FALSE]
    for (nm in value_names) ans[[nm]] <- sum_or_na(z[[nm]])
    ans
  })
  row.names(out) <- NULL
  do.call(rbind, out)
}

add_labels <- function(x, species) {
  x$species <- species
  x$fishery_name <- fishery_name(x$target_stock)
  x <- x[!is.na(x$fishery_name), , drop = FALSE]
  x$catch_type <- catch_type(x$species, x$fishery_name)
  x
}

# Put an integer-mm observation into a 5-mm joint-model bin. Observations below
# or above a species' modeled range are accumulated in the first or last bin.
assign_model_bin <- function(size, species) {
  bins <- model_sizes[[species]]
  raw_bin <- floor((size - 25) / 5) * 5 + 27.5
  pmin(pmax(raw_bin, min(bins)), max(bins))
}

# ---- Annual total catch ---------------------------------------------------

snow_total <- add_labels(
  read_catch_file("2026_final_snow", "total_catch.csv"), "snow"
)
# The joint model's fishery likelihood is male-only.
snow_total <- snow_total[snow_total$group != "female", ]
snow_total <- aggregate_values(
  snow_total,
  c("total_catch_n", "total_catch_wt"),
  c("crab_year", "species", "fishery_name", "catch_type")
)
names(snow_total)[names(snow_total) == "total_catch_wt"] <- "biomass_t"

tanner_total <- add_labels(
  read_catch_file("2026_final_tanner", "total_catch_sex_shell.csv"), "tanner"
)
tanner_total <- tanner_total[tanner_total$sex == "male", ]
tanner_total <- aggregate_values(
  tanner_total,
  c("total_catch_n", "total_catch_t"),
  c("crab_year", "species", "fishery_name", "catch_type")
)
names(tanner_total)[names(tanner_total) == "total_catch_t"] <- "biomass_t"

total_catch <- rbind(snow_total, tanner_total)
names(total_catch)[names(total_catch) == "total_catch_n"] <- "numbers"
total_catch$component <- "total"

# ---- Annual retained catch ----------------------------------------------

# In each retained-catch file, dir_retained_* is catch in that species' own
# fishery and inc_retained_* is catch in the other crab fishery. Before the
# directed/incidental fields were reported separately, retain the available
# tot_retained_* value for the directed fishery and leave incidental as NA.
prepare_retained <- function(x, species, biomass_suffix) {
  dir_n <- x$dir_retained_n
  inc_n <- x$inc_retained_n
  dir_biomass <- x[[paste0("dir_retained_", biomass_suffix)]]
  inc_biomass <- x[[paste0("inc_retained_", biomass_suffix)]]
  total_n <- x$tot_retained_n
  total_biomass <- x[[paste0("tot_retained_", biomass_suffix)]]

  no_split <- is.na(dir_n) & is.na(inc_n)
  dir_n[no_split] <- total_n[no_split]
  dir_biomass[no_split] <- total_biomass[no_split]

  directed <- data.frame(
    crab_year = x$crab_year,
    species = species,
    fishery_name = species,
    catch_type = "directed",
    numbers = dir_n,
    biomass_t = dir_biomass
  )
  incidental <- data.frame(
    crab_year = x$crab_year,
    species = species,
    fishery_name = if (species == "snow") "tanner" else "snow",
    catch_type = "incidental",
    numbers = inc_n,
    biomass_t = inc_biomass
  )
  rbind(directed, incidental)
}

snow_retained <- prepare_retained(
  read_catch_file("2026_final_snow", "retained_catch.csv"),
  species = "snow", biomass_suffix = "wt"
)
tanner_retained <- prepare_retained(
  read_catch_file("2026_final_tanner", "retained_catch.csv"),
  species = "tanner", biomass_suffix = "t"
)

# EBT and WBT rows can occur in the same year. Combine them after assigning
# Tanner retained catch to its directed or incidental fishery.
retained_catch <- aggregate_values(
  rbind(snow_retained, tanner_retained),
  c("numbers", "biomass_t"),
  c("crab_year", "species", "fishery_name", "catch_type")
)
retained_catch$component <- "retained"

annual_catch <- rbind(
  total_catch[, names(retained_catch)],
  retained_catch
)
annual_catch <- annual_catch[order(
  annual_catch$species, annual_catch$fishery_name,
  annual_catch$component, annual_catch$crab_year
), ]
row.names(annual_catch) <- NULL

# ---- Size compositions ---------------------------------------------------

snow_total_directed <- read_catch_file(
  "2026_final_snow", "directed_total_composition.csv"
)
snow_total_incidental <- read_catch_file(
  "2026_final_snow", "crab_bycatch_composition.csv"
)
snow_total_size <- rbind(snow_total_directed, snow_total_incidental)
snow_total_size <- add_labels(snow_total_size, "snow")
# Snow total-composition files include sex; the joint fishery model is male.
snow_total_size <- snow_total_size[snow_total_size$sex == 1, ]

snow_retained_size <- add_labels(read_catch_file(
  "2026_final_snow", "retained_catch_composition.csv"
), "snow")

tanner_total_size <- add_labels(read_catch_file(
  "2026_final_tanner", "observer_total_size_composition.csv"
), "tanner")
# These observations are assumed to be male-only. Shell condition is summed
# below when counts are aggregated by year, fishery, component, and size bin.

tanner_retained_size <- add_labels(read_catch_file(
  "2026_final_tanner", "retained_catch_size_composition.csv"
), "tanner")

standardize_size <- function(x, component) {
  x$component <- component
  x$model_size <- assign_model_bin(x$size, x$species[1])
  x[, c(
    "crab_year", "species", "fishery_name", "catch_type", "component",
    "model_size", "total"
  )]
}

size_composition <- rbind(
  standardize_size(snow_total_size, "total"),
  standardize_size(snow_retained_size, "retained"),
  standardize_size(tanner_total_size, "total"),
  standardize_size(tanner_retained_size, "retained")
)

size_composition <- aggregate_values(
  size_composition, "total",
  c("crab_year", "species", "fishery_name", "catch_type", "component", "model_size")
)
names(size_composition)[names(size_composition) == "total"] <- "sample_n"

composition_id <- interaction(
  size_composition[c("crab_year", "species", "fishery_name", "component")],
  drop = TRUE, lex.order = TRUE
)
annual_sample_n <- ave(size_composition$sample_n, composition_id, FUN = sum)
size_composition$proportion <- size_composition$sample_n / annual_sample_n
size_composition <- size_composition[order(
  size_composition$species, size_composition$fishery_name,
  size_composition$component, size_composition$crab_year,
  size_composition$model_size
), ]
row.names(size_composition) <- NULL

# Convert a long composition table into a year-by-size matrix. Missing years
# remain NA; observed bins with no crab are zero.
composition_matrix <- function(species, fishery, component) {
  z <- size_composition[
    size_composition$species == species &
      size_composition$fishery_name == fishery &
      size_composition$component == component,
  ]
  if (nrow(z) == 0L) return(NULL)
  years <- sort(unique(z$crab_year))
  sizes <- model_sizes[[species]]
  ans <- matrix(0, nrow = length(years), ncol = length(sizes),
                dimnames = list(year = years, size = sizes))
  ans[cbind(match(z$crab_year, years), match(z$model_size, sizes))] <- z$proportion
  ans
}

size_matrices <- list()
for (species in c("snow", "tanner")) {
  for (fishery in c("snow", "tanner")) {
    for (component in c("total", "retained")) {
      nm <- paste(species, fishery, component, sep = "_")
      size_matrices[[nm]] <- composition_matrix(species, fishery, component)
    }
  }
}

# ---- Checks and output ---------------------------------------------------

stopifnot(all(abs(aggregate(
  proportion ~ crab_year + species + fishery_name + component,
  size_composition, sum
)$proportion - 1) < 1e-10))
stopifnot(all(annual_catch$biomass_t >= 0, na.rm = TRUE))
stopifnot(all(annual_catch$numbers >= 0, na.rm = TRUE))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(annual_catch, file.path(output_dir, "annual_catch.csv"), row.names = FALSE)
write.csv(size_composition, file.path(output_dir, "size_composition.csv"), row.names = FALSE)

joint_model_catch <- list(
  annual_catch = annual_catch,
  size_composition = size_composition,
  size_matrices = size_matrices,
  model_sizes = model_sizes,
  notes = c(
    "All biomass is in metric tonnes.",
    "EBT and WBT are combined as the Tanner fishery.",
    "Annual total catch is male-only to match the joint fishery model.",
    "Snow total size compositions are male-only (sex code 1).",
    "Tanner total compositions are assumed male-only and shell conditions are summed.",
    "Retained composition files do not contain sex; retained crab are assumed male.",
    "Incidental retained totals use inc_retained_* from each species file.",
    "When the historical directed/incidental split is absent, total retained catch is assigned to the directed fishery and incidental retained catch remains NA."
  )
)
saveRDS(joint_model_catch, file.path(output_dir, "joint_model_catch.rds"))

cat("Wrote joint-model catch data to", output_dir, "\n")
cat("Annual catch rows:", nrow(annual_catch), "\n")
cat("Size-composition rows:", nrow(size_composition), "\n")
