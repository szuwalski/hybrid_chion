#!/usr/bin/env Rscript

# Paired-station bootstrap of NMFS/BSFRF relative efficiency separately for
# Snow and Tanner crab and for males and females. This intentionally retains
# the empirical, potentially non-monotone size pattern used by the joint model.

input_file <- "data/bsfrf/snow/bsfrf_nmfs_area_numbers_at_size.rds"
output_dir <- "analyses/bsfrf_selectivity"
output_file <- file.path(output_dir,
  "bootstrap_selectivity_results_by_species_sex.rds")
n_boot <- 1000L
bootstrap_seed <- 20260920L
minimum_size <- 27.5
species <- c("SNOW", "TANNER")
sexes <- c("MALE", "FEMALE")
gears <- c("BSFRF_NEPHROPS", "NMFS_STANDARD_83_112")
sizes <- seq(minimum_size, 177.5, by = 5)

x <- readRDS(input_file)
station <- x$station_numbers_at_size
station <- station[station$SPECIES %in% species & station$SEX %in% sexes &
  station$GEAR %in% gears & station$SIZE_BIN_5MM %in% sizes, ]

combo <- expand.grid(SPECIES = species, GEAR = gears, SEX = sexes,
  SIZE_BIN_5MM = sizes, stringsAsFactors = FALSE)
combo$COMBO <- seq_len(nrow(combo))
station <- merge(station, combo,
  by = c("SPECIES", "GEAR", "SEX", "SIZE_BIN_5MM"), all.x = TRUE)

# Pair the two gears within each species-year station set.
stratum_id <- interaction(station$SPECIES, station$YEAR, drop = TRUE)
strata <- lapply(split(station, stratum_id), function(d) {
  stations <- sort(unique(d$STATION_ID))
  mat <- matrix(0, nrow = length(stations), ncol = nrow(combo))
  idx <- match(d$STATION_ID, stations)
  for (i in seq_len(nrow(d)))
    mat[idx[i], d$COMBO[i]] <- mat[idx[i], d$COMBO[i]] + d$ESTIMATED_NUMBERS[i]
  mat
})

totals_to_ratio <- function(total) {
  z <- cbind(combo, TOTAL = total)
  wide <- reshape(z[c("SPECIES", "SEX", "SIZE_BIN_5MM", "GEAR", "TOTAL")],
    idvar = c("SPECIES", "SEX", "SIZE_BIN_5MM"), timevar = "GEAR",
    direction = "wide")
  names(wide) <- sub("^TOTAL\\.", "", names(wide))
  wide$RATIO <- with(wide, ifelse(BSFRF_NEPHROPS > 0,
    NMFS_STANDARD_83_112 / BSFRF_NEPHROPS, NA_real_))
  wide
}

observed_ratio <- totals_to_ratio(Reduce(`+`, lapply(strata, colSums)))
set.seed(bootstrap_seed)
draws <- vector("list", n_boot)
for (b in seq_len(n_boot)) {
  total <- numeric(nrow(combo))
  for (mat in strata) {
    rows <- sample.int(nrow(mat), nrow(mat), replace = TRUE)
    total <- total + colSums(mat[rows, , drop = FALSE])
  }
  z <- totals_to_ratio(total)
  z$REPLICATE <- b
  draws[[b]] <- z[c("REPLICATE", "SPECIES", "SEX", "SIZE_BIN_5MM", "RATIO")]
}
draws <- do.call(rbind, draws)
row.names(draws) <- NULL

summarize_ratio <- function(z) {
  v <- z$RATIO[is.finite(z$RATIO)]
  if (!length(v)) return(data.frame(LOWER = NA, MEDIAN = NA, UPPER = NA,
    N_VALID = 0L))
  data.frame(LOWER = unname(quantile(v, .025)), MEDIAN = median(v),
    UPPER = unname(quantile(v, .975)), N_VALID = length(v))
}
groups <- split(draws, interaction(draws$SPECIES, draws$SEX,
  draws$SIZE_BIN_5MM, drop = TRUE))
summary <- do.call(rbind, lapply(groups, function(z) cbind(
  data.frame(SPECIES = z$SPECIES[1], SEX = z$SEX[1],
    SIZE_BIN_5MM = z$SIZE_BIN_5MM[1]), summarize_ratio(z))))
summary <- merge(summary,
  observed_ratio[c("SPECIES", "SEX", "SIZE_BIN_5MM", "RATIO")],
  by = c("SPECIES", "SEX", "SIZE_BIN_5MM"), all.x = TRUE)
names(summary)[names(summary) == "RATIO"] <- "OBSERVED_RATIO"

# Direct counts of measured specimens at common stations, by species and sex.
count_measured <- function(bio, haul) {
  eligible <- haul$HAUL_ID[haul$COMMON_STATION]
  d <- bio[bio$HAUL_ID %in% eligible & bio$SPECIES %in% species &
    bio$SEX %in% sexes & is.finite(bio$SIZE_BIN_5MM) &
    bio$SIZE_BIN_5MM %in% sizes, ]
  aggregate(HAUL_ID ~ SPECIES + SEX + SIZE_BIN_5MM, d, length)
}
support <- rbind(count_measured(x$bsfrf_bio, x$bsfrf_haul),
  count_measured(x$nmfs_bio, x$nmfs_haul))
support <- aggregate(HAUL_ID ~ SPECIES + SEX + SIZE_BIN_5MM, support, sum)
names(support)[4] <- "N_MEASURED_SUPPORT"

saveRDS(list(
  metadata = list(n_boot = n_boot, seed = bootstrap_seed,
    resampling_unit = "paired station within species-year",
    minimum_size_mm = minimum_size),
  empirical_summary = summary,
  empirical_draws = draws,
  measured_support = support), output_file)
write.csv(summary, file.path(output_dir,
  "bootstrap_empirical_selectivity_by_species_sex.csv"), row.names = FALSE)
write.csv(support, file.path(output_dir,
  "measured_support_by_species_sex.csv"), row.names = FALSE)
cat("Saved", normalizePath(output_file, winslash = "/"), "\n")
