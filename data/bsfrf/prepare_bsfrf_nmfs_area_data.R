#!/usr/bin/env Rscript

# Build common-area, numbers-at-size comparisons for BSFRF Nephrops and
# standard NMFS EBS survey trawls. The output retains haul-level observations,
# station-level densities, and area-year estimates for snow and Tanner crab.

options(stringsAsFactors = FALSE)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
snow_dir <- file.path(root, "data", "bsfrf", "snow")
source_dir <- file.path(snow_dir, "data")
parental_cache <- file.path(root, "data", "cache", "crabpack_parental_raw.rds")
tanner_file <- file.path(root, "data", "bsfrf", "tanner",
                          "rda_SBS.TannerCrab.DataObjects.RData")
output_file <- file.path(snow_dir, "bsfrf_nmfs_area_numbers_at_size.rds")

required <- c(parental_cache, tanner_file,
              file.path(source_dir, paste0(c(2009, 2010, 2016, 2017, 2018), ".csv")),
              file.path(source_dir, "AllStation.csv"),
              file.path(source_dir, "AllStnLoc.csv"))
missing_files <- required[!file.exists(required)]
if (length(missing_files)) {
  stop("Missing required source files:\n", paste(missing_files, collapse = "\n"))
}

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
chr <- function(x) trimws(as.character(x))

normalize_station <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "#N/A")] <- NA_character_
  simple <- grepl("^[A-Z]+[- ]?[0-9]{1,2}$", x)
  ans <- x
  ans[simple] <- sub("^([A-Z]+)[- ]?([0-9]{1,2})$", "\\1-\\2", x[simple])
  pieces <- strsplit(ans[simple], "-", fixed = TRUE)
  ans[simple] <- vapply(pieces, function(z) paste0(z[1], "-", sprintf("%02d", as.integer(z[2]))), "")
  ans
}

haversine_nm <- function(lat1, lon1, lat2, lon2) {
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  3440.065 * 2 * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
}

first_nonmissing <- function(x) {
  ok <- !is.na(x)
  if (is.character(x)) ok <- ok & nzchar(trimws(x))
  if (!any(ok)) return(if (is.character(x)) NA_character_ else NA_real_)
  x[which(ok)[1]]
}

parental <- readRDS(parental_cache)

# The station grid is derived from the saved standard-survey haul pull.
grid_all <- rbind(
  parental$snow$haul[, c("STATION_ID", "MID_LATITUDE", "MID_LONGITUDE")],
  parental$tanner$haul[, c("STATION_ID", "MID_LATITUDE", "MID_LONGITUDE")]
)
grid_all$STATION_ID <- normalize_station(grid_all$STATION_ID)
grid_all <- grid_all[!is.na(grid_all$STATION_ID) &
                       is.finite(grid_all$MID_LATITUDE) &
                       is.finite(grid_all$MID_LONGITUDE), ]
grid <- aggregate(cbind(MID_LATITUDE, MID_LONGITUDE) ~ STATION_ID,
                  data = grid_all, FUN = mean, na.rm = TRUE)

assign_station <- function(station_raw, latitude, longitude, max_distance_nm = 15) {
  raw_norm <- normalize_station(station_raw)
  n <- length(raw_norm)
  assigned <- rep(NA_character_, n)
  method <- rep("unmatched", n)
  distance <- rep(NA_real_, n)

  exact <- match(raw_norm, grid$STATION_ID)
  has_exact <- !is.na(exact)
  assigned[has_exact] <- grid$STATION_ID[exact[has_exact]]
  method[has_exact] <- "station_id"
  distance[has_exact] <- 0

  candidates <- which(!has_exact & is.finite(latitude) & is.finite(longitude))
  for (i in candidates) {
    d <- haversine_nm(latitude[i], longitude[i], grid$MID_LATITUDE, grid$MID_LONGITUDE)
    j <- which.min(d)
    distance[i] <- d[j]
    if (is.finite(d[j]) && d[j] <= max_distance_nm) {
      assigned[i] <- grid$STATION_ID[j]
      method[i] <- "nearest_coordinate"
    } else {
      method[i] <- "coordinate_outside_15nm"
    }
  }
  data.frame(STATION_ID = assigned, STATION_MATCH_METHOD = method,
             STATION_MATCH_DISTANCE_NM = distance)
}

# Legacy lookup supplies coordinates where annual files only contain station labels.
legacy_station <- read.csv(file.path(source_dir, "AllStation.csv"), check.names = FALSE)
legacy_location <- read.csv(file.path(source_dir, "AllStnLoc.csv"), check.names = FALSE)
legacy_lookup <- data.frame(
  STATION_RAW = chr(legacy_station[[1]]),
  LATITUDE = num(legacy_location[[1]]),
  LONGITUDE = num(legacy_location[[2]])
)
legacy_lookup$STATION_KEY <- gsub("[^A-Z0-9]", "", toupper(legacy_lookup$STATION_RAW))

lookup_legacy_coords <- function(station_raw) {
  station_key <- gsub("[^A-Z0-9]", "", toupper(chr(station_raw)))
  idx <- match(station_key, legacy_lookup$STATION_KEY)
  data.frame(LATITUDE = legacy_lookup$LATITUDE[idx],
             LONGITUDE = legacy_lookup$LONGITUDE[idx])
}

fill_retow_station <- function(haul, station) {
  haul <- chr(haul)
  station <- chr(station)
  missing <- is.na(station) | station %in% c("", "NA", "N/A", "#N/A")
  retow <- missing & grepl("[Rr]$", haul)
  if (any(retow)) {
    base_haul <- sub("[Rr]$", "", haul[retow])
    station[retow] <- vapply(base_haul, function(z) {
      first_nonmissing(station[haul == z])
    }, "")
  }
  station
}

make_raw_bsfrf <- function(year, species) {
  d <- read.csv(file.path(source_dir, paste0(year, ".csv")), check.names = FALSE)

  if (year == 2009) {
    d[["NMFS station"]] <- fill_retow_station(d[["Haul#"]], d[["NMFS station"]])
    keep <- !is.na(num(d[["Aswept"]])) & num(d[["Aswept"]]) > 0
    keep <- keep & num(d[["Species"]]) == if (species == "SNOW") 1 else 2
    d <- d[keep, ]
    coords <- lookup_legacy_coords(d[["NMFS station"]])
    haul_id <- paste(year, "PILOT", chr(d[["Haul#"]]), sep = ";")
    out <- data.frame(
      YEAR = year, SPECIES = species, SOURCE_DESIGN = "pilot_area",
      SOURCE_FILE = "2009.csv", HAUL_ID = haul_id,
      STATION_RAW = chr(d[["NMFS station"]]), DATE = chr(d[["Date"]]),
      LATITUDE = coords$LATITUDE, LONGITUDE = coords$LONGITUDE,
      AREA_SWEPT_NM2 = num(d[["Aswept"]]), SEX_CODE = num(d[["Sex"]]),
      SIZE_MM = num(d[["CarapW"]]), SAMPLING_FACTOR = num(d[["Num Crab"]])
    )
  } else if (year == 2010) {
    if (species != "SNOW") return(NULL)
    d <- d[toupper(chr(d[["TowType"]])) == "SELECTIVITYSBS" &
             num(d[["Species"]]) == 1, ]
    coords <- lookup_legacy_coords(normalize_station(d[["NMFSstn"]]))
    haul_id <- paste(year, "SBS", chr(d[["Haul#"]]), sep = ";")
    out <- data.frame(
      YEAR = year, SPECIES = species, SOURCE_DESIGN = "selectivity_sbs",
      SOURCE_FILE = "2010.csv", HAUL_ID = haul_id,
      STATION_RAW = chr(d[["NMFSstn"]]), DATE = chr(d[["Date"]]),
      LATITUDE = coords$LATITUDE, LONGITUDE = coords$LONGITUDE,
      AREA_SWEPT_NM2 = num(d[["AreaSweptnm2"]]), SEX_CODE = num(d[["Sex"]]),
      SIZE_MM = num(d[["CarapW"]]), SAMPLING_FACTOR = num(d[["NumCrab"]])
    )
  } else if (year == 2016) {
    if (species != "SNOW") return(NULL)
    keep <- toupper(chr(d[["Study"]])) %in% c("SBS", "SBSPLUS") &
      toupper(chr(d[["Species"]])) == "OPI"
    d <- d[keep, ]
    out <- data.frame(
      YEAR = year, SPECIES = species,
      SOURCE_DESIGN = ifelse(toupper(chr(d[["Study_New16"]])) == "SBSPLUS" |
                               toupper(chr(d[["Study"]])) == "SBSPLUS", "sbsplus", "sbs"),
      SOURCE_FILE = "2016.csv", HAUL_ID = paste(year, chr(d[["Boat.Tow"]]), sep = ";"),
      STATION_RAW = chr(d[["NMFS_Stn"]]), DATE = chr(d[["Date"]]),
      LATITUDE = num(d[["MidTowLatitude"]]), LONGITUDE = num(d[["MidTowLongitude"]]),
      AREA_SWEPT_NM2 = num(d[["Aswept_nm2"]]), SEX_CODE = num(d[["Sex"]]),
      SIZE_MM = num(d[["Carapace"]]), SAMPLING_FACTOR = num(d[["SampFactor"]])
    )
  } else if (year == 2017) {
    if (species != "SNOW") return(NULL)
    keep <- toupper(chr(d[["IDX-SBS"]])) == "SBS" &
      toupper(chr(d[["Species"]])) == "OPI"
    d <- d[keep, ]
    out <- data.frame(
      YEAR = year, SPECIES = species, SOURCE_DESIGN = "sbs",
      SOURCE_FILE = "2017.csv",
      HAUL_ID = paste(year, chr(d[["Boat"]]), chr(d[["Haul"]]), "SBS", sep = ";"),
      STATION_RAW = chr(d[["Station"]]), DATE = chr(d[["Date"]]),
      LATITUDE = num(d[["Latitude"]]), LONGITUDE = num(d[["Longitude"]]),
      AREA_SWEPT_NM2 = num(d[["Aswept_nm2"]]), SEX_CODE = num(d[["Sex"]]),
      SIZE_MM = num(d[["Carapace"]]), SAMPLING_FACTOR = num(d[["SampFactor"]])
    )
  } else if (year == 2018) {
    if (species != "SNOW") return(NULL)
    keep <- toupper(chr(d[["Study"]])) == "SBS" &
      toupper(chr(d[["Species"]])) == "OPI"
    d <- d[keep, ]
    coords <- lookup_legacy_coords(d[["Station"]])
    out <- data.frame(
      YEAR = year, SPECIES = species, SOURCE_DESIGN = "sbs",
      SOURCE_FILE = "2018.csv",
      HAUL_ID = paste(year, chr(d[["Boat"]]), chr(d[["Haul"]]), "SBS", sep = ";"),
      STATION_RAW = chr(d[["Station"]]), DATE = chr(d[["Date"]]),
      LATITUDE = coords$LATITUDE, LONGITUDE = coords$LONGITUDE,
      AREA_SWEPT_NM2 = num(d[["Aswept_nm2"]]), SEX_CODE = num(d[["Sex"]]),
      SIZE_MM = num(d[["Carapace"]]),
      SAMPLING_FACTOR = num(d[["SUBSAMPLING FACTOR"]])
    )
  } else stop("Unsupported raw BSFRF year: ", year)

  out$SEX <- ifelse(out$SEX_CODE == 1, "MALE",
                    ifelse(out$SEX_CODE == 2, "FEMALE", "UNKNOWN"))
  out
}

make_raw_bsfrf_haul <- function(year, species) {
  d <- read.csv(file.path(source_dir, paste0(year, ".csv")), check.names = FALSE)

  if (year == 2009) {
    d <- d[!is.na(num(d[["Aswept"]])) & num(d[["Aswept"]]) > 0, ]
    d[["NMFS station"]] <- fill_retow_station(d[["Haul#"]], d[["NMFS station"]])
    coords <- lookup_legacy_coords(d[["NMFS station"]])
    out <- data.frame(
      YEAR = year, SPECIES = species, SOURCE_DESIGN = "pilot_area",
      SOURCE_FILE = "2009.csv",
      HAUL_ID = paste(year, "PILOT", chr(d[["Haul#"]]), sep = ";"),
      STATION_RAW = chr(d[["NMFS station"]]), DATE = chr(d[["Date"]]),
      LATITUDE = coords$LATITUDE, LONGITUDE = coords$LONGITUDE,
      AREA_SWEPT_NM2 = num(d[["Aswept"]])
    )
  } else if (year == 2010) {
    if (species != "SNOW") return(NULL)
    d <- d[toupper(chr(d[["TowType"]])) == "SELECTIVITYSBS" &
             is.finite(num(d[["AreaSweptnm2"]])) & num(d[["AreaSweptnm2"]]) > 0, ]
    coords <- lookup_legacy_coords(normalize_station(d[["NMFSstn"]]))
    out <- data.frame(
      YEAR = year, SPECIES = species, SOURCE_DESIGN = "selectivity_sbs",
      SOURCE_FILE = "2010.csv",
      HAUL_ID = paste(year, "SBS", chr(d[["Haul#"]]), sep = ";"),
      STATION_RAW = chr(d[["NMFSstn"]]), DATE = chr(d[["Date"]]),
      LATITUDE = coords$LATITUDE, LONGITUDE = coords$LONGITUDE,
      AREA_SWEPT_NM2 = num(d[["AreaSweptnm2"]])
    )
  } else if (year == 2016) {
    if (species != "SNOW") return(NULL)
    d <- d[toupper(chr(d[["Study"]])) %in% c("SBS", "SBSPLUS") &
             is.finite(num(d[["Aswept_nm2"]])) & num(d[["Aswept_nm2"]]) > 0, ]
    out <- data.frame(
      YEAR = year, SPECIES = species,
      SOURCE_DESIGN = ifelse(toupper(chr(d[["Study_New16"]])) == "SBSPLUS" |
                               toupper(chr(d[["Study"]])) == "SBSPLUS", "sbsplus", "sbs"),
      SOURCE_FILE = "2016.csv",
      HAUL_ID = paste(year, chr(d[["Boat.Tow"]]), sep = ";"),
      STATION_RAW = chr(d[["NMFS_Stn"]]), DATE = chr(d[["Date"]]),
      LATITUDE = num(d[["MidTowLatitude"]]), LONGITUDE = num(d[["MidTowLongitude"]]),
      AREA_SWEPT_NM2 = num(d[["Aswept_nm2"]])
    )
  } else if (year == 2017) {
    if (species != "SNOW") return(NULL)
    d <- d[toupper(chr(d[["IDX-SBS"]])) == "SBS" &
             is.finite(num(d[["Aswept_nm2"]])) & num(d[["Aswept_nm2"]]) > 0, ]
    out <- data.frame(
      YEAR = year, SPECIES = species, SOURCE_DESIGN = "sbs",
      SOURCE_FILE = "2017.csv",
      HAUL_ID = paste(year, chr(d[["Boat"]]), chr(d[["Haul"]]), "SBS", sep = ";"),
      STATION_RAW = chr(d[["Station"]]), DATE = chr(d[["Date"]]),
      LATITUDE = num(d[["Latitude"]]), LONGITUDE = num(d[["Longitude"]]),
      AREA_SWEPT_NM2 = num(d[["Aswept_nm2"]])
    )
  } else if (year == 2018) {
    if (species != "SNOW") return(NULL)
    d <- d[toupper(chr(d[["Study"]])) == "SBS" &
             is.finite(num(d[["Aswept_nm2"]])) & num(d[["Aswept_nm2"]]) > 0, ]
    coords <- lookup_legacy_coords(d[["Station"]])
    out <- data.frame(
      YEAR = year, SPECIES = species, SOURCE_DESIGN = "sbs",
      SOURCE_FILE = "2018.csv",
      HAUL_ID = paste(year, chr(d[["Boat"]]), chr(d[["Haul"]]), "SBS", sep = ";"),
      STATION_RAW = chr(d[["Station"]]), DATE = chr(d[["Date"]]),
      LATITUDE = coords$LATITUDE, LONGITUDE = coords$LONGITUDE,
      AREA_SWEPT_NM2 = num(d[["Aswept_nm2"]])
    )
  } else stop("Unsupported raw BSFRF year: ", year)

  unique(out)
}

# Load the curated Tanner object for Tanner-only work conducted in 2013-2018.
tanner_env <- new.env(parent = emptyenv())
load(tanner_file, envir = tanner_env)
tanner_obj <- tanner_env$obj

tanner_hd <- tanner_obj$dfrHD_BSFRF_SBS
tanner_id <- tanner_obj$dfrID_BSFRF_SBS
tanner_join <- merge(tanner_id, tanner_hd[, c("YEAR", "GIS_STATION", "HAULJOIN",
                                              "START_DATE", "MID_LATITUDE",
                                              "MID_LONGITUDE", "AREA_SWEPT_VARIABLE")],
                     by = "HAULJOIN", all.x = TRUE, sort = FALSE)
tanner_clean <- data.frame(
  YEAR = num(tanner_join$YEAR), SPECIES = "TANNER", SOURCE_DESIGN = "sbs",
  SOURCE_FILE = "rda_SBS.TannerCrab.DataObjects.RData",
  HAUL_ID = paste(num(tanner_join$YEAR), "TANNER_CLEAN", chr(tanner_join$HAULJOIN), sep = ";"),
  STATION_RAW = chr(tanner_join$GIS_STATION), DATE = chr(tanner_join$START_DATE),
  LATITUDE = num(tanner_join$MID_LATITUDE), LONGITUDE = num(tanner_join$MID_LONGITUDE),
  AREA_SWEPT_NM2 = num(tanner_join$AREA_SWEPT_VARIABLE),
  SEX_CODE = num(tanner_join$SEX_CODE), SIZE_MM = num(tanner_join$SIZE),
  SAMPLING_FACTOR = num(tanner_join$SAMPLING_FACTOR), SEX = chr(tanner_join$SEX)
)
tanner_clean$SEX[!tanner_clean$SEX %in% c("MALE", "FEMALE")] <- "UNKNOWN"

bsfrf_bio <- do.call(rbind, Filter(Negate(is.null), list(
  make_raw_bsfrf(2009, "SNOW"),
  make_raw_bsfrf(2009, "TANNER"),
  make_raw_bsfrf(2010, "SNOW"),
  make_raw_bsfrf(2016, "SNOW"),
  make_raw_bsfrf(2017, "SNOW"),
  make_raw_bsfrf(2018, "SNOW"),
  tanner_clean
)))

valid_bio <- is.finite(bsfrf_bio$AREA_SWEPT_NM2) & bsfrf_bio$AREA_SWEPT_NM2 > 0 &
  is.finite(bsfrf_bio$SIZE_MM) & bsfrf_bio$SIZE_MM > 0 &
  is.finite(bsfrf_bio$SAMPLING_FACTOR) & bsfrf_bio$SAMPLING_FACTOR > 0
bsfrf_bio$VALID_BIO_RECORD <- valid_bio

# One haul row is retained even where the target species, sex, or size bin has zero catch.
raw_haul <- do.call(rbind, Filter(Negate(is.null), list(
  make_raw_bsfrf_haul(2009, "SNOW"),
  make_raw_bsfrf_haul(2009, "TANNER"),
  make_raw_bsfrf_haul(2010, "SNOW"),
  make_raw_bsfrf_haul(2016, "SNOW"),
  make_raw_bsfrf_haul(2017, "SNOW"),
  make_raw_bsfrf_haul(2018, "SNOW")
)))
tanner_haul <- data.frame(
  YEAR = num(tanner_hd$YEAR), SPECIES = "TANNER", SOURCE_DESIGN = "sbs",
  SOURCE_FILE = "rda_SBS.TannerCrab.DataObjects.RData",
  HAUL_ID = paste(num(tanner_hd$YEAR), "TANNER_CLEAN", chr(tanner_hd$HAULJOIN), sep = ";"),
  STATION_RAW = chr(tanner_hd$GIS_STATION), DATE = chr(tanner_hd$START_DATE),
  LATITUDE = num(tanner_hd$MID_LATITUDE), LONGITUDE = num(tanner_hd$MID_LONGITUDE),
  AREA_SWEPT_NM2 = num(tanner_hd$AREA_SWEPT_VARIABLE)
)
bsfrf_haul_raw <- unique(rbind(raw_haul, tanner_haul))
haul_key <- paste(bsfrf_haul_raw$SPECIES, bsfrf_haul_raw$HAUL_ID, sep = "|")
haul_conflict_count <- vapply(split(bsfrf_haul_raw$STATION_RAW, haul_key),
                              function(z) length(unique(z[!is.na(z) & nzchar(z)])), 0L)
haul_conflict_keys <- names(haul_conflict_count)[haul_conflict_count > 1]
haul_metadata_conflicts <- bsfrf_haul_raw[haul_key %in% haul_conflict_keys, ]
bsfrf_haul <- do.call(rbind, lapply(split(bsfrf_haul_raw, haul_key), function(z) {
  data.frame(
    YEAR = z$YEAR[1], SPECIES = z$SPECIES[1],
    SOURCE_DESIGN = first_nonmissing(z$SOURCE_DESIGN),
    SOURCE_FILE = first_nonmissing(z$SOURCE_FILE), HAUL_ID = z$HAUL_ID[1],
    STATION_RAW = first_nonmissing(z$STATION_RAW), DATE = first_nonmissing(z$DATE),
    LATITUDE = first_nonmissing(z$LATITUDE), LONGITUDE = first_nonmissing(z$LONGITUDE),
    AREA_SWEPT_NM2 = first_nonmissing(z$AREA_SWEPT_NM2)
  )
}))
row.names(bsfrf_haul) <- NULL
bsfrf_haul <- bsfrf_haul[is.finite(bsfrf_haul$AREA_SWEPT_NM2) &
                           bsfrf_haul$AREA_SWEPT_NM2 > 0, ]
station_assignment <- assign_station(bsfrf_haul$STATION_RAW,
                                     bsfrf_haul$LATITUDE, bsfrf_haul$LONGITUDE)
bsfrf_haul <- cbind(bsfrf_haul, station_assignment)

bsfrf_bio <- merge(bsfrf_bio,
                   bsfrf_haul[, c("HAUL_ID", "SPECIES", "STATION_ID",
                                  "STATION_MATCH_METHOD", "STATION_MATCH_DISTANCE_NM")],
                   by = c("HAUL_ID", "SPECIES"), all.x = TRUE, sort = FALSE)
bsfrf_bio$SIZE_BIN_1MM <- floor(bsfrf_bio$SIZE_MM) + 0.5
bsfrf_bio$SIZE_BIN_5MM <- floor(bsfrf_bio$SIZE_MM / 5) * 5 + 2.5

make_nmfs <- function(species) {
  src <- parental[[tolower(species)]]
  years <- sort(unique(bsfrf_haul$YEAR[bsfrf_haul$SPECIES == species]))
  hd <- src$haul
  hd <- hd[hd$YEAR %in% years & hd$HAUL_TYPE == 3 & hd$PERFORMANCE == 0 &
             is.finite(hd$AREA_SWEPT) & hd$AREA_SWEPT > 0, ]
  hd_out <- data.frame(
    YEAR = hd$YEAR, SPECIES = species, SOURCE_DESIGN = "nmfs_standard",
    SOURCE_FILE = "crabpack_parental_raw.rds",
    HAUL_ID = paste(hd$YEAR, "NMFS", hd$HAULJOIN, sep = ";"),
    HAULJOIN = hd$HAULJOIN, STATION_RAW = hd$STATION_ID,
    STATION_ID = normalize_station(hd$STATION_ID), DATE = hd$START_DATE,
    LATITUDE = hd$MID_LATITUDE, LONGITUDE = hd$MID_LONGITUDE,
    AREA_SWEPT_NM2 = hd$AREA_SWEPT,
    STATION_MATCH_METHOD = "station_id", STATION_MATCH_DISTANCE_NM = 0
  )

  id <- src$specimen
  id <- id[id$YEAR %in% years, ]
  id_out <- data.frame(
    YEAR = id$YEAR, SPECIES = species,
    HAUL_ID = paste(id$YEAR, "NMFS", id$HAULJOIN, sep = ";"),
    HAULJOIN = id$HAULJOIN,
    SEX_CODE = id$SEX, SEX = ifelse(id$SEX == 1, "MALE",
                                    ifelse(id$SEX == 2, "FEMALE", "UNKNOWN")),
    SIZE_MM = id$SIZE, SAMPLING_FACTOR = id$SAMPLING_FACTOR
  )
  id_out <- id_out[id_out$HAUL_ID %in% hd_out$HAUL_ID, ]
  list(haul = hd_out, bio = id_out)
}

nmfs_snow <- make_nmfs("SNOW")
nmfs_tanner <- make_nmfs("TANNER")
nmfs_haul <- rbind(nmfs_snow$haul, nmfs_tanner$haul)
nmfs_bio <- rbind(nmfs_snow$bio, nmfs_tanner$bio)
nmfs_bio$SIZE_BIN_1MM <- floor(nmfs_bio$SIZE_MM) + 0.5
nmfs_bio$SIZE_BIN_5MM <- floor(nmfs_bio$SIZE_MM / 5) * 5 + 2.5

# Restrict both gears to the same station cells within each species-year.
coverage_keys <- unique(bsfrf_haul[!is.na(bsfrf_haul$STATION_ID),
                                   c("YEAR", "SPECIES", "STATION_ID")])
nmfs_keys <- unique(nmfs_haul[, c("YEAR", "SPECIES", "STATION_ID")])
common_station <- merge(coverage_keys, nmfs_keys,
                        by = c("YEAR", "SPECIES", "STATION_ID"))
common_station$COMMON_STATION <- TRUE

bsfrf_haul <- merge(bsfrf_haul, common_station, all.x = TRUE,
                    by = c("YEAR", "SPECIES", "STATION_ID"), sort = FALSE)
nmfs_haul <- merge(nmfs_haul, common_station, all.x = TRUE,
                  by = c("YEAR", "SPECIES", "STATION_ID"), sort = FALSE)
bsfrf_haul$COMMON_STATION[is.na(bsfrf_haul$COMMON_STATION)] <- FALSE
nmfs_haul$COMMON_STATION[is.na(nmfs_haul$COMMON_STATION)] <- FALSE

# Station areas: prefer the Tanner survey-design values, then the station's
# historical median, and finally a documented nominal 400 nm2 grid cell.
station_design <- tanner_obj$dfrSD_SBS
station_design$STATION_ID <- normalize_station(station_design$GIS_STATION)
station_design <- station_design[is.finite(station_design$STATION_AREA) &
                                   station_design$STATION_AREA > 0,
                                 c("YEAR", "STATION_ID", "STATION_AREA")]
station_design <- aggregate(STATION_AREA ~ YEAR + STATION_ID,
                            station_design, mean, na.rm = TRUE)
station_median <- aggregate(STATION_AREA ~ STATION_ID,
                            station_design, median, na.rm = TRUE)
names(station_median)[2] <- "HISTORICAL_STATION_AREA"

station_area <- merge(unique(common_station[, c("YEAR", "STATION_ID")]), station_design,
                      by = c("YEAR", "STATION_ID"), all.x = TRUE)
station_area <- merge(station_area, station_median, by = "STATION_ID", all.x = TRUE)
station_area$STATION_AREA_NM2 <- ifelse(is.finite(station_area$STATION_AREA),
                                        station_area$STATION_AREA,
                                        ifelse(is.finite(station_area$HISTORICAL_STATION_AREA),
                                               station_area$HISTORICAL_STATION_AREA, 400))
station_area$STATION_AREA_SOURCE <- ifelse(is.finite(station_area$STATION_AREA),
                                            "year_specific_tanner_design",
                                            ifelse(is.finite(station_area$HISTORICAL_STATION_AREA),
                                                   "historical_station_median",
                                                   "nominal_400_nm2"))
station_area <- station_area[, c("YEAR", "STATION_ID", "STATION_AREA_NM2",
                                 "STATION_AREA_SOURCE")]

size_bins <- seq(2.5, 177.5, by = 5)
sex_levels <- c("MALE", "FEMALE")

make_haul_numbers <- function(haul, bio, gear) {
  haul <- haul[haul$COMMON_STATION, ]
  bio <- bio[bio$SEX %in% sex_levels & is.finite(bio$SIZE_MM) &
               is.finite(bio$SAMPLING_FACTOR) & bio$SAMPLING_FACTOR > 0, ]
  bio$SIZE_BIN_5MM <- floor(bio$SIZE_MM / 5) * 5 + 2.5
  bio$SIZE_BIN_1MM <- floor(bio$SIZE_MM) + 0.5
  bio <- bio[bio$SIZE_BIN_5MM %in% size_bins, ]
  caught <- aggregate(SAMPLING_FACTOR ~ HAUL_ID + SPECIES + SEX + SIZE_BIN_5MM,
                      bio, sum, na.rm = TRUE)
  names(caught)[5] <- "EXPANDED_CATCH"

  grid_rows <- do.call(rbind, lapply(seq_len(nrow(haul)), function(i) {
    expand.grid(HAUL_ID = haul$HAUL_ID[i], SPECIES = haul$SPECIES[i],
                SEX = sex_levels, SIZE_BIN_5MM = size_bins,
                stringsAsFactors = FALSE)
  }))
  out <- merge(grid_rows, caught,
               by = c("HAUL_ID", "SPECIES", "SEX", "SIZE_BIN_5MM"), all.x = TRUE)
  out$EXPANDED_CATCH[is.na(out$EXPANDED_CATCH)] <- 0
  out <- merge(out, haul[, c("HAUL_ID", "SPECIES", "YEAR", "SOURCE_DESIGN",
                              "STATION_ID", "AREA_SWEPT_NM2")],
               by = c("HAUL_ID", "SPECIES"), all.x = TRUE)
  out$GEAR <- gear
  out$DENSITY_NM2 <- out$EXPANDED_CATCH / out$AREA_SWEPT_NM2
  out
}

bsfrf_numbers <- make_haul_numbers(bsfrf_haul, bsfrf_bio, "BSFRF_NEPHROPS")
nmfs_numbers <- make_haul_numbers(nmfs_haul, nmfs_bio, "NMFS_STANDARD_83_112")
haul_numbers_at_size <- rbind(bsfrf_numbers, nmfs_numbers)

station_numbers_at_size <- aggregate(
  cbind(EXPANDED_CATCH, DENSITY_NM2) ~ YEAR + SPECIES + GEAR + STATION_ID +
    SEX + SIZE_BIN_5MM,
  haul_numbers_at_size, mean, na.rm = TRUE
)
haul_counts <- aggregate(HAUL_ID ~ YEAR + SPECIES + GEAR + STATION_ID,
                         haul_numbers_at_size, function(x) length(unique(x)))
names(haul_counts)[5] <- "N_HAULS"
station_numbers_at_size <- merge(station_numbers_at_size, haul_counts,
                                 by = c("YEAR", "SPECIES", "GEAR", "STATION_ID"),
                                 all.x = TRUE)
station_numbers_at_size <- merge(station_numbers_at_size, station_area,
                                 by = c("YEAR", "STATION_ID"), all.x = TRUE)
station_numbers_at_size$ESTIMATED_NUMBERS <-
  station_numbers_at_size$DENSITY_NM2 * station_numbers_at_size$STATION_AREA_NM2

area_year_numbers_at_size <- aggregate(
  cbind(ESTIMATED_NUMBERS, STATION_AREA_NM2) ~ YEAR + SPECIES + GEAR + SEX + SIZE_BIN_5MM,
  station_numbers_at_size, sum, na.rm = TRUE
)
area_year_numbers_at_size$MEAN_DENSITY_NM2 <-
  area_year_numbers_at_size$ESTIMATED_NUMBERS /
  area_year_numbers_at_size$STATION_AREA_NM2

# Add same-area gear ratios. Infinite ratios are left explicit rather than hidden.
wide <- reshape(area_year_numbers_at_size[, c("YEAR", "SPECIES", "GEAR", "SEX",
                                               "SIZE_BIN_5MM", "MEAN_DENSITY_NM2")],
                idvar = c("YEAR", "SPECIES", "SEX", "SIZE_BIN_5MM"),
                timevar = "GEAR", direction = "wide")
names(wide) <- sub("MEAN_DENSITY_NM2\\.", "", names(wide))
wide$NMFS_TO_BSFRF_RATIO <- wide$NMFS_STANDARD_83_112 / wide$BSFRF_NEPHROPS
area_year_comparison <- wide

coverage <- aggregate(HAUL_ID ~ YEAR + SPECIES + GEAR,
                      haul_numbers_at_size, function(x) length(unique(x)))
names(coverage)[4] <- "N_HAULS"
station_coverage <- aggregate(STATION_ID ~ YEAR + SPECIES + GEAR,
                              haul_numbers_at_size, function(x) length(unique(x)))
names(station_coverage)[4] <- "N_COMMON_STATIONS"
coverage <- merge(coverage, station_coverage, by = c("YEAR", "SPECIES", "GEAR"))

station_crosswalk <- unique(bsfrf_haul[, c("YEAR", "SPECIES", "SOURCE_DESIGN",
                                            "STATION_RAW", "LATITUDE", "LONGITUDE",
                                            "STATION_ID", "STATION_MATCH_METHOD",
                                            "STATION_MATCH_DISTANCE_NM", "COMMON_STATION")])

qa <- list(
  invalid_bsfrf_bio_records = bsfrf_bio[!bsfrf_bio$VALID_BIO_RECORD, ],
  conflicting_bsfrf_haul_metadata = haul_metadata_conflicts,
  unmatched_bsfrf_hauls = bsfrf_haul[is.na(bsfrf_haul$STATION_ID), ],
  bsfrf_hauls_without_nmfs_station_year = bsfrf_haul[!bsfrf_haul$COMMON_STATION, ],
  station_area_sources = as.data.frame(table(station_area$STATION_AREA_SOURCE),
                                      stringsAsFactors = FALSE),
  source_reconciliation = data.frame(
    YEAR = 2009,
    SPECIES = c("SNOW", "TANNER"),
    CLEAN_VALID_EXPANDED_CATCH = c(
      sum(bsfrf_bio$SAMPLING_FACTOR[bsfrf_bio$YEAR == 2009 &
                                      bsfrf_bio$SPECIES == "SNOW" &
                                      bsfrf_bio$VALID_BIO_RECORD], na.rm = TRUE),
      sum(bsfrf_bio$SAMPLING_FACTOR[bsfrf_bio$YEAR == 2009 &
                                      bsfrf_bio$SPECIES == "TANNER" &
                                      bsfrf_bio$VALID_BIO_RECORD], na.rm = TRUE)
    ),
    PUBLISHED_PRELIMINARY_TOTAL = c(21459, 5784),
    NOTE = c(
      "Source CSV totals 21,364 after recorded sampling factors; differs from preliminary memo by -95.",
      "Source CSV totals 5,790; one record has missing size, leaving 5,789 valid; preliminary memo reports 5,784."
    )
  ),
  coverage = coverage
)

metadata <- list(
  description = paste(
    "Area-swept numbers-at-size for snow and Tanner crab caught by BSFRF",
    "Nephrops and standard NMFS 83/112 trawls in common NMFS station cells",
    "within the same year."
  ),
  created = as.character(Sys.time()),
  units = list(size = "mm carapace width", area = "square nautical miles",
               density = "expanded crab per square nautical mile",
               abundance = "expanded crab in common station footprint"),
  primary_filters = list(species = c("SNOW", "TANNER"), sex = sex_levels,
                         nmfs_haul_type = 3, nmfs_performance = 0,
                         common_station_year = TRUE,
                         bsfrf_2016 = c("SBS", "SBSplus")),
  size_bins = "5-mm bins with midpoint labels (e.g., 25-29.999 mm -> 27.5)",
  station_assignment = paste(
    "Normalized exact NMFS station ID where possible; otherwise nearest saved",
    "NMFS station center within 15 nautical miles using tow coordinates."
  ),
  station_area = paste(
    "Year-specific Tanner design area when available, otherwise station historical",
    "median, otherwise nominal 400 nm2. Inspect STATION_AREA_SOURCE before using",
    "absolute abundance; density comparisons do not depend on this assumption."
  ),
  exclusions = c(
    "Hybrids and non-target crab species",
    "Unknown sex from primary numbers-at-size tables",
    "Modified NMFS trawl from primary comparison",
    "2009 24-tow side-by-side workbook from primary area survey; 2009.csv pilot survey retained"
  ),
  sources = required
)

# Structural checks are deliberately fatal so a refreshed source cannot silently
# change keys, common areas, or zero-filled output structure.
if (anyDuplicated(paste(bsfrf_haul$SPECIES, bsfrf_haul$HAUL_ID)))
  stop("Duplicate BSFRF species-haul keys remain after metadata reconciliation")
if (anyDuplicated(paste(nmfs_haul$SPECIES, nmfs_haul$HAUL_ID)))
  stop("Duplicate NMFS species-haul keys")
station_key <- with(station_numbers_at_size,
                    paste(YEAR, SPECIES, GEAR, STATION_ID, SEX, SIZE_BIN_5MM))
if (anyDuplicated(station_key)) stop("Duplicate station-size output keys")
if (any(!is.finite(area_year_numbers_at_size$MEAN_DENSITY_NM2)) ||
    any(area_year_numbers_at_size$MEAN_DENSITY_NM2 < 0))
  stop("Invalid area-year density estimate")
area_check <- aggregate(STATION_AREA_NM2 ~ YEAR + SPECIES + GEAR,
                        area_year_numbers_at_size, unique)
area_wide <- reshape(area_check, idvar = c("YEAR", "SPECIES"),
                     timevar = "GEAR", direction = "wide")
if (any(abs(area_wide$STATION_AREA_NM2.BSFRF_NEPHROPS -
            area_wide$STATION_AREA_NM2.NMFS_STANDARD_83_112) > 1e-8))
  stop("Gear footprints differ within a species-year comparison")

result <- list(
  metadata = metadata,
  bsfrf_haul = bsfrf_haul,
  bsfrf_bio = bsfrf_bio,
  nmfs_haul = nmfs_haul,
  nmfs_bio = nmfs_bio,
  station_crosswalk = station_crosswalk,
  haul_numbers_at_size = haul_numbers_at_size,
  station_numbers_at_size = station_numbers_at_size,
  area_year_numbers_at_size = area_year_numbers_at_size,
  area_year_comparison = area_year_comparison,
  coverage = coverage,
  qa = qa
)

saveRDS(result, output_file, compress = "xz")

cat("Saved:", output_file, "\n")
cat("BSFRF haul rows:", nrow(bsfrf_haul), "\n")
cat("NMFS haul rows:", nrow(nmfs_haul), "\n")
cat("Common-area haul-size rows:", nrow(haul_numbers_at_size), "\n")
cat("Area-year comparison rows:", nrow(area_year_comparison), "\n")
print(coverage)
