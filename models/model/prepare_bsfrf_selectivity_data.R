#!/usr/bin/env Rscript

base_file <- "models/model/data/two_fishery_base_input.rds"
bootstrap_file <- "analyses/bsfrf_selectivity/bootstrap_selectivity_results_by_species_sex.rds"
output_file <- "models/model/data/model_input.rds"

d <- readRDS(base_file)
boot <- readRDS(bootstrap_file)
empirical <- boot$empirical_summary
empirical$SEX <- toupper(empirical$SEX)
support <- boot$measured_support

build_prior <- function(species, sex, sizes, anchor_size) {
  z <- empirical[empirical$SPECIES == species & empirical$SEX == sex &
                   is.finite(empirical$MEDIAN), ]
  z <- merge(z, support[c("SPECIES", "SEX", "SIZE_BIN_5MM",
                           "N_MEASURED_SUPPORT")],
             by = c("SPECIES", "SEX", "SIZE_BIN_5MM"), all.x = TRUE)
  z <- z[order(z$SIZE_BIN_5MM), ]
  lower_ratio <- approx(z$SIZE_BIN_5MM, pmax(z$LOWER, 0.0025),
                        xout = sizes, rule = 2)$y
  upper_ratio <- approx(z$SIZE_BIN_5MM,
                        pmax(z$UPPER, pmax(z$LOWER, 0.0025) * 1.01),
                        xout = sizes, rule = 2)$y
  log_sd <- pmin(1.0, pmax(0.01,
                           (log(upper_ratio) - log(lower_ratio)) / (2 * 1.96)))
  raw_log_sd <- pmin(1.0, pmax(0.01,
    (log(pmax(z$UPPER, 0.0025)) - log(pmax(z$LOWER, 0.0025))) / (2 * 1.96)))
  # Smooth the prior mean using inverse-bootstrap-variance weights. A moderate
  # spar retains the supported juvenile hump while preventing sparse terminal
  # bins from creating a spike in the prior center.
  reliable <- is.finite(z$N_MEASURED_SUPPORT) & z$N_MEASURED_SUPPORT >= 100
  center_smoother <- smooth.spline(
    z$SIZE_BIN_5MM[reliable], log(pmax(z$MEDIAN[reliable], 0.005)),
    w = 1 / raw_log_sd[reliable]^2, spar = 0.45
  )
  prediction_sizes <- pmin(pmax(sizes, min(z$SIZE_BIN_5MM[reliable])),
                           max(z$SIZE_BIN_5MM[reliable]))
  center_ratio <- exp(predict(center_smoother, prediction_sizes)$y)
  anchor <- match(anchor_size, sizes)
  stopifnot(!is.na(anchor))
  log_q_center <- rep(log(center_ratio[anchor]), 3)
  log_shape_center <- log(center_ratio) - log(center_ratio[anchor])
  # A five-dimensional natural-spline basis replaces the former free value at
  # every size bin.  Subtracting the anchor row makes every basis column zero
  # at the anchor, keeping overall scale and size shape separately identified.
  spline_basis <- splines::ns(sizes, df = 5, Boundary.knots = range(sizes))
  spline_basis <- sweep(spline_basis, 2, spline_basis[anchor, ], "-")
  shape_coef_center <- as.numeric(qr.solve(spline_basis, log_shape_center))
  size_z <- (sizes - anchor_size) / 30

  modern_index <- which(sizes %in% z$SIZE_BIN_5MM[reliable])
  modern_sizes <- sizes[modern_index]
  draws <- boot$empirical_draws
  draws$SEX <- toupper(draws$SEX)
  draws <- draws[draws$SPECIES == species & draws$SEX == sex &
                   draws$SIZE_BIN_5MM %in% modern_sizes,
                 c("REPLICATE", "SIZE_BIN_5MM", "RATIO")]
  draws$LOG_RATIO <- ifelse(draws$RATIO > 0, log(draws$RATIO), NA_real_)
  wide <- reshape(draws[c("REPLICATE", "SIZE_BIN_5MM", "LOG_RATIO")],
                  idvar = "REPLICATE", timevar = "SIZE_BIN_5MM",
                  direction = "wide")
  draw_columns <- paste0("LOG_RATIO.", modern_sizes)
  covariance <- cov(wide[, draw_columns, drop = FALSE],
                    use = "pairwise.complete.obs")
  # Small diagonal shrinkage stabilizes inversion without discarding the
  # cross-size dependence produced by paired-station resampling.
  covariance <- 0.95 * covariance + 0.05 * diag(diag(covariance))
  ee <- eigen(covariance, symmetric = TRUE)
  ee$values <- pmax(ee$values, max(ee$values) * 1e-6)
  covariance <- ee$vectors %*% (ee$values * t(ee$vectors))
  precision <- solve(covariance)

  # Bootstrap information is strongest for the modern gear/footprint era.
  era_multiplier <- c(3, 2, 1)
  shape_sd <- outer(log_sd, era_multiplier)
  q_base_sd <- log_sd[anchor]
  list(
    sizes = sizes,
    species = species,
    sex = sex,
    anchor_size = anchor_size,
    anchor_index = anchor,
    log_q_center = log_q_center,
    log_q_sd = q_base_sd * era_multiplier,
    log_shape_center = log_shape_center,
    spline_basis = unclass(spline_basis),
    shape_coef_center = shape_coef_center,
    size_z = size_z,
    historical_tilt_sd = c(0.40, 0.25),
    log_shape_sd = shape_sd,
    modern_index = modern_index,
    modern_sizes = modern_sizes,
    modern_log_center = log(center_ratio[modern_index]),
    modern_covariance = covariance,
    modern_precision = precision,
    era_labels = c("1975-1981", "1982-1988", "1989-present"),
    bootstrap_center_ratio = center_ratio,
    prior_center_smoothing_spar = 0.45,
    prior_center_min_measured_support = 100,
    bootstrap_log_sd = log_sd
  )
}

d$bsfrf_survey_prior <- list(
  snow_male = build_prior("SNOW", "MALE", d$male_sizes, 132.5),
  tanner_male = build_prior("TANNER", "MALE", d$male_sizes, 132.5),
  snow_female = build_prior("SNOW", "FEMALE", d$female_sizes, 102.5),
  tanner_female = build_prior("TANNER", "FEMALE", d$female_sizes, 102.5),
  source = normalizePath(bootstrap_file, winslash = "/", mustWork = TRUE),
  n_boot = boot$metadata$n_boot,
  minimum_size_mm = boot$metadata$minimum_size_mm,
  parameterization = paste(
    "log(selectivity) = era log scale + shared anchored 5-df natural spline",
    "+ historical linear size tilt;",
    "shape estimated separately by Snow/Tanner species and sex;",
    "hybrids use the midpoint of parental curves"
  ),
  era_definition = c("1975-1981", "1982-1988", "1989-present"),
  uncertainty_multiplier = c(3, 2, 1)
)
d$model <- "JMR_two_fishery_bsfrf_species_sex_selectivity"
d$model_note <- paste(
  "Two-fishery model with species- and sex-specific, three-era survey",
  "selectivity; modern bootstrap MVN calibration and hierarchical historical",
  "scale/tilt priors."
)

saveRDS(d, output_file)
cat("Saved", normalizePath(output_file, winslash = "/"), "\n")
for (group in c("snow_male", "tanner_male", "snow_female", "tanner_female")) {
  p <- d$bsfrf_survey_prior[[group]]
  cat(group, "anchor", p$anchor_size, "modern scale center",
      exp(p$log_q_center[3]), "modern log-SD", p$log_q_sd[3], "\n")
}
