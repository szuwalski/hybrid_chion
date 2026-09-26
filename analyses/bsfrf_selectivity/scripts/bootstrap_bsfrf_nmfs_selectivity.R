#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
input_file <- file.path(root, "data", "bsfrf", "snow",
                        "bsfrf_nmfs_area_numbers_at_size.rds")
output_dir <- file.path(root, "analyses", "bsfrf_selectivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_boot <- 1000L
bootstrap_seed <- 20260919L
minimum_size <- 27.5
gears <- c("BSFRF_NEPHROPS", "NMFS_STANDARD_83_112")
sexes <- c("MALE", "FEMALE")
sizes <- seq(minimum_size, 177.5, by = 5)

x <- readRDS(input_file)
station <- x$station_numbers_at_size
station <- station[
  station$SPECIES %in% c("SNOW", "TANNER") &
    station$SEX %in% sexes & station$GEAR %in% gears &
    station$SIZE_BIN_5MM %in% sizes, ]

# Fixed fitting weights are the direct numbers of specimens measured at size.
support_file <- file.path(output_dir,
                          "species_combined_pooled_selectivity_by_sex.csv")
support <- read.csv(support_file, stringsAsFactors = FALSE)
support <- support[support$SIZE_BIN_5MM >= minimum_size, ]
support$SEX <- toupper(support$SEX)

combo <- expand.grid(
  GEAR = gears, SEX = sexes, SIZE_BIN_5MM = sizes,
  stringsAsFactors = FALSE
)
combo$COMBO <- seq_len(nrow(combo))
station <- merge(station, combo,
                 by = c("GEAR", "SEX", "SIZE_BIN_5MM"), all.x = TRUE)

# One station-by-variable matrix per species-year. Sampling rows of these
# matrices keeps NMFS and BSFRF observations paired within a station cell.
stratum_id <- interaction(station$SPECIES, station$YEAR, drop = TRUE)
strata <- lapply(split(station, stratum_id), function(d) {
  stations <- sort(unique(d$STATION_ID))
  mat <- matrix(0, nrow = length(stations), ncol = nrow(combo),
                dimnames = list(stations, as.character(combo$COMBO)))
  idx_row <- match(d$STATION_ID, stations)
  for (i in seq_len(nrow(d))) {
    mat[idx_row[i], d$COMBO[i]] <-
      mat[idx_row[i], d$COMBO[i]] + d$ESTIMATED_NUMBERS[i]
  }
  mat
})

fit_logistic <- function(ratio_data) {
  d <- merge(ratio_data,
             support[c("SEX", "SIZE_BIN_5MM", "N_MEASURED_SUPPORT")],
             by = c("SEX", "SIZE_BIN_5MM"), all.x = TRUE)
  d$N_MEASURED_SUPPORT[is.na(d$N_MEASURED_SUPPORT)] <- 1
  d <- d[is.finite(d$RATIO), ]
  fits <- lapply(split(d, d$SEX), function(z) {
    objective <- function(par) {
      pred <- exp(par[3]) * plogis((z$SIZE_BIN_5MM - par[1]) / exp(par[2]))
      weighted.mean((z$RATIO - pred)^2, z$N_MEASURED_SUPPORT)
    }
    starts <- rbind(
      c(40, log(3), log(0.4)), c(80, log(15), log(0.7)),
      c(100, log(40), log(0.8)), c(160, log(70), log(1.5)),
      c(220, log(100), log(3))
    )
    candidates <- lapply(seq_len(nrow(starts)), function(i) {
      optim(starts[i, ], objective, method = "L-BFGS-B",
            lower = c(0, log(0.5), log(0.01)),
            upper = c(250, log(150), log(5)))
    })
    fit <- candidates[[which.min(vapply(candidates, `[[`, numeric(1), "value"))]]
    data.frame(
      SEX = z$SEX[1], ASYMPTOTE = exp(fit$par[3]), L50_MM = fit$par[1],
      SCALE_MM = exp(fit$par[2]),
      L95_MM = fit$par[1] + qlogis(0.95) * exp(fit$par[2]),
      WEIGHTED_MSE = fit$value, CONVERGENCE = fit$convergence
    )
  })
  do.call(rbind, fits)
}

totals_to_ratio <- function(total_vector) {
  z <- cbind(combo, TOTAL = total_vector)
  wide <- reshape(z[c("SEX", "SIZE_BIN_5MM", "GEAR", "TOTAL")],
                  idvar = c("SEX", "SIZE_BIN_5MM"), timevar = "GEAR",
                  direction = "wide")
  names(wide) <- sub("^TOTAL\\.", "", names(wide))
  wide$RATIO <- with(
    wide,
    ifelse(BSFRF_NEPHROPS > 0,
           NMFS_STANDARD_83_112 / BSFRF_NEPHROPS, NA_real_)
  )
  wide
}

observed_total <- Reduce(`+`, lapply(strata, colSums))
observed_ratio <- totals_to_ratio(observed_total)
observed_fit <- fit_logistic(observed_ratio)

set.seed(bootstrap_seed)
bootstrap_ratio <- vector("list", n_boot)
bootstrap_parameters <- vector("list", n_boot)

for (b in seq_len(n_boot)) {
  total <- numeric(nrow(combo))
  for (mat in strata) {
    sampled_rows <- sample.int(nrow(mat), nrow(mat), replace = TRUE)
    total <- total + colSums(mat[sampled_rows, , drop = FALSE])
  }
  rb <- totals_to_ratio(total)
  rb$REPLICATE <- b
  bootstrap_ratio[[b]] <- rb[c("REPLICATE", "SEX", "SIZE_BIN_5MM", "RATIO")]
  fb <- fit_logistic(rb)
  fb$REPLICATE <- b
  bootstrap_parameters[[b]] <- fb
}

bootstrap_ratio <- do.call(rbind, bootstrap_ratio)
bootstrap_parameters <- do.call(rbind, bootstrap_parameters)
row.names(bootstrap_ratio) <- NULL
row.names(bootstrap_parameters) <- NULL

summarize_quantiles <- function(z) {
  z <- z[is.finite(z)]
  if (!length(z)) return(c(LOWER = NA, MEDIAN = NA, UPPER = NA, N_VALID = 0))
  c(LOWER = unname(quantile(z, 0.025)), MEDIAN = median(z),
    UPPER = unname(quantile(z, 0.975)), N_VALID = length(z))
}

ratio_groups <- split(
  bootstrap_ratio,
  interaction(bootstrap_ratio$SEX, bootstrap_ratio$SIZE_BIN_5MM, drop = TRUE)
)
ratio_summary <- do.call(rbind, lapply(ratio_groups, function(z) {
  q <- summarize_quantiles(z$RATIO)
  data.frame(SEX = z$SEX[1], SIZE_BIN_5MM = z$SIZE_BIN_5MM[1],
             LOWER = q["LOWER"], MEDIAN = q["MEDIAN"],
             UPPER = q["UPPER"], N_VALID = q["N_VALID"])
}))
ratio_summary <- merge(
  ratio_summary,
  observed_ratio[c("SEX", "SIZE_BIN_5MM", "RATIO")],
  by = c("SEX", "SIZE_BIN_5MM"), all.x = TRUE
)
names(ratio_summary)[names(ratio_summary) == "RATIO"] <- "OBSERVED_RATIO"

parameter_names <- c("ASYMPTOTE", "L50_MM", "SCALE_MM", "L95_MM")
parameter_summary <- do.call(rbind, lapply(sexes, function(sex) {
  z <- bootstrap_parameters[bootstrap_parameters$SEX == sex &
                              bootstrap_parameters$CONVERGENCE == 0, ]
  do.call(rbind, lapply(parameter_names, function(parameter) {
    q <- summarize_quantiles(z[[parameter]])
    data.frame(
      SEX = sex, PARAMETER = parameter,
      ESTIMATE = observed_fit[observed_fit$SEX == sex, parameter],
      LOWER = q["LOWER"], MEDIAN = q["MEDIAN"], UPPER = q["UPPER"],
      N_VALID = q["N_VALID"]
    )
  }))
}))

curve_grid <- expand.grid(SEX = sexes, SIZE_BIN_5MM = seq(27.5, 177.5, by = 0.5),
                          stringsAsFactors = FALSE)
curve_replicates <- merge(
  curve_grid, bootstrap_parameters[bootstrap_parameters$CONVERGENCE == 0, ],
  by = "SEX", all.x = TRUE
)
curve_replicates$FITTED_SELECTIVITY <- with(
  curve_replicates,
  ASYMPTOTE * plogis((SIZE_BIN_5MM - L50_MM) / SCALE_MM)
)
curve_groups <- split(
  curve_replicates,
  interaction(curve_replicates$SEX, curve_replicates$SIZE_BIN_5MM, drop = TRUE)
)
curve_summary <- do.call(rbind, lapply(curve_groups, function(z) {
  q <- summarize_quantiles(z$FITTED_SELECTIVITY)
  data.frame(SEX = z$SEX[1], SIZE_BIN_5MM = z$SIZE_BIN_5MM[1],
             LOWER = q["LOWER"], MEDIAN = q["MEDIAN"],
             UPPER = q["UPPER"], N_VALID = q["N_VALID"])
}))
curve_summary <- merge(curve_summary, observed_fit, by = "SEX", all.x = TRUE)
curve_summary$ESTIMATE <- with(
  curve_summary,
  ASYMPTOTE * plogis((SIZE_BIN_5MM - L50_MM) / SCALE_MM)
)

# Quantify whether the recurring juvenile elevation exceeds the immediately
# following size range. Windows were specified from the diagnostic plots rather
# than selected independently in every bootstrap replicate.
hump_windows <- data.frame(
  SEX = sexes,
  HUMP_MIN = c(42.5, 42.5), HUMP_MAX = c(62.5, 67.5),
  POST_MIN = c(67.5, 72.5), POST_MAX = c(82.5, 102.5)
)
hump_draws <- do.call(rbind, lapply(seq_len(nrow(hump_windows)), function(i) {
  w <- hump_windows[i, ]
  z <- bootstrap_ratio[bootstrap_ratio$SEX == w$SEX, ]
  reps <- split(z, z$REPLICATE)
  do.call(rbind, lapply(reps, function(d) {
    hump <- mean(d$RATIO[d$SIZE_BIN_5MM >= w$HUMP_MIN &
                           d$SIZE_BIN_5MM <= w$HUMP_MAX], na.rm = TRUE)
    post <- mean(d$RATIO[d$SIZE_BIN_5MM >= w$POST_MIN &
                           d$SIZE_BIN_5MM <= w$POST_MAX], na.rm = TRUE)
    data.frame(REPLICATE = d$REPLICATE[1], SEX = w$SEX,
               HUMP_MEAN = hump, POST_HUMP_MEAN = post,
               DIFFERENCE = hump - post)
  }))
}))
hump_summary <- do.call(rbind, lapply(split(hump_draws, hump_draws$SEX), function(z) {
  q <- summarize_quantiles(z$DIFFERENCE)
  data.frame(
    SEX = z$SEX[1], DIFFERENCE_LOWER = q["LOWER"],
    DIFFERENCE_MEDIAN = q["MEDIAN"], DIFFERENCE_UPPER = q["UPPER"],
    PROBABILITY_POSITIVE = mean(z$DIFFERENCE > 0, na.rm = TRUE),
    N_VALID = sum(is.finite(z$DIFFERENCE))
  )
}))

bootstrap_ratio$KEY <- paste(bootstrap_ratio$SEX,
                             bootstrap_ratio$SIZE_BIN_5MM, sep = "_")
ratio_wide <- reshape(bootstrap_ratio[c("REPLICATE", "KEY", "RATIO")],
                      idvar = "REPLICATE", timevar = "KEY", direction = "wide")
ratio_covariance <- cov(ratio_wide[, setdiff(names(ratio_wide), "REPLICATE")],
                        use = "pairwise.complete.obs")

write.csv(bootstrap_parameters,
          file.path(output_dir, "bootstrap_logistic_parameter_draws.csv"),
          row.names = FALSE)
write.csv(bootstrap_ratio,
          file.path(output_dir, "bootstrap_empirical_selectivity_draws.csv"),
          row.names = FALSE)
write.csv(parameter_summary,
          file.path(output_dir, "bootstrap_logistic_parameter_summary.csv"),
          row.names = FALSE)
write.csv(ratio_summary,
          file.path(output_dir, "bootstrap_empirical_selectivity_at_size.csv"),
          row.names = FALSE)
write.csv(curve_summary,
          file.path(output_dir, "bootstrap_logistic_curve_at_size.csv"),
          row.names = FALSE)
write.csv(hump_draws,
          file.path(output_dir, "bootstrap_juvenile_hump_contrast_draws.csv"),
          row.names = FALSE)
write.csv(hump_summary,
          file.path(output_dir, "bootstrap_juvenile_hump_contrast_summary.csv"),
          row.names = FALSE)
write.csv(ratio_covariance,
          file.path(output_dir, "bootstrap_empirical_selectivity_covariance.csv"))
saveRDS(
  list(metadata = list(n_boot = n_boot, seed = bootstrap_seed,
                       resampling_unit = "paired station within species-year",
                       minimum_size_mm = minimum_size),
       empirical_summary = ratio_summary,
       empirical_draws = bootstrap_ratio,
       empirical_covariance = ratio_covariance,
       juvenile_hump_summary = hump_summary,
       juvenile_hump_draws = hump_draws,
       logistic_parameter_summary = parameter_summary,
       logistic_parameter_draws = bootstrap_parameters,
       logistic_curve_summary = curve_summary),
  file.path(output_dir, "bootstrap_selectivity_results.rds")
)

ratio_plot <- ratio_summary
ratio_plot$LOWER_PLOT <- pmax(0, ratio_plot$LOWER)
ratio_plot$UPPER_PLOT <- pmin(2, ratio_plot$UPPER)
ratio_plot$OBSERVED_PLOT <- pmin(2, ratio_plot$OBSERVED_RATIO)
ratio_plot_finite <- ratio_plot[is.finite(ratio_plot$OBSERVED_PLOT) &
                                  is.finite(ratio_plot$LOWER_PLOT) &
                                  is.finite(ratio_plot$UPPER_PLOT), ]

p_empirical <- ggplot(ratio_plot_finite, aes(SIZE_BIN_5MM, OBSERVED_PLOT)) +
  geom_hline(yintercept = 1, color = "grey40", linetype = "dashed") +
  geom_ribbon(aes(ymin = LOWER_PLOT, ymax = UPPER_PLOT),
              fill = "#56B4E9", alpha = 0.28) +
  geom_line(color = "#0072B2", linewidth = 0.9) +
  geom_point(color = "#0072B2", size = 1.2) +
  facet_wrap(~ SEX, ncol = 1) +
  scale_x_continuous(breaks = seq(20, 180, 20), limits = c(25, 180)) +
  scale_y_continuous(breaks = seq(0, 2, 0.25), limits = c(0, 2)) +
  labs(
    title = "Station-bootstrap uncertainty in empirical implied selectivity",
    subtitle = "Species pooled; 95% paired-station bootstrap intervals",
    x = "Carapace width (mm)", y = "NMFS / BSFRF",
    caption = "All species-years included; intervals truncated at 2 for display."
  ) + theme_bw(base_size = 10)
ggsave(file.path(output_dir, "bootstrap_empirical_selectivity_at_size.png"),
       p_empirical, width = 9.5, height = 8.5, dpi = 300, bg = "white")

p_logistic <- ggplot(curve_summary, aes(SIZE_BIN_5MM, ESTIMATE)) +
  geom_hline(yintercept = 1, color = "grey40", linetype = "dashed") +
  geom_ribbon(aes(ymin = LOWER, ymax = UPPER), fill = "#E69F00", alpha = 0.3) +
  geom_line(color = "#D55E00", linewidth = 1) +
  geom_point(data = ratio_plot_finite,
             aes(x = SIZE_BIN_5MM, y = OBSERVED_PLOT), inherit.aes = FALSE,
             color = "grey45", alpha = 0.65, size = 1.2) +
  facet_wrap(~ SEX, ncol = 1) +
  scale_x_continuous(breaks = seq(20, 180, 20), limits = c(25, 180)) +
  scale_y_continuous(breaks = seq(0, 2, 0.25), limits = c(0, 2)) +
  labs(
    title = "Bootstrap uncertainty in fitted logistic selectivity",
    subtitle = "Three-parameter logistic over the joint-model size domain",
    x = "Carapace width (mm)", y = "NMFS / BSFRF",
    caption = "Ribbon is the 95% paired-station bootstrap interval."
  ) + theme_bw(base_size = 10)
ggsave(file.path(output_dir, "bootstrap_logistic_selectivity_at_size.png"),
       p_logistic, width = 9.5, height = 8.5, dpi = 300, bg = "white")

cat("Completed", n_boot, "paired-station bootstrap replicates.\n")
print(parameter_summary)
