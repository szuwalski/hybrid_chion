#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_file <- file.path(root, "data", "bsfrf", "snow",
                        "bsfrf_nmfs_area_numbers_at_size.rds")
output_dir <- file.path(root, "analyses", "bsfrf_selectivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Run data/bsfrf/prepare_bsfrf_nmfs_area_data.R first")
x <- readRDS(input_file)

gear_labels <- c(
  BSFRF_NEPHROPS = "BSFRF Nephrops",
  NMFS_STANDARD_83_112 = "NMFS standard 83/112"
)
gear_colors <- c(
  "BSFRF Nephrops" = "#0072B2",
  "NMFS standard 83/112" = "#D55E00"
)
gear_types <- c(
  "BSFRF Nephrops" = "solid",
  "NMFS standard 83/112" = "22"
)

theme_crab <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      strip.background = element_rect(fill = "grey92", color = "grey55"),
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title.position = "plot",
      plot.caption = element_text(hjust = 0, color = "grey30"),
      axis.text.x = element_text(size = 8)
    )
}

numbers <- x$area_year_numbers_at_size
size_bins <- sort(unique(numbers$SIZE_BIN_5MM))
numbers$GEAR_LABEL <- unname(gear_labels[numbers$GEAR])
numbers$SEX <- factor(numbers$SEX, levels = c("MALE", "FEMALE"),
                      labels = c("Male", "Female"))

plot_numbers <- function(species) {
  d <- numbers[numbers$SPECIES == species, ]
  years <- sort(unique(d$YEAR))
  d$YEAR <- factor(d$YEAR, levels = years)
  d$FACET_LABEL <- factor(paste(d$YEAR, d$SEX, sep = " | "),
                          levels = c(paste(years, "Male", sep = " | "),
                                     paste(years, "Female", sep = " | ")))
  ggplot(d, aes(x = SIZE_BIN_5MM, y = ESTIMATED_NUMBERS,
                color = GEAR_LABEL, linetype = GEAR_LABEL,
                group = GEAR_LABEL)) +
    geom_line(linewidth = 0.75, na.rm = TRUE) +
    facet_wrap(~ FACET_LABEL, nrow = 2, scales = "free_y") +
    scale_color_manual(values = gear_colors) +
    scale_linetype_manual(values = gear_types) +
    scale_x_continuous(breaks = seq(0, 180, 20), limits = c(0, 180),
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()),
                       expand = expansion(mult = c(0, 0.06))) +
    labs(
      title = paste0(tools::toTitleCase(tolower(species)),
                     " crab numbers at size in common survey cells"),
      subtitle = "Each panel uses the NMFS station cells sampled by both gears in that year",
      x = "Carapace width (mm; 5-mm bin midpoint)",
      y = "Estimated number of crab",
      caption = paste(
        "BSFRF densities are averaged within station before expansion.",
        "Station-area assumptions are recorded in the source RDS."
      )
    ) +
    theme_crab(9)
}

comparison <- x$area_year_comparison
comparison$SEX <- factor(comparison$SEX, levels = c("MALE", "FEMALE"),
                         labels = c("Male", "Female"))
comparison$RATIO_STATUS <- ifelse(
  comparison$BSFRF_NEPHROPS == 0 & comparison$NMFS_STANDARD_83_112 > 0,
  "BSFRF zero; ratio infinite",
  ifelse(comparison$BSFRF_NEPHROPS == 0 & comparison$NMFS_STANDARD_83_112 == 0,
         "Both gears zero", "Finite")
)
comparison$RATIO_PLOT <- comparison$NMFS_TO_BSFRF_RATIO
comparison$RATIO_PLOT[comparison$RATIO_STATUS == "BSFRF zero; ratio infinite"] <- 4
comparison$RATIO_PLOT[comparison$RATIO_STATUS == "Both gears zero"] <- NA_real_

plot_selectivity <- function(species) {
  d <- comparison[comparison$SPECIES == species, ]
  years <- sort(unique(d$YEAR))
  d$YEAR <- factor(d$YEAR, levels = years)
  finite <- d[d$RATIO_STATUS == "Finite", ]
  infinite <- d[d$RATIO_STATUS == "BSFRF zero; ratio infinite", ]

  ggplot(finite, aes(x = SIZE_BIN_5MM, y = RATIO_PLOT, group = 1)) +
    geom_hline(yintercept = 1, color = "grey35", linetype = "dashed", linewidth = 0.45) +
    geom_line(color = "#009E73", linewidth = 0.75, na.rm = TRUE) +
    geom_point(data = infinite, shape = 24, size = 2.2, fill = "#CC79A7",
               color = "black", stroke = 0.35) +
    facet_grid(SEX ~ YEAR) +
    scale_x_continuous(breaks = seq(0, 180, 20), limits = c(0, 180),
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(breaks = seq(0, 4, 0.5), limits = c(0, 4),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = paste0("Implied NMFS selectivity at size: ",
                     tools::toTitleCase(tolower(species)), " crab"),
      subtitle = "NMFS estimated numbers divided by BSFRF estimated numbers in the same station cells",
      x = "Carapace width (mm; 5-mm bin midpoint)",
      y = "NMFS / BSFRF",
      caption = paste(
        "Dashed line is equal catchability (ratio = 1).",
        "Triangles at 4 mark NMFS-positive bins with zero BSFRF catch (infinite ratio);",
        "bins where both gears caught zero are omitted."
      )
    ) +
    theme_crab(9) + theme(legend.position = "none")
}

for (species in c("SNOW", "TANNER")) {
  width <- if (species == "SNOW") 12 else 15
  ggsave(file.path(output_dir, paste0(tolower(species), "_numbers_at_size_by_gear.png")),
         plot_numbers(species), width = width, height = 6.2, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0(tolower(species), "_implied_selectivity_at_size.png")),
         plot_selectivity(species), width = width, height = 6.2, dpi = 300, bg = "white")
}

# Unweighted biological measurements. These are specimen records actually
# measured, not sampling-factor expansions and not area-swept abundance.
count_observed <- function(bio, haul, gear) {
  eligible <- haul$HAUL_ID[haul$COMMON_STATION]
  d <- bio[bio$HAUL_ID %in% eligible & bio$SEX %in% c("MALE", "FEMALE") &
             is.finite(bio$SIZE_BIN_5MM) & bio$SIZE_BIN_5MM %in% size_bins, ]
  d$GEAR <- gear
  aggregate(HAUL_ID ~ YEAR + SPECIES + GEAR + SEX + SIZE_BIN_5MM,
            d, length)
}

observed <- rbind(
  count_observed(x$bsfrf_bio, x$bsfrf_haul, "BSFRF_NEPHROPS"),
  count_observed(x$nmfs_bio, x$nmfs_haul, "NMFS_STANDARD_83_112")
)
names(observed)[6] <- "N_MEASURED"
observed$GEAR_LABEL <- unname(gear_labels[observed$GEAR])
observed$SEX <- factor(observed$SEX, levels = c("MALE", "FEMALE"),
                       labels = c("Male", "Female"))

plot_observed <- function(species) {
  d <- observed[observed$SPECIES == species, ]
  years <- sort(unique(d$YEAR))
  d$FACET_LABEL <- factor(paste(d$YEAR, d$SEX, sep = " | "),
                          levels = c(paste(years, "Male", sep = " | "),
                                     paste(years, "Female", sep = " | ")))
  ggplot(d, aes(x = SIZE_BIN_5MM, y = N_MEASURED,
                color = GEAR_LABEL, linetype = GEAR_LABEL,
                group = GEAR_LABEL)) +
    geom_line(linewidth = 0.75, na.rm = TRUE) +
    facet_wrap(~ FACET_LABEL, nrow = 2, scales = "free_y") +
    scale_color_manual(values = gear_colors) +
    scale_linetype_manual(values = gear_types) +
    scale_x_continuous(breaks = seq(0, 180, 20), limits = c(0, 180),
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = label_comma(),
                       expand = expansion(mult = c(0, 0.06))) +
    labs(
      title = paste0("Measured ", tools::toTitleCase(tolower(species)),
                     " crab specimens at size"),
      subtitle = "Unweighted biological records from hauls in the common station footprint",
      x = "Carapace width (mm; 5-mm bin midpoint)",
      y = "Number of measured crab",
      caption = paste(
        "Counts are not multiplied by sampling factors or divided by area swept.",
        "They represent direct size measurements, not total catch when subsampling occurred."
      )
    ) +
    theme_crab(9)
}

for (species in c("SNOW", "TANNER")) {
  width <- if (species == "SNOW") 12 else 15
  ggsave(file.path(output_dir, paste0(tolower(species),
                                      "_measured_specimens_at_size_by_gear.png")),
         plot_observed(species), width = width, height = 6.2,
         dpi = 300, bg = "white")
}

sample_totals <- aggregate(N_MEASURED ~ YEAR + SPECIES + GEAR_LABEL + SEX,
                           observed, sum)
sample_totals$SPECIES <- tools::toTitleCase(tolower(sample_totals$SPECIES))
p_sample_totals <- ggplot(sample_totals,
                          aes(x = factor(YEAR), y = N_MEASURED,
                              fill = GEAR_LABEL)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  geom_text(aes(label = comma(N_MEASURED)),
            position = position_dodge(width = 0.78), vjust = -0.25,
            size = 2.6) +
  facet_grid(SEX ~ SPECIES, scales = "free_y") +
  scale_fill_manual(values = gear_colors) +
  scale_y_continuous(labels = label_comma(),
                     expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Measured sample sizes supporting numbers-at-size comparisons",
    subtitle = "Totals across all sizes in common station-year footprints",
    x = "Year", y = "Number of measured crab"
  ) +
  theme_crab(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "measured_sample_sizes_by_year.png"),
       p_sample_totals, width = 11, height = 7.5, dpi = 300, bg = "white")

# Robust cross-year summary: median finite ratio and interquartile range at size.
finite_ratio <- comparison[comparison$RATIO_STATUS == "Finite", ]
pooled_split <- split(finite_ratio,
                      interaction(finite_ratio$SPECIES, finite_ratio$SEX,
                                  finite_ratio$SIZE_BIN_5MM, drop = TRUE))
pooled <- do.call(rbind, lapply(pooled_split, function(z) {
  data.frame(
    SPECIES = z$SPECIES[1], SEX = z$SEX[1], SIZE_BIN_5MM = z$SIZE_BIN_5MM[1],
    MEDIAN_RATIO = median(z$NMFS_TO_BSFRF_RATIO, na.rm = TRUE),
    Q25 = unname(quantile(z$NMFS_TO_BSFRF_RATIO, 0.25, na.rm = TRUE)),
    Q75 = unname(quantile(z$NMFS_TO_BSFRF_RATIO, 0.75, na.rm = TRUE)),
    N_YEARS = length(unique(z$YEAR))
  )
}))
row.names(pooled) <- NULL

p_pooled <- ggplot(pooled, aes(x = SIZE_BIN_5MM, y = MEDIAN_RATIO)) +
  geom_hline(yintercept = 1, color = "grey35", linetype = "dashed", linewidth = 0.45) +
  geom_ribbon(aes(ymin = Q25, ymax = Q75), fill = "#56B4E9", alpha = 0.28) +
  geom_line(color = "#0072B2", linewidth = 0.9) +
  facet_grid(SEX ~ SPECIES, scales = "free") +
  scale_x_continuous(breaks = seq(0, 180, 20), limits = c(0, 180)) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Cross-year summary of implied NMFS selectivity at size",
    subtitle = "Median and interquartile range of finite year-specific NMFS/BSFRF ratios",
    x = "Carapace width (mm; 5-mm bin midpoint)", y = "NMFS / BSFRF",
    caption = "Only years with positive BSFRF estimated numbers in a size bin contribute to its summary."
  ) +
  theme_crab(10) + theme(legend.position = "none")
ggsave(file.path(output_dir, "pooled_implied_selectivity_at_size.png"), p_pooled,
       width = 10, height = 7, dpi = 300, bg = "white")

# Species-combined selectivity by sex using every available species-year. The
# abundance estimates are pooled within gear and size before taking the ratio.
# A three-parameter logistic estimates its asymptote (relative catchability),
# midpoint, and slope scale while remaining monotonic.
combined_year <- aggregate(
  cbind(BSFRF_NEPHROPS, NMFS_STANDARD_83_112) ~ YEAR + SEX + SIZE_BIN_5MM,
  comparison, sum, na.rm = TRUE
)
combined_year$NMFS_TO_BSFRF_RATIO <- with(
  combined_year,
  ifelse(BSFRF_NEPHROPS > 0, NMFS_STANDARD_83_112 / BSFRF_NEPHROPS, NA_real_)
)
combined_pooled <- aggregate(
  cbind(BSFRF_NEPHROPS, NMFS_STANDARD_83_112) ~ SEX + SIZE_BIN_5MM,
  comparison, sum, na.rm = TRUE
)
combined_pooled$NMFS_TO_BSFRF_RATIO <- with(
  combined_pooled,
  ifelse(BSFRF_NEPHROPS > 0, NMFS_STANDARD_83_112 / BSFRF_NEPHROPS, NA_real_)
)
combined_pooled$RATIO_PLOT <- pmin(combined_pooled$NMFS_TO_BSFRF_RATIO, 2)
combined_pooled$ABOVE_PLOT_LIMIT <- combined_pooled$NMFS_TO_BSFRF_RATIO > 2
measured_support <- aggregate(N_MEASURED ~ SEX + SIZE_BIN_5MM, observed, sum)
names(measured_support)[3] <- "N_MEASURED_SUPPORT"
combined_pooled <- merge(combined_pooled, measured_support,
                         by = c("SEX", "SIZE_BIN_5MM"), all.x = TRUE)
combined_pooled$N_MEASURED_SUPPORT[is.na(combined_pooled$N_MEASURED_SUPPORT)] <- 1

fit_logistic <- function(d) {
  d <- d[is.finite(d$NMFS_TO_BSFRF_RATIO), ]
  objective <- function(par) {
    predicted <- exp(par[3]) *
      plogis((d$SIZE_BIN_5MM - par[1]) / exp(par[2]))
    weighted.mean((d$NMFS_TO_BSFRF_RATIO - predicted)^2,
                  w = d$N_MEASURED_SUPPORT)
  }
  fit <- optim(c(80, log(15), log(0.6)), objective, method = "L-BFGS-B",
               lower = c(0, log(0.5), log(0.01)),
               upper = c(250, log(150), log(5)))
  grid <- data.frame(
    SEX = d$SEX[1],
    SIZE_BIN_5MM = seq(min(d$SIZE_BIN_5MM), max(d$SIZE_BIN_5MM), by = 0.5)
  )
  grid$FITTED_SELECTIVITY <- exp(fit$par[3]) * plogis(
    (grid$SIZE_BIN_5MM - fit$par[1]) / exp(fit$par[2])
  )
  list(
    curve = grid,
    coefficients = data.frame(
      SEX = d$SEX[1], ASYMPTOTE = exp(fit$par[3]),
      L50_MM = fit$par[1], SCALE_MM = exp(fit$par[2]),
      L95_MM = fit$par[1] + qlogis(0.95) * exp(fit$par[2]),
      WEIGHTED_MSE = fit$value, N_SIZE_BINS = nrow(d),
      N_MEASURED_SUPPORT = sum(d$N_MEASURED_SUPPORT),
      CONVERGENCE = fit$convergence
    )
  )
}

logistic_fits <- lapply(split(combined_pooled, combined_pooled$SEX), fit_logistic)
logistic_curves <- do.call(rbind, lapply(logistic_fits, `[[`, "curve"))
logistic_coefficients <- do.call(rbind, lapply(logistic_fits, `[[`, "coefficients"))
row.names(logistic_curves) <- NULL
row.names(logistic_coefficients) <- NULL

# Sensitivity fit over the size domain represented in the joint model.
joint_model_min_size <- 27.5
joint_domain_pooled <- combined_pooled[
  combined_pooled$SIZE_BIN_5MM >= joint_model_min_size, ]
joint_domain_fits <- lapply(split(joint_domain_pooled, joint_domain_pooled$SEX),
                            fit_logistic)
joint_domain_curves <- do.call(rbind, lapply(joint_domain_fits, `[[`, "curve"))
joint_domain_coefficients <- do.call(
  rbind, lapply(joint_domain_fits, `[[`, "coefficients"))
row.names(joint_domain_curves) <- NULL
row.names(joint_domain_coefficients) <- NULL

coefficient_comparison <- rbind(
  transform(logistic_coefficients, FIT_DOMAIN = "All observed sizes"),
  transform(joint_domain_coefficients,
            FIT_DOMAIN = paste0("Joint-model domain (>= ", joint_model_min_size, " mm)"))
)

p_species_combined <- ggplot(
  combined_pooled,
  aes(x = SIZE_BIN_5MM, y = RATIO_PLOT)
) +
  geom_hline(yintercept = 1, color = "grey35", linetype = "dashed", linewidth = 0.45) +
  geom_point(data = combined_pooled[!is.na(combined_pooled$ABOVE_PLOT_LIMIT) &
                                      !combined_pooled$ABOVE_PLOT_LIMIT, ],
             color = "grey45", alpha = 0.7, size = 1.6) +
  geom_point(data = combined_pooled[!is.na(combined_pooled$ABOVE_PLOT_LIMIT) &
                                      combined_pooled$ABOVE_PLOT_LIMIT, ],
             shape = 24, fill = "grey65", color = "grey25", size = 2) +
  geom_line(data = logistic_curves, aes(y = FITTED_SELECTIVITY),
            color = "#D55E00", linewidth = 1.1) +
  facet_wrap(~ SEX, ncol = 1, scales = "free_x") +
  scale_x_continuous(breaks = seq(0, 180, 20), limits = c(0, 180)) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Logistic fit to species-combined implied NMFS selectivity",
    subtitle = "Snow + Tanner pooled across every available year; asymptote, midpoint, and slope scale estimated",
    x = "Carapace width (mm; 5-mm bin midpoint)", y = "NMFS / BSFRF",
    caption = "All species-years are included; fits are weighted by measured specimens at size. Triangles mark ratios above 2."
  ) +
  theme_crab(10) + theme(legend.position = "none")

ggsave(file.path(output_dir, "species_combined_implied_selectivity_by_sex.png"),
       p_species_combined, width = 9.5, height = 8.5, dpi = 300, bg = "white")

p_joint_domain <- ggplot(
  joint_domain_pooled,
  aes(x = SIZE_BIN_5MM, y = RATIO_PLOT)
) +
  geom_hline(yintercept = 1, color = "grey35", linetype = "dashed", linewidth = 0.45) +
  geom_point(data = joint_domain_pooled[
    !is.na(joint_domain_pooled$ABOVE_PLOT_LIMIT) &
      !joint_domain_pooled$ABOVE_PLOT_LIMIT, ],
    color = "grey45", alpha = 0.7, size = 1.6) +
  geom_point(data = joint_domain_pooled[
    !is.na(joint_domain_pooled$ABOVE_PLOT_LIMIT) &
      joint_domain_pooled$ABOVE_PLOT_LIMIT, ],
    shape = 24, fill = "grey65", color = "grey25", size = 2) +
  geom_line(data = logistic_curves, aes(y = FITTED_SELECTIVITY),
            color = "grey50", linetype = "dotted", linewidth = 0.9) +
  geom_line(data = joint_domain_curves, aes(y = FITTED_SELECTIVITY),
            color = "#0072B2", linewidth = 1.15) +
  facet_wrap(~ SEX, ncol = 1, scales = "free_x") +
  scale_x_continuous(breaks = seq(20, 180, 20), limits = c(25, 180)) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Logistic sensitivity fit over the joint-model size domain",
    subtitle = "Blue: fit to bins >= 27.5 mm; dotted grey: previous all-size fit",
    x = "Carapace width (mm; 5-mm bin midpoint)", y = "NMFS / BSFRF",
    caption = "All species-years are retained; fits are weighted by measured specimens at size."
  ) +
  theme_crab(10) + theme(legend.position = "none")
ggsave(file.path(output_dir, "species_combined_logistic_joint_model_sizes.png"),
       p_joint_domain, width = 9.5, height = 8.5, dpi = 300, bg = "white")
write.csv(combined_year,
          file.path(output_dir, "species_combined_year_selectivity_by_sex.csv"),
          row.names = FALSE)
write.csv(combined_pooled,
          file.path(output_dir, "species_combined_pooled_selectivity_by_sex.csv"),
          row.names = FALSE)
write.csv(logistic_curves,
          file.path(output_dir, "species_combined_logistic_curve_by_sex.csv"),
          row.names = FALSE)
write.csv(logistic_coefficients,
          file.path(output_dir, "species_combined_logistic_coefficients_by_sex.csv"),
          row.names = FALSE)
write.csv(joint_domain_curves,
          file.path(output_dir, "species_combined_logistic_joint_model_sizes_curve.csv"),
          row.names = FALSE)
write.csv(joint_domain_coefficients,
          file.path(output_dir, "species_combined_logistic_joint_model_sizes_coefficients.csv"),
          row.names = FALSE)
write.csv(coefficient_comparison,
          file.path(output_dir, "species_combined_logistic_domain_comparison.csv"),
          row.names = FALSE)

# Gear-agreement diagnostic. A pseudo-log scale retains real zeros.
p_scatter <- ggplot(comparison,
                    aes(x = BSFRF_NEPHROPS, y = NMFS_STANDARD_83_112,
                        color = factor(YEAR))) +
  geom_abline(slope = 1, intercept = 0, color = "grey25", linetype = "dashed") +
  geom_point(alpha = 0.75, size = 1.35) +
  facet_grid(SEX ~ SPECIES) +
  scale_x_continuous(trans = pseudo_log_trans(base = 10),
                     labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10),
                     labels = label_number(scale_cut = cut_short_scale())) +
  scale_color_viridis_d(option = "D", end = 0.9) +
  labs(
    title = "NMFS versus BSFRF mean density at size",
    subtitle = "Each point is a year-sex-species-size bin; dashed line is 1:1",
    x = "BSFRF crab per square nautical mile",
    y = "NMFS crab per square nautical mile", color = "Year"
  ) +
  theme_crab(10) + theme(legend.position = "right")
ggsave(file.path(output_dir, "nmfs_vs_bsfrf_numbers_at_size.png"), p_scatter,
       width = 10, height = 8, dpi = 300, bg = "white")

# Sampling support by species-year.
coverage <- x$coverage[x$coverage$GEAR == "BSFRF_NEPHROPS", ]
coverage$SPECIES <- tools::toTitleCase(tolower(coverage$SPECIES))
p_coverage <- ggplot(coverage, aes(x = YEAR, y = N_COMMON_STATIONS,
                                  color = SPECIES, group = SPECIES)) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c("Snow" = "#0072B2", "Tanner" = "#D55E00")) +
  scale_x_continuous(breaks = sort(unique(coverage$YEAR))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Common station support for gear comparisons",
    x = "Year", y = "NMFS station cells sampled by both gears", color = "Species",
    caption = "Missing species-years indicate that the BSFRF source did not target that species."
  ) +
  theme_crab(10)
ggsave(file.path(output_dir, "common_station_coverage.png"), p_coverage,
       width = 8.5, height = 5.5, dpi = 300, bg = "white")

cat("Saved plots to", output_dir, "\n")
print(file.info(list.files(output_dir, pattern = "\\.png$", full.names = TRUE))[, "size", drop = FALSE])
