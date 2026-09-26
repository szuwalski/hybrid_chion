library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
argval <- function(key, default) {
  z <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(z)) sub(paste0("^", key, "="), "", z[1]) else default
}

fit <- readRDS(argval("--fit", "models/model/fits/model_fit.rds"))
d <- fit$data
r <- fit$report
out <- argval("--out", "models/model/plots/model_diagnostics")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
year_run <- function(z, cols) {
  key <- interaction(z[cols], drop = TRUE, lex.order = TRUE)
  run <- integer(nrow(z))
  for (ii in split(seq_len(nrow(z)), key)) {
    jj <- ii[order(z$year[ii])]
    run[jj] <- cumsum(c(TRUE, diff(z$year[jj]) != 1))
  }
  interaction(key, run, drop = TRUE)
}

source_info <- data.frame(
  species = c("snow", "snow", "tanner", "tanner"),
  fishery = c("snow", "tanner", "tanner", "snow"),
  block = c("snow_male", "snow_male", "tanner_male", "tanner_male"),
  prefix = c("sm", "sm", "tm", "tm"),
  role = c("directed", "incidental", "directed", "incidental"),
  stringsAsFactors = FALSE
)

prediction <- function(prefix, role, component) {
  stem <- if (component == "total") "caught" else
    if (component == "retained") "ret" else "discard"
  r[[paste(prefix, stem, role, sep = "_")]]
}

catch_rows <- list()
size_rows <- list()
selectivity_rows <- list()
f_rows <- list()
l50_rows <- list()

for (j in seq_len(nrow(source_info))) {
  q <- source_info[j, ]
  b <- d[[q$block]]
  fd <- b$data$fisheries[[q$fishery]]
  years <- b$data$years
  sizes <- b$data$sizes

  for (component in c("total", "retained")) {
    pred <- prediction(q$prefix, q$role, component)
    obs_years <- fd[[paste0(component, "_years")]]
    obs <- fd[[paste0(component, "_numbers")]]
    iy <- match(obs_years, years)
    catch_rows[[length(catch_rows) + 1L]] <- data.frame(
      year = obs_years, observed = obs,
      predicted = rowSums(pred[iy, , drop = FALSE]),
      species = q$species, fishery = q$fishery,
      catch_type = ifelse(q$role == "directed", "directed", "incidental"),
      component = component
    )

    obs_size <- fd[[paste0(component, "_size")]]
    obs_size_years <- fd[[paste0(component, "_size_years")]]
    if (length(obs_size_years)) {
      iy <- match(obs_size_years, years)
      pp <- pred[iy, , drop = FALSE]
      pp <- pp / rowSums(pp)
      size_rows[[length(size_rows) + 1L]] <- data.frame(
        year = rep(obs_size_years, each = length(sizes)),
        size = rep(sizes, times = length(obs_size_years)),
        observed = as.vector(t(obs_size)), predicted = as.vector(t(pp)),
        species = q$species, fishery = q$fishery,
        catch_type = ifelse(q$role == "directed", "directed", "incidental"),
        component = component
      )
    }
  }

  # Historical snow discard compositions retained in the directed fishery.
  if (length(fd$discard_size_years)) {
    pred <- prediction(q$prefix, q$role, "discard")
    iy <- match(fd$discard_size_years, years)
    pp <- pred[iy, , drop = FALSE]
    pp <- pp / rowSums(pp)
    size_rows[[length(size_rows) + 1L]] <- data.frame(
      year = rep(fd$discard_size_years, each = length(sizes)),
      size = rep(sizes, times = length(fd$discard_size_years)),
      observed = as.vector(t(fd$discard_size)), predicted = as.vector(t(pp)),
      species = q$species, fishery = q$fishery,
      catch_type = ifelse(q$role == "directed", "directed", "incidental"),
      component = "discard_historical"
    )
  }

  for (component in c("total", "retention")) {
    stem <- if (component == "total") "fish_total" else "fish_retain"
    mat <- r[[paste(q$prefix, stem, q$role, sep = "_")]]
    active <- years %in% fd$active_years
    z <- mat[active, , drop = FALSE]
    selectivity_rows[[length(selectivity_rows) + 1L]] <- data.frame(
      year = rep(years[active], each = length(sizes)),
      size = rep(sizes, times = sum(active)), value = as.vector(t(z)),
      species = q$species, fishery = q$fishery,
      catch_type = ifelse(q$role == "directed", "directed", "incidental"),
      component = component
    )
  }
  F_name <- paste(q$prefix, "F", q$role, sep = "_")
  F_value <- r[[F_name]]
  f_rows[[length(f_rows) + 1L]] <- data.frame(
    year = years, F = F_value, species = q$species, fishery = q$fishery,
    catch_type = ifelse(q$role == "directed", "directed", "incidental")
  )
  role_name <- if (q$role == "directed") "directed" else "incidental"
  total_l50 <- r[[paste0(q$prefix, "_total_l50_", role_name)]]
  retained_l50 <- r[[paste0(q$prefix, "_ret_l50_", role_name)]]
  l50_rows[[length(l50_rows) + 1L]] <- rbind(
    data.frame(year = years, L50 = total_l50, component = "total"),
    data.frame(year = years, L50 = retained_l50, component = "retained")
  ) |>
    transform(species = q$species, fishery = q$fishery,
              catch_type = role_name)
}

catch_fits <- do.call(rbind, catch_rows)
size_fits <- do.call(rbind, size_rows)
selectivity <- do.call(rbind, selectivity_rows)
fishing_mortality <- do.call(rbind, f_rows)
midpoints <- do.call(rbind, l50_rows)
catch_fits$residual <- log(catch_fits$observed) - log(catch_fits$predicted)
write.csv(catch_fits, file.path(out, "catch_fits.csv"), row.names = FALSE)
write.csv(size_fits, file.path(out, "fishery_size_composition_fits.csv"), row.names = FALSE)
write.csv(selectivity, file.path(out, "annual_fishery_selectivity.csv"), row.names = FALSE)
write.csv(fishing_mortality, file.path(out, "annual_fishing_mortality.csv"), row.names = FALSE)
write.csv(midpoints, file.path(out, "annual_selectivity_midpoints.csv"), row.names = FALSE)

labeller_text <- function(x)
  paste(tools::toTitleCase(x$species), "in", tools::toTitleCase(x$fishery), "fishery")
catch_fits$source <- labeller_text(catch_fits)
size_fits$source <- labeller_text(size_fits)
selectivity$source <- labeller_text(selectivity)
fishing_mortality$source <- labeller_text(fishing_mortality)
midpoints$source <- labeller_text(midpoints)
catch_fits$line_group <- year_run(catch_fits, c("source", "component"))

p <- ggplot(catch_fits, aes(year, group = line_group)) +
  geom_line(aes(y = observed, colour = "Observed"), linewidth = .65) +
  geom_point(aes(y = observed, colour = "Observed"), size = 1.2) +
  geom_line(aes(y = predicted, colour = "Predicted"), linewidth = .75) +
  facet_grid(component ~ source, scales = "free_y") +
  scale_y_log10(labels = scales::label_number()) +
  scale_colour_manual(values = c(Observed = "black", Predicted = "#D55E00")) +
  labs(x = NULL, y = "Catch (numbers; log scale)", colour = NULL,
       title = "Two-fishery catch fits") +
  theme_bw(base_size = 10) + theme(legend.position = "top")
ggsave(file.path(out, "catch_fits.png"), p, width = 14, height = 7, dpi = 180)

for (key in unique(interaction(size_fits$source, size_fits$component, drop = TRUE))) {
  z <- size_fits[interaction(size_fits$source, size_fits$component, drop = TRUE) == key, ]
  nr <- ceiling(length(unique(z$year)) / 5)
  p <- ggplot(z, aes(size)) +
    geom_col(aes(y = observed), fill = "grey65", colour = "grey35",
             width = 4.5) +
    geom_line(aes(y = predicted), colour = "#D55E00", linewidth = .65) +
    facet_wrap(~year, ncol = 5) +
    labs(x = "Carapace width (mm)", y = "Proportion",
         title = paste(unique(z$source), "-", unique(z$component))) +
    theme_bw(base_size = 9)
  filename <- gsub("[^a-z0-9]+", "_", tolower(paste(unique(z$source), unique(z$component))))
  ggsave(file.path(out, paste0("fishery_size_composition_", filename, ".png")),
         p, width = 12, height = max(3, 1.65 * nr), dpi = 170)
}

p <- ggplot(selectivity, aes(size, value, group = year, colour = year)) +
  geom_line(alpha = .45, linewidth = .45) +
  facet_grid(component ~ source, scales = "free_x") +
  scale_colour_viridis_c() +
  labs(x = "Carapace width (mm)", y = "Selectivity", colour = "Year",
       title = "Annual fishery selectivity curves") +
  theme_bw(base_size = 10)
ggsave(file.path(out, "annual_fishery_selectivity.png"), p,
       width = 14, height = 7, dpi = 180)

fm_plot <- subset(fishing_mortality, F > 0)
fm_plot$line_group <- year_run(fm_plot, c("source"))
p <- ggplot(fm_plot, aes(year, F, group = line_group)) +
  geom_line(colour = "#0072B2", linewidth = .75) +
  geom_point(colour = "#0072B2", size = 1) +
  facet_wrap(~source, ncol = 2, scales = "free_y") +
  labs(x = NULL, y = "Fully selected fishing mortality",
       title = "Annual fishing mortality by catch source") +
  theme_bw(base_size = 10)
ggsave(file.path(out, "annual_fishing_mortality.png"), p,
       width = 11, height = 7, dpi = 180)

midpoints$line_group <- year_run(midpoints, c("source", "component"))
p <- ggplot(midpoints, aes(year, L50, colour = component, group = line_group)) +
  geom_line(linewidth = .7) +
  facet_wrap(~source, ncol = 2) +
  scale_colour_manual(values = c(total = "#0072B2", retained = "#D55E00")) +
  labs(x = NULL, y = "Selectivity midpoint (mm)", colour = NULL,
       title = "Annual total and retained selectivity midpoints") +
  theme_bw(base_size = 10) + theme(legend.position = "top")
ggsave(file.path(out, "annual_selectivity_midpoints.png"), p,
       width = 11, height = 7, dpi = 180)

cat("Wrote two-fishery diagnostics to", out, "\n")
