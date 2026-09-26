#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
source("models/model/R/joint_classification_model.R")

args<-commandArgs(trailingOnly=TRUE)
argval<-function(key,default){z<-grep(paste0("^",key,"="),args,value=TRUE);if(length(z))sub(paste0("^",key,"="),"",z[1]) else default}
fit_file <- argval("--fit", "models/model/fits/model_fit.rds")
out_dir <- argval("--out", "models/model/plots/model_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
f <- readRDS(fit_file)
p <- f$par_list
prior <- f$data$bsfrf_survey_prior

era_labels <- c("1975-1981", "1982-1988", "1989-present")
make_curves <- function(group) {
  pr <- prior[[group]]
  parts <- strsplit(group, "_", fixed = TRUE)[[1]]
  species <- tools::toTitleCase(parts[1])
  sex <- tools::toTitleCase(parts[2])
  prefix <- paste0("survey_", parts[1], "_", substr(parts[2], 1, 1))
  offsets<-function(x)if(parts[1]=="snow")c(x[2],x[2],0) else c(x,0)
  if(parts[2]=="male"){
    sp<-paste0(prefix,"_stack")
    fitted<-stacked_logistic_survey_curves(pr$sizes,
      p[[paste0(sp,"_midpoint1")]],p[[paste0(sp,"_log_gap")]],
      p[[paste0(sp,"_log_slope1")]],p[[paste0(sp,"_log_slope2")]],
      p[[paste0(sp,"_logit_weight")]],p[[paste0(sp,"_log_q")]],
      p[[paste0(sp,"_midpoint1_offset")]],p[[paste0(sp,"_log_gap_offset")]],
      p[[paste0(sp,"_log_slope1_offset")]],p[[paste0(sp,"_log_slope2_offset")]],
      p[[paste0(sp,"_logit_weight_offset")]],p[[paste0(sp,"_log_q_offset")]],
      parts[1]=="snow")
  }else{
  midpoint_raw <- p[[paste0(prefix, "_midpoint")]] +
    offsets(p[[paste0(prefix, "_midpoint_offset")]])
  lower_mid<-(min(pr$sizes)+5)/30;upper_mid<-(max(pr$sizes)-5)/30
  midpoint<-lower_mid+(upper_mid-lower_mid)*plogis(midpoint_raw)
  log_slope <- p[[paste0(prefix, "_log_slope")]] +
    offsets(p[[paste0(prefix, "_log_slope_offset")]])
  log_q <- p[[paste0(prefix,"_log_q")]]+
    offsets(p[[paste0(prefix,"_log_q_offset")]])
  fitted <- sapply(seq_len(3), function(era){
    shape<-plogis(exp(log_slope[era]) * (pr$sizes / 30 - midpoint[era]))
    exp(log_q[era])*shape/tail(shape,1)
  })
  }
  if(!is.null(f$data$survey_q_fixed_at_terminal_size)&&
     isTRUE(f$data$survey_q_fixed_at_terminal_size))
    fitted<-sweep(fitted,2,fitted[pr$anchor_index,],"/")
  do.call(rbind, lapply(seq_len(3), function(era) data.frame(
    SPECIES = species, SEX = sex, ERA = era_labels[era],
    SIZE_MM = pr$sizes, FITTED = fitted[, era],
    BSFRF_START_TARGET = pr$bootstrap_center_ratio
  )))
}

groups <- c("snow_male", "tanner_male", "snow_female", "tanner_female")
curves <- do.call(rbind, lapply(groups, make_curves))
hybrid_curves <- do.call(rbind,lapply(c("Male","Female"),function(sex){
  snow<-subset(curves,SPECIES=="Snow"&SEX==sex)
  tanner<-subset(curves,SPECIES=="Tanner"&SEX==sex)
  weights<-rep(plogis(p[[if(sex=="Male")"hm_sel_blend" else "hf_sel_blend"]]),3)
  snow$SPECIES<-"Hybrid"
  snow$FITTED<-weights[match(snow$ERA,era_labels)]*snow$FITTED+
    (1-weights[match(snow$ERA,era_labels)])*tanner$FITTED
  snow$BSFRF_START_TARGET<-NA_real_
  snow
}))
curves<-rbind(curves,hybrid_curves)
curves$ERA <- factor(curves$ERA, levels = era_labels)
scales <- do.call(rbind, lapply(groups, function(group) {
  pr <- prior[[group]]; parts <- strsplit(group, "_", fixed = TRUE)[[1]]
  prefix <- paste0("survey_", parts[1], "_", substr(parts[2], 1, 1))
  offsets<-function(x)if(parts[1]=="snow")c(x[2],x[2],0) else c(x,0)
  if(parts[2]=="male"){
    sp<-paste0(prefix,"_stack")
    m1<-p[[paste0(sp,"_midpoint1")]]+offsets(p[[paste0(sp,"_midpoint1_offset")]])
    gap<-exp(p[[paste0(sp,"_log_gap")]]+offsets(p[[paste0(sp,"_log_gap_offset")]]))
    b1<-exp(p[[paste0(sp,"_log_slope1")]]+offsets(p[[paste0(sp,"_log_slope1_offset")]]))/30
    b2<-exp(p[[paste0(sp,"_log_slope2")]]+offsets(p[[paste0(sp,"_log_slope2_offset")]]))/30
    w<-plogis(p[[paste0(sp,"_logit_weight")]]+offsets(p[[paste0(sp,"_logit_weight_offset")]]))
    q<-exp(p[[paste0(sp,"_log_q")]]+offsets(p[[paste0(sp,"_log_q_offset")]]))
    return(data.frame(SPECIES=tools::toTitleCase(parts[1]),SEX="Male",ERA=era_labels,
      L50_COMPONENT1_MM=30*m1,L50_COMPONENT2_MM=30*(m1+gap),
      SLOPE1_PER_MM=b1,SLOPE2_PER_MM=b2,COMPONENT1_WEIGHT=w,CATCHABILITY=q,
      LOG_Q_OFFSET=offsets(p[[paste0(sp,"_log_q_offset")]])))
  }
  midpoint_raw <- p[[paste0(prefix,"_midpoint")]]+
    offsets(p[[paste0(prefix,"_midpoint_offset")]])
  lower_mid<-(min(pr$sizes)+5)/30;upper_mid<-(max(pr$sizes)-5)/30
  midpoint<-lower_mid+(upper_mid-lower_mid)*plogis(midpoint_raw)
  log_slope <- p[[paste0(prefix,"_log_slope")]]+
    offsets(p[[paste0(prefix,"_log_slope_offset")]])
  log_q<-p[[paste0(prefix,"_log_q")]]+offsets(p[[paste0(prefix,"_log_q_offset")]])
  data.frame(SPECIES = tools::toTitleCase(parts[1]),
    SEX = tools::toTitleCase(parts[2]), ERA = era_labels,
    L50_COMPONENT1_MM = 30*midpoint,L50_COMPONENT2_MM=NA_real_,
    SLOPE1_PER_MM=exp(log_slope)/30,SLOPE2_PER_MM=NA_real_,COMPONENT1_WEIGHT=1,
    CATCHABILITY=exp(log_q),
    LOG_Q_OFFSET=offsets(p[[paste0(prefix,"_log_q_offset")]]))
}))

write.csv(curves, file.path(out_dir, "shared_survey_selectivity_by_era.csv"),
          row.names = FALSE)
write.csv(scales, file.path(out_dir, "shared_survey_scale_by_era.csv"),
          row.names = FALSE)

p_curve <- ggplot(curves, aes(SIZE_MM, FITTED, color = ERA)) +
  geom_line(data=subset(curves,!is.na(BSFRF_START_TARGET)),
            aes(y = BSFRF_START_TARGET), color = "grey60", linetype = "dotted",
            linewidth = 0.8) +
  geom_line(linewidth = 1) +
  facet_grid(SEX ~ SPECIES, scales = "free_x") +
  scale_color_manual(values = c("#CC79A7", "#E69F00", "#0072B2")) +
  labs(
    title = "Survey selectivity by species, sex, and era",
    subtitle = "Dotted grey curves are species-sex BSFRF prior centers",
    x = "Carapace width (mm)", y = "Survey efficiency (q × selectivity)", color = "Survey era",
    caption = "Males use stacked logistics; females use single logistics. Hybrids are estimated parental blends."
  ) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "shared_survey_selectivity_by_era.png"), p_curve,
       width = 9, height = 8, dpi = 300, bg = "white")

cat("Model objective", f$fit$objective, "max gradient", f$max_gradient, "\n")
print(scales, row.names = FALSE)
