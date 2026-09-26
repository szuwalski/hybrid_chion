#!/usr/bin/env Rscript

required <- c("RTMB","ggplot2","patchwork","ggridges","GGally","scales",
  "rlang","mgcv","dplyr","tidyr")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)){
  message("Installing missing packages: ",paste(missing,collapse=", "))
  install.packages(missing,repos="https://cloud.r-project.org")
} else message("All required model and analysis packages are installed.")
