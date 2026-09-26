#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
refit <- "--refit" %in% args
project <- "--project" %in% args
draw_arg <- grep("^--draws=",args,value=TRUE)
draws <- if(length(draw_arg)) as.integer(sub("^--draws=","",draw_arg[1])) else 100L
if(!is.finite(draws)||draws<20)stop("--draws must be at least 20")

rscript <- file.path(R.home("bin"),"Rscript")
run <- function(script,extra=character()){
  message("\n==> ",script)
  status <- system2(rscript,c(shQuote(script),extra))
  if(!identical(status,0L))stop("Command failed: ",script,call.=FALSE)
}

if(refit)run("models/model/run_model.R")
run("models/model/check_hessian.R")
run("models/model/plot_model.R")
for(script in c("mortality_gam.R","terminal_molt_gam.R",
  "female_terminal_molt_gam.R","recruitment_centroid_gam.R",
  "male_abundance_surface_gam.R"))
  run(file.path("analyses/GAM/scripts",script))
if(project)run("models/model/projection/project_model.R",
  c("--years=15",paste0("--draws=",draws),"--seed=20260901"))
message("\nCurrent model diagnostics and GAM analyses completed successfully.")
