library(RTMB)
source("models/model/R/joint_classification_model.R")
args<-commandArgs(trailingOnly=TRUE)
argval<-function(key,default){z<-grep(paste0("^",key,"="),args,value=TRUE);if(length(z))sub(paste0("^",key,"="),"",z[1]) else default}
modern_prior_sd_multiplier<-as.numeric(argval("--modern-prior-sd-mult","1"))
fix_survey_q<-as.logical(as.integer(argval("--fix-survey-q","0")))
legacy_survey<-as.logical(as.integer(argval("--legacy-survey","1")))
fixed_growth<-as.logical(as.integer(argval("--fixed-growth","1")))
bounded_mortality<-as.logical(as.integer(argval("--bounded-mortality","1")))
mortality_link<-argval("--mortality-link","hazard")
if(!mortality_link%in%c("exponential","logit","hazard"))stop("Unknown mortality link")
domed_incidental<-as.logical(as.integer(argval("--domed-incidental","1")))
power_misclass<-as.logical(as.integer(argval("--power-misclass","1")))
closure_2022<-as.logical(as.integer(argval("--closure-2022","1")))
immature_snow_m_sd<-as.numeric(argval("--immature-snow-m-sd","0.7"))
all_mortality_deviation_sd<-as.numeric(argval("--all-mortality-deviation-sd","1"))
snow_male_mature_baseline_sd<-as.numeric(argval("--snow-male-mature-baseline-sd","0.10"))
snow_female_mature_baseline_sd<-as.numeric(argval("--snow-female-mature-baseline-sd","0.10"))
snow_female_immature_baseline_sd<-as.numeric(argval("--snow-female-immature-baseline-sd","0.05"))
snow_female_immature_deviation_sd<-as.numeric(argval("--snow-female-immature-deviation-sd","0.6"))
female_survey_era_midpoint_sd_mm<-as.numeric(argval("--female-survey-era-midpoint-sd-mm","3"))
female_survey_era_log_slope_sd<-as.numeric(argval("--female-survey-era-log-slope-sd","0.10"))
if(!is.finite(modern_prior_sd_multiplier)||modern_prior_sd_multiplier<1)
  stop("--modern-prior-sd-mult must be finite and at least 1")
d<-readRDS("models/model/data/model_input.rds")
if(fix_survey_q){
  # Express each BSFRF curve as shape relative to its largest supported anchor.
  # The anchor itself is fixed at one and is removed from the modern MVN so its
  # zero-variance log contrast does not make the covariance singular.
  normalize_prior_to_anchor<-function(pr){
    anchor_value<-pr$bootstrap_center_ratio[pr$anchor_index]
    pr$bootstrap_center_ratio<-pr$bootstrap_center_ratio/anchor_value
    pr$log_q_center[]<-0
    pr$log_q_sd[]<-0
    pr$log_shape_center<-log(pr$bootstrap_center_ratio)
    pr$shape_coef_center<-as.numeric(qr.solve(pr$spline_basis,pr$log_shape_center))
    aa<-match(pr$anchor_index,pr$modern_index)
    if(!is.na(aa)){
      keep<-setdiff(seq_along(pr$modern_index),aa)
      C<-pr$modern_covariance
      Crel<-C[keep,keep,drop=FALSE]-C[keep,aa]-
        matrix(C[aa,keep],nrow=length(keep),ncol=length(keep),byrow=TRUE)+C[aa,aa]
      ee<-eigen(Crel,symmetric=TRUE)
      ee$values<-pmax(ee$values,max(ee$values)*1e-6)
      Crel<-ee$vectors%*%(ee$values*t(ee$vectors))
      pr$modern_index<-pr$modern_index[keep]
      pr$modern_sizes<-pr$sizes[pr$modern_index]
      pr$modern_covariance<-Crel
      pr$modern_precision<-solve(Crel)
    }
    pr$modern_log_center<-log(pr$bootstrap_center_ratio[pr$modern_index])
    pr
  }
  for(group in c("snow_male","tanner_male","snow_female","tanner_female"))
    d$bsfrf_survey_prior[[group]]<-normalize_prior_to_anchor(d$bsfrf_survey_prior[[group]])
}
if(modern_prior_sd_multiplier!=1){
  for(group in c("snow_male","tanner_male","snow_female","tanner_female")){
    d$bsfrf_survey_prior[[group]]$modern_covariance<-
      d$bsfrf_survey_prior[[group]]$modern_covariance*modern_prior_sd_multiplier^2
    d$bsfrf_survey_prior[[group]]$modern_precision<-
      d$bsfrf_survey_prior[[group]]$modern_precision/modern_prior_sd_multiplier^2
  }
}
d$modern_bsfrf_prior_sd_multiplier<-modern_prior_sd_multiplier
d$survey_q_fixed_at_terminal_size<-fix_survey_q
d$legacy_two_fishery_survey<-as.integer(legacy_survey)
if(xor(is.finite(female_survey_era_midpoint_sd_mm),
       is.finite(female_survey_era_log_slope_sd)))
  stop("Both female survey-era offset SDs must be supplied together")
if(is.finite(female_survey_era_midpoint_sd_mm)){
  if(female_survey_era_midpoint_sd_mm<=0||female_survey_era_log_slope_sd<=0)
    stop("Female survey-era offset SDs must be positive")
  d$female_survey_era_midpoint_sd_mm<-female_survey_era_midpoint_sd_mm
  d$female_survey_era_log_slope_sd<-female_survey_era_log_slope_sd
}
# The directed snow-crab fishery was closed in 2022. The retained-catch record
# contains only two crab and is not a meaningful fishery abundance observation.
# Exclude that scalar retained-catch datum while retaining total-catch and all
# size-composition observations.
snow_fishery<-d$snow_male$data$fisheries$snow
exclude_retained_year<-2022L
drop_retained<-snow_fishery$retained_years==exclude_retained_year
if(closure_2022&&(sum(drop_retained)!=1L||snow_fishery$retained_numbers[drop_retained]!=2))
  stop("Expected the single two-crab retained-catch observation in 2022")
if(closure_2022){
 snow_fishery$retained_years<-snow_fishery$retained_years[!drop_retained]
 snow_fishery$retained_numbers<-snow_fishery$retained_numbers[!drop_retained]
# `active_years` was assembled upstream and still included 2022 because of the
# two-crab record. Explicitly set fishing mortality to zero in the closure year.
snow_fishery$active_years<-setdiff(snow_fishery$active_years,exclude_retained_year)
d$snow_male$data$fisheries$snow<-snow_fishery
# The closure also fixes Tanner bycatch mortality in the Snow fishery to zero.
d$tanner_male$data$fisheries$snow$active_years<-setdiff(
  d$tanner_male$data$fisheries$snow$active_years,exclude_retained_year)
d$excluded_observations<-c(d$excluded_observations,
  "Snow crab in directed snow fishery: 2022 retained catch excluded and F fixed to zero (fishery closed)")
}
# Permit the survey-supported episodic mortality of immature male snow crab.
# This is the annual deviation SD on the bounded-logit scale; the baseline
# mortality prior remains centered at M = 0.30 with logit-scale SD 0.20.
d$snow_male$data$Mdi_prior_sd<-immature_snow_m_sd
# The reviewed growth CSV is a protected source input and overrides the legacy
# observations embedded in the prepared RDS.
growth_path<-"data/growth_fit_observations.csv"
growth_input<-read.csv(growth_path,stringsAsFactors=FALSE)
required_growth_cols<-c("stock","sex","premolt","postmolt")
if(!identical(names(growth_input),required_growth_cols)||
   any(!is.finite(growth_input$premolt))||any(!is.finite(growth_input$postmolt)))
  stop("Invalid growth input: expected finite stock, sex, premolt, postmolt columns")
growth_group<-function(stock,sex){
  z<-growth_input[growth_input$stock==stock&growth_input$sex==sex,]
  if(!nrow(z))stop("No growth observations for ",stock," ",sex)
  z
}
gsm<-growth_group("Snow","Male");gtm<-growth_group("Tanner","Male")
gsf<-growth_group("Snow","Female");gtf<-growth_group("Tanner","Female")
if(fixed_growth){
 d$snow_male$data$growth_premolt<-gsm$premolt
 d$snow_male$data$growth_postmolt<-gsm$postmolt
 d$tanner_male$data$growth_premolt<-gtm$premolt
 d$tanner_male$data$growth_postmolt<-gtm$postmolt
 d$snow_female$data$growth<-list(premolt=gsf$premolt,postmolt=gsf$postmolt)
 d$tanner_female$data$growth<-list(premolt=gtf$premolt,postmolt=gtf$postmolt)
}
d$bounded_mortality<-as.integer(bounded_mortality)
d$mortality_link<-mortality_link
d$domed_incidental_selectivity<-as.integer(domed_incidental)
d$power_misclassification<-as.integer(power_misclass)
for(block in c("snow_male","tanner_male","hybrid_male","snow_female","tanner_female","hybrid_female")){
 d[[block]]$data$bounded_mortality<-as.integer(bounded_mortality)
 d[[block]]$data$mortality_link<-mortality_link
 if(block%in%c("snow_male","tanner_male"))
   d[[block]]$data$domed_incidental_selectivity<-as.integer(domed_incidental)
}
if(is.finite(all_mortality_deviation_sd)){
  if(all_mortality_deviation_sd<=0)stop("Mortality deviation SD must be positive")
  for(block in c("snow_male","tanner_male","snow_female","tanner_female","hybrid_female")){
    d[[block]]$data$Mdi_prior_sd<-all_mortality_deviation_sd
    d[[block]]$data$Mdm_prior_sd<-all_mortality_deviation_sd
  }
  d$hybrid_male$data$Mdev_prior_sd<-all_mortality_deviation_sd
}
if(is.finite(snow_male_mature_baseline_sd)){
  if(snow_male_mature_baseline_sd<=0)
    stop("Mature Snow-male baseline mortality prior SD must be positive")
  d$snow_male$data$logMm_prior_sd<-snow_male_mature_baseline_sd
}
if(is.finite(snow_female_mature_baseline_sd)){
  if(snow_female_mature_baseline_sd<=0)
    stop("Mature Snow-female baseline mortality prior SD must be positive")
  d$snow_female$data$logMm_prior_sd<-snow_female_mature_baseline_sd
}
if(is.finite(snow_female_immature_baseline_sd)){
  if(snow_female_immature_baseline_sd<=0)
    stop("Immature Snow-female baseline mortality prior SD must be positive")
  d$snow_female$data$logMi_prior_sd<-snow_female_immature_baseline_sd
}
if(is.finite(snow_female_immature_deviation_sd)){
  if(snow_female_immature_deviation_sd<=0)
    stop("Immature Snow-female mortality-deviation SD must be positive")
  d$snow_female$data$Mdi_prior_sd<-snow_female_immature_deviation_sd
}
d$tanner_male$data$simple_survey_selectivity<-1
d$hybrid_male$data$simple_survey_selectivity<-1
d$snow_male$data$simple_survey_selectivity<-1
d$snow_female$data$simple_survey_selectivity<-1
d$hybrid_female$data$simple_survey_selectivity<-1
d$tanner_female$data$simple_survey_selectivity<-1
p<-J1_parameters(d)
if(fix_survey_q){
  for(pre in c("snow_m","tanner_m")){
    p[[paste0("survey_",pre,"_stack_log_q")]]<-0
    p[[paste0("survey_",pre,"_stack_log_q_offset")]][]<-0
  }
  for(pre in c("snow_f","tanner_f")){
    p[[paste0("survey_",pre,"_log_q")]]<-0
    p[[paste0("survey_",pre,"_log_q_offset")]][]<-0
  }
}
path<-argval("--fit","models/model/fits/model_fit.rds")
current_path<-"models/model/fits/model_fit.rds"
warm_path<-"models/archive/rtmb_joint_two_fishery_bsfrf/fits/JMR_two_fishery_logistic_MAP.rds"
bridge_warm_path<-argval("--warm","")
if(file.exists(path)||file.exists(bridge_warm_path)||file.exists(current_path)||file.exists(warm_path)){
 old_object<-readRDS(if(file.exists(path))path else if(file.exists(bridge_warm_path))bridge_warm_path else if(file.exists(current_path))current_path else warm_path);old<-old_object$par_list
 for(nm in intersect(names(p),names(old))){
  if(grepl("^survey_(snow|tanner)_[mf]_",nm))next
  if(nm%in%c("survey_tm","survey_sm","survey_sf","survey_tf"))next
  # These names had different meanings under the old logistic curve.
  if(nm%in%c("class_m_intercept","class_m_logslope",
             "class_f_intercept","class_f_logslope"))next
  # Parental baseline-M coordinates changed from log excess mortality to a
  # bounded-logit scale; initialize them under the new parameterization.
  if(grepl("^(sm|tm|sf|tf)_logM[im]$",nm))next
  if(grepl("_(logF|Fdev|fish_|ret_)",nm))next
  if(nm%in%c("sm_logNi","sm_logNm")&&length(p[[nm]])==length(old[[nm]])+1L){
   p[[nm]][-1]<-old[[nm]]
  }else if(nm%in%c("sf_logNi","sf_logNm","tf_logNi","tf_logNm","hf_logNi","hf_logNm")&&
      length(p[[nm]])>length(old[[nm]])){
   # Split only the warm-start pooled tail, using restored initial proportions.
   n<-length(old[[nm]]);tail_idx<-n:length(p[[nm]])
   weights<-exp(p[[nm]][tail_idx]);weights<-weights/sum(weights)
   p[[nm]][seq_len(n-1)]<-old[[nm]][seq_len(n-1)]
   p[[nm]][tail_idx]<-old[[nm]][n]+log(weights)
  }else if(nm%in%c("sf_logNi","sf_logNm","tf_logNi","tf_logNm","hf_logNi","hf_logNm")&&
      length(p[[nm]])<length(old[[nm]])){
   n<-length(p[[nm]])
   p[[nm]]<-c(old[[nm]][seq_len(n-1)],log(sum(exp(old[[nm]][n:length(old[[nm]])]))))
  }else if(length(p[[nm]])==length(old[[nm]]))p[[nm]]<-old[[nm]]
 }
}
if(legacy_survey){
  legacy_fit_path<-"models/archive/rtmb_joint_two_fishery/fits/JMR_two_fishery_MAP.rds"
  legacy_par<-readRDS(legacy_fit_path)$par_list
  for(nm in c("survey_sm","survey_tm","survey_sf","survey_tf"))p[[nm]]<-legacy_par[[nm]]
}
# Condition the assessment on parental growth estimated from the direct molt
# observations.  These MLEs are not allowed to trade off against abundance or
# size-composition likelihoods; hybrid growth remains an estimated blend of
# the resulting parental curves.
direct_growth_fit<-function(premolt,postmolt,center){
  X<-cbind(1,premolt);coef<-as.numeric(solve(crossprod(X),crossprod(X,postmolt)))
  resid<-postmolt-as.vector(X%*%coef)
  c(level=coef[1]+coef[2]*center,slope30=coef[2]*30,
    log_sd=log(sqrt(mean(resid^2))),alpha=coef[1],beta=coef[2],
    sd=sqrt(mean(resid^2)),n=length(premolt))
}
growth_fits<-list(
  sm=direct_growth_fit(d$snow_male$data$growth_premolt,
    d$snow_male$data$growth_postmolt,d$growth_center_male),
  tm=direct_growth_fit(d$tanner_male$data$growth_premolt,
    d$tanner_male$data$growth_postmolt,d$growth_center_male),
  sf=direct_growth_fit(d$snow_female$data$growth$premolt,
    d$snow_female$data$growth$postmolt,d$growth_center_female),
  tf=direct_growth_fit(d$tanner_female$data$growth$premolt,
    d$tanner_female$data$growth$postmolt,d$growth_center_female))
for(pre in names(growth_fits)){
  p[[paste0(pre,"_growth_level")]]<-unname(growth_fits[[pre]]["level"])
  p[[paste0(pre,"_growth_slope30")]]<-unname(growth_fits[[pre]]["slope30"])
  p[[paste0(pre,"_log_growth_sd")]]<-unname(growth_fits[[pre]]["log_sd"])
}
inactive<-c("class_f_intercept","class_f_logslope","common_rec",
  "log_rec_sigma","log_M_sigma","log_rec_loading","log_M_loading",
  "hf_selshift_pre","hf_selt_pre","hf_selshift","hf_selt","sf_mat50","sf_logmatslope",
  "tf_mat50","tf_logmatslope","hf_mat50","hf_logmatslope")
if(!legacy_survey)inactive<-c(inactive,"survey_sm","survey_tm","survey_sf","survey_tf")
inactive<-c(inactive,unlist(lapply(c("snow_m","tanner_m"),function(pre)
  paste0("survey_",pre,c("_midpoint","_log_slope","_log_q",
    "_midpoint_offset","_log_slope_offset","_log_q_offset")))))
inactive<-c(inactive,unlist(lapply(c("snow_m","tanner_m"),function(pre)
  paste0("survey_",pre,"_stack_",c("log_gap_offset","log_slope1_offset",
    "log_slope2_offset","logit_weight_offset")))))
if(fix_survey_q)inactive<-c(inactive,
  "survey_snow_m_stack_log_q","survey_snow_m_stack_log_q_offset",
  "survey_tanner_m_stack_log_q","survey_tanner_m_stack_log_q_offset",
  "survey_snow_f_log_q","survey_snow_f_log_q_offset",
    "survey_tanner_f_log_q","survey_tanner_f_log_q_offset")
if(legacy_survey){
  current_survey_names<-grep("^survey_(snow|tanner)_[mf]_",names(p),value=TRUE)
  inactive<-c(inactive,current_survey_names)
}
if(fixed_growth)inactive<-c(inactive,unlist(lapply(c("sm","tm","sf","tf"),function(pre)
  paste0(pre,c("_growth_level","_growth_slope30","_log_growth_sd")))))
inactive<-c(inactive,unlist(lapply(c('sm','tm','sf','tf'),function(pre)
 paste0(pre,c('_selshift_pre','_selt_pre','_selshift','_selt')))))
inactive<-c(inactive,unlist(lapply(c("sm","tm"),function(pre)
 paste0(pre,c("_fish_bdev","_fish_desc_a","_fish_desc_logb",
   "_ret_bdev","_inc_ret_bdev")))))
if(!domed_incidental)inactive<-c(inactive,unlist(lapply(c("sm","tm"),function(pre)
 paste0(pre,c("_inc_fish_bdev","_inc_fish_desc_a","_inc_fish_desc_logb")))))
if(is.null(d$parental_mature_hybrid_recruitment)||d$parental_mature_hybrid_recruitment!=1){
 p$jr_parent_beta[]<-0
 inactive<-c(inactive,"jr_parent_beta")
}
if(!is.null(d$snow_male$fixed_survey_sel)){
  p$sm_selshift_pre<-0;p$sm_selt_pre<-0;p$sm_selshift<-0;p$sm_selt<-0
  inactive<-c(inactive,"sm_selshift_pre","sm_selt_pre","sm_selshift","sm_selt")
}
if(!is.null(d$tanner_male$fixed_survey_sel)){
  p$tm_selshift_pre<-0;p$tm_selt_pre<-0;p$tm_selshift<-0;p$tm_selt<-0
  inactive<-c(inactive,"tm_selshift_pre","tm_selt_pre","tm_selshift","tm_selt")
}
if(is.null(d$shared_mortality_driver)||d$shared_mortality_driver!=1){
  p$common_M[]<-0;p$M_loading[]<-0
  inactive<-c(inactive,"common_M","M_loading")
}
map<-setNames(lapply(inactive,function(nm)
  factor(rep(NA,length(p[[nm]])))),inactive)
active_map<-function(active){z<-rep(NA_integer_,length(active));z[active]<-seq_len(sum(active));factor(z)}
if(legacy_survey)for(nm in c("survey_sm","survey_tm","survey_sf","survey_tf"))
  map[[nm]]<-active_map(rep(c(TRUE,FALSE,TRUE,FALSE,FALSE),2))
# Snow has no parental survey data in 1975-1981. Its first historical curve is
# structurally the 1982-1988 curve; the unused first offset coordinate is fixed.
if(!legacy_survey){
  for(sex in "f")for(suffix in c("midpoint_offset","log_slope_offset",
    if(!fix_survey_q)"log_q_offset"))
    map[[paste0("survey_snow_",sex,"_",suffix)]]<-active_map(c(FALSE,TRUE))
  for(suffix in c("midpoint1","log_gap","log_slope1","log_slope2","logit_weight",
    if(!fix_survey_q)"log_q"))
    map[[paste0("survey_snow_m_stack_",suffix,"_offset")]]<-active_map(c(FALSE,TRUE))
}
fish_map<-function(block,fishery,field){
  yrs<-head(d[[block]]$data$years,-1)
  active_years<-d[[block]]$data$fisheries[[fishery]][[field]]
  active_map(yrs%in%active_years)
}
for(z in list(c("sm","snow_male","snow","tanner"),
              c("tm","tanner_male","tanner","snow"))){
  pre<-z[1];block<-z[2];directed<-z[3];incidental<-z[4]
  map[[paste0(pre,"_Fdev")]]<-fish_map(block,directed,"active_years")
  map[[paste0(pre,"_fish_adev")]]<-fish_map(block,directed,"total_size_years")
  map[[paste0(pre,"_ret_adev")]]<-fish_map(block,directed,"retained_size_years")
  map[[paste0(pre,"_inc_Fdev")]]<-fish_map(block,incidental,"active_years")
  map[[paste0(pre,"_inc_fish_adev")]]<-fish_map(block,incidental,"total_size_years")
  map[[paste0(pre,"_inc_fish_bdev")]]<-fish_map(block,incidental,"total_size_years")
  map[[paste0(pre,"_inc_ret_adev")]]<-fish_map(block,incidental,"retained_size_years")
}
M_active<-function(block)!head(d[[block]]$data$years,-1)%in%d$fixed_mortality_deviation_years
for(z in list(c("sm_Mdi","sm_Mdm","snow_male"),
  c("tm_Mdi","tm_Mdm","tanner_male"),c("hm_Mdev","hm_Mdev","hybrid_male"),
  c("sf_Mdi","sf_Mdm","snow_female"),c("tf_Mdi","tf_Mdm","tanner_female"),
  c("hf_Mdi","hf_Mdm","hybrid_female"))){
  active<-M_active(z[3])
  p[[z[1]]][!active]<-0;p[[z[2]]][!active]<-0
  map[[z[1]]]<-active_map(active)
  map[[z[2]]]<-active_map(active)
}
common_M_active<-!head(d$years,-1)%in%d$fixed_mortality_deviation_years
p$common_M[!common_M_active]<-0
if(!is.null(d$shared_mortality_driver)&&d$shared_mortality_driver==1)
  map$common_M<-active_map(common_M_active)
# Stage 1 estimates mean curves and annual F with all annual selectivity
# midpoint deviations fixed. Stage 2 releases the supported midpoint changes.
stage_map<-map
for(pre in c("sm","tm"))for(suf in c("fish_adev","ret_adev",
  "inc_fish_adev","inc_fish_bdev","inc_ret_adev")){
  nm<-paste0(pre,"_",suf)
  stage_map[[nm]]<-factor(rep(NA,length(p[[nm]])))
}
obj_stage<-MakeADFun(function(x)do.call(J1_objective(d),x),p,map=stage_map,silent=TRUE)
stage_fit<-nlminb(obj_stage$par,obj_stage$fn,obj_stage$gr,
  control=list(eval.max=6000,iter.max=2500,rel.tol=1e-11,x.tol=1e-8))
stage_par<-obj_stage$env$parList(stage_fit$par)
for(nm in intersect(names(p),names(stage_par)))
  if(length(p[[nm]])==length(stage_par[[nm]]))p[[nm]]<-stage_par[[nm]]
cat("JMR stage-1 objective",stage_fit$objective,"gradient",
  max(abs(obj_stage$gr(stage_fit$par))),"\n")
obj<-MakeADFun(function(x)do.call(J1_objective(d),x),p,map=map,silent=TRUE)
start<-obj$par
if(file.exists(path)){
 oldfit<-readRDS(path)$fit$par
 if(length(oldfit)==length(start)&&identical(names(oldfit),names(start))&&all(is.finite(oldfit)))start<-oldfit
}
cat("JMR initial objective",obj$fn(start),"parameters",length(start),"\n")
if("--no-optimize"%in%commandArgs(trailingOnly=TRUE))quit(save="no",status=0)
fit<-nlminb(start,obj$fn,obj$gr,
  control=list(eval.max=10000,iter.max=5000,rel.tol=1e-12,x.tol=1e-9))
for(i in 2:8){
 fit<-nlminb(fit$par,obj$fn,obj$gr,
   control=list(eval.max=10000,iter.max=5000,rel.tol=1e-13,x.tol=1e-10))
 if(max(abs(obj$gr(fit$par)))<1e-3)break
}
for(cycle in 1:6){
 g<-obj$gr(fit$par);idx<-order(abs(g),decreasing=TRUE)[1:24]
 for(j in idx){
  base<-fit$par;f1<-function(value){z<-base;z[j]<-value;obj$fn(z)}
  z<-optimize(f1,base[j]+c(-.08,.08),tol=1e-12)
  if(z$objective<fit$objective){fit$par[j]<-z$minimum;fit$objective<-z$objective}
 }
 if(max(abs(obj$gr(fit$par)))<1e-3)break
}
# Correlated ascending/descending selectivity coordinates can leave nlminb at
# false convergence even after scalar polishing. A multivariate quasi-Newton
# pass resolves those directions, followed by a final bounded optimizer pass.
if(max(abs(obj$gr(fit$par)))>=1e-3){
  bfgs<-optim(fit$par,obj$fn,obj$gr,method="BFGS",
    control=list(maxit=2500,reltol=1e-12))
  if(is.finite(bfgs$value)&&bfgs$value<=fit$objective){
    fit$par<-bfgs$par;fit$objective<-bfgs$value
  }
  fit<-nlminb(fit$par,obj$fn,obj$gr,
    control=list(eval.max=15000,iter.max=7500,rel.tol=1e-13,x.tol=1e-10))
}
obj$fn(fit$par);report<-obj$report();gradient<-obj$gr(fit$par)
max_gradient<-max(abs(gradient));gradient_order<-order(abs(gradient),decreasing=TRUE)
saveRDS(list(model="JMR_two_fishery_logistic_MAP",parent="JMR_two_fishery_bsfrf_MAP",data=d,fit=fit,
  report=report,max_gradient=max_gradient,par_list=obj$env$parList(fit$par),map=map,
  parental_growth_source="fixed direct Gaussian MLE from premolt-postmolt observations",
  parental_growth_input=growth_path,
  parental_growth_input_md5=unname(tools::md5sum(growth_path)),
  parental_growth_fits=growth_fits,
  fit_type="joint penalized mode / MAP",gradient=gradient,
  mortality_parameterization="floor_plus_scaled_logistic_v1",
  fishery_selectivity_parameterization="directed_logistic_incidental_double_logistic_v3",
  modern_bsfrf_prior_sd_multiplier=modern_prior_sd_multiplier,
  survey_q_fixed_at_terminal_size=fix_survey_q,
  legacy_two_fishery_survey=legacy_survey,
  bridge_settings=list(fixed_growth=fixed_growth,bounded_mortality=bounded_mortality,
    mortality_link=mortality_link,
    all_mortality_deviation_sd=all_mortality_deviation_sd,
    snow_male_mature_baseline_sd=snow_male_mature_baseline_sd,
    snow_female_mature_baseline_sd=snow_female_mature_baseline_sd,
    snow_female_immature_baseline_sd=snow_female_immature_baseline_sd,
    snow_female_immature_deviation_sd=snow_female_immature_deviation_sd,
    female_survey_era_midpoint_sd_mm=female_survey_era_midpoint_sd_mm,
    female_survey_era_log_slope_sd=female_survey_era_log_slope_sd,
    domed_incidental=domed_incidental,power_misclassification=power_misclass,
    closure_2022=closure_2022,immature_snow_m_sd=immature_snow_m_sd),
  selectivity_parameterization=if(legacy_survey)
    "legacy_two_era_single_logistic_original_priors_v1" else
    "bsfrf_regularized_stacked_male_logistic_female_q_v19"),path)
cat("JMR code",fit$convergence,"objective",fit$objective,"gradient",max_gradient,"\n")
print(head(data.frame(parameter=names(fit$par)[gradient_order],
  gradient=gradient[gradient_order]),10),row.names=FALSE)
cat("JMR recruitment coefficients",report$jr_parent_beta,"\n")
cat("JMR mortality loadings\n");print(report$M_loading)
