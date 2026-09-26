# =============================================================================
# Shared utilities
# =============================================================================
source('models/model/R/survey_selectivity.R')

# Three survey eras. Each species-sex group has a modern ascending logistic;
# the two historical curves are offsets from its modern midpoint and log-slope.
survey_era_index <- function(years) {
  ifelse(years <= 1981, 1L, ifelse(years <= 1988, 2L, 3L))
}

logistic_survey_curves <- function(common_sizes, midpoint, log_slope,
                                   log_q, midpoint_offset, log_slope_offset,
                                   log_q_offset, no_early_data=FALSE) {
  offsets<-function(x)if(no_early_data)c(x[2],x[2],0) else c(x,0)
  midpoint_raw <- midpoint + offsets(midpoint_offset)
  lower_mid<-(min(common_sizes)+5)/30
  upper_mid<-(max(common_sizes)-5)/30
  era_midpoint <- lower_mid+(upper_mid-lower_mid)*plogis(midpoint_raw)
  era_log_slope <- log_slope + offsets(log_slope_offset)
  era_log_q <- log_q + offsets(log_q_offset)
  curves <- sapply(seq_len(3), function(era){
    shape<-plogis(exp(era_log_slope[era]) *
      (common_sizes / 30 - era_midpoint[era]))
    # Normalize within the modeled grid so q is efficiency at the terminal
    # size, rather than an unidentified extrapolated asymptote.
    exp(era_log_q[era])*shape/shape[length(shape)]
  })
  curves
}

logistic_survey_matrix <- function(block, common_sizes, midpoint, log_slope,
                                   log_q, midpoint_offset, log_slope_offset,
                                   log_q_offset, no_early_data=FALSE) {
  idx <- match(block$data$sizes, common_sizes)
  stopifnot(!anyNA(idx))
  curves<-logistic_survey_curves(common_sizes,midpoint,log_slope,log_q,
    midpoint_offset,log_slope_offset,log_q_offset,no_early_data)
  t(curves[idx, survey_era_index(block$data$years), drop = FALSE])
}

logistic_survey_prior_nll <- function(prior, curves, no_early_data=FALSE) {
  eps<-1e-12
  # Modern paired-station bootstrap covariance retains dependence among sizes.
  modern_residual<-log(curves[prior$modern_index,3]+eps)-prior$modern_log_center
  nll<-0.5*sum(modern_residual*
    as.vector(prior$modern_precision%*%modern_residual))
  # Historical curves use the same experimental center with uncertainty
  # inflated by era. Snow has no earliest parental observations and shares its
  # earliest curve with 1982-1988, so that curve is penalized only once.
  historical_eras<-if(no_early_data)2L else 1:2
  for(era in historical_eras){
    sd<-sqrt(prior$log_q_sd[era]^2+prior$log_shape_sd[,era]^2)
    residual<-log(curves[,era]+eps)-log(prior$bootstrap_center_ratio+eps)
    nll<-nll+0.5*sum((residual/sd)^2)
  }
  nll
}

stacked_logistic_survey_curves <- function(common_sizes,midpoint1,log_gap,
  log_slope1,log_slope2,logit_weight,log_q,midpoint1_offset,
  log_gap_offset,log_slope1_offset,log_slope2_offset,
  logit_weight_offset,log_q_offset,no_early_data=FALSE){
  offsets<-function(x)if(no_early_data)c(x[2],x[2],0) else c(x,0)
  m1<-midpoint1+offsets(midpoint1_offset)
  gap<-exp(log_gap+offsets(log_gap_offset))
  b1<-exp(log_slope1+offsets(log_slope1_offset))
  b2<-exp(log_slope2+offsets(log_slope2_offset))
  w<-plogis(logit_weight+offsets(logit_weight_offset))
  lq<-log_q+offsets(log_q_offset)
  sapply(seq_len(3),function(era){
    z<-common_sizes/30
    shape<-w[era]*plogis(b1[era]*(z-m1[era]))+
      (1-w[era])*plogis(b2[era]*(z-m1[era]-gap[era]))
    exp(lq[era])*shape/tail(shape,1)
  })
}

stacked_logistic_survey_matrix <- function(block,common_sizes,...){
  curves<-stacked_logistic_survey_curves(common_sizes,...)
  idx<-match(block$data$sizes,common_sizes);stopifnot(!anyNA(idx))
  t(curves[idx,survey_era_index(block$data$years),drop=FALSE])
}

# Adjust a fixed standalone selectivity curve on the logit scale. `shift`
# changes its overall level, while `tilt` changes its size-dependent slope.
joint_sel <- function(base, shift, tilt, z) {
  zz <- (z - mean(z)) / 30
  1 / (1 + exp(-(log(base / (1 - base)) + shift + tilt * zz)))
}

# Fixed linear projection matrix with constant endpoint extrapolation. Multiplying
# this matrix by a parental curve puts it on the hybrid size grid and remains
# compatible with automatic differentiation when the parental curve is fitted.
curve_projection <- function(from, to) {
  W <- matrix(0, length(to), length(from))
  for (i in seq_along(to)) {
    if (to[i] <= from[1]) {
      W[i, 1] <- 1
    } else if (to[i] >= from[length(from)]) {
      W[i, length(from)] <- 1
    } else {
      j <- max(which(from <= to[i]))
      w <- (to[i] - from[j]) / (from[j + 1] - from[j])
      W[i, j] <- 1 - w
      W[i, j + 1] <- w
    }
  }
  W
}

# Reconstruct the omitted final annual deviation and impose a zero-sum
# constraint. The active version applies this only in years with fishing.
center_last <- function(x) {
  z <- c(x, 0)
  z - sum(z) / length(z)
}

center_active_last <- function(x, active) {
  z <- c(x, 0)
  out <- z * 0
  use <- which(active > 0)
  if (length(use)) {
    out[use] <- z[use] - sum(z[use]) / length(use)
  }
  out
}

# Lognormal abundance kernel and multinomial-style composition kernel.
ln_kernel <- function(obs, pred, cv, eps = 1e-12) {
  (log(pred + eps) - log(obs))^2 / (2 * log1p(cv^2))
}

ce <- function(obs, pred, w, eps = 1e-12) {
  oo <- obs / (sum(obs) + eps)
  pp <- pred / (sum(pred) + eps)
  -w * sum(oo * log(pp + eps))
}

# Proper Gaussian regularization for an initial log abundance-at-size vector.
# Second differences discourage unsupported spikes. Weak priors on its mean and
# mean slope remove the constant and linear null spaces of that curvature
# penalty. Centers come from the corresponding accepted standalone fit.
initial_smooth_nll <- function(
  x, center, curve_sd = 1, level_sd = 2, slope_sd = 1
) {
  center <- log(pmax(center, 1))
  -sum(dnorm(diff(x, differences = 2), 0, curve_sd, log = TRUE)) -
    dnorm(mean(x), mean(center), level_sd, log = TRUE) -
    dnorm(mean(diff(x)), mean(diff(center)), slope_sd, log = TRUE)
}

# Gaussian postmolt-on-premolt growth model and its size-transition matrix.
# Destination rows and source columns follow the population projection G %*% N.
# Shrinkage is disallowed, and the final size bin is a plus group.
growth_transition <- function(sizes, alpha, beta, sd) {
  ns <- length(sizes)
  zero <- alpha * 0
  edges <- c(sizes - 2.5, max(sizes) + 2.5)
  G <- zero + matrix(0, ns, ns)
  for (source in seq_len(ns)) {
    if (source == ns) {
      G[ns, source] <- 1
    } else {
      mu <- alpha + beta * sizes[source]
      G[source, source] <- pnorm(edges[source + 1L], mu, sd)
      if (source + 1L < ns) {
        for (destination in seq.int(source + 1L, ns - 1L)) {
          G[destination, source] <-
            pnorm(edges[destination + 1L], mu, sd) -
            pnorm(edges[destination], mu, sd)
        }
      }
      G[ns, source] <- 1 - pnorm(edges[ns], mu, sd)
    }
  }
  G
}

growth_data_nll <- function(premolt, postmolt, alpha, beta, sd) {
  -sum(dnorm(postmolt, alpha + beta * premolt, sd, log = TRUE))
}

# Growth optimizer coordinates use fitted postmolt width at a fixed premolt
# center and change per 30 mm. This is exactly equivalent to alpha and beta,
# while avoiding their strong correlation and very different numerical scales.
growth_coordinates <- function(alpha, beta, center, scale = 30) {
  c(level = alpha + beta * center, slope30 = beta * scale)
}

# Three-bin recruitment allocation with the third bin as the reference. Two
# unconstrained logits produce positive proportions that sum exactly to one.
recruitment_proportions <- function(logits) {
  z <- exp(c(logits, 0))
  z / sum(z)
}


# =============================================================================
# Population-process blocks
# =============================================================================

female_block <- function(
  b, logNi, logNm, logmeanrec, rec_logits, recdev, logMi, logMm, Mdi, Mdm,
  mat50, logmatslope, selshift_pre, selt_pre, selshift, selt, penalize_process = TRUE,
  G = b$G, sel_override = NULL, logMi_sd = .10, logMm_sd = .20, Mdi_sd = .30,
  Mdm_sd = .50, M_override = NULL, rec_offset = NULL, M_offset = NULL
) {
  d <- b$data
  if(!is.null(d$logMi_prior_sd))logMi_sd<-d$logMi_prior_sd
  if(!is.null(d$logMm_prior_sd))logMm_sd<-d$logMm_prior_sd
  if(!is.null(d$Mdi_prior_sd))Mdi_sd<-d$Mdi_prior_sd
  if(!is.null(d$Mdm_prior_sd))Mdm_sd<-d$Mdm_prior_sd
  ny <- length(d$years)
  ns <- length(d$sizes)
  zero <- logNi[1] * 0
  M_floor <- if(is.null(d$mortality_floor)) 0 else d$mortality_floor
  M_link<-if(!is.null(d$mortality_link))d$mortality_link else if(is.null(d$bounded_mortality)||d$bounded_mortality==1)"logit" else "exponential"
  M_center <- switch(M_link,
    logit=qlogis(pmin(pmax((.30-M_floor)/(1-M_floor),1e-6),1-1e-6)),
    hazard=log(-log(pmin(pmax((1-.30)/(1-M_floor),1e-6),1-1e-6))),
    exponential=log(pmax(.30-M_floor,1e-6)))
  M_from_eta<-function(eta)switch(M_link,
    logit=M_floor+(1-M_floor)/(1+exp(-eta)),
    hazard=-log(1-(M_floor+(1-M_floor)*(1-exp(-exp(eta))))),
    exponential=M_floor+exp(eta))
  M_active <- as.numeric(!(head(d$years,-1) %in%
    if(is.null(d$fixed_mortality_deviation_years)) integer(0) else d$fixed_mortality_deviation_years))
  center_M <- function(z) center_active_last(z,M_active)

  sel <- if (is.null(sel_override)) {
    spre<-joint_sel(b$survey_base,selshift_pre,selt_pre,d$sizes)
    spost<-joint_sel(b$survey_base,selshift,selt,d$sizes)
    out<-zero+matrix(0,ny,ns)
    out[d$years<1988,]<-matrix(spre,sum(d$years<1988),ns,byrow=TRUE)
    out[d$years>=1988,]<-matrix(spost,sum(d$years>=1988),ns,byrow=TRUE)
    out
  } else {
    if(is.null(dim(sel_override)))matrix(sel_override,ny,ns,byrow=TRUE) else sel_override
  }
  pm <- d$p_mature_by_year

  # The tanh transforms keep penalized annual deviations in a numerically
  # stable range. They are omitted when shared process drivers are used.
  di <- if (penalize_process) 4 * tanh(center_M(Mdi) / 4) else center_M(Mdi)
  dm <- if (penalize_process) 4 * tanh(center_M(Mdm) / 4) else center_M(Mdm)
  if(!is.null(M_offset)){di<-di+center_M(M_offset$Mi)
    dm<-dm+center_M(M_offset$Mm)}
  dr <- if (penalize_process) 5 * tanh(center_last(recdev) / 5) else center_last(recdev)
  if(!is.null(rec_offset))dr<-dr+center_last(rec_offset)
  Mi <- if (is.null(M_override)) M_from_eta(logMi+di) else M_override$Mi
  Mm <- if (is.null(M_override)) M_from_eta(logMm+dm) else M_override$Mm
  rec <- exp(logmeanrec + dr)

  # Recruitment is allocated to the first three size bins using an estimated
  # three-category simplex.
  rp <- recruitment_proportions(rec_logits)
  Ni <- zero + matrix(0, ny, ns)
  Nm <- zero + matrix(0, ny, ns)
  Ni[1, ] <- exp(logNi)
  Nm[1, ] <- exp(logNm)
  for (y in seq_len(ny - 1)) {
    spring <- isTRUE(d$spring_recruitment)
    before <- if(spring) .59 else 1
    after <- if(spring) .41 else 0
    # Spring event: fishery (where modeled), growth, recruitment, maturation;
    # then survival to the next survey. Legacy fits retain their old timing.
    grown <- as.vector(G %*% (Ni[y, ] * exp(-before*Mi[y])))
    grown[1:3] <- grown[1:3] + rec[y] * rp
    # Ogives are dated by receiving survey year: observation y+1 describes
    # terminal molt during the process interval y -> y+1.
    newly <- grown * pm[y+1L,]
    Ni[y + 1, ] <- (grown - newly)*exp(-after*Mi[y])
    Nm[y + 1, ] <- (Nm[y, ] * exp(-before*Mm[y]) + newly)*exp(-after*Mm[y])
  }
  vi <- Ni*sel
  vm <- Nm*sel
  nll <- if(is.null(M_override))
    -dnorm(logMi,M_center,logMi_sd,log=TRUE)-dnorm(logMm,M_center,logMm_sd,log=TRUE) else
    -dnorm(logMi,0,1.5,log=TRUE)-dnorm(logMm,0,1.5,log=TRUE)
  nll <- nll-
    sum(dnorm(rec_logits,log(c(.75,.20)/.05),1.5,log=TRUE))
  if(is.null(sel_override))nll<-nll-dnorm(selshift_pre,0,.35,log=TRUE)-
    dnorm(selt_pre,0,.35,log=TRUE)-dnorm(selshift,0,.35,log=TRUE)-dnorm(selt,0,.35,log=TRUE)
  nll<-nll+initial_smooth_nll(logNi,b$init_i)+initial_smooth_nll(logNm,b$init_m)
  if(penalize_process)nll<-nll-sum(dnorm(recdev,0,.6,log=TRUE))-
    sum(dnorm(Mdi,0,Mdi_sd,log=TRUE))-sum(dnorm(Mdm,0,Mdm_sd,log=TRUE))
  list(Ni=Ni,Nm=Nm,vi=vi,vm=vm,Mi=Mi,Mm=Mm,rec=rec,rec_prop=rp,sel=sel,pm=pm,nll=nll)
}

# Male parental stocks include a fishery between two natural-mortality periods.
# Immature animals grow and may terminally molt into the mature state.
parent_male_block <- function(b,logNi,logNm,logmeanrec,rec_logits,recdev,logMi,logMm,
  Mdi,Mdm,logF,Fdev,fish_a,fish_logb,fish_adev,fish_bdev,
  fish_desc_a,fish_desc_logb,
  ret_a,ret_logb,ret_adev,ret_bdev,
  inc_logF,inc_Fdev,inc_fish_a,inc_fish_logb,inc_fish_adev,inc_fish_bdev,
  inc_fish_desc_a,inc_fish_desc_logb,
  inc_ret_a,inc_ret_logb,inc_ret_adev,inc_ret_bdev,
  selshift_pre,selt_pre,selshift,selt,
  penalize_process=TRUE,G=b$G,
  M_offset=NULL,sel_override=NULL) {
  d<-b$data;ny<-length(d$years);ns<-length(d$sizes);zero<-logNi[1]*0
  M_floor<-if(is.null(d$mortality_floor))0 else d$mortality_floor
  M_link<-if(!is.null(d$mortality_link))d$mortality_link else if(is.null(d$bounded_mortality)||d$bounded_mortality==1)"logit" else "exponential"
  M_center<-switch(M_link,
    logit=qlogis(pmin(pmax((.30-M_floor)/(1-M_floor),1e-6),1-1e-6)),
    hazard=log(-log(pmin(pmax((1-.30)/(1-M_floor),1e-6),1-1e-6))),
    exponential=log(pmax(.30-M_floor,1e-6)))
  M_from_eta<-function(eta)switch(M_link,
    logit=M_floor+(1-M_floor)/(1+exp(-eta)),
    hazard=-log(1-(M_floor+(1-M_floor)*(1-exp(-exp(eta))))),
    exponential=M_floor+exp(eta))
  logMi_prior_sd<-if(is.null(d$logMi_prior_sd)).10 else d$logMi_prior_sd
  logMm_prior_sd<-if(is.null(d$logMm_prior_sd)).30 else d$logMm_prior_sd
  Mdi_prior_sd<-if(is.null(d$Mdi_prior_sd)).30 else d$Mdi_prior_sd
  Mdm_prior_sd<-if(is.null(d$Mdm_prior_sd)).50 else d$Mdm_prior_sd
  M_active<-as.numeric(!(head(d$years,-1)%in%
    if(is.null(d$fixed_mortality_deviation_years))integer(0) else d$fixed_mortality_deviation_years))
  center_M<-function(z)center_active_last(z,M_active)
  sel<-if(is.null(sel_override)){
    spre<-joint_sel(b$survey_base,selshift_pre,selt_pre,d$sizes)
    spost<-joint_sel(b$survey_base,selshift,selt,d$sizes)
    out<-zero+matrix(0,ny,ns)
    out[d$years<1988,]<-matrix(spre,sum(d$years<1988),ns,byrow=TRUE)
    out[d$years>=1988,]<-matrix(spost,sum(d$years>=1988),ns,byrow=TRUE)
    out
  }else if(is.null(dim(sel_override)))matrix(sel_override,ny,ns,byrow=TRUE) else sel_override
  di<-if(penalize_process)4*tanh(center_M(Mdi)/4) else center_M(Mdi)
  dm<-if(penalize_process)4*tanh(center_M(Mdm)/4) else center_M(Mdm)
  if(!is.null(M_offset)){di<-di+center_M(M_offset$Mi)
    dm<-dm+center_M(M_offset$Mm)}
  dr<-if(penalize_process)5*tanh(center_last(recdev)/5) else center_last(recdev)
  Mi<-M_from_eta(logMi+di)
  Mm<-M_from_eta(logMm+dm)
  rec<-exp(logmeanrec+dr)
  fishery_names<-names(d$fisheries)
  directed_name<-b$stock
  incidental_name<-setdiff(fishery_names,directed_name)
  softplus<-function(x)log(1+exp(x))
  lower_l50<-max(50,min(d$sizes)+10)
  upper_l50<-max(d$sizes)-10
  bounded_l50<-function(x)lower_l50+(upper_l50-lower_l50)/(1+exp(-x))
  build_total<-function(fd,logF,Fdev,l50_raw,log_slope,l50_dev,
    desc_offset_raw,desc_log_slope,desc_dev,dome=FALSE,dev_sd=7.5){
    F_active<-as.numeric(d$years%in%fd$active_years)
    total_dev_active<-as.numeric(d$years%in%fd$total_size_years)
    fdev<-5*tanh(center_active_last(Fdev,F_active)/5)
    mean_l50<-bounded_l50(l50_raw)
    annual_l50<-mean_l50+dev_sd*center_active_last(l50_dev,total_dev_active)
    annual_desc_l50<-annual_l50+softplus(desc_offset_raw+
      dev_sd*center_active_last(desc_dev,total_dev_active))
    slope<-exp(log_slope)/30;desc_slope<-exp(desc_log_slope)/30
    total_sel<-zero+matrix(0,ny,ns)
    for(y in seq_len(ny)){
      ascending<-1/(1+exp(-slope*(d$sizes-annual_l50[y])))
      if(dome){
        descending<-1/(1+exp(desc_slope*(d$sizes-annual_desc_l50[y])))
        raw_total<-ascending*descending
        smooth_peak<-sum(raw_total^20)^(1/20)
        total_sel[y,]<-F_active[y]*raw_total/smooth_peak
      }else{
        # Directed-fishery encounter selectivity is ascending logistic and is
        # normalized at the terminal modeled size so F remains fully selected.
        total_sel[y,]<-F_active[y]*ascending/ascending[ns]
      }
    }
    list(F=F_active*exp(logF+fdev),total_sel=total_sel,F_active=F_active,
      total_dev_active=total_dev_active,mean_l50=mean_l50,
      annual_l50=annual_l50,annual_desc_l50=annual_desc_l50,
      slope=slope,desc_slope=desc_slope)
  }
  directed<-build_total(d$fisheries[[directed_name]],logF,Fdev,
    fish_a,fish_logb,fish_adev,fish_desc_a,fish_desc_logb,fish_bdev,FALSE)
  incidental<-build_total(d$fisheries[[incidental_name]],inc_logF,inc_Fdev,
    inc_fish_a,inc_fish_logb,inc_fish_adev,
    inc_fish_desc_a,inc_fish_desc_logb,inc_fish_bdev,
    !is.null(d$domed_incidental_selectivity)&&d$domed_incidental_selectivity==1)
  directed$ret_dev_active<-as.numeric(d$years%in%
    d$fisheries[[directed_name]]$retained_size_years)
  incidental$ret_dev_active<-as.numeric(d$years%in%
    d$fisheries[[incidental_name]]$retained_size_years)
  direct_ret_dev<-center_active_last(ret_adev,directed$ret_dev_active)
  inc_ret_dev<-center_active_last(inc_ret_adev,incidental$ret_dev_active)
  # Retention midpoints are guaranteed to lie to the right of their matching
  # total-catch midpoint. Incidental retention is additionally shrunk toward
  # directed retention below. Only midpoint deviations vary annually.
  directed$ret_l50<-directed$annual_l50+softplus(ret_a+5*direct_ret_dev)
  incidental$ret_l50<-incidental$annual_l50+softplus(inc_ret_a+3*inc_ret_dev)
  directed$ret_slope<-exp(ret_logb)/30
  incidental$ret_slope<-exp(ret_logb+inc_ret_logb)/30
  directed$retain_sel<-zero+matrix(0,ny,ns)
  incidental$retain_sel<-zero+matrix(0,ny,ns)
  for(y in seq_len(ny)){
    directed$retain_sel[y,]<-directed$F_active[y]/
      (1+exp(-directed$ret_slope*(d$sizes-directed$ret_l50[y])))
    incidental$retain_sel[y,]<-incidental$F_active[y]/
      (1+exp(-incidental$ret_slope*(d$sizes-incidental$ret_l50[y])))
  }
  rp<-recruitment_proportions(rec_logits)
  Ni<-zero+matrix(0,ny,ns);Nm<-zero+matrix(0,ny,ns)
  make_zero<-function()zero+matrix(0,ny,ns)
  caught_directed<-ret_directed<-discard_directed<-make_zero()
  caught_incidental<-ret_incidental<-discard_incidental<-make_zero()
  Ni[1,]<-exp(logNi);Nm[1,]<-exp(logNm)
  for(y in seq_len(ny-1)) {
    ti<-Ni[y,]*exp(-.59*Mi[y]);tm<-Nm[y,]*exp(-.59*Mm[y])
    hd<-directed$F[y]*directed$total_sel[y,]
    hi<-incidental$F[y]*incidental$total_sel[y,]
    h<-hd+hi
    removed_fraction<-1-exp(-h)
    share_d<-hd/(h+1e-12);share_i<-hi/(h+1e-12)
    cid<-ti*removed_fraction*share_d;cmd<-tm*removed_fraction*share_d
    cii<-ti*removed_fraction*share_i;cmi<-tm*removed_fraction*share_i
    caught_directed[y,]<-cid+cmd
    caught_incidental[y,]<-cii+cmi
    ret_directed[y,]<-caught_directed[y,]*directed$retain_sel[y,]
    ret_incidental[y,]<-caught_incidental[y,]*incidental$retain_sel[y,]
    discard_directed[y,]<-caught_directed[y,]-ret_directed[y,]
    discard_incidental[y,]<-caught_incidental[y,]-ret_incidental[y,]
    ti<-ti*exp(-h)+cid*(1-directed$retain_sel[y,])*.70+
      cii*(1-incidental$retain_sel[y,])*.70
    tm<-tm*exp(-h)+cmd*(1-directed$retain_sel[y,])*.70+
      cmi*(1-incidental$retain_sel[y,])*.70
    grown<-as.vector(G%*%ti);grown[1:3]<-grown[1:3]+rec[y]*rp
    # Observation-year y+1 ogive applies to the y -> y+1 transition.
    ni<-grown*(1-d$prop_term_molt[y+1L,]);nm<-tm+grown*d$prop_term_molt[y+1L,]
    Ni[y+1,]<-ni*exp(-.41*Mi[y]);Nm[y+1,]<-nm*exp(-.41*Mm[y])
  }
  vi<-Ni*sel;vm<-Nm*sel
  da<-directed$F_active[seq_along(Fdev)]>0
  dta<-directed$total_dev_active[seq_along(fish_adev)]>0
  dra<-directed$ret_dev_active[seq_along(ret_adev)]>0
  ia<-incidental$F_active[seq_along(inc_Fdev)]>0
  ita<-incidental$total_dev_active[seq_along(inc_fish_adev)]>0
  ira<-incidental$ret_dev_active[seq_along(inc_ret_adev)]>0
  direct_ret_mean<-directed$mean_l50+softplus(ret_a)
  inc_ret_mean<-incidental$mean_l50+softplus(inc_ret_a)
  nll <- -sum(dnorm(Fdev[da],0,.50,log=TRUE))-
    sum(dnorm(fish_adev[dta],0,1,log=TRUE))-
    sum(dnorm(ret_adev[dra],0,1,log=TRUE))-
    sum(dnorm(inc_Fdev[ia],0,.50,log=TRUE))-
    sum(dnorm(inc_fish_adev[ita],0,1,log=TRUE))-
    sum(dnorm(inc_fish_bdev[ita],0,1,log=TRUE))-
    sum(dnorm(inc_ret_adev[ira],0,1,log=TRUE))-
    dnorm(fish_a,0,1.5,log=TRUE)-dnorm(fish_logb,log(4),.7,log=TRUE)-
    dnorm(softplus(ret_a),10,7.5,log=TRUE)-dnorm(ret_logb,log(5),.7,log=TRUE)-
    dnorm(inc_fish_a,0,1.5,log=TRUE)-dnorm(inc_fish_logb,log(4),.7,log=TRUE)-
    dnorm(softplus(inc_fish_desc_a),40,20,log=TRUE)-
    dnorm(inc_fish_desc_logb,log(4),.7,log=TRUE)-
    dnorm(inc_ret_mean,direct_ret_mean,10,log=TRUE)-
    dnorm(inc_ret_logb,0,.25,log=TRUE)-
    dnorm(logMi,M_center,logMi_prior_sd,log=TRUE)-
    dnorm(logMm,M_center,logMm_prior_sd,log=TRUE)-
    sum(dnorm(rec_logits,log(c(.75,.20)/.05),1.5,log=TRUE))
  if(is.null(sel_override)){
    if(!is.null(d$survey_sel_prior)&&!is.null(d$survey_sel_prior_sd)){
      use<-is.finite(d$survey_sel_prior_sd)&d$survey_sel_prior_sd>0
      nll<-nll-sum(dnorm(spre[use],d$survey_sel_prior[use],
        d$survey_sel_prior_sd[use],log=TRUE))-
        sum(dnorm(spost[use],d$survey_sel_prior[use],
        d$survey_sel_prior_sd[use],log=TRUE))
    }else nll<-nll-dnorm(selshift_pre,0,.35,log=TRUE)-dnorm(selt_pre,0,.35,log=TRUE)-
      dnorm(selshift,0,.35,log=TRUE)-dnorm(selt,0,.35,log=TRUE)
  }
  nll<-nll+initial_smooth_nll(logNi,b$init_i)+initial_smooth_nll(logNm,b$init_m)
  if(penalize_process)nll<-nll-sum(dnorm(recdev,0,.6,log=TRUE))-
    sum(dnorm(Mdi[M_active>0],0,Mdi_prior_sd,log=TRUE))-
    sum(dnorm(Mdm[M_active>0],0,Mdm_prior_sd,log=TRUE))
  add_fishery_likelihood<-function(fd,pred_total,pred_ret,pred_discard){
    val<-zero+0
    iy<-match(fd$total_years,d$years)
    for(k in seq_along(iy))val<-val+ln_kernel(fd$total_numbers[k],sum(pred_total[iy[k],]),fd$sigma_total)
    iy<-match(fd$retained_years,d$years)
    for(k in seq_along(iy))val<-val+ln_kernel(fd$retained_numbers[k],sum(pred_ret[iy[k],]),fd$sigma_retained)
    iy<-match(fd$total_size_years,d$years)
    for(k in seq_along(iy))val<-val+ce(fd$total_size[k,],pred_total[iy[k],],fd$total_comp_weight[k])
    iy<-match(fd$retained_size_years,d$years)
    for(k in seq_along(iy))val<-val+ce(fd$retained_size[k,],pred_ret[iy[k],],fd$retained_comp_weight[k])
    iy<-match(fd$discard_size_years,d$years)
    for(k in seq_along(iy))val<-val+ce(fd$discard_size[k,],pred_discard[iy[k],],fd$discard_comp_weight[k])
    val
  }
  nll<-nll+add_fishery_likelihood(d$fisheries[[directed_name]],
    caught_directed,ret_directed,discard_directed)+
    add_fishery_likelihood(d$fisheries[[incidental_name]],
      caught_incidental,ret_incidental,discard_incidental)
  list(Ni=Ni,Nm=Nm,vi=vi,vm=vm,Mi=Mi,Mm=Mm,rec=rec,rec_prop=rp,sel=sel,
    F_directed=directed$F,F_incidental=incidental$F,
    fish_total_directed=directed$total_sel,fish_retain_directed=directed$retain_sel,
    fish_total_incidental=incidental$total_sel,fish_retain_incidental=incidental$retain_sel,
    caught_directed=caught_directed,ret_directed=ret_directed,
    discard_directed=discard_directed,caught_incidental=caught_incidental,
    ret_incidental=ret_incidental,discard_incidental=discard_incidental,
    total_l50_directed=directed$annual_l50,ret_l50_directed=directed$ret_l50,
    total_l50_incidental=incidental$annual_l50,ret_l50_incidental=incidental$ret_l50,
    desc_l50_directed=directed$annual_desc_l50,
    desc_l50_incidental=incidental$annual_desc_l50,
    nll=nll)
}

# Hybrid males use the same growth/recruitment logic, but their terminal-molt,
# selectivity, and mortality curves are constrained by the parental estimates.
hybrid_male_block <- function(b,logN,logmeanrec,rec_logits,recdev,logM,Mdev,ptm,sel,
  tm_blend,sel_blend,penalize_process=TRUE,M_override=NULL,G=b$G,
  rec_offset=NULL) {
  d<-b$data;ny<-length(d$years);ns<-length(d$sizes);zero<-logN[1]*0
  M_active<-as.numeric(!(head(d$years,-1)%in%
    if(is.null(d$fixed_mortality_deviation_years))integer(0) else d$fixed_mortality_deviation_years))
  centered_Mdev<-center_active_last(Mdev,M_active)
  dm<-if(penalize_process)4*tanh(centered_Mdev/4) else centered_Mdev
  dr<-if(penalize_process)5*tanh(center_last(recdev)/5) else center_last(recdev)
  if(!is.null(rec_offset))dr<-dr+center_last(rec_offset)
  M_floor<-if(is.null(d$mortality_floor))0 else d$mortality_floor
  M_link<-if(!is.null(d$mortality_link))d$mortality_link else if(is.null(d$bounded_mortality)||d$bounded_mortality==1)"logit" else "exponential"
  M_center<-switch(M_link,
    logit=qlogis(pmin(pmax((.30-M_floor)/(1-M_floor),1e-6),1-1e-6)),
    hazard=log(-log(pmin(pmax((1-.30)/(1-M_floor),1e-6),1-1e-6))),
    exponential=log(.30))
  M<-if(is.null(M_override))switch(M_link,
    logit=M_floor+(1-M_floor)/(1+exp(-(logM+dm))),
    hazard=-log(1-(M_floor+(1-M_floor)*(1-exp(-exp(logM+dm))))),
    exponential=exp(logM+dm)) else M_override
  rec<-exp(logmeanrec+dr)
  rp<-recruitment_proportions(rec_logits);U<-zero+matrix(0,ny,ns);T<-zero+matrix(0,ny,ns)
  N0<-exp(logN);U[1,]<-N0*(1-ptm);T[1,]<-N0*ptm
  for(y in seq_len(ny-1)) {
    before<-if(isTRUE(d$spring_recruitment)) .59 else 1
    after<-1-before
    grown<-as.vector(G%*%(U[y,]*exp(-before*M[y])));grown[1:3]<-grown[1:3]+rec[y]*rp
    nt<-grown*ptm;U[y+1,]<-(grown-nt)*exp(-after*M[y])
    T[y+1,]<-(T[y,]*exp(-before*M[y])+nt)*exp(-after*M[y])
  }
  N<-U+T;v<-if(is.null(dim(sel)))sweep(N,2,sel,"*") else N*sel
  nll<-if(is.null(M_override))-dnorm(logM,M_center,.10,log=TRUE) else -dnorm(logM,0,1.5,log=TRUE)
  nll<-nll-dnorm(tm_blend,0,1.5,log=TRUE)-
    dnorm(sel_blend,0,1.5,log=TRUE)-
    sum(dnorm(rec_logits,log(c(.75,.20)/.05),1.5,log=TRUE))
  nll<-nll+initial_smooth_nll(logN,b$init_N)
  if(penalize_process)nll<-nll-sum(dnorm(recdev,0,.6,log=TRUE))-
    sum(dnorm(Mdev,0,if(is.null(d$Mdev_prior_sd)).50 else d$Mdev_prior_sd,log=TRUE))
  list(U=U,T=T,N=N,v=v,M=M,rec=rec,rec_prop=rp,sel=sel,ptm=ptm,nll=nll)
}

# =============================================================================
# Parameter initialization
# =============================================================================

# Prefixes identify the six modeled groups:
#   sm/tm/hm = snow, Tanner, and hybrid males
#   sf/tf/hf = snow, Tanner, and hybrid females
J1_parameters <- function(d) {
  M_eta_start<-function(value,b){
    floor<-if(is.null(b$data$mortality_floor))0 else b$data$mortality_floor
    link<-if(!is.null(b$data$mortality_link))b$data$mortality_link else if(!is.null(d$bounded_mortality)&&d$bounded_mortality==0)"exponential" else "logit"
    switch(link,
      logit=qlogis(pmin(pmax((value-floor)/(1-floor),1e-6),1-1e-6)),
      hazard=log(-log(pmin(pmax((1-value)/(1-floor),1e-6),1-1e-6))),
      exponential=log(value))
  }
  logistic_start<-function(curves,sizes){
    active<-rowSums(curves)>0
    q<-apply(curves[active,,drop=FALSE],2,median)
    z<-(sizes-mean(sizes))/30
    fit<-lm(qlogis(pmin(pmax(q,1e-4),1-1e-4))~z)
    c(a=unname(coef(fit)[1]),logb=log(pmax(unname(coef(fit)[2]),.05)))
  }
  fp<-function(b,p)setNames(
    {
    list(log(pmax(b$init_i,1)),log(pmax(b$init_m,1)),log(b$init_rec),log(c(.75,.20)/.05),
      rep(0,length(b$data$years)-1),M_eta_start(b$init_Mi,b),M_eta_start(b$init_Mm,b),
      rep(0,length(b$data$years)-1),rep(0,length(b$data$years)-1),
      log(median(b$init_F)),rep(0,length(b$data$years)-1),0,log(4),
      rep(0,length(b$data$years)-1),rep(0,length(b$data$years)-1),40,log(4),10,log(5),
      rep(0,length(b$data$years)-1),rep(0,length(b$data$years)-1),
      log(.05),rep(0,length(b$data$years)-1),0,log(4),
      rep(0,length(b$data$years)-1),rep(0,length(b$data$years)-1),40,log(4),10,0,
      rep(0,length(b$data$years)-1),rep(0,length(b$data$years)-1),0,0,0,0)},
      paste0(p,c("_logNi","_logNm","_logrec","_rec_logits","_recdev","_logMi","_logMm",
        "_Mdi","_Mdm","_logF","_Fdev","_fish_a","_fish_logb","_fish_adev","_fish_bdev","_fish_desc_a","_fish_desc_logb",
        "_ret_a","_ret_logb","_ret_adev","_ret_bdev",
        "_inc_logF","_inc_Fdev","_inc_fish_a","_inc_fish_logb",
        "_inc_fish_adev","_inc_fish_bdev","_inc_fish_desc_a","_inc_fish_desc_logb","_inc_ret_a","_inc_ret_logb",
        "_inc_ret_adev","_inc_ret_bdev",
        "_selshift_pre","_selt_pre","_selshift","_selt")))
  ff<-function(b,p)setNames(list(log(pmax(b$init_i,1)),log(pmax(b$init_m,1)),
    log(b$init_rec),log(c(.75,.20)/.05),rep(0,length(b$data$years)-1),M_eta_start(b$init_Mi,b),M_eta_start(b$init_Mm,b),
    rep(0,length(b$data$years)-1),rep(0,length(b$data$years)-1),
    unname(b$mat_start["mid"]),unname(b$mat_start["logslope"]),0,0,0,0),
    paste0(p,c("_logNi","_logNm","_logrec","_rec_logits","_recdev","_logMi","_logMm",
      "_Mdi","_Mdm","_mat50","_logmatslope","_selshift_pre","_selt_pre","_selshift","_selt")))
  hm<-d$hybrid_male
  out<-c(fp(d$snow_male,"sm"),fp(d$tanner_male,"tm"),
    list(hm_logN=log(pmax(hm$init_N,1)),hm_logrec=log(hm$init_rec),
      hm_rec_logits=log(c(.75,.20)/.05),
      hm_recdev=rep(0,length(hm$data$years)-1),hm_logM=log(hm$init_M),
      hm_Mdev=rep(0,length(hm$data$years)-1),hm_tm_blend=0,hm_sel_blend=0),
    ff(d$snow_female,"sf"),ff(d$tanner_female,"tf"),ff(d$hybrid_female,"hf"),
    list(hf_sel_blend=0),
    list(
      sm_growth_level=unname(d$snow_male$growth_start["alpha"] + d$snow_male$growth_start["beta"] * d$growth_center_male),
      sm_growth_slope30=unname(d$snow_male$growth_start["beta"] * 30),
      sm_log_growth_sd=log(unname(d$snow_male$growth_start["sd"])),
      tm_growth_level=unname(d$tanner_male$growth_start["alpha"] + d$tanner_male$growth_start["beta"] * d$growth_center_male),
      tm_growth_slope30=unname(d$tanner_male$growth_start["beta"] * 30),
      tm_log_growth_sd=log(unname(d$tanner_male$growth_start["sd"])),
      hm_growth_mean_blend=0,hm_growth_sd_blend=0,
      sf_growth_level=unname(d$snow_female$growth_start["alpha"] + d$snow_female$growth_start["beta"] * d$growth_center_female),
      sf_growth_slope30=unname(d$snow_female$growth_start["beta"] * 30),
      sf_log_growth_sd=log(unname(d$snow_female$growth_start["sd"])),
      tf_growth_level=unname(d$tanner_female$growth_start["alpha"] + d$tanner_female$growth_start["beta"] * d$growth_center_female),
      tf_growth_slope30=unname(d$tanner_female$growth_start["beta"] * 30),
      tf_log_growth_sd=log(unname(d$tanner_female$growth_start["sd"])),
      hf_growth_mean_blend=0,hf_growth_sd_blend=0),
    # Classification parameters are logit(error at 27.5 mm) and log(power).
    # The curve is structurally zero at and above 90 mm.
    list(class_f_intercept=if(!is.null(d$power_misclassification)&&d$power_misclassification==0)-2 else 0,
      class_f_logslope=0,
      class_m_intercept=if(!is.null(d$power_misclassification)&&d$power_misclassification==0)-2 else 0,
      class_m_logslope=0,
      common_rec=rep(0,length(d$years)-1),common_M=rep(0,length(d$years)-1),
      log_rec_sigma=rep(log(.5),6),log_M_sigma=rep(log(.35),11),
      log_rec_loading=rep(log(.1),6),log_M_loading=rep(log(.1),11),
      jr_parent_beta=c(0,0),M_loading=rep(0,11)))
  out$survey_sm<-stacked_survey_start(d$snow_male)
  out$survey_tm<-stacked_survey_start(d$tanner_male)
  out$survey_sf<-stacked_survey_start(d$snow_female)
  out$survey_tf<-stacked_survey_start(d$tanner_female)
  for(group in c("snow_male","tanner_male","snow_female","tanner_female")){
    pr<-d$bsfrf_survey_prior[[group]]
    prefix<-paste0("survey_",sub("male$","m",sub("female$","f",group)))
    target<-pmax(pr$bootstrap_center_ratio,1e-5)
    lower_mid<-(min(pr$sizes)+5)/30;upper_mid<-(max(pr$sizes)-5)/30
    midfun<-function(raw)lower_mid+(upper_mid-lower_mid)*plogis(raw)
    loss<-function(x){shape<-plogis(exp(x[2])*(pr$sizes/30-midfun(x[1])));
      sum((exp(x[3])*shape/tail(shape,1)-target)^2)}
    start_mid<-qlogis((pr$anchor_size/30-lower_mid)/(upper_mid-lower_mid))
    init<-optim(c(start_mid,log(2),log(tail(target,1))),loss,method="BFGS")$par
    out[[paste0(prefix,"_midpoint")]]<-init[1]
    out[[paste0(prefix,"_log_slope")]]<-init[2]
    out[[paste0(prefix,"_log_q")]]<-init[3]
    out[[paste0(prefix,"_midpoint_offset")]]<-rep(0,2)
    out[[paste0(prefix,"_log_slope_offset")]]<-rep(0,2)
    out[[paste0(prefix,"_log_q_offset")]]<-rep(0,2)
  }
  for(group in c("snow_male","tanner_male")){
    pr<-d$bsfrf_survey_prior[[group]]
    prefix<-paste0("survey_",sub("male$","m",group),"_stack")
    target<-pmax(pr$bootstrap_center_ratio,1e-5);z<-pr$sizes/30
    loss<-function(x){m1<-x[1];m2<-m1+exp(x[2]);w<-plogis(x[5])
      shape<-w*plogis(exp(x[3])*(z-m1))+(1-w)*plogis(exp(x[4])*(z-m2))
      sum((exp(x[6])*shape/tail(shape,1)-target)^2)}
    starts<-lapply(c(.25,.5,.75),function(w)c(45/30,log(45/30),log(3),log(2),qlogis(w),log(tail(target,1))))
    fits<-lapply(starts,function(x)optim(x,loss,method="BFGS",control=list(maxit=2000)))
    init<-fits[[which.min(vapply(fits,function(x)x$value,numeric(1)))]]$par
    for(i in seq_along(c("midpoint1","log_gap","log_slope1","log_slope2","logit_weight","log_q"))){
      nm<-c("midpoint1","log_gap","log_slope1","log_slope2","logit_weight","log_q")[i]
      out[[paste0(prefix,"_",nm)]]<-init[i]
      out[[paste0(prefix,"_",nm,"_offset")]]<-rep(0,2)
    }
  }
  out
}

# =============================================================================
# J1 joint objective
# =============================================================================

# This outer function precomputes fixed interpolation matrices. The returned
# inner function is the objective passed to RTMB.
J1_objective <- function(d) {
 Wsm_m<-curve_projection(d$snow_male$data$sizes,d$hybrid_male$data$sizes)
 Wtm_m<-curve_projection(d$tanner_male$data$sizes,d$hybrid_male$data$sizes)
 Wsf_f<-curve_projection(d$snow_female$data$sizes,d$hybrid_female$data$sizes)
 Wtf_f<-curve_projection(d$tanner_female$data$sizes,d$hybrid_female$data$sizes)
 sm_tm_h<-as.vector(Wsm_m%*%apply(d$snow_male$data$prop_term_molt,2,median))
 tm_tm_h<-as.vector(Wtm_m%*%apply(d$tanner_male$data$prop_term_molt,2,median))
 function(
 sm_logNi,sm_logNm,sm_logrec,sm_rec_logits,sm_recdev,sm_logMi,sm_logMm,sm_Mdi,sm_Mdm,sm_logF,sm_Fdev,sm_fish_a,sm_fish_logb,sm_fish_adev,sm_fish_bdev,sm_fish_desc_a,sm_fish_desc_logb,sm_ret_a,sm_ret_logb,sm_ret_adev,sm_ret_bdev,sm_inc_logF,sm_inc_Fdev,sm_inc_fish_a,sm_inc_fish_logb,sm_inc_fish_adev,sm_inc_fish_bdev,sm_inc_fish_desc_a,sm_inc_fish_desc_logb,sm_inc_ret_a,sm_inc_ret_logb,sm_inc_ret_adev,sm_inc_ret_bdev,sm_selshift_pre,sm_selt_pre,sm_selshift,sm_selt,
 tm_logNi,tm_logNm,tm_logrec,tm_rec_logits,tm_recdev,tm_logMi,tm_logMm,tm_Mdi,tm_Mdm,tm_logF,tm_Fdev,tm_fish_a,tm_fish_logb,tm_fish_adev,tm_fish_bdev,tm_fish_desc_a,tm_fish_desc_logb,tm_ret_a,tm_ret_logb,tm_ret_adev,tm_ret_bdev,tm_inc_logF,tm_inc_Fdev,tm_inc_fish_a,tm_inc_fish_logb,tm_inc_fish_adev,tm_inc_fish_bdev,tm_inc_fish_desc_a,tm_inc_fish_desc_logb,tm_inc_ret_a,tm_inc_ret_logb,tm_inc_ret_adev,tm_inc_ret_bdev,tm_selshift_pre,tm_selt_pre,tm_selshift,tm_selt,
 hm_logN,hm_logrec,hm_rec_logits,hm_recdev,hm_logM,hm_Mdev,hm_tm_blend,hm_sel_blend,
 sf_logNi,sf_logNm,sf_logrec,sf_rec_logits,sf_recdev,sf_logMi,sf_logMm,sf_Mdi,sf_Mdm,sf_mat50,sf_logmatslope,sf_selshift_pre,sf_selt_pre,sf_selshift,sf_selt,
 tf_logNi,tf_logNm,tf_logrec,tf_rec_logits,tf_recdev,tf_logMi,tf_logMm,tf_Mdi,tf_Mdm,tf_mat50,tf_logmatslope,tf_selshift_pre,tf_selt_pre,tf_selshift,tf_selt,
 hf_logNi,hf_logNm,hf_logrec,hf_rec_logits,hf_recdev,hf_logMi,hf_logMm,hf_Mdi,hf_Mdm,hf_mat50,hf_logmatslope,hf_selshift_pre,hf_selt_pre,hf_selshift,hf_selt,hf_sel_blend,
 sm_growth_level,sm_growth_slope30,sm_log_growth_sd,
 tm_growth_level,tm_growth_slope30,tm_log_growth_sd,
 hm_growth_mean_blend,hm_growth_sd_blend,
 sf_growth_level,sf_growth_slope30,sf_log_growth_sd,
 tf_growth_level,tf_growth_slope30,tf_log_growth_sd,
 hf_growth_mean_blend,hf_growth_sd_blend,
 class_f_intercept,class_f_logslope,class_m_intercept,class_m_logslope,
 common_rec,common_M,log_rec_sigma,log_M_sigma,log_rec_loading,log_M_loading,
  jr_parent_beta,M_loading,survey_sm,survey_tm,survey_sf,survey_tf,
  survey_snow_m_midpoint,survey_tanner_m_midpoint,
  survey_snow_f_midpoint,survey_tanner_f_midpoint,
  survey_snow_m_log_slope,survey_tanner_m_log_slope,
  survey_snow_f_log_slope,survey_tanner_f_log_slope,
  survey_snow_m_log_q,survey_tanner_m_log_q,
  survey_snow_f_log_q,survey_tanner_f_log_q,
  survey_snow_m_midpoint_offset,survey_tanner_m_midpoint_offset,
  survey_snow_f_midpoint_offset,survey_tanner_f_midpoint_offset,
  survey_snow_m_log_slope_offset,survey_tanner_m_log_slope_offset,
  survey_snow_f_log_slope_offset,survey_tanner_f_log_slope_offset,
  survey_snow_m_log_q_offset,survey_tanner_m_log_q_offset,
  survey_snow_f_log_q_offset,survey_tanner_f_log_q_offset,
  survey_snow_m_stack_midpoint1,survey_snow_m_stack_log_gap,
  survey_snow_m_stack_log_slope1,survey_snow_m_stack_log_slope2,
  survey_snow_m_stack_logit_weight,survey_snow_m_stack_log_q,
  survey_snow_m_stack_midpoint1_offset,survey_snow_m_stack_log_gap_offset,
  survey_snow_m_stack_log_slope1_offset,survey_snow_m_stack_log_slope2_offset,
  survey_snow_m_stack_logit_weight_offset,survey_snow_m_stack_log_q_offset,
  survey_tanner_m_stack_midpoint1,survey_tanner_m_stack_log_gap,
  survey_tanner_m_stack_log_slope1,survey_tanner_m_stack_log_slope2,
  survey_tanner_m_stack_logit_weight,survey_tanner_m_stack_log_q,
  survey_tanner_m_stack_midpoint1_offset,survey_tanner_m_stack_log_gap_offset,
  survey_tanner_m_stack_log_slope1_offset,survey_tanner_m_stack_log_slope2_offset,
  survey_tanner_m_stack_logit_weight_offset,survey_tanner_m_stack_log_q_offset) {
  # J1 normally keeps annual processes stock-specific. This switch retains the
  # machinery needed by nested models that add common recruitment/M drivers.
  shared<-FALSE
  shared_M<-!is.null(d$shared_mortality_driver)&&d$shared_mortality_driver==1
  eff<-function(raw,common,b,ls,ll)exp(ls)*raw+exp(ll)*common[b$year_map[-length(b$year_map)]]
  if(shared){
    re<-list(eff(sm_recdev,common_rec,d$snow_male,log_rec_sigma[1],log_rec_loading[1]),eff(tm_recdev,common_rec,d$tanner_male,log_rec_sigma[2],log_rec_loading[2]),eff(hm_recdev,common_rec,d$hybrid_male,log_rec_sigma[3],log_rec_loading[3]),eff(sf_recdev,common_rec,d$snow_female,log_rec_sigma[4],log_rec_loading[4]),eff(tf_recdev,common_rec,d$tanner_female,log_rec_sigma[5],log_rec_loading[5]),eff(hf_recdev,common_rec,d$hybrid_female,log_rec_sigma[6],log_rec_loading[6]))
    me<-list(eff(sm_Mdi,common_M,d$snow_male,log_M_sigma[1],log_M_loading[1]),eff(sm_Mdm,common_M,d$snow_male,log_M_sigma[2],log_M_loading[2]),eff(tm_Mdi,common_M,d$tanner_male,log_M_sigma[3],log_M_loading[3]),eff(tm_Mdm,common_M,d$tanner_male,log_M_sigma[4],log_M_loading[4]),eff(hm_Mdev,common_M,d$hybrid_male,log_M_sigma[5],log_M_loading[5]),eff(sf_Mdi,common_M,d$snow_female,log_M_sigma[6],log_M_loading[6]),eff(sf_Mdm,common_M,d$snow_female,log_M_sigma[7],log_M_loading[7]),eff(tf_Mdi,common_M,d$tanner_female,log_M_sigma[8],log_M_loading[8]),eff(tf_Mdm,common_M,d$tanner_female,log_M_sigma[9],log_M_loading[9]),eff(hf_Mdi,common_M,d$hybrid_female,log_M_sigma[10],log_M_loading[10]),eff(hf_Mdm,common_M,d$hybrid_female,log_M_sigma[11],log_M_loading[11]))
  }else{re<-list(sm_recdev,tm_recdev,hm_recdev,sf_recdev,tf_recdev,hf_recdev);me<-list(sm_Mdi,sm_Mdm,tm_Mdi,tm_Mdm,hm_Mdev,sf_Mdi,sf_Mdm,tf_Mdi,tf_Mdm,hf_Mdi,hf_Mdm)}
  common_M_active<-as.numeric(!(head(d$years,-1)%in%
    if(is.null(d$fixed_mortality_deviation_years))integer(0) else d$fixed_mortality_deviation_years))
  common_M_centered<-center_active_last(common_M,common_M_active)
  Moff<-function(b,k)M_loading[k]*common_M_centered[b$year_map[-length(b$year_map)]]
  center_block_M<-function(z,b)center_active_last(z,as.numeric(!(head(b$data$years,-1)%in%
    if(is.null(b$data$fixed_mortality_deviation_years))integer(0) else b$data$fixed_mortality_deviation_years)))
  mo<-if(shared_M)list(
    sm=list(Mi=Moff(d$snow_male,1),Mm=Moff(d$snow_male,2)),
    tm=list(Mi=Moff(d$tanner_male,3),Mm=Moff(d$tanner_male,4)),
    hm=Moff(d$hybrid_male,5),
    sf=list(Mi=Moff(d$snow_female,6),Mm=Moff(d$snow_female,7)),
    tf=list(Mi=Moff(d$tanner_female,8),Mm=Moff(d$tanner_female,9)),
    hf=list(Mi=Moff(d$hybrid_female,10),Mm=Moff(d$hybrid_female,11))) else NULL
  # Parental growth regressions are fitted directly to increment observations.
  # Hybrid regression parameters are separate estimated convex blends, bounded
  # componentwise by the simultaneously fitted parental estimates for each sex.
  blend<-function(x)1/(1+exp(-x))
  sm_growth_beta<-sm_growth_slope30/30;tm_growth_beta<-tm_growth_slope30/30
  sf_growth_beta<-sf_growth_slope30/30;tf_growth_beta<-tf_growth_slope30/30
  sm_growth_alpha<-sm_growth_level-sm_growth_beta*d$growth_center_male
  tm_growth_alpha<-tm_growth_level-tm_growth_beta*d$growth_center_male
  sf_growth_alpha<-sf_growth_level-sf_growth_beta*d$growth_center_female
  tf_growth_alpha<-tf_growth_level-tf_growth_beta*d$growth_center_female
  sm_growth_sd<-exp(sm_log_growth_sd);tm_growth_sd<-exp(tm_log_growth_sd)
  sf_growth_sd<-exp(sf_log_growth_sd);tf_growth_sd<-exp(tf_log_growth_sd)
  hmw<-blend(hm_growth_mean_blend);hms<-blend(hm_growth_sd_blend)
  hfw<-blend(hf_growth_mean_blend);hfs<-blend(hf_growth_sd_blend)
  hm_growth_alpha<-hmw*sm_growth_alpha+(1-hmw)*tm_growth_alpha
  hm_growth_beta<-hmw*sm_growth_beta+(1-hmw)*tm_growth_beta
  hm_growth_sd<-hms*sm_growth_sd+(1-hms)*tm_growth_sd
  hf_growth_alpha<-hfw*sf_growth_alpha+(1-hfw)*tf_growth_alpha
  hf_growth_beta<-hfw*sf_growth_beta+(1-hfw)*tf_growth_beta
  hf_growth_sd<-hfs*sf_growth_sd+(1-hfs)*tf_growth_sd
  sm_G<-growth_transition(d$snow_male$data$sizes,sm_growth_alpha,sm_growth_beta,sm_growth_sd)
  tm_G<-growth_transition(d$tanner_male$data$sizes,tm_growth_alpha,tm_growth_beta,tm_growth_sd)
  hm_G<-growth_transition(d$hybrid_male$data$sizes,hm_growth_alpha,hm_growth_beta,hm_growth_sd)
  sf_G<-growth_transition(d$snow_female$data$sizes,sf_growth_alpha,sf_growth_beta,sf_growth_sd)
  tf_G<-growth_transition(d$tanner_female$data$sizes,tf_growth_alpha,tf_growth_beta,tf_growth_sd)
  hf_G<-growth_transition(d$hybrid_female$data$sizes,hf_growth_alpha,hf_growth_beta,hf_growth_sd)
  # Evaluate all six stock-by-sex population processes independently before
  # they enter the shared classification and observation layers.
  legacy_survey<-!is.null(d$legacy_two_fishery_survey)&&d$legacy_two_fishery_survey==1
  curve_sm<-stacked_logistic_survey_curves(d$male_sizes,survey_snow_m_stack_midpoint1,survey_snow_m_stack_log_gap,survey_snow_m_stack_log_slope1,survey_snow_m_stack_log_slope2,survey_snow_m_stack_logit_weight,survey_snow_m_stack_log_q,survey_snow_m_stack_midpoint1_offset,survey_snow_m_stack_log_gap_offset,survey_snow_m_stack_log_slope1_offset,survey_snow_m_stack_log_slope2_offset,survey_snow_m_stack_logit_weight_offset,survey_snow_m_stack_log_q_offset,TRUE)
  curve_tm<-stacked_logistic_survey_curves(d$male_sizes,survey_tanner_m_stack_midpoint1,survey_tanner_m_stack_log_gap,survey_tanner_m_stack_log_slope1,survey_tanner_m_stack_log_slope2,survey_tanner_m_stack_logit_weight,survey_tanner_m_stack_log_q,survey_tanner_m_stack_midpoint1_offset,survey_tanner_m_stack_log_gap_offset,survey_tanner_m_stack_log_slope1_offset,survey_tanner_m_stack_log_slope2_offset,survey_tanner_m_stack_logit_weight_offset,survey_tanner_m_stack_log_q_offset,FALSE)
  curve_sf<-logistic_survey_curves(d$female_sizes,survey_snow_f_midpoint,survey_snow_f_log_slope,survey_snow_f_log_q,survey_snow_f_midpoint_offset,survey_snow_f_log_slope_offset,survey_snow_f_log_q_offset,TRUE)
  curve_tf<-logistic_survey_curves(d$female_sizes,survey_tanner_f_midpoint,survey_tanner_f_log_slope,survey_tanner_f_log_q,survey_tanner_f_midpoint_offset,survey_tanner_f_log_slope_offset,survey_tanner_f_log_q_offset)
  if(!is.null(d$survey_q_fixed_at_terminal_size)&&d$survey_q_fixed_at_terminal_size){
    normalize_at_anchor<-function(curve,sizes,prior){
      ii<-match(prior$anchor_size,sizes)
      sweep(curve,2,curve[ii,],"/")
    }
    curve_sm<-normalize_at_anchor(curve_sm,d$male_sizes,d$bsfrf_survey_prior$snow_male)
    curve_tm<-normalize_at_anchor(curve_tm,d$male_sizes,d$bsfrf_survey_prior$tanner_male)
    curve_sf<-normalize_at_anchor(curve_sf,d$female_sizes,d$bsfrf_survey_prior$snow_female)
    curve_tf<-normalize_at_anchor(curve_tf,d$female_sizes,d$bsfrf_survey_prior$tanner_female)
  }
  select_curve<-function(block,common_sizes,curves){idx<-match(block$data$sizes,common_sizes);t(curves[idx,survey_era_index(block$data$years),drop=FALSE])}
  if(legacy_survey){
    sel_sm<-stacked_survey_matrix(d$snow_male,survey_sm)
    sel_tm<-stacked_survey_matrix(d$tanner_male,survey_tm)
    sel_sf<-stacked_survey_matrix(d$snow_female,survey_sf)
    sel_tf<-stacked_survey_matrix(d$tanner_female,survey_tf)
  }else{
    sel_sm<-select_curve(d$snow_male,d$male_sizes,curve_sm)
    sel_tm<-select_curve(d$tanner_male,d$male_sizes,curve_tm)
    sel_sf<-select_curve(d$snow_female,d$female_sizes,curve_sf)
    sel_tf<-select_curve(d$tanner_female,d$female_sizes,curve_tf)
  }
  sm<-parent_male_block(d$snow_male,sm_logNi,sm_logNm,sm_logrec,sm_rec_logits,re[[1]],sm_logMi,sm_logMm,me[[1]],me[[2]],sm_logF,sm_Fdev,sm_fish_a,sm_fish_logb,sm_fish_adev,sm_fish_bdev,sm_fish_desc_a,sm_fish_desc_logb,sm_ret_a,sm_ret_logb,sm_ret_adev,sm_ret_bdev,sm_inc_logF,sm_inc_Fdev,sm_inc_fish_a,sm_inc_fish_logb,sm_inc_fish_adev,sm_inc_fish_bdev,sm_inc_fish_desc_a,sm_inc_fish_desc_logb,sm_inc_ret_a,sm_inc_ret_logb,sm_inc_ret_adev,sm_inc_ret_bdev,sm_selshift_pre,sm_selt_pre,sm_selshift,sm_selt,!shared,G=sm_G,M_offset=if(shared_M)mo$sm else NULL,sel_override=sel_sm)
  tm<-parent_male_block(d$tanner_male,tm_logNi,tm_logNm,tm_logrec,tm_rec_logits,re[[2]],tm_logMi,tm_logMm,me[[3]],me[[4]],tm_logF,tm_Fdev,tm_fish_a,tm_fish_logb,tm_fish_adev,tm_fish_bdev,tm_fish_desc_a,tm_fish_desc_logb,tm_ret_a,tm_ret_logb,tm_ret_adev,tm_ret_bdev,tm_inc_logF,tm_inc_Fdev,tm_inc_fish_a,tm_inc_fish_logb,tm_inc_fish_adev,tm_inc_fish_bdev,tm_inc_fish_desc_a,tm_inc_fish_desc_logb,tm_inc_ret_a,tm_inc_ret_logb,tm_inc_ret_adev,tm_inc_ret_bdev,tm_selshift_pre,tm_selt_pre,tm_selshift,tm_selt,!shared,G=tm_G,M_offset=if(shared_M)mo$tm else NULL,sel_override=sel_tm)
  sf<-female_block(d$snow_female,sf_logNi,sf_logNm,sf_logrec,sf_rec_logits,re[[4]],sf_logMi,sf_logMm,me[[6]],me[[7]],sf_mat50,sf_logmatslope,sf_selshift_pre,sf_selt_pre,sf_selshift,sf_selt,!shared,G=sf_G,M_offset=if(shared_M)mo$sf else NULL,sel_override=sel_sf)
  tf<-female_block(d$tanner_female,tf_logNi,tf_logNm,tf_logrec,tf_rec_logits,re[[5]],tf_logMi,tf_logMm,me[[8]],me[[9]],tf_mat50,tf_logmatslope,tf_selshift_pre,tf_selt_pre,tf_selshift,tf_selt,!shared,G=tf_G,M_offset=if(shared_M)mo$tf else NULL,sel_override=sel_tf)
  # JR links hybrid recruitment to the lag-5 snow-female/Tanner-male pathway.
  # Mature snow females are multiplied by mate sufficiency and the Tanner share
  # among >=95-mm mature males. Fixed J1-based centering/scaling preserves a
  # meaningful nested comparison while the covariate itself changes with the
  # simultaneously estimated latent parental states. A quadratic basis is the
  # parsimonious joint analogue of the selected exploratory smooth.
  jr_on<-!is.null(d$parental_mature_hybrid_recruitment)&&
    d$parental_mature_hybrid_recruitment==1
  jr_driver<-function(hybrid){
    out<-sm_logNi[1]*0+rep(0,length(hybrid$data$years)-1)
    for(y in seq_along(out)){
      py<-hybrid$data$years[y]-d$jr_lag
      isf<-match(py,d$snow_female$data$years)
      itm<-match(py,d$tanner_male$data$years)
      ism<-match(py,d$snow_male$data$years)
      if(!is.na(isf)&&!is.na(itm)&&!is.na(ism)){
        snow_female<-sum(sf$Nm[isf,])
        tanner_male<-sum(tm$Nm[itm,d$tanner_male$data$sizes>=d$jr_large_male_min])
        snow_male<-sum(sm$Nm[ism,d$snow_male$data$sizes>=d$jr_large_male_min])
        male_total<-tanner_male+snow_male
        sufficiency<-1-exp(-male_total/d$jr_male_capacity_K)
        pathway<-snow_female*sufficiency*(tanner_male+1)/(male_total+2)
        out[y]<-(log(pathway+1)-d$jr_driver_center)/d$jr_driver_scale
      }
    }
    out
  }
  if(jr_on){
    hm_jr_driver<-jr_driver(d$hybrid_male)
    hf_jr_driver<-jr_driver(d$hybrid_female)
  }else{
    hm_jr_driver<-sm_logNi[1]*0+rep(0,length(d$hybrid_male$data$years)-1)
    hf_jr_driver<-sm_logNi[1]*0+rep(0,length(d$hybrid_female$data$years)-1)
  }
  jr_effect<-function(x)jr_parent_beta[1]*x+jr_parent_beta[2]*x^2
  hm_jr_effect<-if(jr_on)jr_effect(hm_jr_driver) else NULL
  hf_jr_effect<-if(jr_on)jr_effect(hf_jr_driver) else NULL
  # Logistic blend weights constrain hybrid terminal molt and selectivity to
  # lie between the corresponding snow and Tanner curves.
  wtm<-1/(1+exp(-hm_tm_blend));wms<-1/(1+exp(-hm_sel_blend));wfs<-1/(1+exp(-hf_sel_blend))
  hm_ptm_parent_snow<-sm_tm_h;hm_ptm_parent_tanner<-tm_tm_h
  sel_on_hybrid<-function(x,parent,hybrid,W){
    yy<-pmin(pmax(hybrid$data$years,min(parent$data$years)),max(parent$data$years))
    x[match(yy,parent$data$years),,drop=FALSE]%*%t(W)
  }
  # Evaluate the continuous parental curves on hybrid sizes, rather than
  # holding the last parental bin constant. This preserves the asymptote of 1.
  if(legacy_survey){
    hm_sel_parent_snow<-stacked_survey_matrix(d$hybrid_male,survey_sm)
    hm_sel_parent_tanner<-stacked_survey_matrix(d$hybrid_male,survey_tm)
    hf_sel_parent_snow<-stacked_survey_matrix(d$hybrid_female,survey_sf)
    hf_sel_parent_tanner<-stacked_survey_matrix(d$hybrid_female,survey_tf)
  }else{
    hm_sel_parent_snow<-select_curve(d$hybrid_male,d$male_sizes,curve_sm)
    hm_sel_parent_tanner<-select_curve(d$hybrid_male,d$male_sizes,curve_tm)
    hf_sel_parent_snow<-select_curve(d$hybrid_female,d$female_sizes,curve_sf)
    hf_sel_parent_tanner<-select_curve(d$hybrid_female,d$female_sizes,curve_tf)
  }
  hm_ptm<-wtm*hm_ptm_parent_snow+(1-wtm)*hm_ptm_parent_tanner
  hm_sel_constrained<-wms*hm_sel_parent_snow+(1-wms)*hm_sel_parent_tanner
  hf_sel_constrained<-wfs*hf_sel_parent_snow+(1-wfs)*hf_sel_parent_tanner
  # Hybrid mortality is a differentiable convex blend of annual parental
  # mortality. This guarantees it remains inside the parental envelope without
  # introducing the nondifferentiability of pmin/pmax clipping. Female bounds
  # match maturity state. The current hybrid-male model has one M for both
  # states, so its two stock endpoints average parental immature and mature M.
  annual_on_hybrid<-function(x,parent,hybrid){
    # Constant endpoint extrapolation covers the hybrid terminal reporting year
    # when a parental series ends one year earlier.
    yy<-pmin(pmax(hybrid$data$years,min(parent$data$years)),max(parent$data$years))
    x[match(yy,parent$data$years)]
  }
  hm_w<-1/(1+exp(-(hm_logM+(if(!shared)4*tanh(center_block_M(me[[5]],d$hybrid_male)/4) else center_block_M(me[[5]],d$hybrid_male))+
    (if(shared_M)center_block_M(mo$hm,d$hybrid_male) else 0))))
  hm_snow_endpoint<-.5*(annual_on_hybrid(sm$Mi,d$snow_male,d$hybrid_male)+annual_on_hybrid(sm$Mm,d$snow_male,d$hybrid_male))
  hm_tanner_endpoint<-.5*(annual_on_hybrid(tm$Mi,d$tanner_male,d$hybrid_male)+annual_on_hybrid(tm$Mm,d$tanner_male,d$hybrid_male))
  hm_M_bound<-hm_w*hm_snow_endpoint+(1-hm_w)*hm_tanner_endpoint
  hf_wi<-1/(1+exp(-(hf_logMi+(if(!shared)4*tanh(center_block_M(me[[10]],d$hybrid_female)/4) else center_block_M(me[[10]],d$hybrid_female))+
    (if(shared_M)center_block_M(mo$hf$Mi,d$hybrid_female) else 0))))
  hf_wm<-1/(1+exp(-(hf_logMm+(if(!shared)4*tanh(center_block_M(me[[11]],d$hybrid_female)/4) else center_block_M(me[[11]],d$hybrid_female))+
    (if(shared_M)center_block_M(mo$hf$Mm,d$hybrid_female) else 0))))
  hf_M_bound<-list(
    Mi=hf_wi*annual_on_hybrid(sf$Mi,d$snow_female,d$hybrid_female)+(1-hf_wi)*annual_on_hybrid(tf$Mi,d$tanner_female,d$hybrid_female),
    Mm=hf_wm*annual_on_hybrid(sf$Mm,d$snow_female,d$hybrid_female)+(1-hf_wm)*annual_on_hybrid(tf$Mm,d$tanner_female,d$hybrid_female))
  hm<-hybrid_male_block(d$hybrid_male,hm_logN,hm_logrec,hm_rec_logits,re[[3]],hm_logM,me[[5]],hm_ptm,hm_sel_constrained,hm_tm_blend,hm_sel_blend,!shared,M_override=hm_M_bound,G=hm_G,rec_offset=hm_jr_effect)
  hf<-female_block(d$hybrid_female,hf_logNi,hf_logNm,hf_logrec,hf_rec_logits,re[[6]],hf_logMi,hf_logMm,me[[10]],me[[11]],hf_mat50,hf_logmatslope,hf_selshift_pre,hf_selt_pre,hf_selshift,hf_selt,!shared,sel_override=hf_sel_constrained,
    logMi_sd=.08,logMm_sd=.12,Mdi_sd=.50,Mdm_sd=.50,M_override=hf_M_bound,G=hf_G,rec_offset=hf_jr_effect)
  hf$nll<-hf$nll-dnorm(hf_sel_blend,0,1.5,log=TRUE)
  if(legacy_survey){
    shared_survey_prior<-
      stacked_survey_prior_nll(d$snow_male,survey_sm)+
      stacked_survey_prior_nll(d$tanner_male,survey_tm)+
      stacked_survey_prior_nll(d$snow_female,survey_sf)+
      stacked_survey_prior_nll(d$tanner_female,survey_tf)
    # Historical female logistic curves are offsets from the modern curve in
    # interpretable coordinates. The midpoint parameter is expressed in units
    # of 30 mm; slope is already on the log scale. These penalties directly
    # partially pool each species' early curve toward its current curve while
    # retaining the existing curve-level BSFRF priors for both eras.
    if(!is.null(d$female_survey_era_midpoint_sd_mm)&&
       !is.null(d$female_survey_era_log_slope_sd)){
      female_era_offset_nll<-function(theta){
        -dnorm(30*(theta[1,1]-theta[1,2]),0,
          d$female_survey_era_midpoint_sd_mm,log=TRUE)-
        dnorm(theta[3,1]-theta[3,2],0,
          d$female_survey_era_log_slope_sd,log=TRUE)
      }
      shared_survey_prior<-shared_survey_prior+
        female_era_offset_nll(survey_sf)+
        female_era_offset_nll(survey_tf)
    }
  }else shared_survey_prior<-
      logistic_survey_prior_nll(d$bsfrf_survey_prior$snow_male,curve_sm,TRUE)+
      logistic_survey_prior_nll(d$bsfrf_survey_prior$tanner_male,curve_tm)+
      logistic_survey_prior_nll(d$bsfrf_survey_prior$snow_female,curve_sf,TRUE)+
      logistic_survey_prior_nll(d$bsfrf_survey_prior$tanner_female,curve_tf)
  # Weak internal mixture regularization prevents a component with zero weight
  # from leaving its gap and slopes unidentified. Curve-level calibration is
  # still supplied by the BSFRF priors above.
  if(!legacy_survey)for(x in list(
    c(survey_snow_m_stack_log_gap,survey_snow_m_stack_log_slope1,
      survey_snow_m_stack_log_slope2,survey_snow_m_stack_logit_weight),
    c(survey_tanner_m_stack_log_gap,survey_tanner_m_stack_log_slope1,
      survey_tanner_m_stack_log_slope2,survey_tanner_m_stack_logit_weight))){
    shared_survey_prior<-shared_survey_prior-
      dnorm(x[1],log(1.5),1,log=TRUE)-
      dnorm(x[2],log(3),1,log=TRUE)-dnorm(x[3],log(2),1,log=TRUE)-
      dnorm(x[4],0,1.5,log=TRUE)
  }
  # Place every prediction on common joint year/size grids. The zero inherits
  # the automatic-differentiation type from a model parameter.
  zero<-sm_logNi[1]*0;eps<-1e-12;ny<-length(d$years)
  align2<-function(x,b,sizes){z<-zero+matrix(0,ny,length(sizes));z[b$year_map,b$size_map]<-x;z}
  align3<-function(i,m,b,sizes){z<-zero+array(0,c(ny,length(sizes),2));z[b$year_map,b$size_map,1]<-i;z[b$year_map,b$size_map,2]<-m;z}
  SM<-align3(sm$vi,sm$vm,d$snow_male,d$male_sizes);TM<-align3(tm$vi,tm$vm,d$tanner_male,d$male_sizes)
  HM<-align2(hm$v,d$hybrid_male,d$male_sizes)
  SF<-align3(sf$vi,sf$vm,d$snow_female,d$female_sizes);TF<-align3(tf$vi,tf$vm,d$tanner_female,d$female_sizes)
  HF<-align3(hf$vi,hf$vm,d$hybrid_female,d$female_sizes)
  # Probability that a true hybrid is recorded under a parental label.  Error
  # is estimated at the smallest modeled bin and declines as a power function,
  # reaching exactly zero at 90 mm.  Large hybrids are therefore treated as
  # correctly identified by construction, rather than through a soft prior.
  classification_curve<-function(sizes,logit_p_small,log_power){
    active<-as.numeric(sizes<90)
    # A positive inactive-bin base avoids the undefined AD derivative of
    # 0^gamma; multiplying by the fixed mask still gives an exact zero.
    base<-ifelse(sizes<90,(90-sizes)/(90-27.5),1)
    active*d$classification_on/(1+exp(-logit_p_small))*base^exp(log_power)
  }
  power_misclass<-is.null(d$power_misclassification)||d$power_misclassification==1
  if(power_misclass){
    pm<-classification_curve(d$male_sizes,class_m_intercept,class_m_logslope)
    if(!is.null(d$classification_shared_sex)&&d$classification_shared_sex==1)
      pf<-classification_curve(d$female_sizes,class_m_intercept,class_m_logslope) else
      pf<-classification_curve(d$female_sizes,class_f_intercept,class_f_logslope)
  }else{
    pm<-d$classification_on/(1+exp(-(class_m_intercept-
      exp(class_m_logslope)*(d$male_sizes-60)/30)))
    if(!is.null(d$classification_shared_sex)&&d$classification_shared_sex==1)
      pf<-d$classification_on/(1+exp(-(class_m_intercept-
        exp(class_m_logslope)*(d$female_sizes-60)/30))) else
      pf<-d$classification_on/(1+exp(-(class_f_intercept-
        exp(class_f_logslope)*(d$female_sizes-60)/30)))
  }
  # Misclassified hybrids are removed from the observed-hybrid prediction and
  # allocated between parental labels in proportion to vulnerable abundance.
  obsSM<-SM;obsTM<-TM;obsHM<-HM
  Pm<-SM[,,1]+SM[,,2]+TM[,,1]+TM[,,2]
  leakm<-sweep(HM,2,pm,"*")
  obsHM<-HM-leakm
  for(k in 1:2){obsSM[,,k]<-SM[,,k]+leakm*SM[,,k]/(Pm+eps);obsTM[,,k]<-TM[,,k]+leakm*TM[,,k]/(Pm+eps)}
  obsSF<-SF;obsTF<-TF;obsHF<-HF
  for(k in 1:2){P<-SF[,,k]+TF[,,k];L<-sweep(HF[,,k],2,pf,"*");obsHF[,,k]<-HF[,,k]-L
    obsSF[,,k]<-SF[,,k]+L*SF[,,k]/(P+eps);obsTF[,,k]<-TF[,,k]+L*TF[,,k]/(P+eps)}
  # Add abundance and size-composition likelihoods after classification error
  # has transformed the vulnerable population predictions.
  nll<-sm$nll+tm$nll+hm$nll+sf$nll+tf$nll+hf$nll
  abundance_comp<-zero+0;size_comp<-zero+0
  male_like<-function(b,pred,hybrid=FALSE){a<-zero+0;s<-zero+0;dd<-b$data
    if(hybrid){for(y in which(dd$survey_observed)){yy<-b$year_map[y];pr<-sum(pred[yy,]);oo<-dd$survey_total[y];cv<-dd$survey_cv[y]
      if(is.finite(oo)&&oo>0&&is.finite(cv)&&cv>0)a<-a+ln_kernel(oo,pr,cv)
      if(sum(dd$survey_at_size[y,])>0)s<-s+ce(dd$survey_at_size[y,],pred[yy,b$size_map],b$survey_comp_weight[y])}} else for(k in 1:2)for(y in seq_along(dd$years)){
      yy<-b$year_map[y];oo<-if(k==1)dd$imm_n_obs[y] else dd$mat_n_obs[y];cv<-if(k==1)dd$sigma_numbers_imm[y] else dd$sigma_numbers_mat[y]
      if(dd$years[y]!=2020&&is.finite(oo)&&oo>0&&is.finite(cv)&&cv>0)a<-a+ln_kernel(oo,sum(pred[yy,b$size_map,k]),cv)
      os<-if(k==1)dd$imm_n_size_obs[y,] else dd$mat_n_size_obs[y,];w<-if(k==1)b$survey_comp_weight_imm[y] else b$survey_comp_weight_mat[y]
      if(sum(os)>0)s<-s+ce(os,pred[yy,b$size_map,k],w)};c(a,s)}
  female_like<-function(b,pred){a<-zero+0;s<-zero+0;dd<-b$data
    for(y in which(dd$survey_observed)){yy<-b$year_map[y]
      for(k in 1:2){oo<-dd$survey_total[y,k];cv<-dd$survey_cv[y,k]
        if(is.finite(oo)&&oo>0&&is.finite(cv)&&cv>0)a<-a+ln_kernel(oo,sum(pred[yy,b$size_map,k]),cv)}
      oi<-dd$survey_at_size[y,,1];om<-dd$survey_at_size[y,,2]
      pi<-pred[yy,b$size_map,1];pmat<-pred[yy,b$size_map,2]
      if(sum(oi)>0)s<-s+ce(oi,pi,b$survey_comp_weight_imm[y])
      if(sum(om)>0)s<-s+ce(om,pmat,b$survey_comp_weight_mat[y])};c(a,s)}
  q<-male_like(d$snow_male,obsSM);abundance_comp<-abundance_comp+q[1];size_comp<-size_comp+q[2]
  q<-male_like(d$tanner_male,obsTM);abundance_comp<-abundance_comp+q[1];size_comp<-size_comp+q[2]
  q<-male_like(d$hybrid_male,obsHM,TRUE);abundance_comp<-abundance_comp+q[1];size_comp<-size_comp+q[2]
  q<-female_like(d$snow_female,obsSF);abundance_comp<-abundance_comp+q[1];size_comp<-size_comp+q[2]
  q<-female_like(d$tanner_female,obsTF);abundance_comp<-abundance_comp+q[1];size_comp<-size_comp+q[2]
  q<-female_like(d$hybrid_female,obsHF);abundance_comp<-abundance_comp+q[1];size_comp<-size_comp+q[2]
  # No observation-error prior: identification is supplied by the likelihood
  # below 90 mm and by the structural zero-error rule at larger sizes.
  class_prior<-if(power_misclass)zero+0 else
    -dnorm(class_m_intercept,-2,1.5,log=TRUE)-dnorm(class_m_logslope,0,1,log=TRUE)
  if(!power_misclass&&(is.null(d$classification_shared_sex)||d$classification_shared_sex!=1))
    class_prior<-class_prior-dnorm(class_f_intercept,-2,1.5,log=TRUE)-
      dnorm(class_f_logslope,0,1,log=TRUE)
  driver_prior<-zero+0
  if(shared){
    driver_prior<- -sum(dnorm(common_rec,0,1,log=TRUE))-sum(dnorm(common_M,0,1,log=TRUE))-
      sum(dnorm(sm_recdev,0,1,log=TRUE))-sum(dnorm(tm_recdev,0,1,log=TRUE))-sum(dnorm(hm_recdev,0,1,log=TRUE))-sum(dnorm(sf_recdev,0,1,log=TRUE))-sum(dnorm(tf_recdev,0,1,log=TRUE))-sum(dnorm(hf_recdev,0,1,log=TRUE))-
      sum(dnorm(sm_Mdi,0,1,log=TRUE))-sum(dnorm(sm_Mdm,0,1,log=TRUE))-sum(dnorm(tm_Mdi,0,1,log=TRUE))-sum(dnorm(tm_Mdm,0,1,log=TRUE))-sum(dnorm(hm_Mdev,0,1,log=TRUE))-sum(dnorm(sf_Mdi,0,1,log=TRUE))-sum(dnorm(sf_Mdm,0,1,log=TRUE))-sum(dnorm(tf_Mdi,0,1,log=TRUE))-sum(dnorm(tf_Mdm,0,1,log=TRUE))-sum(dnorm(hf_Mdi,0,1,log=TRUE))-sum(dnorm(hf_Mdm,0,1,log=TRUE))-
      sum(dnorm(log_rec_sigma,log(.5),.7,log=TRUE))-sum(dnorm(log_M_sigma,log(.35),.35,log=TRUE))-sum(dnorm(log_rec_loading,log(.1),1,log=TRUE))-sum(dnorm(log_M_loading,log(.1),1,log=TRUE))
  }
  growth_nll<-
    growth_data_nll(d$snow_male$data$growth_premolt,d$snow_male$data$growth_postmolt,sm_growth_alpha,sm_growth_beta,sm_growth_sd)+
    growth_data_nll(d$tanner_male$data$growth_premolt,d$tanner_male$data$growth_postmolt,tm_growth_alpha,tm_growth_beta,tm_growth_sd)+
    growth_data_nll(d$snow_female$data$growth$premolt,d$snow_female$data$growth$postmolt,sf_growth_alpha,sf_growth_beta,sf_growth_sd)+
    growth_data_nll(d$tanner_female$data$growth$premolt,d$tanner_female$data$growth$postmolt,tf_growth_alpha,tf_growth_beta,tf_growth_sd)
  hybrid_growth_prior<--sum(dnorm(c(hm_growth_mean_blend,hm_growth_sd_blend,
    hf_growth_mean_blend,hf_growth_sd_blend),0,1.5,log=TRUE))
  jr_prior<-if(jr_on).5*sum(jr_parent_beta^2) else zero+0
  shared_M_prior<-if(shared_M).5*sum(common_M^2)+.5*sum((M_loading/.35)^2) else zero+0
  components<-c(sm$nll,tm$nll,hm$nll,sf$nll,tf$nll,hf$nll,abundance_comp,size_comp,class_prior,driver_prior,growth_nll,hybrid_growth_prior,jr_prior,shared_M_prior,shared_survey_prior)

  # REPORT values are saved for plotting and diagnostics. ADREPORT values also
  # receive delta-method standard errors from the fitted Hessian.
  sm_Ni<-sm$Ni;sm_Nm<-sm$Nm;tm_Ni<-tm$Ni;tm_Nm<-tm$Nm
  hm_U<-hm$U;hm_T<-hm$T;hm_N<-hm$N
  sf_Ni<-sf$Ni;sf_Nm<-sf$Nm;tf_Ni<-tf$Ni;tf_Nm<-tf$Nm;hf_Ni<-hf$Ni;hf_Nm<-hf$Nm
  REPORT(sm_Ni);REPORT(sm_Nm);REPORT(tm_Ni);REPORT(tm_Nm)
  REPORT(hm_U);REPORT(hm_T);REPORT(hm_N)
  REPORT(sf_Ni);REPORT(sf_Nm);REPORT(tf_Ni);REPORT(tf_Nm);REPORT(hf_Ni);REPORT(hf_Nm)
  REPORT(obsSM);REPORT(obsTM);REPORT(obsHM);REPORT(obsSF);REPORT(obsTF);REPORT(obsHF)
  REPORT(SM);REPORT(TM);REPORT(HM);REPORT(SF);REPORT(TF);REPORT(HF);REPORT(pf);REPORT(pm)
  sm_Mi<-sm$Mi;sm_Mm<-sm$Mm;sm_rec<-sm$rec;sm_sel<-sm$sel
  tm_Mi<-tm$Mi;tm_Mm<-tm$Mm;tm_rec<-tm$rec;tm_sel<-tm$sel
  sm_F_directed<-sm$F_directed;sm_F_incidental<-sm$F_incidental
  tm_F_directed<-tm$F_directed;tm_F_incidental<-tm$F_incidental
  sm_fish_total_directed<-sm$fish_total_directed
  sm_fish_retain_directed<-sm$fish_retain_directed
  sm_fish_total_incidental<-sm$fish_total_incidental
  sm_fish_retain_incidental<-sm$fish_retain_incidental
  tm_fish_total_directed<-tm$fish_total_directed
  tm_fish_retain_directed<-tm$fish_retain_directed
  tm_fish_total_incidental<-tm$fish_total_incidental
  tm_fish_retain_incidental<-tm$fish_retain_incidental
  sm_caught_directed<-sm$caught_directed;sm_ret_directed<-sm$ret_directed
  sm_discard_directed<-sm$discard_directed
  sm_caught_incidental<-sm$caught_incidental;sm_ret_incidental<-sm$ret_incidental
  sm_discard_incidental<-sm$discard_incidental
  tm_caught_directed<-tm$caught_directed;tm_ret_directed<-tm$ret_directed
  tm_discard_directed<-tm$discard_directed
  tm_caught_incidental<-tm$caught_incidental;tm_ret_incidental<-tm$ret_incidental
  tm_discard_incidental<-tm$discard_incidental
  sm_total_l50_directed<-sm$total_l50_directed;sm_ret_l50_directed<-sm$ret_l50_directed
  sm_total_l50_incidental<-sm$total_l50_incidental;sm_ret_l50_incidental<-sm$ret_l50_incidental
  tm_total_l50_directed<-tm$total_l50_directed;tm_ret_l50_directed<-tm$ret_l50_directed
  tm_total_l50_incidental<-tm$total_l50_incidental;tm_ret_l50_incidental<-tm$ret_l50_incidental
  hm_M<-hm$M;hm_rec<-hm$rec;hm_sel<-hm$sel;hm_ptm<-hm$ptm
  sf_Mi<-sf$Mi;sf_Mm<-sf$Mm;sf_rec<-sf$rec;sf_sel<-sf$sel;sf_pm<-sf$pm
  tf_Mi<-tf$Mi;tf_Mm<-tf$Mm;tf_rec<-tf$rec;tf_sel<-tf$sel;tf_pm<-tf$pm
  hf_Mi<-hf$Mi;hf_Mm<-hf$Mm;hf_rec<-hf$rec;hf_sel<-hf$sel;hf_pm<-hf$pm
  REPORT(sm_Mi);REPORT(sm_Mm);REPORT(sm_rec);REPORT(sm_sel)
  REPORT(tm_Mi);REPORT(tm_Mm);REPORT(tm_rec);REPORT(tm_sel)
  REPORT(sm_F_directed);REPORT(sm_F_incidental)
  REPORT(tm_F_directed);REPORT(tm_F_incidental)
  REPORT(sm_fish_total_directed);REPORT(sm_fish_retain_directed)
  REPORT(sm_fish_total_incidental);REPORT(sm_fish_retain_incidental)
  REPORT(tm_fish_total_directed);REPORT(tm_fish_retain_directed)
  REPORT(tm_fish_total_incidental);REPORT(tm_fish_retain_incidental)
  REPORT(sm_caught_directed);REPORT(sm_ret_directed);REPORT(sm_discard_directed)
  REPORT(sm_caught_incidental);REPORT(sm_ret_incidental);REPORT(sm_discard_incidental)
  REPORT(tm_caught_directed);REPORT(tm_ret_directed);REPORT(tm_discard_directed)
  REPORT(tm_caught_incidental);REPORT(tm_ret_incidental);REPORT(tm_discard_incidental)
  REPORT(sm_total_l50_directed);REPORT(sm_ret_l50_directed)
  REPORT(sm_total_l50_incidental);REPORT(sm_ret_l50_incidental)
  REPORT(tm_total_l50_directed);REPORT(tm_ret_l50_directed)
  REPORT(tm_total_l50_incidental);REPORT(tm_ret_l50_incidental)
  REPORT(hm_M);REPORT(hm_rec);REPORT(hm_sel);REPORT(hm_ptm)
  REPORT(sf_Mi);REPORT(sf_Mm);REPORT(sf_rec);REPORT(sf_sel);REPORT(sf_pm)
  REPORT(tf_Mi);REPORT(tf_Mm);REPORT(tf_rec);REPORT(tf_sel);REPORT(tf_pm)
  REPORT(hf_Mi);REPORT(hf_Mm);REPORT(hf_rec);REPORT(hf_sel);REPORT(hf_pm)
  REPORT(sel_sm);REPORT(sel_tm);REPORT(sel_sf);REPORT(sel_tf)
  REPORT(shared_survey_prior)
  sm_rec_prop<-sm$rec_prop;tm_rec_prop<-tm$rec_prop;hm_rec_prop<-hm$rec_prop
  sf_rec_prop<-sf$rec_prop;tf_rec_prop<-tf$rec_prop;hf_rec_prop<-hf$rec_prop
  REPORT(sm_rec_prop);REPORT(tm_rec_prop);REPORT(hm_rec_prop)
  REPORT(sf_rec_prop);REPORT(tf_rec_prop);REPORT(hf_rec_prop)
  ADREPORT(sm_rec_prop);ADREPORT(tm_rec_prop);ADREPORT(hm_rec_prop)
  ADREPORT(sf_rec_prop);ADREPORT(tf_rec_prop);ADREPORT(hf_rec_prop)
  REPORT(hm_w);REPORT(hf_wi);REPORT(hf_wm)
  REPORT(hm_snow_endpoint);REPORT(hm_tanner_endpoint)
  REPORT(hf_M_bound)
  REPORT(wtm);REPORT(wms);REPORT(wfs)
  REPORT(sm_G);REPORT(tm_G);REPORT(hm_G);REPORT(sf_G);REPORT(tf_G);REPORT(hf_G)
  REPORT(sm_growth_alpha);REPORT(sm_growth_beta);REPORT(sm_growth_sd)
  REPORT(tm_growth_alpha);REPORT(tm_growth_beta);REPORT(tm_growth_sd)
  REPORT(hm_growth_alpha);REPORT(hm_growth_beta);REPORT(hm_growth_sd)
  REPORT(sf_growth_alpha);REPORT(sf_growth_beta);REPORT(sf_growth_sd)
  REPORT(tf_growth_alpha);REPORT(tf_growth_beta);REPORT(tf_growth_sd)
  REPORT(hf_growth_alpha);REPORT(hf_growth_beta);REPORT(hf_growth_sd)
  REPORT(hmw);REPORT(hms);REPORT(hfw);REPORT(hfs)
  REPORT(hm_ptm_parent_snow);REPORT(hm_ptm_parent_tanner)
  REPORT(hm_sel_parent_snow);REPORT(hm_sel_parent_tanner)
  REPORT(hf_sel_parent_snow);REPORT(hf_sel_parent_tanner)
  ADREPORT(wtm);ADREPORT(wms);ADREPORT(wfs)
  ADREPORT(sm_growth_alpha);ADREPORT(sm_growth_beta);ADREPORT(sm_growth_sd)
  ADREPORT(tm_growth_alpha);ADREPORT(tm_growth_beta);ADREPORT(tm_growth_sd)
  ADREPORT(hm_growth_alpha);ADREPORT(hm_growth_beta);ADREPORT(hm_growth_sd)
  ADREPORT(sf_growth_alpha);ADREPORT(sf_growth_beta);ADREPORT(sf_growth_sd)
  ADREPORT(tf_growth_alpha);ADREPORT(tf_growth_beta);ADREPORT(tf_growth_sd)
  ADREPORT(hf_growth_alpha);ADREPORT(hf_growth_beta);ADREPORT(hf_growth_sd)
  ADREPORT(pf);ADREPORT(pm)
  # Delta-method targets used for uncertainty ribbons. Annual abundance totals
  # are reported instead of every abundance-at-size cell to keep sdreport small.
  sm_Ni_total<-rowSums(sm_Ni);sm_Nm_total<-rowSums(sm_Nm)
  tm_Ni_total<-rowSums(tm_Ni);tm_Nm_total<-rowSums(tm_Nm);hm_N_total<-rowSums(hm_N)
  sm_male_gt101<-rowSums((sm_Ni+sm_Nm)[,d$snow_male$data$sizes>101,drop=FALSE])
  tm_male_gt101<-rowSums((tm_Ni+tm_Nm)[,d$tanner_male$data$sizes>101,drop=FALSE])
  hm_male_gt101<-rowSums(hm_N[,d$hybrid_male$data$sizes>101,drop=FALSE])
  sf_Ni_total<-rowSums(sf_Ni);sf_Nm_total<-rowSums(sf_Nm)
  tf_Ni_total<-rowSums(tf_Ni);tf_Nm_total<-rowSums(tf_Nm)
  hf_Ni_total<-rowSums(hf_Ni);hf_Nm_total<-rowSums(hf_Nm)
  ADREPORT(sm_Ni_total);ADREPORT(sm_Nm_total);ADREPORT(tm_Ni_total);ADREPORT(tm_Nm_total)
  ADREPORT(hm_N_total);ADREPORT(sf_Ni_total);ADREPORT(sf_Nm_total)
  ADREPORT(sm_male_gt101);ADREPORT(tm_male_gt101);ADREPORT(hm_male_gt101)
  ADREPORT(tf_Ni_total);ADREPORT(tf_Nm_total);ADREPORT(hf_Ni_total);ADREPORT(hf_Nm_total)
  ADREPORT(sm_Mi);ADREPORT(sm_Mm);ADREPORT(sm_rec)
  ADREPORT(sm_F_directed);ADREPORT(sm_F_incidental)
  ADREPORT(tm_Mi);ADREPORT(tm_Mm);ADREPORT(tm_rec)
  ADREPORT(tm_F_directed);ADREPORT(tm_F_incidental)
  ADREPORT(hm_M);ADREPORT(hm_rec);ADREPORT(sf_Mi);ADREPORT(sf_Mm);ADREPORT(sf_rec)
  ADREPORT(tf_Mi);ADREPORT(tf_Mm);ADREPORT(tf_rec);ADREPORT(hf_Mi);ADREPORT(hf_Mm);ADREPORT(hf_rec)
  rec_variance_explained<-exp(2*log_rec_loading)/(exp(2*log_rec_loading)+exp(2*log_rec_sigma))
  M_variance_explained<-exp(2*log_M_loading)/(exp(2*log_M_loading)+exp(2*log_M_sigma))
  REPORT(common_rec);REPORT(common_M);REPORT(log_rec_sigma);REPORT(log_M_sigma)
  REPORT(log_rec_loading);REPORT(log_M_loading);REPORT(rec_variance_explained);REPORT(M_variance_explained)
  REPORT(jr_parent_beta);REPORT(hm_jr_driver);REPORT(hf_jr_driver)
  if(jr_on)ADREPORT(jr_parent_beta)
  REPORT(M_loading);REPORT(common_M_centered)
  if(shared_M){ADREPORT(M_loading);ADREPORT(common_M_centered)}
  REPORT(components)
  sum(components)
 }
}

