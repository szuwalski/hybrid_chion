suppressPackageStartupMessages({library(mgcv);library(ggplot2);library(dplyr)})

args<-commandArgs(trailingOnly=TRUE)
argval<-function(key,default){z<-grep(paste0("^",key,"="),args,value=TRUE);if(length(z))sub(paste0("^",key,"="),"",z[1]) else default}
fit_path<-argval("--fit","models/model/fits/model_fit.rds")
x <- readRDS(fit_path)
d <- x$data;r <- x$report
out <- argval("--out","analyses/GAM")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
opath <- function(name) file.path(out,paste0("hybrid_recruitment_",name))

large <- function(N,sizes,years,name){
  z <- data.frame(parent_year=years,
    value=rowSums(N[,sizes>=d$jr_large_male_min,drop=FALSE]))
  names(z)[2] <- name;z
}
parents <- merge(
  large(r$tm_Nm,d$tanner_male$data$sizes,d$tanner_male$data$years,"tanner_male"),
  large(r$sm_Nm,d$snow_male$data$sizes,d$snow_male$data$years,"snow_male"),
  by="parent_year") %>%
  left_join(data.frame(parent_year=d$snow_female$data$years,
    snow_female=as.numeric(rowSums(r$sf_Nm))),by="parent_year") %>%
  mutate(year=parent_year+d$jr_lag,tanner_male=tanner_male/1e9,
    snow_male=snow_male/1e9,snow_female=snow_female/1e9)
responses <- bind_rows(
  data.frame(year=d$hybrid_female$data$years,sex="Female",recruitment=as.numeric(r$hf_rec)),
  data.frame(year=d$hybrid_male$data$years,sex="Male",recruitment=as.numeric(r$hm_rec)))
dat <- inner_join(responses,parents,by="year") %>%
  filter(is.finite(tanner_male),is.finite(snow_male),is.finite(snow_female),
    is.finite(recruitment)) %>%
  mutate(sex=factor(sex,c("Female","Male")),log_recruitment=log1p(recruitment),
    log_snow_female=log(snow_female))

# Strongly constrain the female smooth and moderately constrain the bivariate
# male surface to avoid fitting fine-scale structure in sparse areas.
mod <- gam(log_recruitment ~ sex+s(log_snow_female,k=3)+
    s(tanner_male,snow_male,k=12),data=dat,
  select=TRUE,method="REML")
sm <- as.data.frame(summary(mod)$s.table);sm$term<-rownames(sm);rownames(sm)<-NULL
stats <- data.frame(n=nrow(dat),first_year=min(dat$year),last_year=max(dat$year),
  deviance_explained=summary(mod)$dev.expl,adjusted_r_squared=summary(mod)$r.sq,
  REML=mod$gcv.ubre,AIC=AIC(mod),
  female_smooth_basis_dimension=mod$smooth[[1]]$bs.dim,
  male_surface_basis_dimension=mod$smooth[[2]]$bs.dim)

pr <- predict(mod,newdata=dat,se.fit=TRUE)
fits <- dat %>% transmute(year,sex,observed=recruitment,
  fitted=pmax(0,expm1(as.numeric(pr$fit))),
  lower=pmax(0,expm1(as.numeric(pr$fit-1.96*pr$se.fit))),
  upper=pmax(0,expm1(as.numeric(pr$fit+1.96*pr$se.fit))))
write.csv(dat,opath("gam_data.csv"),row.names=FALSE)
write.csv(fits,opath("observed_fitted_recruitment.csv"),row.names=FALSE)
write.csv(stats,opath("model_performance.csv"),row.names=FALSE)
write.csv(sm,opath("smooth_significance.csv"),row.names=FALSE)
capture.output(summary(mod),file=opath("gam_summary.txt"))

pfit <- ggplot(fits,aes(year,observed))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#80B1D3",alpha=.25)+
  geom_line(aes(y=fitted),colour="#1F78B4",linewidth=.9)+geom_point(size=1.6)+
  facet_wrap(~sex,ncol=1,scales="free_y")+theme_bw()+
  labs(title="Hybrid recruitment fitted by parental abundances",
    subtitle=paste0("s(log snow females, k=3) + s(Tanner males, snow males, k=12), parental abundance lagged ",
      d$jr_lag," years"),
    x="Recruitment year",y="Recruitment")
ggsave(opath("observed_fitted_recruitment.png"),pfit,width=11,height=8,dpi=220)

gx <- seq(min(dat$tanner_male),max(dat$tanner_male),length.out=80)
gy <- seq(min(dat$snow_male),max(dat$snow_male),length.out=80)
surface <- expand.grid(tanner_male=gx,snow_male=gy,sex=factor("Female",levels(dat$sex)),
  snow_female=median(dat$snow_female))
surface$log_snow_female <- log(surface$snow_female)
term <- predict(mod,newdata=surface,type="terms",se.fit=TRUE)
j <- grep("tanner_male,snow_male",colnames(term$fit),fixed=TRUE)
surface$effect <- term$fit[,j];surface$se <- term$se.fit[,j]
write.csv(surface,opath("male_abundance_smooth_surface.csv"),row.names=FALSE)
psurface <- ggplot(surface,aes(tanner_male,snow_male,fill=effect))+
  geom_raster()+geom_contour(aes(z=effect),colour="white",alpha=.7)+
  geom_point(data=distinct(dat,tanner_male,snow_male),inherit.aes=FALSE,
    aes(tanner_male,snow_male),size=.8,alpha=.45)+
  scale_fill_viridis_c()+coord_cartesian(expand=FALSE)+theme_bw()+
  labs(title="Joint Tanner- and snow-male abundance smooth",
    subtitle="Color is the centered partial effect on log(recruitment + 1); points show observed predictor combinations",
    x="Large Tanner males (billions)",y="Large snow males (billions)",fill="Partial effect")
ggsave(opath("male_abundance_smooth_surface.png"),psurface,width=10,height=8,dpi=220)

female_grid <- data.frame(snow_female=seq(min(dat$snow_female),max(dat$snow_female),
  length.out=160),tanner_male=median(dat$tanner_male),snow_male=median(dat$snow_male),
  sex=factor("Female",levels(dat$sex)))
female_grid$log_snow_female <- log(female_grid$snow_female)
female_term <- predict(mod,newdata=female_grid,type="terms",se.fit=TRUE)
fj <- match("s(log_snow_female)",colnames(female_term$fit))
female_grid$effect <- female_term$fit[,fj]
female_grid$lower <- female_grid$effect-1.96*female_term$se.fit[,fj]
female_grid$upper <- female_grid$effect+1.96*female_term$se.fit[,fj]
write.csv(female_grid,opath("snow_female_smooth.csv"),row.names=FALSE)
pfemale <- ggplot(female_grid,aes(snow_female,effect))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#FB9A99",alpha=.3)+
  geom_line(colour="#E31A1C",linewidth=1)+theme_bw()+
  labs(title="Mature snow-female abundance smooth",
    subtitle="Log-biomass smooth; centered partial effect on log(recruitment + 1), with approximately 95% CI",
    x="Mature snow females (billions)",y="Partial effect")
ggsave(opath("snow_female_smooth.png"),pfemale,width=9,height=6,dpi=220)

# Projection smooth figures use flat 2D surfaces; no perspective views.

# Compare alternative representations of mature snow-female abundance while
# holding the joint male surface constant.
common_rhs <- "sex+s(tanner_male,snow_male,k=12)"
female_models <- list(
  `Raw smooth`=gam(as.formula(paste("log_recruitment~s(snow_female,k=3)+",common_rhs)),
    data=dat,select=TRUE,method="REML"),
  `Log smooth`=mod,
  `Female exposure offset`=gam(as.formula(paste("log_recruitment~offset(log_snow_female)+",common_rhs)),
    data=dat,select=TRUE,method="REML"))
female_comparison <- bind_rows(lapply(names(female_models),function(nm){
  m <- female_models[[nm]];s <- summary(m)
  data.frame(model=nm,n=nrow(dat),deviance_explained=s$dev.expl,
    adjusted_r_squared=s$r.sq,REML=unname(m$gcv.ubre),AIC=AIC(m))
}))
write.csv(female_comparison,opath("female_formulation_comparison.csv"),row.names=FALSE)

female_x <- seq(min(dat$snow_female),max(dat$snow_female),length.out=160)
base_new <- data.frame(snow_female=female_x,log_snow_female=log(female_x),
  tanner_male=median(dat$tanner_male),snow_male=median(dat$snow_male),
  sex=factor("Female",levels(dat$sex)))
female_effects <- bind_rows(lapply(c("Raw smooth","Log smooth"),function(nm){
  m <- female_models[[nm]];pt <- predict(m,newdata=base_new,type="terms",se.fit=TRUE)
  target <- if(nm=="Raw smooth") "s(snow_female)" else "s(log_snow_female)"
  jv <- match(target,colnames(pt$fit))
  data.frame(snow_female=female_x,effect=pt$fit[,jv],
    lower=pt$fit[,jv]-1.96*pt$se.fit[,jv],upper=pt$fit[,jv]+1.96*pt$se.fit[,jv],model=nm)
}),data.frame(snow_female=female_x,effect=log(female_x)-mean(log(dat$snow_female)),
  lower=NA_real_,upper=NA_real_,model="Female exposure offset"))
write.csv(female_effects,opath("female_formulation_effects.csv"),row.names=FALSE)
pcompare <- ggplot(female_effects,aes(snow_female,effect,colour=model,fill=model))+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.13,colour=NA)+geom_line(linewidth=1)+
  theme_bw()+labs(title="Alternative mature snow-female formulations",
    subtitle="Centered effects on log(recruitment + 1); offset fixes the female coefficient at one",
    x="Mature snow females (billions)",y="Female contribution",colour=NULL,fill=NULL)
ggsave(opath("female_formulation_comparison.png"),pcompare,width=10,height=7,dpi=220)

comparison_fits <- bind_rows(lapply(names(female_models),function(nm){
  pp <- predict(female_models[[nm]],newdata=dat,se.fit=TRUE)
  dat %>% transmute(year,sex,model=nm,observed=recruitment,
    fitted=pmax(0,expm1(as.numeric(pp$fit))),
    lower=pmax(0,expm1(as.numeric(pp$fit-1.96*pp$se.fit))),
    upper=pmax(0,expm1(as.numeric(pp$fit+1.96*pp$se.fit))))
}))
write.csv(comparison_fits,opath("female_formulation_fitted_recruitment.csv"),row.names=FALSE)
pcomparison_fit <- ggplot(comparison_fits,aes(year,observed))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#80B1D3",alpha=.22)+
  geom_line(aes(y=fitted),colour="#1F78B4",linewidth=.8)+geom_point(size=1.2)+
  facet_grid(sex~model,scales="free_y")+theme_bw()+
  labs(title="Recruitment fits under alternative female formulations",
    x="Recruitment year",y="Recruitment")
ggsave(opath("female_formulation_fitted_recruitment.png"),pcomparison_fit,
  width=15,height=8,dpi=220)

# Alternative biological parameterization: replace the bivariate male-abundance
# surface with the Tanner proportion among all large males.
dat <- dat %>% mutate(tanner_male_fraction=tanner_male/(tanner_male+snow_male))
composition_mod <- gam(log_recruitment ~ sex+s(log_snow_female,k=3)+
    s(tanner_male_fraction,k=4,bs="ts"),data=dat,select=TRUE,method="REML")
composition_summary <- summary(composition_mod)
composition_stats <- data.frame(model="Log female plus male composition",n=nrow(dat),
  deviance_explained=composition_summary$dev.expl,
  adjusted_r_squared=composition_summary$r.sq,REML=unname(composition_mod$gcv.ubre),
  AIC=AIC(composition_mod))
composition_smooths <- as.data.frame(composition_summary$s.table)
composition_smooths$term <- rownames(composition_smooths);rownames(composition_smooths)<-NULL
write.csv(composition_stats,opath("male_composition_model_performance.csv"),row.names=FALSE)
write.csv(composition_smooths,opath("male_composition_smooth_significance.csv"),row.names=FALSE)
capture.output(composition_summary,file=opath("male_composition_gam_summary.txt"))

composition_grid <- data.frame(tanner_male_fraction=seq(min(dat$tanner_male_fraction),
  max(dat$tanner_male_fraction),length.out=160),
  log_snow_female=median(dat$log_snow_female),
  sex=factor("Female",levels(dat$sex)))
cpt <- predict(composition_mod,newdata=composition_grid,type="terms",se.fit=TRUE)
cj <- match("s(tanner_male_fraction)",colnames(cpt$fit))
composition_grid$effect <- cpt$fit[,cj]
composition_grid$lower <- composition_grid$effect-1.96*cpt$se.fit[,cj]
composition_grid$upper <- composition_grid$effect+1.96*cpt$se.fit[,cj]
write.csv(composition_grid,opath("male_composition_smooth.csv"),row.names=FALSE)
pcomposition <- ggplot(composition_grid,aes(tanner_male_fraction,effect))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#CAB2D6",alpha=.35)+
  geom_line(colour="#6A3D9A",linewidth=1)+
  geom_rug(data=distinct(dat,tanner_male_fraction),aes(y=NULL),inherit.aes=FALSE,
    sides="b",alpha=.45)+theme_bw()+
  labs(title="Large-male composition smooth",
    subtitle="Tanner proportion among large Tanner and snow males; ribbon is approximately 95% CI",
    x="Proportion of large males that are Tanner crab",y="Partial effect on log(recruitment + 1)")
ggsave(opath("male_composition_smooth.png"),pcomposition,width=9,height=6,dpi=220)

cp <- predict(composition_mod,newdata=dat,se.fit=TRUE)
composition_fits <- dat %>% transmute(year,sex,observed=recruitment,
  fitted=pmax(0,expm1(as.numeric(cp$fit))),
  lower=pmax(0,expm1(as.numeric(cp$fit-1.96*cp$se.fit))),
  upper=pmax(0,expm1(as.numeric(cp$fit+1.96*cp$se.fit))))
write.csv(composition_fits,opath("male_composition_fitted_recruitment.csv"),row.names=FALSE)
pcomposition_fit <- ggplot(composition_fits,aes(year,observed))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#CAB2D6",alpha=.3)+
  geom_line(aes(y=fitted),colour="#6A3D9A",linewidth=.9)+geom_point(size=1.5)+
  facet_wrap(~sex,ncol=1,scales="free_y")+theme_bw()+
  labs(title="Hybrid recruitment fitted using large-male composition",
    x="Recruitment year",y="Recruitment")
ggsave(opath("male_composition_fitted_recruitment.png"),pcomposition_fit,
  width=11,height=8,dpi=220)

print(stats);print(sm)
print(female_comparison)
print(composition_stats);print(composition_smooths)
