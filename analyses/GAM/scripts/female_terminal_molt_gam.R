suppressPackageStartupMessages({
  library(mgcv)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args<-commandArgs(trailingOnly=TRUE)
argval<-function(key,default){z<-grep(paste0("^",key,"="),args,value=TRUE);if(length(z))sub(paste0("^",key,"="),"",z[1]) else default}
fit_path <- argval("--fit","models/model/fits/model_fit.rds")
parent_cache <- "data/cache/crabpack_parental_raw.rds"
hybrid_cache <- "data/cache/crabpack_hybrid_raw.rds"
out <- argval("--out","analyses/GAM")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
opath <- function(name) file.path(out,paste0("female_terminal_molt_",name))

x <- readRDS(fit_path);d <- x$data;r <- x$report
parental <- readRDS(parent_cache);hybrid <- readRDS(hybrid_cache)

temperature_occupied <- function(raw,stock){
  z <- raw$specimen
  h <- raw$haul %>% transmute(HAULJOIN,temperature=GEAR_TEMPERATURE) %>%
    distinct(HAULJOIN,.keep_all=TRUE)
  z$temperature <- h$temperature[match(z$HAULJOIN,h$HAULJOIN)]
  z <- z[z$SEX==2 & is.finite(z$CLUTCH_SIZE) & z$CLUTCH_SIZE==0 &
    is.finite(z$YEAR) & is.finite(z$temperature) &
    is.finite(z$SAMPLING_FACTOR) & z$SAMPLING_FACTOR>0 &
    is.finite(z$AREA_SWEPT) & z$AREA_SWEPT>0,]
  z$cpue <- z$SAMPLING_FACTOR/z$AREA_SWEPT
  z %>% group_by(process_year=YEAR) %>%
    summarise(temperature=sum(temperature*cpue)/sum(cpue),
      temperature_cpue=sum(cpue),temperature_records=n(),.groups="drop") %>%
    mutate(species=stock)
}

abundance_covariates <- function(Ni,Nm,years,stock){
  data.frame(species=stock,process_year=years,
    immature_females=rowSums(Ni)/1e9,mature_females=rowSums(Nm)/1e9)
}

observed_response <- function(block,stock){
  p <- as.data.frame(block$data$raw_maturity_probability)
  n <- as.data.frame(block$data$maturity_specimen_n)
  names(p) <- names(n) <- as.character(block$data$sizes)
  p$observation_year <- n$observation_year <- block$data$years
  pp <- p %>% pivot_longer(-observation_year,names_to="size",values_to="p_terminal_molt")
  nn <- n %>% pivot_longer(-observation_year,names_to="size",values_to="specimen_n")
  left_join(pp,nn,by=c("observation_year","size")) %>%
    filter(is.finite(p_terminal_molt),is.finite(specimen_n),specimen_n>0) %>%
    mutate(species=stock,size=as.numeric(size),process_year=observation_year-1L)
}

temperatures <- bind_rows(
  temperature_occupied(parental$snow,"Snow"),
  temperature_occupied(parental$tanner,"Tanner"),
  temperature_occupied(hybrid,"Hybrid"))
abundances <- bind_rows(
  abundance_covariates(r$sf_Ni,r$sf_Nm,d$snow_female$data$years,"Snow"),
  abundance_covariates(r$tf_Ni,r$tf_Nm,d$tanner_female$data$years,"Tanner"),
  abundance_covariates(r$hf_Ni,r$hf_Nm,d$hybrid_female$data$years,"Hybrid"))
responses <- bind_rows(
  observed_response(d$snow_female,"Snow"),
  observed_response(d$tanner_female,"Tanner"),
  observed_response(d$hybrid_female,"Hybrid"))

dat <- responses %>%
  inner_join(temperatures,by=c("species","process_year")) %>%
  inner_join(abundances,by=c("species","process_year")) %>%
  filter(is.finite(temperature),is.finite(immature_females),is.finite(mature_females)) %>%
  group_by(species) %>%
  mutate(p_model=pmin(.999,pmax(.001,p_terminal_molt)),
    log_immature_females=log1p(immature_females),
    log_mature_females=log1p(mature_females),
    obs_weight=sqrt(specimen_n)/mean(sqrt(specimen_n))) %>% ungroup()
write.csv(dat,opath("female_terminal_molt_gam_data.csv"),row.names=FALSE)
write.csv(temperatures,opath("temperature_occupied_immature_females.csv"),row.names=FALSE)

fit_one <- function(z){
  gam(p_model ~ s(size,k=8,bs="ts")+s(temperature,k=4,bs="ts")+
      s(log_immature_females,log_mature_females,k=8),
    family=betar(link="logit"),weights=obs_weight,data=z,
    select=TRUE,method="REML")
}
models <- lapply(split(dat,dat$species),fit_one)
saveRDS(models,opath("female_terminal_molt_gams.rds"))

performance <- bind_rows(lapply(names(models),function(sp){
  z <- dat[dat$species==sp,];m <- models[[sp]];s <- summary(m)
  data.frame(species=sp,n=nrow(z),years=n_distinct(z$process_year),
    first_process_year=min(z$process_year),last_process_year=max(z$process_year),
    median_specimen_n=median(z$specimen_n),deviance_explained=s$dev.expl,
    adjusted_r_squared=s$r.sq,REML=unname(m$gcv.ubre),AIC=AIC(m),
    boundary_fraction=mean(z$p_terminal_molt<=0|z$p_terminal_molt>=1))
}))
smooth_significance <- bind_rows(lapply(names(models),function(sp){
  s <- as.data.frame(summary(models[[sp]])$s.table)
  s$term <- rownames(s);s$species <- sp;rownames(s)<-NULL;s
})) %>% select(species,term,everything())
write.csv(performance,opath("model_performance.csv"),row.names=FALSE)
write.csv(smooth_significance,opath("smooth_significance.csv"),row.names=FALSE)
capture.output(lapply(models,summary),file=opath("gam_summaries.txt"))

fitted <- bind_rows(lapply(names(models),function(sp){
  z <- dat[dat$species==sp,];p <- predict(models[[sp]],newdata=z,type="link",se.fit=TRUE)
  z %>% transmute(species,observation_year,process_year,size,specimen_n,
    p_terminal_molt,fitted=plogis(as.numeric(p$fit)),
    lower=plogis(as.numeric(p$fit-1.96*p$se.fit)),upper=plogis(as.numeric(p$fit+1.96*p$se.fit)))
}))
write.csv(fitted,opath("observed_fitted_female_terminal_molt.csv"),row.names=FALSE)

heat <- fitted %>% select(species,process_year,size,Observed=p_terminal_molt,Fitted=fitted) %>%
  pivot_longer(c(Observed,Fitted),names_to="series",values_to="probability")
heat$species <- factor(heat$species,c("Snow","Tanner","Hybrid"))
pheat <- ggplot(heat,aes(process_year,size,fill=probability))+
  geom_tile()+facet_grid(species~series,scales="free_y")+
  scale_fill_viridis_c(limits=c(0,1))+theme_bw()+
  labs(title="Observed and GAM-fitted female terminal-molt probability",
    subtitle="Only observed new-shell proportions; process year is one year before survey year",
    x="Terminal-molt process year",y="Carapace width (mm)",fill="Probability")
ggsave(opath("observed_fitted_female_terminal_molt_heatmaps.png"),pheat,
  width=12,height=11,dpi=220)

smooth_1d <- bind_rows(lapply(names(models),function(sp){
  z <- dat[dat$species==sp,];m <- models[[sp]]
  bind_rows(lapply(c("size","temperature"),function(v){
    nd <- data.frame(size=median(z$size),temperature=median(z$temperature),
      log_immature_females=median(z$log_immature_females),
      log_mature_females=median(z$log_mature_females))
    nd <- nd[rep(1,160),];nd[[v]] <- seq(min(z[[v]]),max(z[[v]]),length.out=160)
    pt <- predict(m,newdata=nd,type="terms",se.fit=TRUE)
    j <- match(paste0("s(",v,")"),colnames(pt$fit))
    data.frame(species=sp,covariate=v,value=nd[[v]],effect=pt$fit[,j],
      lower=pt$fit[,j]-1.96*pt$se.fit[,j],upper=pt$fit[,j]+1.96*pt$se.fit[,j])
  }))
}))
write.csv(smooth_1d,opath("size_temperature_smooths.csv"),row.names=FALSE)
smooth_1d$covariate <- factor(smooth_1d$covariate,c("size","temperature"),
  c("Carapace width (mm)","CPUE-weighted bottom temperature (C)"))
psmooth <- ggplot(smooth_1d,aes(value,effect,colour=species,fill=species))+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.13,colour=NA)+geom_line(linewidth=.9)+
  facet_wrap(~covariate,scales="free_x")+theme_bw()+
  labs(title="Female terminal-molt GAM size and temperature smooths",
    subtitle="Centered partial effects on the logit scale; ribbons are approximately 95% CIs",
    x=NULL,y="Partial effect",colour=NULL,fill=NULL)
ggsave(opath("size_temperature_smooths.png"),psmooth,width=12,height=6,dpi=220)

surfaces <- bind_rows(lapply(names(models),function(sp){
  z <- dat[dat$species==sp,];m <- models[[sp]]
  nd <- expand.grid(immature_females=seq(min(z$immature_females),max(z$immature_females),length.out=80),
    mature_females=seq(min(z$mature_females),max(z$mature_females),length.out=80))
  nd$log_immature_females <- log1p(nd$immature_females)
  nd$log_mature_females <- log1p(nd$mature_females)
  nd$size <- median(z$size);nd$temperature <- median(z$temperature)
  pt <- predict(m,newdata=nd,type="terms",se.fit=TRUE)
  j <- grep("log_immature_females,log_mature_females",colnames(pt$fit),fixed=TRUE)
  nd$effect <- pt$fit[,j];nd$se <- pt$se.fit[,j];nd$species <- sp;nd
}))
write.csv(surfaces,opath("female_abundance_joint_smooth_surfaces.csv"),row.names=FALSE)
points <- dat %>% distinct(species,process_year,immature_females,mature_females)
psurface <- ggplot(surfaces,aes(immature_females,mature_females,fill=effect))+
  geom_tile()+geom_contour(aes(x=immature_females,y=mature_females,z=effect),
    data=surfaces,inherit.aes=FALSE,colour="white",alpha=.65)+
  geom_point(data=points,aes(x=immature_females,y=mature_females),inherit.aes=FALSE,
    size=.7,alpha=.4)+facet_wrap(~species,scales="free")+
  scale_fill_viridis_c()+theme_bw()+
  labs(title="Joint female-abundance effects on terminal molt",
    subtitle="Immature and mature latent abundance; centered logit-scale effects, with observed annual combinations",
    x="Immature females (billions)",y="Mature females (billions)",fill="Partial effect")
ggsave(opath("female_abundance_joint_smooth_surfaces.png"),psurface,
  width=15,height=5.5,dpi=220)

temp_plot <- temperatures %>% mutate(species=factor(species,c("Snow","Tanner","Hybrid")))
ptemp <- ggplot(temp_plot,aes(process_year,temperature,colour=species))+
  geom_line(linewidth=.8)+geom_point(size=1.1)+facet_wrap(~species,ncol=1)+theme_bw()+
  labs(title="Temperature occupied by immature females",
    subtitle="CPUE-weighted gear temperature for females with clutch size zero",
    x="Year",y="Bottom temperature (C)",colour=NULL)
ggsave(opath("temperature_occupied_immature_females.png"),ptemp,
  width=10,height=9,dpi=220)

print(performance)
print(smooth_significance)
