suppressPackageStartupMessages({
  library(mgcv)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args<-commandArgs(trailingOnly=TRUE)
argval<-function(key,default){z<-grep(paste0("^",key,"="),args,value=TRUE);if(length(z))sub(paste0("^",key,"="),"",z[1]) else default}
fit_path <- argval("--fit","models/model/fits/model_fit.rds")
cache_path <- "data/cache/crabpack_parental_raw.rds"
out <- argval("--out","analyses/GAM")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
opath <- function(name) file.path(out,paste0("male_terminal_molt_",name))

x <- readRDS(fit_path)
d <- x$data
r <- x$report
crab <- readRDS(cache_path)

# Terminal-molt decisions represented by the survey-year ogive occurred during
# the preceding annual interval. Consequently observation year y is assigned to
# process year y-1 before environmental and abundance covariates are joined.
temperature_occupied <- function(raw,stock,cutoff,inclusive_small){
  z <- raw$specimen
  h <- raw$haul %>%
    transmute(HAULJOIN,temperature=GEAR_TEMPERATURE) %>%
    distinct(HAULJOIN,.keep_all=TRUE)
  z$temperature <- h$temperature[match(z$HAULJOIN,h$HAULJOIN)]
  keep_size <- if(inclusive_small) z$SIZE_1MM<=cutoff else z$SIZE_1MM<cutoff
  z <- z[z$SEX==1 & keep_size & is.finite(z$YEAR) &
    is.finite(z$temperature) & is.finite(z$SAMPLING_FACTOR) &
    z$SAMPLING_FACTOR>0 & is.finite(z$AREA_SWEPT) & z$AREA_SWEPT>0,]
  z$cpue <- z$SAMPLING_FACTOR/z$AREA_SWEPT
  z %>% group_by(process_year=YEAR) %>%
    summarise(temperature=sum(temperature*cpue)/sum(cpue),
      temperature_cpue=sum(cpue),temperature_records=n(),.groups="drop") %>%
    mutate(species=stock)
}

abundance_covariates <- function(Ni,Nm,sizes,years,stock,cutoff,inclusive_small){
  small <- if(inclusive_small) sizes<=cutoff else sizes<cutoff
  data.frame(species=stock,process_year=years,
    small_males=rowSums(Ni[,small,drop=FALSE]+Nm[,small,drop=FALSE])/1e9,
    large_males=rowSums(Ni[,!small,drop=FALSE]+Nm[,!small,drop=FALSE])/1e9)
}

response_long <- function(block,stock){
  z <- as.data.frame(block$data$prop_term_molt)
  names(z) <- as.character(block$data$sizes)
  z$observation_year <- block$data$years
  z %>% pivot_longer(-observation_year,names_to="size",values_to="p_terminal_molt") %>%
    mutate(species=stock,size=as.numeric(size),process_year=observation_year-1L)
}

temperatures <- bind_rows(
  temperature_occupied(crab$snow,"Snow",76,TRUE),
  temperature_occupied(crab$tanner,"Tanner",101,FALSE))
abundances <- bind_rows(
  abundance_covariates(r$sm_Ni,r$sm_Nm,d$snow_male$data$sizes,
    d$snow_male$data$years,"Snow",76,TRUE),
  abundance_covariates(r$tm_Ni,r$tm_Nm,d$tanner_male$data$sizes,
    d$tanner_male$data$years,"Tanner",101,FALSE))
responses <- bind_rows(
  response_long(d$snow_male,"Snow"),response_long(d$tanner_male,"Tanner"))

dat <- responses %>%
  inner_join(temperatures,by=c("species","process_year")) %>%
  inner_join(abundances,by=c("species","process_year")) %>%
  filter(is.finite(p_terminal_molt),is.finite(size),is.finite(temperature),
    is.finite(small_males),is.finite(large_males)) %>%
  mutate(p_model=pmin(.999,pmax(.001,p_terminal_molt)),
    log_small_males=log1p(small_males),log_large_males=log1p(large_males))
write.csv(dat,opath("terminal_molt_gam_data.csv"),row.names=FALSE)
write.csv(temperatures,opath("temperature_occupied_preterminal_males.csv"),row.names=FALSE)

fit_one <- function(z){
  gam(p_model ~ s(size,k=8,bs="ts")+s(temperature,k=4,bs="ts")+
      s(log_small_males,log_large_males,k=8),
    family=betar(link="logit"),data=z,select=TRUE,method="REML")
}
models <- lapply(split(dat,dat$species),fit_one)
saveRDS(models,opath("terminal_molt_gams.rds"))

performance <- bind_rows(lapply(names(models),function(sp){
  m <- models[[sp]];z <- dat[dat$species==sp,];s <- summary(m)
  data.frame(species=sp,n=nrow(z),years=n_distinct(z$process_year),
    first_process_year=min(z$process_year),last_process_year=max(z$process_year),
    deviance_explained=s$dev.expl,adjusted_r_squared=s$r.sq,
    REML=unname(m$gcv.ubre),AIC=AIC(m),
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
  z %>% transmute(species,observation_year,process_year,size,p_terminal_molt,
    fitted=plogis(as.numeric(p$fit)),lower=plogis(as.numeric(p$fit-1.96*p$se.fit)),
    upper=plogis(as.numeric(p$fit+1.96*p$se.fit)))
}))
write.csv(fitted,opath("observed_fitted_terminal_molt.csv"),row.names=FALSE)

heat <- fitted %>% select(species,process_year,size,Observed=p_terminal_molt,Fitted=fitted) %>%
  pivot_longer(c(Observed,Fitted),names_to="series",values_to="probability")
pheat <- ggplot(heat,aes(process_year,size,fill=probability))+
  geom_tile()+facet_grid(species~series,scales="free_y")+
  scale_fill_viridis_c(limits=c(0,1))+theme_bw()+
  labs(title="Observed and GAM-fitted terminal-molt probability",
    subtitle="Process year is one year earlier than the survey year attached to the observed ogive",
    x="Terminal-molt process year",y="Carapace width (mm)",fill="Probability")
ggsave(opath("observed_fitted_terminal_molt_heatmaps.png"),pheat,
  width=12,height=8,dpi=220)

pscatter <- ggplot(fitted,aes(p_terminal_molt,fitted,colour=size))+
  geom_abline(slope=1,intercept=0,colour="grey55",linetype=2)+geom_point(alpha=.45,size=1)+
  facet_wrap(~species)+scale_colour_viridis_c()+coord_equal(xlim=c(0,1),ylim=c(0,1))+
  theme_bw()+labs(title="Observed versus fitted terminal-molt probability",
    x="Observed probability",y="Fitted probability",colour="Size (mm)")
ggsave(opath("observed_vs_fitted_terminal_molt.png"),pscatter,
  width=10,height=5,dpi=220)

smooth_1d <- bind_rows(lapply(names(models),function(sp){
  z <- dat[dat$species==sp,];m <- models[[sp]]
  bind_rows(lapply(c("size","temperature"),function(v){
    nd <- data.frame(size=median(z$size),temperature=median(z$temperature),
      log_small_males=median(z$log_small_males),log_large_males=median(z$log_large_males))
    nd <- nd[rep(1,160),]
    nd[[v]] <- seq(min(z[[v]]),max(z[[v]]),length.out=160)
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
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.15,colour=NA)+geom_line(linewidth=.9)+
  facet_wrap(~covariate,scales="free_x")+theme_bw()+
  labs(title="Terminal-molt GAM one-dimensional smooths",
    subtitle="Centered partial effects on the logit scale; ribbons are approximately 95% CIs",
    x=NULL,y="Partial effect",colour=NULL,fill=NULL)
ggsave(opath("size_temperature_smooths.png"),psmooth,width=12,height=6,dpi=220)

surfaces <- bind_rows(lapply(names(models),function(sp){
  z <- dat[dat$species==sp,];m <- models[[sp]]
  nd <- expand.grid(small_males=seq(min(z$small_males),max(z$small_males),length.out=80),
    large_males=seq(min(z$large_males),max(z$large_males),length.out=80))
  nd$log_small_males <- log1p(nd$small_males);nd$log_large_males <- log1p(nd$large_males)
  nd$size <- median(z$size);nd$temperature <- median(z$temperature)
  pt <- predict(m,newdata=nd,type="terms",se.fit=TRUE)
  j <- grep("log_small_males,log_large_males",colnames(pt$fit),fixed=TRUE)
  nd$effect <- pt$fit[,j];nd$se <- pt$se.fit[,j];nd$species <- sp;nd
}))
write.csv(surfaces,opath("abundance_joint_smooth_surfaces.csv"),row.names=FALSE)
points <- dat %>% distinct(species,process_year,small_males,large_males)
psurface <- ggplot(surfaces,aes(small_males,large_males,fill=effect))+
  geom_tile()+
  geom_contour(aes(x=small_males,y=large_males,z=effect),data=surfaces,
    inherit.aes=FALSE,colour="white",alpha=.65)+
  geom_point(data=points,aes(x=small_males,y=large_males),inherit.aes=FALSE,
    size=.7,alpha=.4)+facet_wrap(~species,scales="free")+
  scale_fill_viridis_c()+theme_bw()+
  labs(title="Joint effects of small- and large-male abundance on terminal molt",
    subtitle="Centered partial effects on the logit scale; points are observed annual combinations",
    x="Small males (billions)",y="Large males (billions)",fill="Partial effect")
ggsave(opath("abundance_joint_smooth_surfaces.png"),psurface,
  width=12,height=6,dpi=220)

temp_plot <- temperatures %>% mutate(species=factor(species,c("Snow","Tanner")))
ptemp <- ggplot(temp_plot,aes(process_year,temperature,colour=species))+
  geom_line(linewidth=.8)+geom_point(size=1.2)+facet_wrap(~species,ncol=1)+
  theme_bw()+labs(title="Temperature occupied by pre-terminal males",
    subtitle="CPUE-weighted gear temperature; snow <=76 mm and Tanner <101 mm",
    x="Year",y="Bottom temperature (C)",colour=NULL)
ggsave(opath("temperature_occupied_preterminal_males.png"),ptemp,
  width=10,height=7,dpi=220)

print(performance)
print(smooth_significance)
