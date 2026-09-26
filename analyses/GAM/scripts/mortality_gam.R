suppressPackageStartupMessages({
  library(RTMB)
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
sea_ice_path <- "data/sea_ice_yearly.csv"
out <- argval("--out","analyses/GAM")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
opath <- function(name) file.path(out,paste0("mortality_",name))

x <- readRDS(fit_path)
if(is.null(x$sdr) || !isTRUE(x$sdr$pdHess)) stop("The saved JMR fit needs a positive-definite Hessian")
d <- x$data
r <- x$report
parental <- readRDS(parent_cache)
hybrid <- readRDS(hybrid_cache)
sea_ice <- read.csv(sea_ice_path) %>%
  transmute(year=as.integer(year),sea_ice_fraction=as.numeric(sea_ice_fraction))
if(anyDuplicated(sea_ice$year)) stop("sea_ice_yearly.csv contains duplicate years")

# Male maturity is not observed consistently in crabpack, so size partitions
# represent the JMR maturity states: 76 mm for Snow and 101 mm for Tanner.
# Hybrid males have one mortality component and therefore use all sizes.
component_temperature <- function(raw,stock,sex,maturity){
  z <- raw$specimen
  h <- raw$haul %>% transmute(HAULJOIN,temperature=GEAR_TEMPERATURE) %>%
    distinct(HAULJOIN,.keep_all=TRUE)
  z$temperature <- h$temperature[match(z$HAULJOIN,h$HAULJOIN)]
  keep <- z$SEX == ifelse(sex=="Male",1,2)
  if(sex=="Male"){
    if(stock=="Snow"){
      keep <- keep & if(maturity=="Immature") z$SIZE<76 else z$SIZE>=76
    } else if(stock=="Tanner") {
      keep <- keep & if(maturity=="Immature") z$SIZE<101 else z$SIZE>=101
    }
  } else if(maturity=="Immature") {
    keep <- keep & is.finite(z$CLUTCH_SIZE) & z$CLUTCH_SIZE==0
  } else {
    keep <- keep & is.finite(z$CLUTCH_SIZE) & z$CLUTCH_SIZE>0
  }
  z <- z[keep & is.finite(z$YEAR) & is.finite(z$temperature) &
    is.finite(z$SAMPLING_FACTOR) & z$SAMPLING_FACTOR>0 &
    is.finite(z$AREA_SWEPT) & z$AREA_SWEPT>0,]
  z$cpue <- z$SAMPLING_FACTOR/z$AREA_SWEPT
  z %>% group_by(year=YEAR) %>% summarise(
    temperature=sum(temperature*cpue)/sum(cpue),
    temperature_cpue=sum(cpue),temperature_records=n(),.groups="drop") %>%
    mutate(stock=stock,sex=sex,maturity=maturity)
}

spec <- data.frame(
  prefix=c("sm","sm","tm","tm","hm","sf","sf","tf","tf","hf","hf"),
  M_name=c("sm_Mi","sm_Mm","tm_Mi","tm_Mm","hm_M","sf_Mi","sf_Mm",
    "tf_Mi","tf_Mm","hf_Mi","hf_Mm"),
  N_name=c("sm_Ni","sm_Nm","tm_Ni","tm_Nm","hm_N","sf_Ni","sf_Nm",
    "tf_Ni","tf_Nm","hf_Ni","hf_Nm"),
  block=c("snow_male","snow_male","tanner_male","tanner_male","hybrid_male",
    "snow_female","snow_female","tanner_female","tanner_female",
    "hybrid_female","hybrid_female"),
  stock=c("Snow","Snow","Tanner","Tanner","Hybrid","Snow","Snow","Tanner",
    "Tanner","Hybrid","Hybrid"),
  sex=c(rep("Male",5),rep("Female",6)),
  maturity=c("Immature","Mature","Immature","Mature","Total","Immature",
    "Mature","Immature","Mature","Immature","Mature"),stringsAsFactors=FALSE)

raw_for <- function(stock) {
  if(stock=="Snow") parental$snow else if(stock=="Tanner") parental$tanner else hybrid
}

temperatures <- bind_rows(lapply(seq_len(nrow(spec)),function(i){
  maturity <- spec$maturity[i]
  # Hybrid male is a single total component and uses all observed male sizes.
  if(maturity=="Total") maturity <- "Total"
  component_temperature(raw_for(spec$stock[i]),spec$stock[i],spec$sex[i],maturity)
})) %>% distinct(stock,sex,maturity,year,.keep_all=TRUE)
write.csv(temperatures,opath("temperature_occupied_by_component.csv"),row.names=FALSE)

get_se <- function(name,n){
  idx <- which(names(x$sdr$value)==name)
  if(length(idx)!=n) stop("ADREPORT mismatch for ",name)
  as.numeric(x$sdr$sd[idx])
}

rows <- lapply(seq_len(nrow(spec)),function(i){
  s <- spec[i,]
  years <- d[[s$block]]$data$years
  N <- as.matrix(r[[s$N_name]])
  M <- as.numeric(r[[s$M_name]])
  if(nrow(N)!=length(years) || length(M)!=length(years)) stop("Year mismatch for ",s$M_name)
  abundance <- rowSums(N)/1e9
  sizes <- d[[s$block]]$data$sizes
  average_size <- as.numeric(N %*% sizes)/pmax(rowSums(N),.Machine$double.eps)
  M_se <- get_se(s$M_name,length(M))
  data.frame(year=years,stock=s$stock,sex=s$sex,maturity=s$maturity,
    component=paste(s$stock,s$sex,s$maturity,sep=" / "),M,M_se,
    mortality_probability=1-exp(-M),
    mortality_probability_se=exp(-M)*M_se,
    abundance,average_size)
})

mortality <- bind_rows(rows) %>%
  # Restrict inference to 1988 onward, when both survey gear and spatial
  # footprint were stable and mortality estimates became consistently comparable.
  # There was no 2020 survey, so its deviation was fixed. The terminal-year M
  # also has no subsequent survey transition. Neither is a GAM response datum.
  group_by(component) %>% filter(year>=1988,year<max(year),year!=2020) %>% ungroup() %>%
  inner_join(temperatures,by=c("year","stock","sex","maturity")) %>%
  inner_join(sea_ice,by="year") %>%
  filter(is.finite(mortality_probability),is.finite(mortality_probability_se),
    mortality_probability_se>0,is.finite(temperature),is.finite(average_size),
    is.finite(abundance),is.finite(sea_ice_fraction)) %>% group_by(component) %>%
  mutate(scaled_abundance=as.numeric(scale(abundance)),
    response=pmin(.999,pmax(.001,mortality_probability)),
    precision_weight=(1/mortality_probability_se)/mean(1/mortality_probability_se)) %>%
  ungroup()
write.csv(mortality,opath("mortality_gam_data.csv"),row.names=FALSE)

weight_diagnostics <- mortality %>% group_by(component) %>% summarise(
  n=n(),unweighted_mean=mean(mortality_probability),
  inverse_se_weighted_mean=weighted.mean(mortality_probability,precision_weight),
  minimum_weight=min(precision_weight),maximum_weight=max(precision_weight),
  near_floor_years=sum(M<0.12),
  weight_on_near_floor_years=sum(precision_weight[M<0.12])/sum(precision_weight),
  .groups="drop")
write.csv(weight_diagnostics,opath("inverse_se_weight_diagnostics.csv"),row.names=FALSE)

fit_one <- function(z){
  gam(response ~ s(scaled_abundance,k=4,bs="ts")+
      s(temperature,k=4,bs="ts")+s(average_size,k=4,bs="ts")+
      s(sea_ice_fraction,k=4,bs="ts"),
    family=betar(link="logit"),data=z,select=TRUE,method="REML")
}
models <- lapply(split(mortality,mortality$component),fit_one)
saveRDS(models,opath("mortality_gams.rds"))

performance <- bind_rows(lapply(names(models),function(nm){
  z <- mortality[mortality$component==nm,];m <- models[[nm]];sm <- summary(m)
  data.frame(component=nm,n=nrow(z),first_year=min(z$year),last_year=max(z$year),
    deviance_explained=sm$dev.expl,adjusted_r_squared=sm$r.sq,
    REML=unname(m$gcv.ubre),AIC=AIC(m))
}))
significance <- bind_rows(lapply(names(models),function(nm){
  tab <- as.data.frame(summary(models[[nm]])$s.table)
  tab$term <- rownames(tab);tab$component <- nm;rownames(tab)<-NULL;tab
})) %>% select(component,term,everything())
write.csv(performance,opath("model_performance.csv"),row.names=FALSE)
write.csv(significance,opath("smooth_significance.csv"),row.names=FALSE)
capture.output(lapply(models,summary),file=opath("gam_summaries.txt"))

fitted <- bind_rows(lapply(names(models),function(nm){
  z <- mortality[mortality$component==nm,]
  pr <- predict(models[[nm]],newdata=z,type="link",se.fit=TRUE)
  z %>% transmute(year,stock,sex,maturity,component,M,M_se,
    observed=mortality_probability,
    observed_lower=1-exp(-pmax(0,M-1.96*M_se)),
    observed_upper=1-exp(-(M+1.96*M_se)),
    fitted=plogis(as.numeric(pr$fit)),
    fitted_lower=plogis(as.numeric(pr$fit-1.96*pr$se.fit)),
    fitted_upper=plogis(as.numeric(pr$fit+1.96*pr$se.fit)))
}))
write.csv(fitted,opath("observed_fitted_mortality.csv"),row.names=FALSE)

component_levels <- with(spec,paste(stock,sex,maturity,sep=" / "))
fitted$component <- factor(fitted$component,component_levels)
panel_labels <- performance %>% mutate(
  component=factor(component,component_levels),
  label=sprintf("Deviance explained: %.1f%%",100*deviance_explained))
pfit <- ggplot(fitted,aes(year))+geom_ribbon(aes(ymin=fitted_lower,ymax=fitted_upper),
    fill="#2C7FB8",alpha=.22)+geom_line(aes(y=fitted),colour="#2C7FB8",linewidth=.8)+
  geom_errorbar(aes(ymin=observed_lower,ymax=observed_upper),width=0,
    colour="#4D4D4D",alpha=.65)+geom_point(aes(y=observed),size=1.15)+
  geom_text(data=panel_labels,aes(x=Inf,y=Inf,label=label),inherit.aes=FALSE,
    hjust=1.08,vjust=1.35,size=3.1)+
  facet_wrap(~component,ncol=3)+coord_cartesian(ylim=c(0,1))+theme_bw()+
  labs(title="JMR mortality estimates and GAM fits",
    subtitle="1988 onward; unweighted GAM fit; bars retain JMR uncertainty for interpretation",
    x="Mortality process year",y="Annual mortality probability (0-1)")
ggsave(opath("observed_fitted_mortality.png"),pfit,width=15,height=13,dpi=220)

smooths <- bind_rows(lapply(names(models),function(nm){
  z <- mortality[mortality$component==nm,];m <- models[[nm]]
  bind_rows(lapply(c("scaled_abundance","temperature","average_size","sea_ice_fraction"),function(v){
    nd <- data.frame(scaled_abundance=median(z$scaled_abundance),
      temperature=median(z$temperature),average_size=median(z$average_size),
      sea_ice_fraction=median(z$sea_ice_fraction))
    nd <- nd[rep(1,160),]
    nd[[v]] <- seq(min(z[[v]]),max(z[[v]]),length.out=160)
    pt <- predict(m,newdata=nd,type="terms",se.fit=TRUE)
    j <- match(paste0("s(",v,")"),colnames(pt$fit))
    data.frame(component=nm,covariate=v,value=nd[[v]],effect=pt$fit[,j],
      lower=pt$fit[,j]-1.96*pt$se.fit[,j],upper=pt$fit[,j]+1.96*pt$se.fit[,j])
  }))
}))
write.csv(smooths,opath("covariate_smooths.csv"),row.names=FALSE)
smooths <- smooths %>% separate(component,c("stock","sex","maturity"),sep=" / ",remove=FALSE)
smooths$covariate <- factor(smooths$covariate,
  c("scaled_abundance","temperature","average_size","sea_ice_fraction"),
  c("Scaled abundance","CPUE-weighted bottom temperature (C)",
    "Average carapace width (mm)","Annual sea-ice fraction"))
smooths$group <- paste(smooths$sex,smooths$maturity,sep=" / ")
psmooth <- ggplot(smooths,aes(value,effect,colour=group,fill=group))+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.10,colour=NA)+geom_line(linewidth=.8)+
  facet_grid(stock~covariate,scales="free_x")+theme_bw()+
  labs(title="Mortality GAM covariate smooths",
    subtitle="1988 onward; centered partial effects on the logit scale with approximate 95% intervals",
    x=NULL,y="Partial effect",colour=NULL,fill=NULL)
ggsave(opath("covariate_smooths.png"),psmooth,width=15,height=11,dpi=220)

ptemp <- ggplot(filter(temperatures,year>=1988),
    aes(year,temperature,colour=interaction(sex,maturity)))+
  geom_line(linewidth=.7)+geom_point(size=.9)+facet_wrap(~stock,ncol=1)+theme_bw()+
  labs(title="Temperature occupied by mortality-model components",
    subtitle="CPUE-weighted bottom temperature; Snow split at 76 mm, Tanner at 101 mm, Hybrid uses all male sizes",
    x="Year",y="Bottom temperature (C)",colour="Component")
ggsave(opath("temperature_occupied_by_component.png"),ptemp,width=10,height=9,dpi=220)

temperature_mortality <- mortality %>% mutate(
  lower=1-exp(-pmax(0,M-1.96*M_se)),upper=1-exp(-(M+1.96*M_se)),
  highlighted=year %in% c(2018,2019),
  component=factor(component,component_levels))
write.csv(temperature_mortality,opath("mortality_vs_temperature.csv"),row.names=FALSE)
ptm <- ggplot(temperature_mortality,aes(temperature,mortality_probability))+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0,colour="grey75",alpha=.55)+
  geom_text(aes(label=year,colour=highlighted),size=2.35,check_overlap=FALSE)+
  facet_wrap(~component,ncol=3)+coord_cartesian(ylim=c(0,1))+
  scale_colour_manual(values=c(`FALSE`="#333333",`TRUE`="#D7301F"),guide="none")+
  theme_bw()+labs(title="JMR mortality probability versus temperature occupied",
    subtitle="1988 onward; labels are mortality process years; 2018-2019 highlighted in red; 2020 excluded",
    x="CPUE-weighted bottom temperature occupied (C)",
    y="Annual mortality probability (0-1)")
ggsave(opath("mortality_vs_temperature.png"),ptm,width=15,height=13,dpi=240)

print(performance)
print(significance)
