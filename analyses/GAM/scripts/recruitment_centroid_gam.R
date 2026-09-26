suppressPackageStartupMessages({library(mgcv);library(ggplot2);library(dplyr);library(tidyr)})

args<-commandArgs(trailingOnly=TRUE)
argval<-function(key,default){z<-grep(paste0("^",key,"="),args,value=TRUE);if(length(z))sub(paste0("^",key,"="),"",z[1]) else default}
fit_path <- argval("--fit","models/model/fits/model_fit.rds")
centroid_path <- "analyses/GAM/data/annual_component_centroids.csv"
out <- argval("--out","analyses/GAM")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
opath <- function(name) file.path(out,paste0("recruitment_centroid_",name))

x <- readRDS(fit_path); d <- x$data; r <- x$report
centroids <- read.csv(centroid_path)
lag_years <- d$jr_lag

annual_total <- function(N,years,name)
  data.frame(parent_year=years,value=rowSums(N),series=name)
large_total <- function(N,sizes,years,name)
  data.frame(parent_year=years,value=rowSums(N[,sizes>=d$jr_large_male_min,drop=FALSE]),series=name)

abundance <- bind_rows(
  annual_total(r$sf_Nm,d$snow_female$data$years,"snow_female_mature"),
  large_total(r$tm_Nm,d$tanner_male$data$sizes,d$tanner_male$data$years,"tanner_male_large"),
  large_total(r$sm_Nm,d$snow_male$data$sizes,d$snow_male$data$years,"snow_male_large")) %>%
  pivot_wider(names_from=series,values_from=value)

get_centroid <- function(target_stock,target_sex,target_component,prefix){
  z <- centroids %>% filter(.data$stock==.env$target_stock,
    .data$sex==.env$target_sex,.data$component==.env$target_component) %>%
    transmute(parent_year=year,latitude=.data$latitude,longitude=.data$longitude)
  names(z)[2:3] <- paste0(prefix,c("_lat","_lon"));z
}
spatial <- Reduce(function(a,b)full_join(a,b,by="parent_year"),list(
  get_centroid("Snow","Female","Mature","sf"),
  get_centroid("Tanner","Male","Mature","tm"),
  get_centroid("Snow","Male","Mature","sm")))

haversine <- function(lat1,lon1,lat2,lon2){
  rad <- pi/180;dlon <- (lon2-lon1)*rad;dlat <- (lat2-lat1)*rad
  a <- sin(dlat/2)^2+cos(lat1*rad)*cos(lat2*rad)*sin(dlon/2)^2
  6371*2*atan2(sqrt(a),sqrt(pmax(0,1-a)))
}
spatial <- spatial %>% mutate(
  sf_tm_distance_km=haversine(sf_lat,sf_lon,tm_lat,tm_lon),
  sf_sm_distance_km=haversine(sf_lat,sf_lon,sm_lat,sm_lon)) %>%
  select(parent_year,sf_tm_distance_km,sf_sm_distance_km)

predictors <- abundance %>% left_join(spatial,by="parent_year") %>%
  mutate(year=parent_year+lag_years,
    log_snow_female_mature=log1p(snow_female_mature),
    log_tanner_male_large=log1p(tanner_male_large),
    log_snow_male_large=log1p(snow_male_large))
responses <- bind_rows(
  data.frame(year=d$hybrid_female$data$years,sex="Female",recruitment=as.numeric(r$hf_rec)),
  data.frame(year=d$hybrid_male$data$years,sex="Male",recruitment=as.numeric(r$hm_rec)))
dat <- responses %>% inner_join(predictors,by="year") %>%
  filter(if_all(c(log_snow_female_mature,log_tanner_male_large,
    log_snow_male_large,sf_tm_distance_km,sf_sm_distance_km),is.finite)) %>%
  mutate(sex=factor(sex,c("Female","Male")),log_recruitment=log1p(recruitment))

base_formula <- log_recruitment ~ sex+
  s(log_snow_female_mature,k=4,bs="ts")+
  s(log_tanner_male_large,k=4,bs="ts")+
  s(log_snow_male_large,k=4,bs="ts")
spatial_formula <- update(base_formula,.~.+
  s(sf_tm_distance_km,k=4,bs="ts")+
  s(sf_sm_distance_km,k=4,bs="ts"))

base_model <- gam(base_formula,data=dat,select=TRUE,method="REML")
spatial_model <- gam(spatial_formula,data=dat,select=TRUE,method="REML")

model_stats <- bind_rows(
  data.frame(model="Abundance only",deviance_explained=summary(base_model)$dev.expl,
    adjusted_r_squared=summary(base_model)$r.sq,REML=base_model$gcv.ubre,AIC=AIC(base_model)),
  data.frame(model="Abundance plus centroid distance",deviance_explained=summary(spatial_model)$dev.expl,
    adjusted_r_squared=summary(spatial_model)$r.sq,REML=spatial_model$gcv.ubre,AIC=AIC(spatial_model)))
smooth_stats <- as.data.frame(summary(spatial_model)$s.table)
smooth_stats$covariate <- rownames(smooth_stats);rownames(smooth_stats)<-NULL
smooth_stats <- smooth_stats[,c("covariate",setdiff(names(smooth_stats),"covariate"))]

pred <- predict(spatial_model,newdata=dat,se.fit=TRUE)
fits <- dat %>% transmute(year,sex,observed=recruitment,
  fitted=pmax(0,expm1(as.numeric(pred$fit))),
  lower=pmax(0,expm1(as.numeric(pred$fit-1.96*pred$se.fit))),
  upper=pmax(0,expm1(as.numeric(pred$fit+1.96*pred$se.fit))))
write.csv(dat,opath("recruitment_gam_data.csv"),row.names=FALSE)
write.csv(fits,opath("observed_fitted_recruitment.csv"),row.names=FALSE)
write.csv(model_stats,opath("model_performance.csv"),row.names=FALSE)
write.csv(smooth_stats,opath("smooth_significance.csv"),row.names=FALSE)
capture.output(summary(spatial_model),file=opath("gam_summary.txt"))

p_fit <- ggplot(fits,aes(year,observed))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#2C7FB8",alpha=.2)+
  geom_line(aes(y=fitted),colour="#2C7FB8",linewidth=.9)+
  geom_point(size=1.7,colour="#222222")+facet_wrap(~sex,ncol=1,scales="free_y")+
  labs(title="JMR-estimated and centroid-GAM-fitted hybrid recruitment",
    subtitle=paste0("Parental abundance and centroid distances lagged ",lag_years,
      " years; points are JMR estimates, ribbon is GAM 95% confidence interval"),
    x="Recruitment year",y="Recruitment")+theme_bw()
ggsave(opath("observed_fitted_recruitment.png"),p_fit,width=11,height=8,dpi=220)

covariates <- c("log_snow_female_mature","log_tanner_male_large",
  "log_snow_male_large","sf_tm_distance_km","sf_sm_distance_km")
labels <- c(log_snow_female_mature="Log mature snow-female abundance",
  log_tanner_male_large="Log large Tanner-male abundance",
  log_snow_male_large="Log large snow-male abundance",
  sf_tm_distance_km="Snow female to Tanner male centroid distance (km)",
  sf_sm_distance_km="Snow female to snow male centroid distance (km)")
smooths <- lapply(covariates,function(v){
  nd <- dat[rep(1,120),];for(q in covariates)nd[[q]]<-median(dat[[q]])
  nd[[v]]<-seq(min(dat[[v]]),max(dat[[v]]),length.out=120)
  trm <- predict(spatial_model,newdata=nd,type="terms",se.fit=TRUE)
  col <- grep(paste0("^s\\(",v,"\\)$"),colnames(trm$fit))
  data.frame(covariate=labels[[v]],x=nd[[v]],effect=trm$fit[,col],
    lower=trm$fit[,col]-1.96*trm$se.fit[,col],
    upper=trm$fit[,col]+1.96*trm$se.fit[,col])
}) %>% bind_rows()
smooths$covariate<-factor(smooths$covariate,levels=unname(labels))
write.csv(smooths,opath("smooth_estimates.csv"),row.names=FALSE)
p_smooth <- ggplot(smooths,aes(x,effect))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#7FCDBB",alpha=.35)+
  geom_line(colour="#0868AC",linewidth=.9)+geom_hline(yintercept=0,linetype=3)+
  facet_wrap(~covariate,scales="free",ncol=2)+
  labs(title="Recruitment GAM smooths",subtitle="Centered partial effects on log(recruitment + 1) scale; k = 4, shrinkage smooths, select = TRUE",
    x=NULL,y="Partial effect")+theme_bw()
ggsave(opath("covariate_smooths.png"),p_smooth,width=12,height=10,dpi=220)

print(model_stats)
print(smooth_stats)
cat("Analysis years",min(dat$year),"to",max(dat$year),"with",nrow(dat),"sex-year observations\n")
