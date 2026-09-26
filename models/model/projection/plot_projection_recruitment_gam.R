suppressPackageStartupMessages({library(mgcv);library(ggplot2);library(patchwork)})
plot_projection_recruitment_gam<-function(out){
 mod<-readRDS(file.path(out,"hybrid_recruitment_gam.rds"))
 dat<-read.csv(file.path(out,"hybrid_recruitment_gam_training.csv"))
 dat$sex<-factor(dat$sex,levels=mod$xlevels$sex)
 gx<-seq(min(dat$tanner_male),max(dat$tanner_male),length.out=100)
 gy<-seq(min(dat$snow_male),max(dat$snow_male),length.out=100)
 surface<-expand.grid(tanner_male=gx,snow_male=gy)
 surface$sex<-factor("Female",levels=levels(dat$sex))
 surface$log_snow_female<-median(dat$log_snow_female)
 term<-predict(mod,surface,type="terms",se.fit=TRUE,unconditional=TRUE)
 j<-match("s(tanner_male,snow_male)",colnames(term$fit));stopifnot(!is.na(j))
 surface$effect<-term$fit[,j];surface$se<-term$se.fit[,j]
 pts<-unique(dat[c("tanner_male","snow_male")])
 hull<-pts[chull(pts),];hull<-rbind(hull,hull[1,])
 p<-ggplot(surface,aes(tanner_male,snow_male,fill=effect))+
  geom_raster()+geom_contour(aes(tanner_male,snow_male,z=effect),inherit.aes=FALSE,colour="white",alpha=.6)+
  geom_path(data=hull,aes(tanner_male,snow_male),inherit.aes=FALSE,linetype=2)+
  geom_point(data=pts,aes(tanner_male,snow_male),inherit.aes=FALSE,size=1.3)+
  scale_fill_viridis_c()+theme_bw()+
  labs(title="Hybrid recruitment: joint Tanner and snow male smooth",
   subtitle="Exact GAM used in spring-timing projections; dashed boundary encloses training combinations",
   caption="Surface outside the dashed boundary is extrapolation. Points are fitted parental abundance estimates.",
   x="Large mature Tanner males (billions)",y="Large mature snow males (billions)",
   fill="Partial effect on\nlog(recruitment + 1)")
 male_plot<-p+labs(title="Joint male abundance",subtitle=NULL,caption=NULL)
 write.csv(surface,file.path(out,"male_abundance_smooth_surface.csv"),row.names=FALSE)
 fg<-data.frame(snow_female=seq(min(dat$snow_female),max(dat$snow_female),length.out=200),
  tanner_male=median(dat$tanner_male),snow_male=median(dat$snow_male),
  sex=factor("Female",levels=levels(dat$sex)))
 fg$log_snow_female<-log(fg$snow_female)
 term<-predict(mod,fg,type="terms",se.fit=TRUE,unconditional=TRUE)
 j<-match("s(log_snow_female)",colnames(term$fit));stopifnot(!is.na(j))
 fg$effect<-term$fit[,j];fg$lower<-fg$effect-1.96*term$se.fit[,j]
 fg$upper<-fg$effect+1.96*term$se.fit[,j]
 p<-ggplot(fg,aes(snow_female,effect))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#E69F00",alpha=.22)+
  geom_line(linewidth=.9)+geom_rug(data=dat,aes(x=snow_female),inherit.aes=FALSE,alpha=.4)+theme_bw()+
  labs(title="Hybrid recruitment: mature snow-female smooth",
   subtitle="Exact projection GAM; smooth of log abundance, displayed against abundance",
   x="Mature snow females (billions)",y="Partial effect on log(recruitment + 1)",
   caption="Approximate 95% smooth confidence interval; does not propagate first-stage population-fit uncertainty.")
 female_plot<-p+labs(title="Mature snow-female abundance",subtitle=NULL,caption=NULL)
 combined<-(female_plot|male_plot)+plot_annotation(
  title="Hybrid recruitment GAM: all smooth terms",
  subtitle="Current spring-timing projection model; parental abundance lagged five years",
  caption="Centered effects on log(recruitment + 1). Female ribbon: approximate 95% confidence interval.\nMale panel: training points and dashed convex hull; outside the hull is extrapolation. First-stage fit uncertainty is not included.")
 ggsave(file.path(out,"hybrid_recruitment_gam_smooths.png"),combined,width=15,height=7,dpi=220)
 write.csv(fg,file.path(out,"snow_female_smooth.csv"),row.names=FALSE)
 fit<-predict(mod,dat,se.fit=TRUE,unconditional=TRUE)
 dat$fitted<-pmax(0,expm1(fit$fit));dat$lower<-pmax(0,expm1(fit$fit-1.96*fit$se.fit))
 dat$upper<-pmax(0,expm1(fit$fit+1.96*fit$se.fit))
 p<-ggplot(dat,aes(year,recruitment/1e6))+
  geom_ribbon(aes(ymin=lower/1e6,ymax=upper/1e6),fill="#2166AC",alpha=.2)+
  geom_line(aes(y=fitted/1e6),colour="#2166AC",linewidth=.8)+geom_point()+
  facet_wrap(~sex,ncol=1,scales="free_y")+theme_bw()+
  labs(title="Hybrid recruitment: GAM fitted values",
   subtitle="Points: population-model recruitment estimates; line: back-transformed fitted log response",
   x="Recruitment process year",y="Recruitment (millions)",
   caption="Ribbon: approximate confidence interval for the fitted log response, back-transformed; not a prediction interval.")
 ggsave(file.path(out,"gam_fitted_recruitment.png"),p,width=10,height=8,dpi=220)
 write.csv(dat,file.path(out,"gam_fitted_recruitment.csv"),row.names=FALSE)
 stopifnot(all(is.finite(surface$effect)),all(is.finite(fg$effect)),
  max(abs(as.numeric(predict(mod,dat))-as.numeric(fit$fit)))<1e-10)
 cat("Saved smooths from exact projection GAM to",out,"\n")
}
if(sys.nframe()==0L){
 args<-commandArgs(TRUE)
 out<-if(length(args))args[1] else "models/model/projections/projection_15yr"
 plot_projection_recruitment_gam(out)
}
