suppressPackageStartupMessages({library(RTMB);library(ggplot2);library(mgcv)})
source("models/model/projection/R/projection_maturation.R")
args<-commandArgs(trailingOnly=TRUE)
argval<-function(key,default){z<-grep(paste0("^",key,"="),args,value=TRUE);if(length(z))sub(paste0("^",key,"="),"",z[1]) else default}
horizon<-as.integer(argval("--years","15"));ndraw<-as.integer(argval("--draws","100"))
seed<-as.integer(argval("--seed","20260901"))
if(!is.finite(horizon)||horizon<1)stop("--years must be a positive integer")
if(!is.finite(ndraw)||ndraw<20)stop("--draws must be at least 20")
fit_path<-argval("--fit","models/model/fits/model_fit.rds")
out<-argval("--out",paste0("models/model/projection/output/projection_",horizon,"yr"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-readRDS(fit_path);d<-x$data
source("models/model/R/joint_classification_model.R")
stopifnot(!isTRUE(d$parental_mature_hybrid_recruitment==1))
if(is.null(x$sdr)||!isTRUE(x$sdr$pdHess))stop("Projection requires a positive-definite saved Hessian")
mortality_for_output<-function(z){
 if(!is.null(x$bridge_settings$mortality_link)&&x$bridge_settings$mortality_link=="hazard")
  1-exp(-z) else z
}

tailmed<-function(z,n=10)median(tail(as.numeric(z),min(n,length(z))),na.rm=TRUE)
rowmed<-function(z,n=10)apply(tail(z,min(n,nrow(z))),2,median,na.rm=TRUE)

project_parent_male<-function(b,Ni0,Nm0,G,Mi,Mm,F,rec,rp,fs,fr,ptm,n){
 ns<-length(Ni0);Ni<-Nm<-matrix(0,n+1,ns);Ni[1,]<-Ni0;Nm[1,]<-Nm0
 for(y in seq_len(n)){
  ti<-Ni[y,]*exp(-.59*Mi[y]);tm<-Nm[y,]*exp(-.59*Mm[y])
  ci<-ti*(1-exp(-F[y]*fs[y,]));cm<-tm*(1-exp(-F[y]*fs[y,]))
  ti<-ti*exp(-F[y]*fs[y,])+ci*(1-fr[y,])*b$data$catch$discard_survival
  tm<-tm*exp(-F[y]*fs[y,])+cm*(1-fr[y,])*b$data$catch$discard_survival
  grown<-as.vector(G%*%ti);grown[1:3]<-grown[1:3]+rec[y]*rp
  ni<-grown*(1-ptm[y,]);nm<-tm+grown*ptm[y,]
  Ni[y+1,]<-ni*exp(-.41*Mi[y]);Nm[y+1,]<-nm*exp(-.41*Mm[y])
 }
 list(i=Ni,m=Nm)
}
project_parent_male_two_fishery<-function(Ni0,Nm0,G,Mi,Mm,Fd,Fi,rec,rp,
 seld,retd,seli,reti,ptm,n){
 ns<-length(Ni0);Ni<-Nm<-matrix(0,n+1,ns);Ni[1,]<-Ni0;Nm[1,]<-Nm0
 for(y in seq_len(n)){
  ti<-Ni[y,]*exp(-.59*Mi[y]);tm<-Nm[y,]*exp(-.59*Mm[y])
  hd<-Fd[y]*seld[y,];hi<-Fi[y]*seli[y,];h<-hd+hi
  removed<-1-exp(-h);share_d<-hd/(h+1e-12);share_i<-hi/(h+1e-12)
  cid<-ti*removed*share_d;cmd<-tm*removed*share_d
  cii<-ti*removed*share_i;cmi<-tm*removed*share_i
  ti<-ti*exp(-h)+cid*(1-retd[y,])*.70+cii*(1-reti[y,])*.70
  tm<-tm*exp(-h)+cmd*(1-retd[y,])*.70+cmi*(1-reti[y,])*.70
  grown<-as.vector(G%*%ti);grown[1:3]<-grown[1:3]+rec[y]*rp
  ni<-grown*(1-ptm[y,]);nm<-tm+grown*ptm[y,]
  Ni[y+1,]<-ni*exp(-.41*Mi[y]);Nm[y+1,]<-nm*exp(-.41*Mm[y])
 }
 list(i=Ni,m=Nm)
}
project_female<-function(Ni0,Nm0,G,Mi,Mm,rec,rp,pm,n,
 spring=isTRUE(d$snow_female$data$spring_recruitment)){
 ns<-length(Ni0);Ni<-Nm<-matrix(0,n+1,ns);Ni[1,]<-Ni0;Nm[1,]<-Nm0
 before<-if(spring).59 else 1;after<-1-before
 for(y in seq_len(n)){
  grown<-as.vector(G%*%(Ni[y,]*exp(-before*Mi[y])));grown[1:3]<-grown[1:3]+rec[y]*rp
  newly<-grown*pm[y,];Ni[y+1,]<-(grown-newly)*exp(-after*Mi[y])
  Nm[y+1,]<-(Nm[y,]*exp(-before*Mm[y])+newly)*exp(-after*Mm[y])
 }
 list(i=Ni,m=Nm)
}
project_hybrid_male<-function(U0,T0,G,M,rec,rp,ptm,n,
 spring=isTRUE(d$hybrid_male$data$spring_recruitment)){
 ns<-length(U0);U<-T<-matrix(0,n+1,ns);U[1,]<-U0;T[1,]<-T0
 before<-if(spring).59 else 1;after<-1-before
 for(y in seq_len(n)){
  grown<-as.vector(G%*%(U[y,]*exp(-before*M[y])));grown[1:3]<-grown[1:3]+rec[y]*rp
  nt<-grown*ptm;U[y+1,]<-(grown-nt)*exp(-after*M[y])
  T[y+1,]<-(T[y,]*exp(-before*M[y])+nt)*exp(-after*M[y])
 }
 list(i=U,m=T)
}
inside_convex_hull<-function(x,y,hull){
 h2<-rbind(hull,hull[1,]);cross<-diff(h2[,1])*(y-h2[-nrow(h2),2])-
   diff(h2[,2])*(x-h2[-nrow(h2),1])
 all(cross>=-1e-10)||all(cross<=1e-10)
}

project_all<-function(r,par,gam_beta){
 rec_years<-Reduce(intersect,lapply(c("snow_male","tanner_male","hybrid_male","snow_female","tanner_female","hybrid_female"),function(nm)head(d[[nm]]$data$years,-1)))
 M_years<-setdiff(rec_years,2020L)
 F_years<-tail(intersect(head(d$snow_male$data$years,-1),
   head(d$tanner_male$data$years,-1)),10)
 ry<-sample(rec_years,horizon,replace=TRUE);my<-sample(M_years,horizon,replace=TRUE)
 fy<-sample(F_years,horizon,replace=TRUE)
 # Sample complete observation-year ogives paired across stocks and sexes.
 maturity_years<-sample(rec_years+1L,horizon,replace=TRUE)
 mat<-function(z,b)z[match(maturity_years,b$data$years),,drop=FALSE]
 at<-function(z,years,pick)as.numeric(z[match(pick,years)])
 fishmat<-function(prefix,suffix,b,pick)
   r[[paste0(prefix,"_",suffix)]][match(pick,b$data$years),,drop=FALSE]
 if(!is.null(r$sm_F_directed)){
  twofish<-function(prefix,b,Ni,Nm,G,Mi,Mm,rec,rp){
   project_parent_male_two_fishery(tail(Ni,1),tail(Nm,1),G,
    at(Mi,b$data$years,my),at(Mm,b$data$years,my),
    at(r[[paste0(prefix,"_F_directed")]],b$data$years,fy),
    at(r[[paste0(prefix,"_F_incidental")]],b$data$years,fy),
    at(rec,b$data$years,ry),rp,
    fishmat(prefix,"fish_total_directed",b,fy),fishmat(prefix,"fish_retain_directed",b,fy),
    fishmat(prefix,"fish_total_incidental",b,fy),fishmat(prefix,"fish_retain_incidental",b,fy),
    mat(b$data$prop_term_molt,b),horizon)
  }
  sm<-twofish("sm",d$snow_male,r$sm_Ni,r$sm_Nm,r$sm_G,r$sm_Mi,r$sm_Mm,r$sm_rec,r$sm_rec_prop)
  tm<-twofish("tm",d$tanner_male,r$tm_Ni,r$tm_Nm,r$tm_G,r$tm_Mi,r$tm_Mm,r$tm_rec,r$tm_rec_prop)
 }else{
  fs<-function(prefix,b,pick)r[[paste0(prefix,"_fish_total")]][match(pick,b$data$years),,drop=FALSE]
  fr<-function(prefix,b,pick)r[[paste0(prefix,"_fish_retain")]][match(pick,b$data$years),,drop=FALSE]
  sm<-project_parent_male(d$snow_male,tail(r$sm_Ni,1),tail(r$sm_Nm,1),r$sm_G,
    at(r$sm_Mi,d$snow_male$data$years,my),at(r$sm_Mm,d$snow_male$data$years,my),
    at(r$sm_F,d$snow_male$data$years,fy),at(r$sm_rec,d$snow_male$data$years,ry),r$sm_rec_prop,
    fs("sm",d$snow_male,fy),fr("sm",d$snow_male,fy),mat(d$snow_male$data$prop_term_molt,d$snow_male),horizon)
  tm<-project_parent_male(d$tanner_male,tail(r$tm_Ni,1),tail(r$tm_Nm,1),r$tm_G,
    at(r$tm_Mi,d$tanner_male$data$years,my),at(r$tm_Mm,d$tanner_male$data$years,my),
    at(r$tm_F,d$tanner_male$data$years,fy),at(r$tm_rec,d$tanner_male$data$years,ry),r$tm_rec_prop,
    fs("tm",d$tanner_male,fy),fr("tm",d$tanner_male,fy),mat(d$tanner_male$data$prop_term_molt,d$tanner_male),horizon)
 }
 sf<-project_female(tail(r$sf_Ni,1),tail(r$sf_Nm,1),r$sf_G,
   at(r$sf_Mi,d$snow_female$data$years,my),at(r$sf_Mm,d$snow_female$data$years,my),
   at(r$sf_rec,d$snow_female$data$years,ry),r$sf_rec_prop,mat(r$sf_pm,d$snow_female),horizon)
 tf<-project_female(tail(r$tf_Ni,1),tail(r$tf_Nm,1),r$tf_G,
   at(r$tf_Mi,d$tanner_female$data$years,my),at(r$tf_Mm,d$tanner_female$data$years,my),
   at(r$tf_rec,d$tanner_female$data$years,ry),r$tf_rec_prop,mat(r$tf_pm,d$tanner_female),horizon)
 parent_state<-function(hist,proj,years,year){
  if(year<=max(years))hist[match(year,years),] else proj[year-max(years)+1,]
 }
 # rec[y] belongs to process year y, and is first present in survey state y+1.
 future_recruitment_years<-max(d$hybrid_male$data$years)+seq_len(horizon)-1L
 parental<-t(vapply(future_recruitment_years,function(y){
  py<-y-d$jr_lag
  SF<-sum(parent_state(r$sf_Nm,sf$m,d$snow_female$data$years,py))
  TM<-sum(parent_state(r$tm_Nm,tm$m,d$tanner_male$data$years,py)[d$tanner_male$data$sizes>=d$jr_large_male_min])
  SM<-sum(parent_state(r$sm_Nm,sm$m,d$snow_male$data$years,py)[d$snow_male$data$sizes>=d$jr_large_male_min])
  c(parent_year=py,snow_female=SF/1e9,tanner_male=TM/1e9,snow_male=SM/1e9)
 },numeric(4)))
 gam_new<-rbind(
  data.frame(sex=factor("Female",levels=levels(rec_gam_data$sex)),
    log_snow_female=log(parental[,"snow_female"]),tanner_male=parental[,"tanner_male"],
    snow_male=parental[,"snow_male"]),
  data.frame(sex=factor("Male",levels=levels(rec_gam_data$sex)),
    log_snow_female=log(parental[,"snow_female"]),tanner_male=parental[,"tanner_male"],
    snow_male=parental[,"snow_male"]))
 Xg<-predict(rec_gam,newdata=gam_new,type="lpmatrix")
 eta<-as.numeric(Xg%*%gam_beta)
 residual_year<-sample(rec_gam_residual$year,horizon,replace=TRUE)
 female_res<-rec_gam_residual$female[match(residual_year,rec_gam_residual$year)]
 male_res<-rec_gam_residual$male[match(residual_year,rec_gam_residual$year)]
 if(any(!is.finite(c(female_res,male_res))))stop("Missing paired GAM residual for sampled year")
 hf_rec<-pmax(0,expm1(eta[seq_len(horizon)]+female_res))
 hm_rec<-pmax(0,expm1(eta[horizon+seq_len(horizon)]+male_res))
 male_scaled<-cbind((parental[,"tanner_male"]-rec_gam_scale["tm_mean"])/rec_gam_scale["tm_sd"],
   (parental[,"snow_male"]-rec_gam_scale["sm_mean"])/rec_gam_scale["sm_sd"])
 nearest<-apply(male_scaled,1,function(z)min(sqrt(rowSums((rec_gam_training_scaled-
   matrix(z,nrow(rec_gam_training_scaled),2,byrow=TRUE))^2))))
 outside_hull<-vapply(seq_len(horizon),function(i)!inside_convex_hull(
   parental[i,"tanner_male"],parental[i,"snow_male"],rec_gam_hull),logical(1))
 outside_range<-parental[,"tanner_male"]<min(rec_gam_data$tanner_male)|
   parental[,"tanner_male"]>max(rec_gam_data$tanner_male)|
   parental[,"snow_male"]<min(rec_gam_data$snow_male)|
   parental[,"snow_male"]>max(rec_gam_data$snow_male)|
   parental[,"snow_female"]<min(rec_gam_data$snow_female)|
   parental[,"snow_female"]>max(rec_gam_data$snow_female)
 hm<-project_hybrid_male(tail(r$hm_U,1),tail(r$hm_T,1),r$hm_G,
   at(r$hm_M,d$hybrid_male$data$years,my),hm_rec,r$hm_rec_prop,r$hm_ptm,horizon)
 hf<-project_female(tail(r$hf_Ni,1),tail(r$hf_Nm,1),r$hf_G,
   at(r$hf_Mi,d$hybrid_female$data$years,my),at(r$hf_Mm,d$hybrid_female$data$years,my),
   hf_rec,r$hf_rec_prop,mat(r$hf_pm,d$hybrid_female),horizon)
 list(sm=sm,tm=tm,hm=hm,sf=sf,tf=tf,hf=hf,
   sampled_recruitment_years=ry,sampled_mortality_years=my,sampled_fishing_years=fy,
   sampled_maturity_years=maturity_years,
   hybrid_parental=parental,hybrid_residual_years=residual_year,
   hybrid_outside_range=outside_range,hybrid_outside_hull=outside_hull,
   hybrid_nearest_scaled_distance=nearest,hm_rec=hm_rec,hf_rec=hf_rec,
   projected_rec=list(sm=at(r$sm_rec,d$snow_male$data$years,ry),
     tm=at(r$tm_rec,d$tanner_male$data$years,ry),hm=hm_rec,
     sf=at(r$sf_rec,d$snow_female$data$years,ry),
     tf=at(r$tf_rec,d$tanner_female$data$years,ry),hf=hf_rec),
   projected_M=list(sm_i=at(r$sm_Mi,d$snow_male$data$years,my),
     sm_m=at(r$sm_Mm,d$snow_male$data$years,my),
     tm_i=at(r$tm_Mi,d$tanner_male$data$years,my),
     tm_m=at(r$tm_Mm,d$tanner_male$data$years,my),
     hm=at(r$hm_M,d$hybrid_male$data$years,my),
     sf_i=at(r$sf_Mi,d$snow_female$data$years,my),
     sf_m=at(r$sf_Mm,d$snow_female$data$years,my),
     tf_i=at(r$tf_Mi,d$tanner_female$data$years,my),
     tf_m=at(r$tf_Mm,d$tanner_female$data$years,my),
     hf_i=at(r$hf_Mi,d$hybrid_female$data$years,my),
     hf_m=at(r$hf_Mm,d$hybrid_female$data$years,my)))
}

# Median crabpack calculated weight (grams), pre-interpolated to model bins.
# Bundling this small table avoids a dependency on large raw specimen caches.
weight_table<-read.csv("models/model/data/component_weight_at_size.csv")
weight_curve<-function(component,sizes){
 z<-weight_table[weight_table$component==component,]
 out<-z$weight_g[match(sizes,z$size)]
 if(any(!is.finite(out)))stop("Missing weight-at-size values for ",component)
 out
}
weights<-list(
 sm=weight_curve("sm",d$snow_male$data$sizes),
 tm=weight_curve("tm",d$tanner_male$data$sizes),
 hm=weight_curve("hm",d$hybrid_male$data$sizes),
 sf=weight_curve("sf",d$snow_female$data$sizes),
 tf=weight_curve("tf",d$tanner_female$data$sizes),
 hf=weight_curve("hf",d$hybrid_female$data$sizes))

spec<-list(
 sm=list(stock="Snow",sex="Male",ci="Immature",cm="Mature",years=d$snow_male$data$years,sizes=d$snow_male$data$sizes,hi="sm_Ni",hm="sm_Nm"),
 tm=list(stock="Tanner",sex="Male",ci="Immature",cm="Mature",years=d$tanner_male$data$years,sizes=d$tanner_male$data$sizes,hi="tm_Ni",hm="tm_Nm"),
 hm=list(stock="Hybrid",sex="Male",ci="Nonterminal",cm="Terminal",years=d$hybrid_male$data$years,sizes=d$hybrid_male$data$sizes,hi="hm_U",hm="hm_T"),
 sf=list(stock="Snow",sex="Female",ci="Immature",cm="Mature",years=d$snow_female$data$years,sizes=d$snow_female$data$sizes,hi="sf_Ni",hm="sf_Nm"),
 tf=list(stock="Tanner",sex="Female",ci="Immature",cm="Mature",years=d$tanner_female$data$years,sizes=d$tanner_female$data$sizes,hi="tf_Ni",hm="tf_Nm"),
 hf=list(stock="Hybrid",sex="Female",ci="Immature",cm="Mature",years=d$hybrid_female$data$years,sizes=d$hybrid_female$data$sizes,hi="hf_Ni",hm="hf_Nm"))

summarize_run<-function(r,par,draw,gam_beta){
 pr<-project_all(r,par,gam_beta);states<-bio<-process_long<-size_states<-list();ii<-ib<-ip<-0
 for(nm in names(spec)){
  s<-spec[[nm]];hist_i<-r[[s$hi]];hist_m<-r[[s$hm]]
  future_years<-max(s$years)+seq_len(horizon)
  for(band in c("<76 mm","76-101 mm",">101 mm")){
   keep_size<-switch(band,"<76 mm"=s$sizes<76,"76-101 mm"=s$sizes>=76&s$sizes<=101,">101 mm"=s$sizes>101)
   for(k in 1:3){
    hh<-switch(k,hist_i,hist_m,hist_i+hist_m)
    pp<-switch(k,pr[[nm]]$i,pr[[nm]]$m,pr[[nm]]$i+pr[[nm]]$m)
    comp<-c(s$ci,s$cm,"Total")[k]
    size_states[[length(size_states)+1]]<-data.frame(draw=draw,year=c(s$years,future_years),
      stock=s$stock,sex=s$sex,component=comp,size_group=band,
      period=c(rep("Historical",length(s$years)),rep("Projected",horizon)),
      value=c(rowSums(hh[,keep_size,drop=FALSE]),rowSums(pp[-1,keep_size,drop=FALSE])))
   }
  }
  for(k in 1:2){
   ii<-ii+1;states[[ii]]<-data.frame(draw=draw,year=s$years,stock=s$stock,sex=s$sex,
     component=if(k==1)s$ci else s$cm,period="Historical",
     value=rowSums(if(k==1)hist_i else hist_m))
   ii<-ii+1;states[[ii]]<-data.frame(draw=draw,year=future_years,stock=s$stock,sex=s$sex,
     component=if(k==1)s$ci else s$cm,period="Projected",
     value=rowSums((if(k==1)pr[[nm]]$i else pr[[nm]]$m)[-1,,drop=FALSE]))
  }
  keep<-s$sizes>101;w<-weights[[nm]][keep]/1e6
  ib<-ib+1;bio[[ib]]<-data.frame(draw=draw,year=s$years,stock=s$stock,sex=s$sex,
    period="Historical",value=as.vector((hist_i[,keep,drop=FALSE]+hist_m[,keep,drop=FALSE])%*%w))
  ib<-ib+1;bio[[ib]]<-data.frame(draw=draw,year=future_years,stock=s$stock,sex=s$sex,
    period="Projected",value=as.vector((pr[[nm]]$i[-1,keep,drop=FALSE]+pr[[nm]]$m[-1,keep,drop=FALSE])%*%w))
 }
 add_process<-function(values,years,stock,sex,component,process,period){
  if(process=="Natural mortality")values<-mortality_for_output(values)
   ip<<-ip+1;process_long[[ip]]<<-data.frame(draw=draw,year=years,stock=stock,
     sex=sex,component=component,process=process,period=period,value=as.numeric(values))
 }
 future_years<-max(d$years)+seq_len(horizon)-1L
 curves<-maturation_curves(r,d)
 for(nm in names(spec)){
  s<-spec[[nm]];selection<-maturation_selection[maturation_selection$key==nm,]
  j<-match(selection$size_mm,s$sizes);v<-curves[[nm]][,j]
  add_process(v[-1],s$years[-1]-1L,s$stock,s$sex,paste0(selection$size_mm," mm"),
    "Terminal molt probability","Historical")
  add_process(v[match(pr$sampled_maturity_years,s$years)],future_years,s$stock,s$sex,
    paste0(selection$size_mm," mm"),"Terminal molt probability","Projected")
 }
 rec_meta<-list(sm=c("Snow","Male"),tm=c("Tanner","Male"),hm=c("Hybrid","Male"),
   sf=c("Snow","Female"),tf=c("Tanner","Female"),hf=c("Hybrid","Female"))
 rec_hist<-list(sm=r$sm_rec,tm=r$tm_rec,hm=r$hm_rec,sf=r$sf_rec,tf=r$tf_rec,hf=r$hf_rec)
 rec_year<-list(sm=d$snow_male$data$years,tm=d$tanner_male$data$years,
   hm=d$hybrid_male$data$years,sf=d$snow_female$data$years,
   tf=d$tanner_female$data$years,hf=d$hybrid_female$data$years)
 for(nm in names(rec_meta)){
   y<-head(rec_year[[nm]],-1L)
   add_process(head(rec_hist[[nm]],length(y)),y,rec_meta[[nm]][1],rec_meta[[nm]][2],"Total","Recruitment","Historical")
   add_process(pr$projected_rec[[nm]],future_years,rec_meta[[nm]][1],rec_meta[[nm]][2],"Total","Recruitment","Projected")
 }
 M_meta<-list(sm_i=c("Snow","Male","Immature"),sm_m=c("Snow","Male","Mature"),
   tm_i=c("Tanner","Male","Immature"),tm_m=c("Tanner","Male","Mature"),hm=c("Hybrid","Male","Total"),
   sf_i=c("Snow","Female","Immature"),sf_m=c("Snow","Female","Mature"),
   tf_i=c("Tanner","Female","Immature"),tf_m=c("Tanner","Female","Mature"),
   hf_i=c("Hybrid","Female","Immature"),hf_m=c("Hybrid","Female","Mature"))
 M_hist<-list(sm_i=r$sm_Mi,sm_m=r$sm_Mm,tm_i=r$tm_Mi,tm_m=r$tm_Mm,hm=r$hm_M,
   sf_i=r$sf_Mi,sf_m=r$sf_Mm,tf_i=r$tf_Mi,tf_m=r$tf_Mm,hf_i=r$hf_Mi,hf_m=r$hf_Mm)
 M_year<-list(sm_i=d$snow_male$data$years,sm_m=d$snow_male$data$years,
   tm_i=d$tanner_male$data$years,tm_m=d$tanner_male$data$years,hm=d$hybrid_male$data$years,
   sf_i=d$snow_female$data$years,sf_m=d$snow_female$data$years,
   tf_i=d$tanner_female$data$years,tf_m=d$tanner_female$data$years,
   hf_i=d$hybrid_female$data$years,hf_m=d$hybrid_female$data$years)
 for(nm in names(M_meta)){
   y<-head(M_year[[nm]],-1L)
   add_process(head(M_hist[[nm]],length(y)),y,M_meta[[nm]][1],M_meta[[nm]][2],M_meta[[nm]][3],"Natural mortality","Historical")
   add_process(pr$projected_M[[nm]],future_years,M_meta[[nm]][1],M_meta[[nm]][2],M_meta[[nm]][3],"Natural mortality","Projected")
 }
 processes<-data.frame(draw=draw,projection_year=future_years,
   recruitment_year=future_years,receiving_survey_year=future_years+1L,
   sampled_recruitment_year=pr$sampled_recruitment_years,
   sampled_mortality_year=pr$sampled_mortality_years,
   sampled_fishing_year=pr$sampled_fishing_years,
   sampled_maturity_observation_year=pr$sampled_maturity_years,
   hybrid_parent_year=pr$hybrid_parental[,"parent_year"],
   mature_snow_female_billions=pr$hybrid_parental[,"snow_female"],
   large_tanner_male_billions=pr$hybrid_parental[,"tanner_male"],
   large_snow_male_billions=pr$hybrid_parental[,"snow_male"],
   hybrid_gam_residual_year=pr$hybrid_residual_years,
   outside_training_range=pr$hybrid_outside_range,
   outside_male_convex_hull=pr$hybrid_outside_hull,
   nearest_training_point_scaled_distance=pr$hybrid_nearest_scaled_distance,
   hybrid_male_recruitment=pr$hm_rec,hybrid_female_recruitment=pr$hf_rec)
 list(states=do.call(rbind,states),biomass=do.call(rbind,bio),processes=processes,
   size_states=do.call(rbind,size_states),
   process_long=do.call(rbind,process_long))
}

# Reconstruct the fitted objective so Hessian draws can be converted to reports.
p<-J1_parameters(d);for(nm in intersect(names(p),names(x$par_list)))p[[nm]]<-x$par_list[[nm]]
obj<-MakeADFun(function(z)do.call(J1_objective(d),z),p,map=x$map,silent=TRUE)
theta<-x$fit$par;r0<-obj$report(theta);par0<-obj$env$parList(theta)
maturation_audit<-select_maturation_sizes(r0,d,spec)
maturation_selection<-subset(maturation_audit,selected)
write.csv(maturation_audit,file.path(out,"terminal_molt_size_selection_audit.csv"),row.names=FALSE)
write.csv(maturation_selection,file.path(out,"terminal_molt_selected_sizes.csv"),row.names=FALSE)
stopifnot(identical(names(obj$par),names(theta)),
 abs(obj$fn(theta)-x$fit$objective)<1e-5)
source("models/model/projection/R/check_projection_transitions.R")
transition_checks<-check_projection_transitions(r0,d)
write.csv(transition_checks,file.path(out,"transition_verification.csv"),row.names=FALSE)
write.csv(check_spring_recruitment_pulse(),file.path(out,"spring_recruitment_pulse_verification.csv"),row.names=FALSE)
large_parent<-function(N,sizes,years,name){
 z<-data.frame(parent_year=years,value=rowSums(N[,sizes>=d$jr_large_male_min,drop=FALSE]))
 names(z)[2]<-name;z
}
rec_parents<-merge(
 large_parent(r0$tm_Nm,d$tanner_male$data$sizes,d$tanner_male$data$years,"tanner_male"),
 large_parent(r0$sm_Nm,d$snow_male$data$sizes,d$snow_male$data$years,"snow_male"),
 by="parent_year")
rec_parents<-merge(rec_parents,data.frame(parent_year=d$snow_female$data$years,
 snow_female=rowSums(r0$sf_Nm)),by="parent_year")
rec_parents$year<-rec_parents$parent_year+d$jr_lag
rec_parents$tanner_male<-rec_parents$tanner_male/1e9
rec_parents$snow_male<-rec_parents$snow_male/1e9
rec_parents$snow_female<-rec_parents$snow_female/1e9
rec_response<-rbind(
 data.frame(year=head(d$hybrid_female$data$years,-1),sex="Female",recruitment=head(as.numeric(r0$hf_rec),-1)),
 data.frame(year=head(d$hybrid_male$data$years,-1),sex="Male",recruitment=head(as.numeric(r0$hm_rec),-1)))
rec_gam_data<-merge(rec_response,rec_parents,by="year")
rec_gam_data<-rec_gam_data[is.finite(rec_gam_data$recruitment)&rec_gam_data$snow_female>0,]
rec_gam_data$sex<-factor(rec_gam_data$sex,c("Female","Male"))
rec_gam_data$log_recruitment<-log1p(rec_gam_data$recruitment)
rec_gam_data$log_snow_female<-log(rec_gam_data$snow_female)
rec_gam<-gam(log_recruitment~sex+s(log_snow_female,k=3)+
 s(tanner_male,snow_male,k=12),data=rec_gam_data,select=TRUE,method="REML")
saveRDS(rec_gam,file.path(out,"hybrid_recruitment_gam.rds"))
write.csv(rec_gam_data,file.path(out,"hybrid_recruitment_gam_training.csv"),row.names=FALSE)
source("models/model/projection/plot_projection_recruitment_gam.R")
plot_projection_recruitment_gam(out)
capture.output(summary(rec_gam),file=file.path(out,"hybrid_recruitment_gam_summary.txt"))
rr<-data.frame(year=rec_gam_data$year,sex=rec_gam_data$sex,
 residual=as.numeric(residuals(rec_gam,type="response")))
rec_gam_residual<-merge(
 setNames(rr[rr$sex=="Female",c("year","residual")],c("year","female")),
 setNames(rr[rr$sex=="Male",c("year","residual")],c("year","male")),by="year")
rec_gam_scale<-c(tm_mean=mean(rec_gam_data$tanner_male),tm_sd=sd(rec_gam_data$tanner_male),
 sm_mean=mean(rec_gam_data$snow_male),sm_sd=sd(rec_gam_data$snow_male))
rec_gam_training_scaled<-cbind(
 (rec_gam_data$tanner_male-rec_gam_scale["tm_mean"])/rec_gam_scale["tm_sd"],
 (rec_gam_data$snow_male-rec_gam_scale["sm_mean"])/rec_gam_scale["sm_sd"])
male_unique<-unique(rec_gam_data[c("tanner_male","snow_male")])
hh<-chull(male_unique$tanner_male,male_unique$snow_male)
rec_gam_hull<-as.matrix(male_unique[hh,c("tanner_male","snow_male")])
Vg<-vcov(rec_gam,unconditional=TRUE);Vg<-(Vg+t(Vg))/2
Cg<-chol(Vg+diag(1e-12,nrow(Vg)))
set.seed(seed)
gam_draw_beta<-sweep(matrix(rnorm(ndraw*length(coef(rec_gam))),ndraw,length(coef(rec_gam)))%*%Cg,
 2,coef(rec_gam),"+")
central<-summarize_run(r0,par0,0,coef(rec_gam))
process_values<-rbind(
 data.frame(stock=c("Snow","Tanner","Hybrid","Snow","Tanner","Hybrid"),sex=rep(c("Male","Female"),each=3),component="Total",process="Recruitment",value=c(tailmed(r0$sm_rec),tailmed(r0$tm_rec),tailmed(r0$hm_rec),tailmed(r0$sf_rec),tailmed(r0$tf_rec),tailmed(r0$hf_rec))),
 data.frame(stock=c("Snow","Snow","Tanner","Tanner","Hybrid","Snow","Snow","Tanner","Tanner","Hybrid","Hybrid"),sex=c(rep("Male",5),rep("Female",6)),component=c("Immature","Mature","Immature","Mature","Total","Immature","Mature","Immature","Mature","Immature","Mature"),process="Natural mortality",value=mortality_for_output(c(tailmed(r0$sm_Mi),tailmed(r0$sm_Mm),tailmed(r0$tm_Mi),tailmed(r0$tm_Mm),tailmed(r0$hm_M),tailmed(r0$sf_Mi),tailmed(r0$sf_Mm),tailmed(r0$tf_Mi),tailmed(r0$tf_Mm),tailmed(r0$hf_Mi),tailmed(r0$hf_Mm)))),
 data.frame(stock=rep(c("Snow","Tanner"),each=2),sex="Male",component=c("Directed","Incidental"),
  process="Fishing mortality",value=if(!is.null(r0$sm_F_directed))
   c(tailmed(r0$sm_F_directed),tailmed(r0$sm_F_incidental),
     tailmed(r0$tm_F_directed),tailmed(r0$tm_F_incidental)) else
   c(tailmed(r0$sm_F),NA,tailmed(r0$tm_F),NA)))
write.csv(process_values,file.path(out,"projection_process_values.csv"),row.names=FALSE)
V<-x$sdr$cov.fixed;V<-(V+t(V))/2
C<-chol(V+diag(1e-12,nrow(V)))
Z<-matrix(rnorm(ndraw*length(theta)),ndraw,length(theta))
draw_theta<-sweep(Z%*%C,2,theta,"+")
state_draws<-bio_draws<-process_draws<-process_long_draws<-list();accepted<-0
size_draws<-list()
draw_status<-data.frame(draw=seq_len(ndraw),status="accepted")
for(i in seq_len(ndraw)){
 z<-try(obj$report(draw_theta[i,]),silent=TRUE)
 if(inherits(z,"try-error")){draw_status$status[i]<-as.character(z);next}
 a<-try(summarize_run(z,obj$env$parList(draw_theta[i,]),i,gam_draw_beta[i,]),silent=TRUE)
 if(inherits(a,"try-error")){draw_status$status[i]<-as.character(a);next}
 if(!all(is.finite(a$states$value))||!all(is.finite(a$biomass$value))||
    !all(is.finite(a$process_long$value))){draw_status$status[i]<-"nonfinite trajectory";next}
 accepted<-accepted+1;state_draws[[accepted]]<-a$states;bio_draws[[accepted]]<-a$biomass
 process_draws[[accepted]]<-a$processes
 process_long_draws[[accepted]]<-a$process_long
 size_draws[[accepted]]<-a$size_states
}
write.csv(draw_status,file.path(out,"draw_status.csv"),row.names=FALSE)
if(accepted<max(20,.5*ndraw))stop("Too few finite uncertainty draws: ",accepted)

intervals<-function(point,draws,keys){
 qfun<-function(p,nm){z<-aggregate(value~.,draws[c(keys,"value")],function(x)quantile(x,p,na.rm=TRUE));names(z)[ncol(z)]<-nm;z}
 out<-Reduce(function(a,b)merge(a,b,by=keys),list(qfun(.025,"lower"),qfun(.5,"draw_median"),qfun(.975,"upper")))
 names(point)[names(point)=="value"]<-"estimate";merge(point[c(keys,"estimate")],out,by=keys)
}
state_out<-intervals(central$states,do.call(rbind,state_draws),c("year","stock","sex","component","period"))
bio_out<-intervals(central$biomass,do.call(rbind,bio_draws),c("year","stock","sex","period"))
process_out<-intervals(central$process_long,do.call(rbind,process_long_draws),
  c("year","stock","sex","component","process","period"))
size_out<-intervals(central$size_states,do.call(rbind,size_draws),
 c("year","stock","sex","component","size_group","period"))
size_out$estimate[size_out$period=="Projected"]<-size_out$draw_median[size_out$period=="Projected"]
write.csv(size_out,file.path(out,"latent_abundance_by_size_group.csv"),row.names=FALSE)
source("models/model/projection/R/plot_projection_sizes.R")
state_out$estimate[state_out$period=="Projected"]<-state_out$draw_median[state_out$period=="Projected"]
bio_out$estimate[bio_out$period=="Projected"]<-bio_out$draw_median[bio_out$period=="Projected"]
process_out$estimate[process_out$period=="Projected"]<-process_out$draw_median[process_out$period=="Projected"]
plot_projection_maturation(process_out,maturation_selection,out,accepted,maturation_audit)
write.csv(state_out,file.path(out,"population_component_trajectories.csv"),row.names=FALSE)
write.csv(bio_out,file.path(out,"biomass_gt101mm_metric_tons.csv"),row.names=FALSE)
process_draw_all<-do.call(rbind,process_draws)
write.csv(process_draw_all,file.path(out,"projection_process_draws.csv"),row.names=FALSE)
write.csv(process_out,file.path(out,"recruitment_mortality_trajectories.csv"),row.names=FALSE)

diag_median<-aggregate(cbind(mature_snow_female_billions,large_tanner_male_billions,
  large_snow_male_billions,nearest_training_point_scaled_distance)~projection_year,
  process_draw_all,median,na.rm=TRUE)
diag_flags<-aggregate(cbind(outside_training_range,outside_male_convex_hull)~projection_year,
  transform(process_draw_all,outside_training_range=as.numeric(outside_training_range),
    outside_male_convex_hull=as.numeric(outside_male_convex_hull)),mean,na.rm=TRUE)
gam_diagnostics<-merge(diag_median,diag_flags,by="projection_year")
names(gam_diagnostics)[names(gam_diagnostics)=="outside_training_range"]<-"proportion_outside_training_range"
names(gam_diagnostics)[names(gam_diagnostics)=="outside_male_convex_hull"]<-"proportion_outside_male_convex_hull"
write.csv(gam_diagnostics,file.path(out,"hybrid_recruitment_gam_extrapolation.csv"),row.names=FALSE)
p_extrap<-ggplot(gam_diagnostics,aes(projection_year,proportion_outside_male_convex_hull))+
  geom_line(colour="#D95F02",linewidth=.9)+geom_point(colour="#D95F02",size=1.8)+
  geom_line(aes(y=proportion_outside_training_range),colour="#1B9E77",linewidth=.9,linetype=2)+
  geom_point(aes(y=proportion_outside_training_range),colour="#1B9E77",size=1.8,shape=1)+
  coord_cartesian(ylim=c(0,1))+theme_bw()+
  labs(title="Hybrid-recruitment GAM extrapolation diagnostics",
    subtitle="Orange: outside observed male-abundance convex hull; green dashed: any predictor outside its observed range",
    x="Projected recruitment year",y="Proportion of projection draws")
ggsave(file.path(out,"hybrid_recruitment_gam_extrapolation.png"),p_extrap,width=10,height=6,dpi=220)

state_out$panel<-paste(state_out$stock,state_out$sex,state_out$component,sep=" / ")
state_out<-state_out[order(match(state_out$stock,c("Snow","Tanner","Hybrid")),
  match(state_out$sex,c("Female","Male")),state_out$component),]
state_out$panel<-factor(state_out$panel,levels=unique(state_out$panel))
p1<-ggplot(state_out,aes(year,estimate,group=interaction(panel,period)))+geom_ribbon(aes(ymin=lower,ymax=upper),fill="#7B3294",alpha=.16)+
  geom_line(aes(linetype=period),colour="#333333",linewidth=.65)+facet_wrap(~panel,scales="free_y",ncol=4)+
  geom_vline(xintercept=max(d$years),linetype=3,colour="grey50")+theme_bw()+theme(legend.position="bottom")+
  labs(title=paste0("JMR population trajectories and ",horizon,"-year projections"),
    subtitle=paste0("Historical process resampling plus parameter uncertainty; 95% intervals from ",accepted," projections"),
    x="Year",y="Abundance",linetype=NULL)
unlink(file.path(out,"population_component_trajectories.png"))
ggsave(file.path(out,"latent_abundance_by_species_sex_maturity.png"),p1,width=16,height=11,dpi=220)

process_out$panel<-paste(process_out$stock,process_out$sex,process_out$component,sep=" / ")
recruitment_out<-subset(process_out,process=="Recruitment")
recruitment_out<-recruitment_out[order(match(recruitment_out$sex,c("Female","Male")),
  match(recruitment_out$stock,c("Snow","Tanner","Hybrid"))),]
recruitment_out$panel<-factor(paste(recruitment_out$stock,recruitment_out$sex,sep=" / "),
  levels=unique(paste(recruitment_out$stock,recruitment_out$sex,sep=" / ")))
p_recruitment<-ggplot(recruitment_out,aes(year,estimate,group=interaction(panel,period)))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#008837",alpha=.16)+
  geom_line(aes(linetype=period),colour="#222222",linewidth=.65)+
  facet_wrap(~panel,scales="free_y",ncol=3)+
  geom_vline(xintercept=max(d$years),linetype=3,colour="grey50")+theme_bw()+
  theme(legend.position="bottom")+
  labs(title="Historical and projected recruitment",
    subtitle=paste0("Projected medians and 95% intervals from ",accepted," stochastic projections"),
    x="Recruitment year (first enters next survey state)",y="Recruitment",linetype=NULL)
ggsave(file.path(out,"recruitment_trajectories.png"),p_recruitment,width=14,height=8,dpi=220)

mortality_out<-subset(process_out,process=="Natural mortality")
mortality_out<-mortality_out[order(match(mortality_out$stock,c("Snow","Tanner","Hybrid")),
  match(mortality_out$sex,c("Female","Male")),mortality_out$component),]
mortality_out$panel<-factor(mortality_out$panel,levels=unique(mortality_out$panel))
p_mortality<-ggplot(mortality_out,aes(year,estimate,group=interaction(panel,period)))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#008837",alpha=.16)+
  geom_line(aes(linetype=period),colour="#222222",linewidth=.65)+
  facet_wrap(~panel,scales="fixed",ncol=4)+
  geom_vline(xintercept=max(d$years),linetype=3,colour="grey50")+
  coord_cartesian(ylim=c(0,1))+theme_bw()+
  theme(legend.position="bottom")+
  labs(title="Historical and projected natural mortality",
    subtitle=paste0("Species occupy rows; projected medians and 95% intervals from ",accepted," stochastic projections"),
    x="Process year (transition to next survey)",y="Annual natural mortality probability",linetype=NULL)
ggsave(file.path(out,"mortality_trajectories.png"),p_mortality,width=16,height=11,dpi=220)
unlink(file.path(out,"recruitment_mortality_trajectories.png"))

hybrid_rec<-subset(process_out,stock=="Hybrid" & process=="Recruitment")
p_hybrid_rec<-ggplot(hybrid_rec,aes(year,estimate,group=period))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#762A83",alpha=.18)+
  geom_line(aes(linetype=period),colour="#222222",linewidth=.8)+
  facet_wrap(~sex,ncol=1,scales="free_y")+
  geom_vline(xintercept=max(d$years),linetype=3,colour="grey50")+theme_bw()+
  theme(legend.position="bottom")+
  labs(title="Historical and projected hybrid recruitment",
    subtitle=paste0("Lag-",d$jr_lag," parental GAM with 2D male smooth; ",accepted," trial draws and 95% intervals"),
    x="Recruitment year (first enters next survey state)",y="Recruitment",linetype=NULL)
ggsave(file.path(out,"hybrid_recruitment_trajectory.png"),p_hybrid_rec,width=10,height=7,dpi=220)

bio_plot<-subset(bio_out,sex=="Male")
biomass_caps<-aggregate(estimate~stock,bio_plot,max,na.rm=TRUE)
names(biomass_caps)[2]<-"panel_cap";biomass_caps$panel_cap<-1.5*biomass_caps$panel_cap
bio_plot<-merge(bio_plot,biomass_caps,by="stock")
bio_plot$plot_lower<-bio_plot$lower
bio_plot$plot_upper<-bio_plot$upper
p2<-ggplot(bio_plot,aes(year,estimate,group=interaction(stock,period)))+
  geom_ribbon(aes(ymin=plot_lower,ymax=plot_upper),fill="#2166AC",alpha=.18)+
  geom_line(aes(linetype=period),colour="#333333",linewidth=.7)+facet_wrap(~stock,scales="free_y",ncol=3)+
  geom_blank(data=biomass_caps,aes(x=max(bio_plot$year),y=panel_cap),inherit.aes=FALSE)+
  geom_vline(xintercept=max(d$years),linetype=3,colour="grey50")+theme_bw()+theme(legend.position="bottom")+
  labs(title=paste0("Historical and projected biomass above 101 mm CW"),
    subtitle="Metric tons from empirical weight-at-size; full 95% intervals",
    x="Year",y="Biomass (metric tons)",linetype=NULL)
ggsave(file.path(out,"biomass_gt101mm_metric_tons.png"),p2,width=13,height=5,dpi=220)

assumptions<-data.frame(setting=c("projection_years","trailing_years","uncertainty_draws_requested",
  "uncertainty_draws_accepted","fishing_mortality","recruitment","natural_mortality",
  "hybrid_recruitment_model","hybrid_recruitment_uncertainty","fishery_selectivity",
  "other_size_processes","biomass_threshold","weight_source"),
  value=c(horizon,10,ndraw,accepted,"sampled common calendar year with replacement from final 10 fitted years",
    "sampled common calendar year with replacement from full fitted history; hybrids use projected parental trajectories",
    "sampled common transition year; terminal year and baseline-fixed 2020 excluded",
    "sex + s(log mature snow females, k=3) + s(large Tanner males, large snow males, k=12); lag 5 years",
    "one unconditional GAM coefficient draw per trajectory plus paired male/female historical residual year sampled with replacement",
    "row corresponding to sampled fishing year","paired historical observation-year maturation ogives; hybrid male terminal molt and hybrid female supported-year median remain fixed; fitted growth per parameter draw",
    ">101 mm carapace width","crabpack CALCULATED_WEIGHT_1MM; grams converted to metric tons"))
write.csv(assumptions,file.path(out,"projection_assumptions.csv"),row.names=FALSE)
writeLines(c(paste("Fit:",fit_path),paste("MD5:",tools::md5sum(fit_path)),
 paste("Seed:",seed),paste("Objective reconstructed:",obj$fn(theta)),
 "Recruitment and mortality are labeled by transition start year y; abundance by receiving survey year y+1.",
 paste("First forecast: recruitment",max(d$years),"parents",max(d$years)-d$jr_lag,"receiving survey",max(d$years)+1L),
 "Maturation rows use observation years; sampled jointly across components, independently of recruitment and mortality.",
 "GAM fitted to MAP states once: first-stage uncertainty in GAM fitting is not propagated.",
 paste("Common spring recruitment timing:",isTRUE(d$hybrid_male$data$spring_recruitment)),
 "Spring order: fishery where modeled, growth, recruitment, maturation, remaining 41% mortality."),
 file.path(out,"run_metadata.txt"))
cat("JMR projection written to",out,"using",accepted,"finite uncertainty draws\n")
