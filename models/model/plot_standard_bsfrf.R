library(ggplot2);library(patchwork);library(ggridges);library(GGally)
args<-commandArgs(trailingOnly=TRUE)
argval<-function(key,default){z<-grep(paste0("^",key,"="),args,value=TRUE);if(length(z))sub(paste0("^",key,"="),"",z[1]) else default}
fit_path<-argval("--fit","models/model/fits/model_fit.rds")
out<-argval("--out","models/model/plots/model_diagnostics")
x<-readRDS(fit_path);d<-x$data;r<-x$report
has_uncertainty<-!is.null(x$sdr)&&isTRUE(x$sdr$pdHess)
unc<-if(has_uncertainty)
  data.frame(name=names(x$sdr$value),value=as.numeric(x$sdr$value),se=as.numeric(x$sdr$sd)) else
  data.frame(name=character(),value=numeric(),se=numeric())

# Compatibility summaries for the standard plots. The full four-source
# fishery diagnostics are produced by plot_two_fishery_fits.R.
r$sm_F<-r$sm_F_directed+r$sm_F_incidental
r$tm_F<-r$tm_F_directed+r$tm_F_incidental
r$sm_fish_total<-r$sm_fish_total_directed
r$tm_fish_total<-r$tm_fish_total_directed
r$sm_fish_retain<-r$sm_fish_retain_directed
r$tm_fish_retain<-r$tm_fish_retain_directed

point_value<-function(name){
 if(!is.null(r[[name]]))return(as.numeric(r[[name]]))
 if(name=="sm_Ni_total")return(rowSums(r$sm_Ni))
 if(name=="sm_Nm_total")return(rowSums(r$sm_Nm))
 if(name=="tm_Ni_total")return(rowSums(r$tm_Ni))
 if(name=="tm_Nm_total")return(rowSums(r$tm_Nm))
 if(name=="hm_N_total")return(rowSums(r$hm_N))
 if(name=="sf_Ni_total")return(rowSums(r$sf_Ni))
 if(name=="sf_Nm_total")return(rowSums(r$sf_Nm))
 if(name=="tf_Ni_total")return(rowSums(r$tf_Ni))
 if(name=="tf_Nm_total")return(rowSums(r$tf_Nm))
 if(name=="hf_Ni_total")return(rowSums(r$hf_Ni))
 if(name=="hf_Nm_total")return(rowSums(r$hf_Nm))
 if(name=="sm_male_gt101")return(rowSums((r$sm_Ni+r$sm_Nm)[,d$snow_male$data$sizes>101,drop=FALSE]))
 if(name=="tm_male_gt101")return(rowSums((r$tm_Ni+r$tm_Nm)[,d$tanner_male$data$sizes>101,drop=FALSE]))
 if(name=="hm_male_gt101")return(rowSums(r$hm_N[,d$hybrid_male$data$sizes>101,drop=FALSE]))
 stop("No point-estimate source for ",name)
}
urow<-function(name,year,stock,sex,maturity,process){
 z<-unc[unc$name==name,]
 value<-if(nrow(z))z$value else point_value(name)
 se<-if(nrow(z))z$se else rep(NA_real_,length(value))
 if(length(value)!=length(year))stop("Value length mismatch for ",name)
 data.frame(year=year,stock=stock,sex=sex,maturity=maturity,process=process,
   value=value,se=se,lower=ifelse(is.na(se),NA,pmax(0,value-1.96*se)),
   upper=ifelse(is.na(se),NA,value+1.96*se))}
dir.create(out,recursive=TRUE,showWarnings=FALSE)
stocks<-c("Snow","Tanner","Hybrid");mats<-c("Immature","Mature")
sexcols<-c(Female="#CC79A7",Male="#0072B2");stockcols<-c(Snow="#2C7BB6",Tanner="#D95F02",Hybrid="#7B3294")
# Return a plotting group that changes whenever adjacent observations are not
# consecutive calendar years. This prevents ggplot from bridging missing years.
year_run<-function(z,cols){
 key<-interaction(z[cols],drop=TRUE,lex.order=TRUE);run<-integer(nrow(z))
 for(ii in split(seq_len(nrow(z)),key)){jj<-ii[order(z$year[ii])]
   run[jj]<-cumsum(c(TRUE,diff(z$year[jj])!=1))}
 interaction(key,run,drop=TRUE)
}

streams<-list(
 list(b=d$snow_male,p=r$obsSM,stock="Snow",sex="Male",states=1:2),
 list(b=d$tanner_male,p=r$obsTM,stock="Tanner",sex="Male",states=1:2),
 list(b=d$hybrid_male,p=r$obsHM,stock="Hybrid",sex="Male",states=1),
 list(b=d$snow_female,p=r$obsSF,stock="Snow",sex="Female",states=1:2),
 list(b=d$tanner_female,p=r$obsTF,stock="Tanner",sex="Female",states=1:2),
 list(b=d$hybrid_female,p=r$obsHF,stock="Hybrid",sex="Female",states=1:2))

# Record and visualize the exact annually varying effective sample sizes used
# by the survey size-composition likelihood. Missing values are years with no
# corresponding composition and never contribute to the objective.
weight_rows<-list();iw<-0
add_weights<-function(b,stock,sex,maturity,field){
 data.frame(year=b$data$years,stock=stock,sex=sex,maturity=maturity,
   component="Survey",sample_size_used=b[[field]],cap=d$size_comp_sample_size_cap)
}
for(st in c("Snow","Tanner")){
 b<-d[[paste0(tolower(st),"_male")]]
 iw<-iw+1;weight_rows[[iw]]<-add_weights(b,st,"Male","Immature","survey_comp_weight_imm")
 iw<-iw+1;weight_rows[[iw]]<-add_weights(b,st,"Male","Mature","survey_comp_weight_mat")
}
iw<-iw+1;weight_rows[[iw]]<-add_weights(d$hybrid_male,"Hybrid","Male","Total","survey_comp_weight")
for(st in c("Snow","Tanner","Hybrid")){
 b<-d[[paste0(tolower(st),"_female")]]
 iw<-iw+1;weight_rows[[iw]]<-add_weights(b,st,"Female","Immature","survey_comp_weight_imm")
 iw<-iw+1;weight_rows[[iw]]<-add_weights(b,st,"Female","Mature","survey_comp_weight_mat")
}
weights_used<-do.call(rbind,weight_rows)
write.csv(weights_used,file.path(out,"annual_size_composition_weights_used.csv"),row.names=FALSE)
weights_plot<-subset(weights_used,is.finite(sample_size_used));weights_plot$line_group<-year_run(weights_plot,c("stock","sex","maturity"))
pw<-ggplot(weights_plot,aes(year,sample_size_used,colour=stock,group=line_group))+
 geom_hline(yintercept=d$size_comp_sample_size_cap,linetype="dashed",colour="grey45")+
 geom_line(linewidth=.7)+geom_point(size=.8)+facet_grid(sex+maturity~stock,scales="free_x")+
 scale_colour_manual(values=stockcols)+labs(title="Annual survey size-composition sample sizes used by JMR",
   subtitle=paste0("Actual annual sample size capped at ",d$size_comp_sample_size_cap,"; dashed line is the cap"),
   x=NULL,y="Effective sample size")+theme_bw()+theme(legend.position="none")
ggsave(file.path(out,"annual_size_composition_weights_used.png"),pw,width=14,height=10,dpi=220)

survey<-sizeobs<-sizepred<-list();ii<-0;is<-0
for(q in streams){b<-q$b;dd<-b$data
 for(k in q$states){mat<-if(q$stock=="Hybrid"&&q$sex=="Male")"Mature" else mats[k]
  if(q$sex=="Male"&&q$stock!="Hybrid"){
   oo<-if(k==1)dd$imm_n_obs else dd$mat_n_obs;cv<-if(k==1)dd$sigma_numbers_imm else dd$sigma_numbers_mat
   os<-if(k==1)dd$imm_n_size_obs else dd$mat_n_size_obs
  }else if(q$sex=="Male"){oo<-dd$survey_total;cv<-dd$survey_cv;os<-dd$survey_at_size
  }else{oo<-dd$survey_total[,k];cv<-dd$survey_cv[,k];os<-dd$survey_at_size[,,k]}
  pp<-if(length(dim(q$p))==2)q$p[b$year_map,b$size_map] else q$p[b$year_map,b$size_map,k]
  ii<-ii+1;survey[[ii]]<-data.frame(year=dd$years,stock=q$stock,sex=q$sex,maturity=mat,
    observed=oo,cv=cv,predicted=rowSums(pp))
  use<-which(rowSums(os)>0)
  if(length(use)){is<-is+1;sizeobs[[is]]<-data.frame(year=rep(dd$years[use],each=ncol(os)),size=rep(dd$sizes,length(use)),
    proportion=as.vector(t(sweep(os[use,,drop=FALSE],1,rowSums(os[use,,drop=FALSE]),"/"))),stock=q$stock,sex=q$sex,maturity=mat)
   sizepred[[is]]<-data.frame(year=rep(dd$years[use],each=ncol(os)),size=rep(dd$sizes,length(use)),
    proportion=as.vector(t(sweep(pp[use,,drop=FALSE],1,rowSums(pp[use,,drop=FALSE]),"/"))),stock=q$stock,sex=q$sex,maturity=mat)}
 }}
survey<-do.call(rbind,survey);survey<-subset(survey,year!=2020&is.finite(observed)&observed>0&is.finite(cv)&cv>0)
survey$line_group<-year_run(survey,c("stock","sex","maturity"))
survey$sdlog<-sqrt(log1p(survey$cv^2));survey$lower<-survey$observed/exp(1.96*survey$sdlog);survey$upper<-survey$observed*exp(1.96*survey$sdlog)
write.csv(survey,file.path(out,"survey_abundance_fits.csv"),row.names=FALSE)
spanel<-function(st,ma,title=FALSE,row=FALSE){z<-subset(survey,stock==st&maturity==ma)
 ggplot(z,aes(year,colour=sex,group=line_group))+geom_errorbar(aes(ymin=lower,ymax=upper),width=.12,alpha=.4,linewidth=.3)+
 geom_point(aes(y=observed,shape=sex),size=1.1)+geom_line(aes(y=predicted),linewidth=.8)+scale_colour_manual(values=sexcols)+
 scale_shape_manual(values=c(Female=16,Male=17))+labs(title=if(title)ma else NULL,x=if(st=="Hybrid")"Year" else NULL,y=if(row)st else NULL)+theme_bw()+theme(legend.position="none")}
main<-((spanel("Snow","Immature",TRUE,TRUE)|spanel("Snow","Mature",TRUE))/
 (spanel("Tanner","Immature",FALSE,TRUE)|spanel("Tanner","Mature"))/
 (spanel("Hybrid","Immature",FALSE,TRUE)|spanel("Hybrid","Mature")))+
 plot_annotation(title="Survey-abundance fits: JMR",subtitle="Points and lognormal 95% CV intervals are observations; lines are joint predictions")
sexlegdat<-data.frame(x=rep(1:2,2),y=rep(1:2,each=2),sex=rep(c("Female","Male"),each=2))
sexleg<-ggplot(sexlegdat,aes(x,y,colour=sex,group=sex))+geom_line()+scale_colour_manual(values=sexcols)+coord_cartesian(xlim=c(-1,0))+theme_void()+theme(legend.position="bottom")
ggsave(file.path(out,"survey_abundance_fits.png"),wrap_plots(main,sexleg,ncol=1,heights=c(1,.05)),width=9,height=11,dpi=300)

so<-do.call(rbind,sizeobs);sp<-do.call(rbind,sizepred);med<-aggregate(proportion~stock+sex+maturity+size,sp,median)
allsize<-sort(unique(c(so$size,med$size)));so$sf<-factor(so$size,allsize);med$sf<-factor(med$size,allsize)
commercial<-c("Snow commercial size"=100.6819,"Tanner commercial size"=134.3525);refcols<-c("Snow commercial size"="#D73027","Tanner commercial size"="#7B3294")
zpanel<-function(st,ma,title=FALSE,row=FALSE){a<-subset(so,stock==st&maturity==ma);m<-subset(med,stock==st&maturity==ma)
 p<-ggplot(a,aes(sf,proportion,fill=sex,group=interaction(sf,sex)))+geom_boxplot(position=position_dodge(.72),width=.7,outlier.size=.25,linewidth=.25,alpha=.7)+
 geom_line(data=m,aes(sf,proportion,colour=sex,group=sex),inherit.aes=FALSE,linewidth=.8)+scale_fill_manual(values=c(Female="#E8B4D2",Male="#8CC6E7"))+
 scale_colour_manual(values=c(sexcols,refcols),limits=c("Female","Male",names(refcols)),drop=FALSE)+
 scale_x_discrete(breaks=function(z)z[seq(1,length(z),4)])+labs(title=if(title)ma else NULL,x=if(st=="Hybrid")"Carapace width (mm)" else NULL,y=if(row)st else NULL)+theme_bw()+theme(legend.position="none")
 if(ma=="Mature")p<-p+geom_vline(data=data.frame(x=approx(allsize,seq_along(allsize),xout=commercial,rule=2)$y,reference=names(commercial)),aes(xintercept=x,colour=reference),inherit.aes=FALSE,linetype="dashed",linewidth=.8);p}
main<-((zpanel("Snow","Immature",TRUE,TRUE)|zpanel("Snow","Mature",TRUE))/
 (zpanel("Tanner","Immature",FALSE,TRUE)|zpanel("Tanner","Mature"))/
 (zpanel("Hybrid","Immature",FALSE,TRUE)|zpanel("Hybrid","Mature")))+
 plot_annotation(title="Size-composition fits: JMR",subtitle="Boxplots are annual observations; curves are median joint predictions")
ld<-data.frame(x=rep(1:2,4),y=rep(1:4,each=2),series=rep(c("Female","Male",names(refcols)),each=2));leg<-ggplot(ld,aes(x,y,colour=series,group=series))+geom_line(aes(linetype=grepl("commercial",series)))+scale_colour_manual(values=c(sexcols,refcols))+scale_linetype_manual(values=c(`FALSE`="solid",`TRUE`="dashed"),guide="none")+coord_cartesian(xlim=c(-1,0))+theme_void()+theme(legend.position="bottom")
ggsave(file.path(out,"size_composition_fits.png"),wrap_plots(main,leg,ncol=1,heights=c(1,.055))+
 plot_annotation(caption=if(!is.null(d$hybrid_female$data$plus_group))
   paste0("Female ",max(d$female_sizes),"-mm class is the final plus group; original upper-size observations retained.") else NULL),width=9,height=11,dpi=300)

# Detailed annual size-composition diagnostics: one figure for each of the 11
# modeled species/sex/maturity observation streams, with years in facets.
annual_keep<-c("year","size","proportion","stock","sex","maturity")
annual_size<-rbind(transform(so[,annual_keep],series="Observed"),transform(sp[,annual_keep],series="Predicted"))
annual_size$display_maturity<-as.character(annual_size$maturity)
annual_size$display_maturity[annual_size$stock=="Hybrid"&annual_size$sex=="Male"]<-"Total"
write.csv(annual_size,file.path(out,"size_composition_by_year.csv"),row.names=FALSE)
annual_size_dir<-out
combos<-unique(annual_size[,c("stock","sex","maturity","display_maturity")])
combos<-combos[order(match(combos$stock,stocks),match(combos$sex,c("Male","Female")),match(combos$display_maturity,c("Immature","Mature","Total"))),]
for(i in seq_len(nrow(combos))){q<-combos[i,];st<-as.character(q$stock);sx<-as.character(q$sex);ma<-as.character(q$maturity);dm<-as.character(q$display_maturity);z<-subset(annual_size,stock==st&sex==sx&maturity==ma)
 ny<-length(unique(z$year));nr<-ceiling(ny/5);linecol<-unname(stockcols[st])
 p<-ggplot()+geom_col(data=subset(z,series=="Observed"),aes(size,proportion),fill="grey82",colour="grey68",linewidth=.18,width=4.5)+
  geom_line(data=subset(z,series=="Predicted"),aes(size,proportion),colour=linecol,linewidth=.8)+
  facet_wrap(~year,ncol=5)+scale_y_continuous(expand=expansion(mult=c(0,.08)))+
  labs(title=paste(st,sx,dm,"size compositions by year"),subtitle="Light-grey bars: observed proportions; coloured lines: JMR predictions",x="Carapace width (mm)",y="Proportion")+theme_bw()+theme(strip.text=element_text(size=9))
 fn<-paste0("size_composition_",tolower(st),"_",tolower(sx),"_",tolower(dm),".png")
 ggsave(file.path(annual_size_dir,fn),p,width=15,height=max(7,2.15*nr+2),dpi=220)
}

# Fishery size-composition diagnostics use the same observed-bar and fitted-line
# convention as the annual survey plots. Tanner's second series is total catch;
# snow crab's is the unretained/discard component used by its likelihood.
if(FALSE){ # Replaced by the four-source diagnostics in plot_two_fishery_fits.R
fish_size<-list();ifs<-0
add_fish_size<-function(b,pred,obs,years,stock,catch_type){
  iy<-match(as.integer(years),as.integer(b$data$years));use<-which(rowSums(obs)>0&!is.na(iy))
  if(!length(use))return(NULL)
  oo<-sweep(obs[use,,drop=FALSE],1,rowSums(obs[use,,drop=FALSE]),"/")
  pp<-pred[iy[use],,drop=FALSE]
  pp<-sweep(pp,1,rowSums(pp),"/")
  base<-data.frame(year=rep(as.integer(years[use]),each=ncol(obs)),
    size=rep(b$data$sizes,length(use)),stock=stock,catch_type=catch_type)
  rbind(transform(base,series="Observed",proportion=as.vector(t(oo))),
    transform(base,series="Predicted",proportion=as.vector(t(pp))))
}
for(q in list(c("Snow","snow_male","sm"),c("Tanner","tanner_male","tm"))){
  b<-d[[q[2]]];cd<-b$data$catch
  ry<-if(is.null(cd$ret_size_years))cd$ret_years else cd$ret_size_years
  sy<-if(is.null(cd$second_size_years))cd$ret_years else cd$second_size_years
  ifs<-ifs+1;fish_size[[ifs]]<-add_fish_size(b,r[[paste0(q[3],"_ret")]],
    cd$ret_size,ry,q[1],"Retained")
  ifs<-ifs+1;fish_size[[ifs]]<-add_fish_size(b,r[[paste0(q[3],"_sec")]],
    cd$second_size,sy,q[1],if(q[1]=="Snow")"Discard / unretained" else "Total catch")
}
fish_size<-do.call(rbind,fish_size)
write.csv(fish_size,file.path(out,"fishery_size_composition_by_year.csv"),row.names=FALSE)

fish_obs<-subset(fish_size,series=="Observed");fish_pred<-subset(fish_size,series=="Predicted")
fish_med<-aggregate(proportion~stock+catch_type+size,fish_pred,median)
p_fsc<-ggplot(fish_obs,aes(factor(size),proportion))+
  geom_boxplot(fill="grey82",colour="grey55",outlier.size=.25,linewidth=.25,width=.7)+
  geom_line(data=fish_med,aes(factor(size),proportion,group=1),
    inherit.aes=FALSE,colour="#0072B2",linewidth=.85)+
  facet_wrap(vars(stock,catch_type),ncol=2,scales="free_x")+
  scale_x_discrete(breaks=function(z)z[seq(1,length(z),4)])+
  labs(title="Fishery size-composition fits: JMR",
    subtitle="Boxplots are annual observations; blue curves are median JMR predictions",
    x="Carapace width (mm)",y="Proportion")+theme_bw()
ggsave(file.path(out,"fishery_size_composition_fits.png"),p_fsc,width=13,height=8,dpi=220)

fish_dir<-file.path(out,"fishery_size_composition_by_year")
dir.create(fish_dir,recursive=TRUE,showWarnings=FALSE)
fish_combos<-unique(fish_size[c("stock","catch_type")])
for(i in seq_len(nrow(fish_combos))){st<-fish_combos$stock[i];ct<-fish_combos$catch_type[i]
  z<-subset(fish_size,stock==st&catch_type==ct);ny<-length(unique(z$year));nr<-ceiling(ny/5)
  p<-ggplot()+geom_col(data=subset(z,series=="Observed"),aes(size,proportion),
    fill="grey82",colour="grey68",linewidth=.18,width=4.5)+
    geom_line(data=subset(z,series=="Predicted"),aes(size,proportion),
      colour=unname(stockcols[st]),linewidth=.8)+facet_wrap(~year,ncol=5)+
    scale_y_continuous(expand=expansion(mult=c(0,.08)))+
    labs(title=paste(st,ct,"fishery size compositions by year"),
      subtitle="Light-grey bars: observed proportions; coloured lines: JMR predictions",
      x="Carapace width (mm)",y="Proportion")+theme_bw()+
    theme(strip.text=element_text(size=9))
  fn<-paste0("fishery_size_composition_",tolower(st),"_",
    gsub("[^a-z]+","_",tolower(ct)),".png")
  ggsave(file.path(fish_dir,fn),p,width=15,height=max(7,2.15*nr+2),dpi=220)
}
}

latent<-rbind(
 urow("sm_Ni_total",d$snow_male$data$years,"Snow","Male","Immature","Latent abundance"),urow("sm_Nm_total",d$snow_male$data$years,"Snow","Male","Mature","Latent abundance"),
 urow("tm_Ni_total",d$tanner_male$data$years,"Tanner","Male","Immature","Latent abundance"),urow("tm_Nm_total",d$tanner_male$data$years,"Tanner","Male","Mature","Latent abundance"),
 urow("hm_N_total",d$hybrid_male$data$years,"Hybrid","Male","Mature","Latent abundance"),
 urow("sf_Ni_total",d$snow_female$data$years,"Snow","Female","Immature","Latent abundance"),urow("sf_Nm_total",d$snow_female$data$years,"Snow","Female","Mature","Latent abundance"),
 urow("tf_Ni_total",d$tanner_female$data$years,"Tanner","Female","Immature","Latent abundance"),urow("tf_Nm_total",d$tanner_female$data$years,"Tanner","Female","Mature","Latent abundance"),
 urow("hf_Ni_total",d$hybrid_female$data$years,"Hybrid","Female","Immature","Latent abundance"),urow("hf_Nm_total",d$hybrid_female$data$years,"Hybrid","Female","Mature","Latent abundance"))
write.csv(latent,file.path(out,"latent_population_trajectories.csv"),row.names=FALSE)
latent$line_group<-year_run(latent,c("stock","sex","maturity"))
lpanel<-function(st,ma,title=FALSE,row=FALSE){z<-subset(latent,stock==st&maturity==ma);ggplot(z,aes(year,value,colour=sex,fill=sex,group=line_group))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.18,colour=NA)+geom_line(linewidth=.85)+scale_colour_manual(values=sexcols)+scale_fill_manual(values=sexcols)+labs(title=if(title)ma else NULL,x=if(st=="Hybrid")"Year" else NULL,y=if(row)st else NULL)+theme_bw()+theme(legend.position="none")}
main<-((lpanel("Snow","Immature",TRUE,TRUE)|lpanel("Snow","Mature",TRUE))/
 (lpanel("Tanner","Immature",FALSE,TRUE)|lpanel("Tanner","Mature"))/
 (lpanel("Hybrid","Immature",FALSE,TRUE)|lpanel("Hybrid","Mature")))+
 plot_annotation(title="Latent populations: JMR",subtitle="Lines are fitted totals; ribbons are approximate 95% inverse-Hessian intervals")
ggsave(file.path(out,"latent_population_fits.png"),wrap_plots(main,sexleg,ncol=1,heights=c(1,.05)),width=9,height=11,dpi=300)

# Four-panel comparison of true latent and observed survey abundance of males
# above 101 mm. Model intervals are delta-method Hessian intervals. Survey
# intervals propagate the maturity-specific total-index CVs to the >101 subset.
large_model<-rbind(
 urow("sm_male_gt101",d$snow_male$data$years,"Snow","Male","Total","Abundance >101 mm"),
 urow("tm_male_gt101",d$tanner_male$data$years,"Tanner","Male","Total","Abundance >101 mm"),
 urow("hm_male_gt101",d$hybrid_male$data$years,"Hybrid","Male","Total","Abundance >101 mm"))
names(large_model)[names(large_model)=="value"]<-"abundance"
parent_survey_gt101<-function(b,stock){dd<-b$data;use<-dd$sizes>101
 ai<-rowSums(dd$imm_n_size_obs[,use,drop=FALSE]);am<-rowSums(dd$mat_n_size_obs[,use,drop=FALSE])
 abundance<-ai+am;sd_abs<-sqrt((ai*dd$sigma_numbers_imm)^2+(am*dd$sigma_numbers_mat)^2)
 cv<-sd_abs/pmax(abundance,1);sdlog<-sqrt(log1p(cv^2))
 data.frame(year=dd$years,stock=stock,abundance=abundance,
   lower=abundance/exp(1.96*sdlog),upper=abundance*exp(1.96*sdlog))}
hybrid_survey_gt101<-function(b){dd<-b$data;abundance<-rowSums(dd$survey_at_size[,dd$sizes>101,drop=FALSE])
 sdlog<-sqrt(log1p(dd$survey_cv^2));data.frame(year=dd$years,stock="Hybrid",abundance=abundance,
   lower=abundance/exp(1.96*sdlog),upper=abundance*exp(1.96*sdlog))}
large_survey<-rbind(parent_survey_gt101(d$snow_male,"Snow"),
 parent_survey_gt101(d$tanner_male,"Tanner"),hybrid_survey_gt101(d$hybrid_male))
large_survey<-subset(large_survey,year!=2020&is.finite(abundance)&abundance>0)
large_model$line_group<-year_run(large_model,c("stock"))
large_survey$line_group<-year_run(large_survey,c("stock"))
add_share<-function(z){common<-Reduce(intersect,lapply(split(z$year,z$stock),unique));z<-subset(z,year%in%common)
 z$proportion<-z$abundance/ave(z$abundance,z$year,FUN=sum);z}
large_model$stock<-factor(large_model$stock,levels=stocks)
large_survey$stock<-factor(large_survey$stock,levels=stocks)
model_share<-add_share(large_model);survey_share<-add_share(large_survey)
large_output<-rbind(transform(large_model,source="Model estimate",proportion=NA_real_),
 transform(large_survey,source="Survey observation",se=NA_real_,process=NA_character_,
   sex="Male",maturity="Total",proportion=NA_real_))
large_output$proportion[large_output$source=="Model estimate"]<-
 model_share$proportion[match(paste(large_output$year[large_output$source=="Model estimate"],
   large_output$stock[large_output$source=="Model estimate"]),paste(model_share$year,model_share$stock))]
large_output$proportion[large_output$source=="Survey observation"]<-
 survey_share$proportion[match(paste(large_output$year[large_output$source=="Survey observation"],
   large_output$stock[large_output$source=="Survey observation"]),paste(survey_share$year,survey_share$stock))]
write.csv(large_output,file.path(out,"male_abundance_gt101mm.csv"),row.names=FALSE)
ymax<-max(c(large_model$upper,large_survey$upper),na.rm=TRUE)/1e6*1.03
p_model_abs<-ggplot(large_model,aes(year,abundance/1e6,colour=stock,fill=stock,group=line_group))+
 geom_ribbon(aes(ymin=lower/1e6,ymax=upper/1e6),alpha=.16,colour=NA)+geom_line(linewidth=.85)+
 scale_colour_manual(values=stockcols)+scale_fill_manual(values=stockcols)+coord_cartesian(ylim=c(0,ymax))+
 labs(title="Model estimates",x=NULL,y="Abundance (millions)",colour="Stock",fill="Stock")+
 theme_bw()+theme(legend.position="bottom")+guides(fill="none")
p_survey_abs<-ggplot(large_survey,aes(year,abundance/1e6,colour=stock,group=line_group))+
 geom_errorbar(aes(ymin=lower/1e6,ymax=upper/1e6),width=.15,alpha=.45,linewidth=.3)+
 geom_point(size=1.25)+geom_line(linewidth=.55,alpha=.7)+scale_colour_manual(values=stockcols)+
 coord_cartesian(ylim=c(0,ymax))+labs(title="Survey observations",x=NULL,y=NULL,colour="Stock")+
 theme_bw()+theme(legend.position="none")
share_panel<-function(z,title_y=FALSE)ggplot(z,aes(year,proportion,fill=stock))+
 geom_col(width=1,colour=NA)+scale_fill_manual(values=stockcols)+
 scale_y_continuous(breaks=seq(0,1,.2),expand=c(0,0))+coord_cartesian(ylim=c(0,1))+
 scale_x_continuous(expand=expansion(mult=c(0,.002)))+
 labs(x="Year",y=if(title_y)"Proportion of total" else NULL,fill="Stock")+
 theme_bw()+theme(legend.position="none")
p_model_share<-share_panel(model_share,TRUE);p_survey_share<-share_panel(survey_share,FALSE)
four_large<-(p_model_abs|p_survey_abs)/(p_model_share|p_survey_share)+
 plot_annotation(title="Male crab above 101 mm carapace width",
   subtitle="Top: abundance with approximate 95% intervals; bottom: annual species composition")+
 plot_layout(heights=c(1.05,1))
ggsave(file.path(out,"male_abundance_and_species_proportion_gt101mm.png"),
  four_large,width=15,height=9,dpi=220)

proc<-rbind(
 urow("sm_Mi",d$snow_male$data$years,"Snow","Male","Immature","Natural mortality"),urow("sm_Mm",d$snow_male$data$years,"Snow","Male","Mature","Natural mortality"),urow("sm_rec",d$snow_male$data$years,"Snow","Male","Total","Recruitment"),urow("sm_F",d$snow_male$data$years,"Snow","Male","Total","Fishing mortality"),
 urow("tm_Mi",d$tanner_male$data$years,"Tanner","Male","Immature","Natural mortality"),urow("tm_Mm",d$tanner_male$data$years,"Tanner","Male","Mature","Natural mortality"),urow("tm_rec",d$tanner_male$data$years,"Tanner","Male","Total","Recruitment"),urow("tm_F",d$tanner_male$data$years,"Tanner","Male","Total","Fishing mortality"),
 urow("hm_M",d$hybrid_male$data$years,"Hybrid","Male","Total","Natural mortality"),urow("hm_rec",d$hybrid_male$data$years,"Hybrid","Male","Total","Recruitment"))
for(a in list(c("Snow","sf",d$snow_female$data$years),c("Tanner","tf",d$tanner_female$data$years),c("Hybrid","hf",d$hybrid_female$data$years))){st<-a[1];pre<-a[2];yy<-as.numeric(a[-c(1,2)]);proc<-rbind(proc,urow(paste0(pre,"_Mi"),yy,st,"Female","Immature","Natural mortality"),urow(paste0(pre,"_Mm"),yy,st,"Female","Mature","Natural mortality"),urow(paste0(pre,"_rec"),yy,st,"Female","Total","Recruitment"))}
# Population dynamics use instantaneous hazards. For the cumulative-hazard
# parameterization, report the corresponding annual death probability so the
# requested mortality scale remains interpretable and bounded by zero and one.
if(!is.null(d$mortality_link)&&d$mortality_link=="hazard"){
  use_M<-proc$process=="Natural mortality"
  proc$value[use_M]<-1-exp(-proc$value[use_M])
  proc$lower[use_M]<-1-exp(-proc$lower[use_M])
  proc$upper[use_M]<-1-exp(-proc$upper[use_M])
}
proc$series<-paste(proc$sex,proc$maturity,sep=" / ");pal<-setNames(hcl.colors(length(unique(proc$series)),"Dark 3"),unique(proc$series))
# Normalize recruitment by the maximum fitted value within each species panel.
# Raw values and intervals remain in value/lower/upper; plot_* columns contain
# the normalized recruitment scale and unchanged mortality scales.
proc$panel_max<-NA_real_
for(st in stocks){use<-proc$stock==st&proc$process=="Recruitment";proc$panel_max[use]<-max(proc$value[use],na.rm=TRUE)}
proc$plot_value<-proc$value;proc$plot_lower<-proc$lower;proc$plot_upper<-proc$upper
use<-proc$process=="Recruitment";proc$plot_value[use]<-proc$value[use]/proc$panel_max[use];proc$plot_lower[use]<-proc$lower[use]/proc$panel_max[use];proc$plot_upper[use]<-proc$upper[use]/proc$panel_max[use]
write.csv(proc,file.path(out,"population_processes.csv"),row.names=FALSE)
proc$stock<-factor(proc$stock,levels=stocks);proc$process<-factor(proc$process,levels=c("Recruitment","Natural mortality","Fishing mortality"))
proc$line_group<-year_run(proc,c("stock","process","series"))
process_panel<-function(label,ylab,ylim=NULL,xlab=NULL){
  z<-subset(proc,process==label)
  q<-ggplot(z,aes(year,plot_value,colour=series,fill=series,group=line_group))+
    geom_ribbon(aes(ymin=plot_lower,ymax=plot_upper),alpha=.16,colour=NA)+
    geom_line(linewidth=.75)+facet_grid(.~stock,drop=FALSE)+
    scale_colour_manual(values=pal)+scale_fill_manual(values=pal)+
    labs(x=xlab,y=ylab,colour="Sex / maturity",fill="Sex / maturity")+
    theme_bw()+guides(colour=guide_legend(nrow=2),fill="none")
  if(!is.null(ylim))q<-q+coord_cartesian(ylim=ylim)
  q
}
p<-(process_panel("Recruitment","Relative recruitment",c(0,1.05))/
    process_panel("Natural mortality","Other mortality",c(0,1))/
    process_panel("Fishing mortality","Fishing mortality",xlab="Year"))+
  plot_layout(guides="collect")+
  plot_annotation(title="Population processes: JMR",
    subtitle="Other mortality is always shown on a 0-1 scale; recruitment is normalized within species; ribbons are approximate 95% inverse-Hessian intervals",
    theme=theme(legend.position="bottom"))
ggsave(file.path(out,"population_processes.png"),p,width=15,height=11,dpi=220)

# Pairwise comparison of every estimated natural-mortality trajectory. Values
# are annual death probabilities (rather than instantaneous hazards), matching
# the scale used in population_processes.png. Correlations use overlapping
# calendar years only, which is important because the component time series do
# not all begin in the same year.
mortality_pairs_long<-subset(proc,process=="Natural mortality",
  select=c(year,stock,sex,maturity,value))
mortality_pairs_long$component<-paste(mortality_pairs_long$stock,
  mortality_pairs_long$sex,mortality_pairs_long$maturity,sep=" / ")
mortality_pairs_long<-mortality_pairs_long[,c("year","component","value")]
mortality_pairs_long$component<-factor(mortality_pairs_long$component,
  levels=c("Snow / Male / Immature","Snow / Male / Mature",
    "Tanner / Male / Immature","Tanner / Male / Mature",
    "Hybrid / Male / Total","Snow / Female / Immature",
    "Snow / Female / Mature","Tanner / Female / Immature",
    "Tanner / Female / Mature","Hybrid / Female / Immature",
    "Hybrid / Female / Mature"))
write.csv(mortality_pairs_long,file.path(out,"mortality_pairs_data.csv"),row.names=FALSE)
mortality_wide<-reshape(mortality_pairs_long,idvar="year",timevar="component",
  direction="wide")
names(mortality_wide)<-sub("^value[.]","",names(mortality_wide))
component_names<-levels(mortality_pairs_long$component)
component_names<-component_names[component_names%in%names(mortality_wide)]
mortality_wide<-mortality_wide[,c("year",component_names),drop=FALSE]

cor_rows<-list();ic<-0
for(i in seq_along(component_names))for(j in seq_along(component_names)){
  xi<-mortality_wide[[component_names[i]]];yj<-mortality_wide[[component_names[j]]]
  ok<-is.finite(xi)&is.finite(yj);rr<-if(sum(ok)>2)cor(xi[ok],yj[ok]) else NA_real_
  pp<-if(sum(ok)>3&&sd(xi[ok])>0&&sd(yj[ok])>0)
    cor.test(xi[ok],yj[ok])$p.value else NA_real_
  ic<-ic+1;cor_rows[[ic]]<-data.frame(x=component_names[i],y=component_names[j],
    n=sum(ok),correlation=rr,p_value=pp)
}
write.csv(do.call(rbind,cor_rows),file.path(out,"mortality_pairs_correlations.csv"),row.names=FALSE)

pair_lower<-function(data,mapping,...){
  ggplot(data,mapping)+geom_point(size=.65,alpha=.8,na.rm=TRUE)+
    geom_smooth(method="lm",formula=y~x,colour="#2036B5",fill="#5368FF",
      linewidth=.55,se=TRUE,alpha=.32,na.rm=TRUE)+
    geom_smooth(method="loess",formula=y~x,colour="#D7191C",fill="#FF4D4D",
      linewidth=.55,se=TRUE,alpha=.30,span=.75,na.rm=TRUE)+
    theme_bw(base_size=6)+theme(panel.grid.minor=element_blank())
}
pair_diag<-function(data,mapping,...){
  ggplot(data,mapping)+geom_density(linewidth=.65,colour="black",na.rm=TRUE)+
    theme_bw(base_size=6)+theme(panel.grid.minor=element_blank())
}
pair_upper<-function(data,mapping,...){
  xx<-rlang::eval_tidy(mapping$x,data);yy<-rlang::eval_tidy(mapping$y,data)
  ok<-is.finite(xx)&is.finite(yy);rr<-if(sum(ok)>2)cor(xx[ok],yy[ok]) else NA_real_
  pp<-if(sum(ok)>3&&sd(xx[ok])>0&&sd(yy[ok])>0)
    cor.test(xx[ok],yy[ok])$p.value else NA_real_
  stars<-if(!is.finite(pp))"" else if(pp<.001)"***" else if(pp<.01)"**" else if(pp<.05)"*" else if(pp<.1)"." else ""
  lab<-if(is.finite(rr))sprintf("r = %.2f%s",rr,stars) else "r = NA"
  ggplot(data.frame(x=1,y=1,r=rr),aes(x,y,fill=r))+geom_tile()+
    annotate("text",x=1,y=1,label=lab,size=10.5,fontface="bold",
      colour="grey15",lineheight=.95)+
    scale_fill_gradient2(low="#D98C8C",mid="#F7F7F7",high="#5B9BE6",
      midpoint=0,limits=c(-1,1),na.value="grey90",guide="none")+
    theme_void()+coord_cartesian(expand=FALSE)
}
pair_labels<-c("Snow M-Imm","Snow M-Mat","Tanner M-Imm","Tanner M-Mat",
  "Hybrid M","Snow F-Imm","Snow F-Mat","Tanner F-Imm","Tanner F-Mat",
  "Hybrid F-Imm","Hybrid F-Mat")
pair_labels<-pair_labels[seq_along(component_names)]
p_pairs<-GGally::ggpairs(mortality_wide,columns=match(component_names,names(mortality_wide)),
  columnLabels=pair_labels,lower=list(continuous=pair_lower),
  diag=list(continuous=pair_diag),upper=list(continuous=pair_upper),
  axisLabels="show",progress=FALSE)+
  theme(strip.text=element_text(size=13,face="bold"),axis.text=element_text(size=5),
    axis.title=element_text(size=5))
ggsave(file.path(out,"mortality_pairs.png"),p_pairs,width=28,height=28,dpi=220,
  limitsize=FALSE)

selmed<-function(x)if(is.null(dim(x)))x else apply(x,2,median)
curves<-rbind(data.frame(size=d$snow_male$data$sizes,stock="Snow",sex="Male",process="Survey selectivity",value=selmed(r$sm_sel)),data.frame(size=d$tanner_male$data$sizes,stock="Tanner",sex="Male",process="Survey selectivity",value=selmed(r$tm_sel)),data.frame(size=d$hybrid_male$data$sizes,stock="Hybrid",sex="Male",process="Survey selectivity",value=selmed(r$hm_sel)),data.frame(size=d$hybrid_male$data$sizes,stock="Hybrid",sex="Male",process="Maturation / terminal molt",value=r$hm_ptm),data.frame(size=d$snow_female$data$sizes,stock="Snow",sex="Female",process="Survey selectivity",value=selmed(r$sf_sel)),data.frame(size=d$snow_female$data$sizes,stock="Snow",sex="Female",process="Median maturation / terminal molt",value=apply(r$sf_pm,2,median)),data.frame(size=d$tanner_female$data$sizes,stock="Tanner",sex="Female",process="Survey selectivity",value=selmed(r$tf_sel)),data.frame(size=d$tanner_female$data$sizes,stock="Tanner",sex="Female",process="Median maturation / terminal molt",value=apply(r$tf_pm,2,median)),data.frame(size=d$hybrid_female$data$sizes,stock="Hybrid",sex="Female",process="Survey selectivity",value=selmed(r$hf_sel)),data.frame(size=d$hybrid_female$data$sizes,stock="Hybrid",sex="Female",process="Median maturation / terminal molt",value=apply(r$hf_pm,2,median)))
for(q in list(c("Snow","snow_male","sm"),c("Tanner","tanner_male","tm"))){b<-d[[q[2]]];curves<-rbind(curves,data.frame(size=b$data$sizes,stock=q[1],sex="Male",process="Maturation / terminal molt",value=apply(b$data$prop_term_molt,2,median)),data.frame(size=b$data$sizes,stock=q[1],sex="Male",process="Total fishery selectivity",value=apply(r[[paste0(q[3],"_fish_total")]],2,median)),data.frame(size=b$data$sizes,stock=q[1],sex="Male",process="Retention selectivity",value=apply(r[[paste0(q[3],"_fish_retain")]],2,median)))}
curves<-subset(curves,process!="Survey selectivity");curves$era<-"Other process"
survey_curve_rows<-do.call(rbind,lapply(list(
  c("Snow","Male","snow_male","sm_sel"),c("Tanner","Male","tanner_male","tm_sel"),
  c("Hybrid","Male","hybrid_male","hm_sel"),c("Snow","Female","snow_female","sf_sel"),
  c("Tanner","Female","tanner_female","tf_sel"),c("Hybrid","Female","hybrid_female","hf_sel")),function(q){
  b<-d[[q[3]]];mat<-r[[q[4]]]
  do.call(rbind,lapply(c("1975-1981","1982-1988","1989-present"),function(er){
    use<-switch(er,"1975-1981"=b$data$years<=1981,
      "1982-1988"=b$data$years>=1982&b$data$years<=1988,
      "1989-present"=b$data$years>=1989)
    data.frame(size=b$data$sizes,stock=q[1],sex=q[2],process="Survey selectivity",
      value=apply(mat[use,,drop=FALSE],2,median),era=er)}))}))
curves<-rbind(curves,survey_curve_rows)
curves$era<-factor(curves$era,levels=c("Other process","1989-present","1982-1988","1975-1981"))
p<-ggplot(curves,aes(size,value,colour=stock,linetype=era,
  group=interaction(stock,era)))+geom_line(linewidth=.9)+facet_grid(sex~process)+
  scale_colour_manual(values=stockcols)+
  scale_linetype_manual(values=c("1975-1981"="dotdash","1982-1988"="dashed",
    "1989-present"="solid","Other process"="solid"),
    breaks=c("1975-1981","1982-1988","1989-present"))+ 
  coord_cartesian(ylim=c(0,1))+labs(title="Selectivity and maturation: JMR",
    x="Carapace width (mm)",y="Probability",colour="Stock",linetype="Survey era")+
  theme_bw()+theme(legend.position="bottom")+
  guides(linetype=guide_legend(order=2),colour=guide_legend(order=1))
write.csv(curves,file.path(out,"selectivity_and_maturation.csv"),row.names=FALSE)
ggsave(file.path(out,"selectivity_and_maturation.png"),p,width=13,height=8,dpi=220)

fish_curves<-do.call(rbind,lapply(list(c("Snow","snow_male","sm"),
  c("Tanner","tanner_male","tm")),function(q){
  b<-d[[q[2]]];active<-rowSums(b$fish_total)>0
  bind<-function(x,process){z<-expand.grid(year=b$data$years[active],size=b$data$sizes)
    z$value<-as.vector(r[[x]][active,,drop=FALSE]);z$stock<-q[1];z$process<-process;z}
  rbind(bind(paste0(q[3],"_fish_total"),"Total fishery selectivity"),
    bind(paste0(q[3],"_fish_retain"),"Retention selectivity"))}))
write.csv(fish_curves,file.path(out,"annual_fishery_selectivity.csv"),row.names=FALSE)
p_fish<-ggplot(fish_curves,aes(size,value,group=year,colour=year))+
  geom_line(alpha=.55,linewidth=.55)+facet_grid(stock~process,scales="free_x")+
  scale_colour_viridis_c()+coord_cartesian(ylim=c(0,1))+
  labs(title="Estimated annual fishery selectivity",x="Carapace width (mm)",
    y="Probability",colour="Year")+theme_bw()
ggsave(file.path(out,"annual_directed_fishery_selectivity.png"),p_fish,width=12,height=7,dpi=220)

survey_sel_era<-do.call(rbind,lapply(list(
  c("Snow","Male","snow_male","sm_sel"),c("Tanner","Male","tanner_male","tm_sel"),
  c("Hybrid","Male","hybrid_male","hm_sel"),c("Snow","Female","snow_female","sf_sel"),
  c("Tanner","Female","tanner_female","tf_sel"),c("Hybrid","Female","hybrid_female","hf_sel")),function(q){
  b<-d[[q[3]]];mat<-r[[q[4]]]
  do.call(rbind,lapply(c("1975-1981","1982-1988","1989-present"),function(er){
    use<-switch(er,"1975-1981"=b$data$years<=1981,
      "1982-1988"=b$data$years>=1982&b$data$years<=1988,
      "1989-present"=b$data$years>=1989)
    data.frame(size=b$data$sizes,stock=q[1],sex=q[2],era=er,
      value=apply(mat[use,,drop=FALSE],2,median))}))}))
write.csv(survey_sel_era,file.path(out,"survey_selectivity_by_era.csv"),row.names=FALSE)
p_sse<-ggplot(survey_sel_era,aes(size,value,colour=stock,linetype=era))+
  geom_line(linewidth=.9)+facet_wrap(~sex,ncol=1,scales="free_x")+
  scale_colour_manual(values=stockcols)+coord_cartesian(ylim=c(0,1))+
  labs(title="Survey selectivity by footprint era",x="Carapace width (mm)",
    y="Probability",colour="Stock",linetype="Survey era")+theme_bw()+
  theme(legend.position="bottom")
ggsave(file.path(out,"survey_selectivity_by_era.png"),p_sse,width=10,height=8,dpi=220)

annual_maturity<-do.call(rbind,lapply(list(Snow=d$snow_female$data,
  Tanner=d$tanner_female$data,Hybrid=d$hybrid_female$data),function(dd){
  z<-expand.grid(observation_year=dd$years,size=dd$sizes)
  z$year<-z$observation_year-1L
  z$probability<-as.vector(dd$p_mature_by_year)
  z$midpoint<-rep(dd$maturity_midpoint_by_year,times=length(dd$sizes))
  src<-dd$maturity_ogive_source
  if(ncol(src)!=length(dd$sizes))src<-src[,match(dd$sizes,as.numeric(colnames(src))),drop=FALSE]
  z$source<-as.vector(src);z}))
annual_maturity$stock<-rep(names(list(Snow=d$snow_female$data,
  Tanner=d$tanner_female$data,Hybrid=d$hybrid_female$data)),
  sapply(list(d$snow_female$data,d$tanner_female$data,d$hybrid_female$data),
    function(dd)length(dd$years)*length(dd$sizes)))
write.csv(annual_maturity,file.path(out,"annual_female_maturation_probability.csv"),row.names=FALSE)
midline<-unique(annual_maturity[c("year","stock","midpoint")])
p_mat<-ggplot(annual_maturity,aes(size,year,fill=probability))+geom_tile()+
  geom_line(data=midline,aes(midpoint,year),inherit.aes=FALSE,colour="white",linewidth=.65)+
  facet_wrap(~stock,nrow=1,scales="free_x")+scale_fill_viridis_c(limits=c(0,1))+
  labs(title="Annual female maturation / terminal-molt probability at size",
    subtitle="Parental annual ogives; hybrid fixed median (n > 10). White line: first 50% crossing",
    x="Carapace width (mm)",y="Process year (observation year minus one)",fill="Probability")+theme_bw()
ggsave(file.path(out,"annual_female_maturation_probability.png"),p_mat,width=14,height=7,dpi=220)

# Direct verification of the three parental-envelope constraints.
constraints<-rbind(
 data.frame(size=d$hybrid_male$data$sizes,sex="Male",process="Terminal molt",curve="Snow parent",value=r$hm_ptm_parent_snow),
 data.frame(size=d$hybrid_male$data$sizes,sex="Male",process="Terminal molt",curve="Tanner parent",value=r$hm_ptm_parent_tanner),
 data.frame(size=d$hybrid_male$data$sizes,sex="Male",process="Terminal molt",curve="Hybrid",value=r$hm_ptm),
 data.frame(size=d$hybrid_male$data$sizes,sex="Male",process="Survey selectivity",curve="Snow parent",value=selmed(r$hm_sel_parent_snow)),
 data.frame(size=d$hybrid_male$data$sizes,sex="Male",process="Survey selectivity",curve="Tanner parent",value=selmed(r$hm_sel_parent_tanner)),
 data.frame(size=d$hybrid_male$data$sizes,sex="Male",process="Survey selectivity",curve="Hybrid",value=selmed(r$hm_sel)),
 data.frame(size=d$hybrid_female$data$sizes,sex="Female",process="Survey selectivity",curve="Snow parent",value=selmed(r$hf_sel_parent_snow)),
 data.frame(size=d$hybrid_female$data$sizes,sex="Female",process="Survey selectivity",curve="Tanner parent",value=selmed(r$hf_sel_parent_tanner)),
 data.frame(size=d$hybrid_female$data$sizes,sex="Female",process="Survey selectivity",curve="Hybrid",value=selmed(r$hf_sel)))
write.csv(constraints,file.path(out,"hybrid_parental_curve_constraints.csv"),row.names=FALSE)
p<-ggplot(constraints,aes(size,value,colour=curve,linetype=curve))+geom_line(linewidth=.95)+facet_grid(sex~process,scales="free_x")+scale_colour_manual(values=c("Snow parent"=unname(stockcols["Snow"]),"Tanner parent"=unname(stockcols["Tanner"]),"Hybrid"=unname(stockcols["Hybrid"])))+scale_linetype_manual(values=c("Snow parent"="dashed","Tanner parent"="dashed","Hybrid"="solid"))+coord_cartesian(ylim=c(0,1))+labs(title="Hybrid curves constrained by parental estimates",subtitle="Survey: continuous parental stacked logistics; terminal molt: projected parental medians",x="Carapace width (mm)",y="Probability",colour=NULL,linetype=NULL)+theme_bw()+theme(legend.position="bottom")
ggsave(file.path(out,"hybrid_parental_curve_constraints.png"),p,width=12,height=7,dpi=220)
weights<-rbind(urow("wtm",NA,"Hybrid","Male","Total","Terminal-molt snow weight"),urow("wms",NA,"Hybrid","Male","Total","Survey-selectivity snow weight"),urow("wfs",NA,"Hybrid","Female","Total","Survey-selectivity snow weight"))
write.csv(weights,file.path(out,"hybrid_parental_blend_weights.csv"),row.names=FALSE)

class_unc<-function(name,size,sex){z<-unc[unc$name==name,]
 value<-if(nrow(z))z$value else point_value(name)
 se<-if(nrow(z))z$se else rep(NA_real_,length(value))
 if(length(value)!=length(size))stop("Classification value length mismatch for ",name)
 data.frame(size=size,sex=sex,probability=value,se=se,
   lower=ifelse(is.na(se),NA,pmax(0,value-1.96*se)),
   upper=ifelse(is.na(se),NA,pmin(1,value+1.96*se)))}
classdat<-rbind(class_unc("pf",d$female_sizes,"Female"),class_unc("pm",d$male_sizes,"Male"))
write.csv(classdat,file.path(out,"classification_ogives.csv"),row.names=FALSE)
p<-ggplot(classdat,aes(size,probability,colour=sex,fill=sex))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.18,colour=NA)+geom_line(linewidth=1)+geom_vline(xintercept=90,linetype="dashed",colour="grey45")+facet_wrap(~sex,ncol=1)+scale_colour_manual(values=sexcols)+scale_fill_manual(values=sexcols)+coord_cartesian(ylim=c(0,1))+labs(title="Probability a true hybrid is assigned a parental label",subtitle="Shared-sex power curve; misclassification is structurally zero at and above 90 mm",x="Carapace width (mm)",y="Probability",colour="Sex",fill="Sex")+theme_bw()+theme(legend.position="none")
ggsave(file.path(out,"classification_ogives.png"),p,width=9,height=8,dpi=220)
component_names <- c("Snow male block","Tanner male block","Hybrid male block","Snow female block","Tanner female block","Hybrid female block","Survey abundance","Survey size composition","Classification constraint (no prior)","Shared-process priors","Parental growth data","Hybrid growth-blend priors","Parental-recruitment coefficient prior","Shared-mortality driver prior","BSFRF survey-selectivity priors")
component_names<-head(component_names,length(r$components))
write.csv(data.frame(component=component_names,value=r$components),file.path(out,"objective_components.csv"),row.names=FALSE)
oc<-read.csv(file.path(out,"objective_components.csv"));p<-ggplot(oc,aes(reorder(component,value),value))+geom_col(fill="#4C78A8")+coord_flip()+labs(title="Objective contributions: JMR",x=NULL,y="Negative log likelihood / penalty")+theme_bw();ggsave(file.path(out,"objective_components.png"),p,width=10,height=7,dpi=220)

# Annual transfers implied by the single sex-specific classification ogives.
reallocation_rows<-function(observed,parent,sex,destination){
 number<-rowSums(observed-parent)
 parent_label_abundance<-rowSums(observed)
 data.frame(year=d$years,sex=sex,destination=destination,number=number,
   parent_label_abundance=parent_label_abundance,
   proportion=ifelse(parent_label_abundance>0,number/parent_label_abundance,NA_real_))
}
# Denominator includes hybrids carrying the parental label before reassignment.
realloc<-rbind(reallocation_rows(r$obsSM,r$SM,"Male","Snow"),
 reallocation_rows(r$obsTM,r$TM,"Male","Tanner"),
 reallocation_rows(r$obsSF,r$SF,"Female","Snow"),
 reallocation_rows(r$obsTF,r$TF,"Female","Tanner"))
realloc$line_group<-year_run(realloc,c("sex","destination"))
write.csv(realloc,file.path(out,"annual_reallocation.csv"),row.names=FALSE)
reallocation_panel<-function(metric,title,ylab){
 ggplot(realloc,aes(year,.data[[metric]],colour=destination,group=line_group))+
   geom_line(linewidth=.8)+facet_wrap(~sex,ncol=1,scales="free_y")+
   scale_colour_manual(values=stockcols[c("Snow","Tanner")])+
   labs(title=title,x="Year",y=ylab,colour="Observed as")+theme_bw()+
   theme(legend.position="bottom")
}
p<-(reallocation_panel("number","Number reallocated","Number")|
 reallocation_panel("proportion","Proportion of parent stock reallocated","Proportion")+
   scale_y_continuous(labels=scales::label_percent(accuracy=.1)))+
 plot_layout(guides="collect")+
 plot_annotation(title="True hybrids assigned parental labels",
   subtitle="Proportions use modeled parent-label abundance before reassignment",
   theme=theme(legend.position="bottom"))
ggsave(file.path(out,"annual_reallocation.png"),p,width=13,height=7,dpi=220)

# Sum maturity states before computing the fraction at each year and size.
reallocation_size_rows<-function(observed,parent,sizes,sex,stock){
 abundance<-apply(observed,c(1,2),sum)
 number<-apply(observed-parent,c(1,2),sum)
 stopifnot(identical(dim(abundance),c(length(d$years),length(sizes))))
 z<-expand.grid(year=d$years,size=sizes)
 z$sex<-sex;z$stock<-stock
 z$number<-as.vector(number)
 z$parent_label_abundance<-as.vector(abundance)
 z$proportion<-ifelse(z$parent_label_abundance>0,
   z$number/z$parent_label_abundance,NA_real_)
 z
}
realloc_size<-rbind(
 reallocation_size_rows(r$obsSM,r$SM,d$male_sizes,"Male","Snow"),
 reallocation_size_rows(r$obsTM,r$TM,d$male_sizes,"Male","Tanner"),
 reallocation_size_rows(r$obsSF,r$SF,d$female_sizes,"Female","Snow"),
 reallocation_size_rows(r$obsTF,r$TF,d$female_sizes,"Female","Tanner"))
write.csv(realloc_size,file.path(out,"reallocation_at_size.csv"),row.names=FALSE)
p<-ggplot(realloc_size,aes(size,year,fill=proportion))+
 geom_tile()+facet_grid(sex~stock)+
 scale_fill_viridis_c(labels=scales::label_percent(),na.value="grey90")+
 scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+
 labs(title="Proportion of parent stock reallocated to hybrids by size and year",
   subtitle="Reallocated number / modeled parent-label abundance; maturity states combined",
   x="Carapace width (mm)",y="Year",fill="Proportion",
   caption=paste("Gray: no modeled parent-label abundance.",
     if(!is.null(d$hybrid_female$data$plus_group))paste0("Female ",max(d$female_sizes),"-mm class is the final plus group; original upper-size observations retained.") else ""))+
 theme_bw()+theme(panel.grid=element_blank(),legend.position="right")
ggsave(file.path(out,"reallocation_at_size.png"),p,width=12,height=9,dpi=220)

growth<-function(G,z,stock,sex)data.frame(premolt=rep(z,each=length(z)),postmolt=rep(z,length(z)),probability=as.vector(G),stock=stock,sex=sex)
gd<-rbind(growth(r$sm_G,d$snow_male$data$sizes,"Snow","Male"),growth(r$tm_G,d$tanner_male$data$sizes,"Tanner","Male"),growth(r$hm_G,d$hybrid_male$data$sizes,"Hybrid","Male"),growth(r$sf_G,d$snow_female$data$sizes,"Snow","Female"),growth(r$tf_G,d$tanner_female$data$sizes,"Tanner","Female"),growth(r$hf_G,d$hybrid_female$data$sizes,"Hybrid","Female"));write.csv(gd,file.path(out,"estimated_growth_matrices.csv"),row.names=FALSE);p<-ggplot(gd,aes(premolt,postmolt,fill=probability))+geom_raster()+facet_grid(sex~stock,scales="free")+scale_fill_viridis_c()+labs(title="Jointly estimated growth-transition matrices: JMR",subtitle="Hybrid regression parameters are bounded between sex-specific parental estimates",x="Premolt width (mm)",y="Postmolt width (mm)",fill="Probability")+theme_bw()+theme(legend.position="bottom");ggsave(file.path(out,"estimated_growth_matrices.png"),p,width=13,height=8,dpi=220)
growth_parameters<-rbind(data.frame(stock=c("Snow","Tanner","Hybrid"),sex="Male",alpha=c(r$sm_growth_alpha,r$tm_growth_alpha,r$hm_growth_alpha),beta=c(r$sm_growth_beta,r$tm_growth_beta,r$hm_growth_beta),sd=c(r$sm_growth_sd,r$tm_growth_sd,r$hm_growth_sd)),data.frame(stock=c("Snow","Tanner","Hybrid"),sex="Female",alpha=c(r$sf_growth_alpha,r$tf_growth_alpha,r$hf_growth_alpha),beta=c(r$sf_growth_beta,r$tf_growth_beta,r$hf_growth_beta),sd=c(r$sf_growth_sd,r$tf_growth_sd,r$hf_growth_sd)));write.csv(growth_parameters,file.path(out,"growth_parameters.csv"),row.names=FALSE)
recruitment_allocation<-rbind(data.frame(stock=c("Snow","Tanner","Hybrid"),sex="Male",bin1=c(r$sm_rec_prop[1],r$tm_rec_prop[1],r$hm_rec_prop[1]),bin2=c(r$sm_rec_prop[2],r$tm_rec_prop[2],r$hm_rec_prop[2]),bin3=c(r$sm_rec_prop[3],r$tm_rec_prop[3],r$hm_rec_prop[3])),data.frame(stock=c("Snow","Tanner","Hybrid"),sex="Female",bin1=c(r$sf_rec_prop[1],r$tf_rec_prop[1],r$hf_rec_prop[1]),bin2=c(r$sf_rec_prop[2],r$tf_rec_prop[2],r$hf_rec_prop[2]),bin3=c(r$sf_rec_prop[3],r$tf_rec_prop[3],r$hf_rec_prop[3])));write.csv(recruitment_allocation,file.path(out,"recruitment_allocation.csv"),row.names=FALSE)

# Growth-data fits (four parental stock-by-sex panels) beside the sex-specific
# parental envelopes and the estimated hybrid mean relationship. A shared
# alpha/beta blend weight guarantees the hybrid line lies inside the envelope.
growth_obs<-rbind(
 data.frame(stock="Snow",sex="Male",premolt=d$snow_male$data$growth_premolt,postmolt=d$snow_male$data$growth_postmolt),
 data.frame(stock="Tanner",sex="Male",premolt=d$tanner_male$data$growth_premolt,postmolt=d$tanner_male$data$growth_postmolt),
 data.frame(stock="Snow",sex="Female",premolt=d$snow_female$data$growth$premolt,postmolt=d$snow_female$data$growth$postmolt),
 data.frame(stock="Tanner",sex="Female",premolt=d$tanner_female$data$growth$premolt,postmolt=d$tanner_female$data$growth$postmolt))
parent_growth<-growth_parameters[growth_parameters$stock!="Hybrid",]
fit_curves<-do.call(rbind,lapply(seq_len(nrow(parent_growth)),function(i){q<-parent_growth[i,];oo<-growth_obs[growth_obs$stock==q$stock&growth_obs$sex==q$sex,];z<-seq(min(oo$premolt),max(oo$premolt),length.out=150);mu<-q$alpha+q$beta*z;data.frame(stock=q$stock,sex=q$sex,premolt=z,mean=mu,lower=mu-1.96*q$sd,upper=mu+1.96*q$sd)}))
p_growth_fits<-ggplot(growth_obs,aes(premolt,postmolt))+geom_point(alpha=.35,size=1.2,colour="#444444")+geom_ribbon(data=fit_curves,aes(x=premolt,ymin=lower,ymax=upper,fill=stock),inherit.aes=FALSE,alpha=.18,colour=NA)+geom_line(data=fit_curves,aes(x=premolt,y=mean,colour=stock),inherit.aes=FALSE,linewidth=.9)+geom_abline(slope=1,intercept=0,linetype=3,colour="grey55")+facet_grid(sex~stock,scales="free")+scale_colour_manual(values=stockcols[c("Snow","Tanner")])+scale_fill_manual(values=stockcols[c("Snow","Tanner")])+labs(title="Parental growth-data fits",subtitle="Lines: fitted means; ribbons: +/- 1.96 residual SD",x="Premolt width (mm)",y="Postmolt width (mm)")+theme_bw()+theme(legend.position="none")
bound_curves<-lapply(c("Male","Female"),function(sx){q<-growth_parameters[growth_parameters$sex==sx,];zlim<-if(sx=="Male")range(d$hybrid_male$data$sizes) else range(d$hybrid_female$data$sizes);z<-seq(zlim[1],zlim[2],length.out=200);curves<-do.call(rbind,lapply(seq_len(nrow(q)),function(i)data.frame(sex=sx,stock=q$stock[i],premolt=z,postmolt=q$alpha[i]+q$beta[i]*z)));wide<-reshape(curves,idvar=c("sex","premolt"),timevar="stock",direction="wide");wide$lower<-pmin(wide$postmolt.Snow,wide$postmolt.Tanner);wide$upper<-pmax(wide$postmolt.Snow,wide$postmolt.Tanner);attr(wide,"reshapeWide")<-NULL;list(curves=curves,bounds=wide)})
growth_lines<-do.call(rbind,lapply(bound_curves,`[[`,"curves"));growth_bounds<-do.call(rbind,lapply(bound_curves,`[[`,"bounds"))
p_growth_bounds<-ggplot()+geom_ribbon(data=growth_bounds,aes(premolt,ymin=lower,ymax=upper),fill="grey75",alpha=.55)+geom_line(data=growth_lines,aes(premolt,postmolt,colour=stock,linetype=stock),linewidth=1)+facet_wrap(~sex,ncol=1,scales="free")+scale_colour_manual(values=stockcols)+scale_linetype_manual(values=c(Snow="solid",Tanner="solid",Hybrid="longdash"))+labs(title="Parentally bounded hybrid growth",subtitle="Grey band is the fitted parental mean-growth envelope",x="Premolt width (mm)",y="Expected postmolt width (mm)",colour=NULL,linetype=NULL)+theme_bw()+theme(legend.position="bottom")
write.csv(growth_obs,file.path(out,"growth_fit_observations.csv"),row.names=FALSE);write.csv(fit_curves,file.path(out,"parental_growth_fit_curves.csv"),row.names=FALSE);write.csv(growth_lines,file.path(out,"parental_and_hybrid_growth_curves.csv"),row.names=FALSE)
p_growth_combined<-p_growth_fits+p_growth_bounds+plot_layout(widths=c(1.45,1));ggsave(file.path(out,"growth_data_fits_and_hybrid_bounds.png"),p_growth_combined,width=16,height=9,dpi=220)

# Annual latent abundance-at-size ridgelines. Absolute abundance is retained in
# the CSV; ridge height is normalized within each stock/sex/maturity panel so
# shapes remain visible despite orders-of-magnitude differences among stocks.
ridge_long<-function(mat,years,sizes,stock,sex,maturity,source="Latent")
 data.frame(year=rep(years,each=length(sizes)),size=rep(sizes,length(years)),
   number=as.vector(t(mat)),stock=stock,sex=sex,maturity=maturity,source=source)
latent_size<-rbind(
 ridge_long(r$sm_Ni,d$snow_male$data$years,d$snow_male$data$sizes,"Snow","Male","Immature"),ridge_long(r$sm_Nm,d$snow_male$data$years,d$snow_male$data$sizes,"Snow","Male","Mature"),
 ridge_long(r$tm_Ni,d$tanner_male$data$years,d$tanner_male$data$sizes,"Tanner","Male","Immature"),ridge_long(r$tm_Nm,d$tanner_male$data$years,d$tanner_male$data$sizes,"Tanner","Male","Mature"),
 ridge_long(r$hm_N,d$hybrid_male$data$years,d$hybrid_male$data$sizes,"Hybrid","Male","Total"),
 ridge_long(r$sf_Ni,d$snow_female$data$years,d$snow_female$data$sizes,"Snow","Female","Immature"),ridge_long(r$sf_Nm,d$snow_female$data$years,d$snow_female$data$sizes,"Snow","Female","Mature"),
 ridge_long(r$tf_Ni,d$tanner_female$data$years,d$tanner_female$data$sizes,"Tanner","Female","Immature"),ridge_long(r$tf_Nm,d$tanner_female$data$years,d$tanner_female$data$sizes,"Tanner","Female","Mature"),
 ridge_long(r$hf_Ni,d$hybrid_female$data$years,d$hybrid_female$data$sizes,"Hybrid","Female","Immature"),ridge_long(r$hf_Nm,d$hybrid_female$data$years,d$hybrid_female$data$sizes,"Hybrid","Female","Mature"))
latent_size$panel<-interaction(latent_size$stock,latent_size$sex,latent_size$maturity,drop=TRUE)
mx<-ave(latent_size$number,latent_size$panel,FUN=max);latent_size$ridge_height<-latent_size$number/pmax(mx,1)
latent_size$stock<-factor(latent_size$stock,levels=stocks)
write.csv(latent_size,file.path(out,"latent_numbers_at_size_by_year.csv"),row.names=FALSE)
p<-ggplot(latent_size,aes(size,year,height=ridge_height,group=year,fill=sex))+geom_ridgeline(stat="identity",scale=4,alpha=.6,colour="white",linewidth=.18)+facet_grid(sex+maturity~stock,scales="free",space="free_y")+scale_fill_manual(values=sexcols)+scale_y_continuous(breaks=seq(1980,2025,5))+labs(title="Latent abundance at size by year",subtitle="Ridge heights are normalized within panel; CSV contains absolute abundance",x="Carapace width (mm)",y="Year",fill="Sex")+theme_bw()+theme(legend.position="bottom",panel.spacing.y=unit(.55,"lines"))
ggsave(file.path(out,"latent_numbers_at_size_ridges.png"),p,width=15,height=16,dpi=220)

# Hybrid observation-space comparison. HM/HF are true-hybrid abundance after
# survey selectivity but before the classification layer removes animals into
# observed snow or Tanner labels. This is neither raw latent abundance nor the
# fitted observed-hybrid stream.
hm_obs<-d$hybrid_male$data;hm_use<-which(hm_obs$survey_observed)
hybrid_obs<-rbind(
 ridge_long(hm_obs$survey_at_size[hm_use,,drop=FALSE],hm_obs$years[hm_use],hm_obs$sizes,"Hybrid","Male","Total","Observed hybrid label"),
 ridge_long(r$HM[d$hybrid_male$year_map[hm_use],d$hybrid_male$size_map,drop=FALSE],hm_obs$years[hm_use],hm_obs$sizes,"Hybrid","Male","Total","Classification-corrected expected"))
hf_obs<-d$hybrid_female$data;hf_use<-which(hf_obs$survey_observed)
for(k in 1:2){ma<-mats[k];hybrid_obs<-rbind(hybrid_obs,
 ridge_long(hf_obs$survey_at_size[hf_use,,k,drop=FALSE][,,1],hf_obs$years[hf_use],hf_obs$sizes,"Hybrid","Female",ma,"Observed hybrid label"),
 ridge_long(r$HF[d$hybrid_female$year_map[hf_use],d$hybrid_female$size_map,k,drop=FALSE][,,1],hf_obs$years[hf_use],hf_obs$sizes,"Hybrid","Female",ma,"Classification-corrected expected"))}
hybrid_obs$panel<-interaction(hybrid_obs$sex,hybrid_obs$maturity,drop=TRUE)
mx<-ave(hybrid_obs$number,hybrid_obs$panel,FUN=max);hybrid_obs$ridge_height<-hybrid_obs$number/pmax(mx,1)
write.csv(hybrid_obs,file.path(out,"hybrid_observed_vs_classification_corrected_at_size.csv"),row.names=FALSE)
hybrid_ridge_plot<-function(sx){z<-subset(hybrid_obs,sex==sx)
 ggplot()+geom_ridgeline(data=subset(z,source=="Observed hybrid label"),aes(size,year,height=ridge_height,group=year),stat="identity",scale=7,fill="grey75",colour="grey40",alpha=.75,linewidth=.3)+geom_ridgeline(data=subset(z,source=="Classification-corrected expected"),aes(size,year,height=ridge_height,group=year),stat="identity",scale=7,fill=NA,colour=unname(stockcols["Hybrid"]),linewidth=.85)+facet_wrap(~maturity,nrow=1,scales="free_x")+scale_y_continuous(breaks=seq(1980,2025,5))+labs(title=paste("Observed and classification-corrected hybrid abundance at size:",tolower(sx)),subtitle="Grey: observed hybrid label; purple: true-hybrid expectation after survey selectivity, before misclassification; ridge overlap scale = 7",x="Carapace width (mm)",y="Year")+theme_bw()+theme(panel.spacing.x=unit(1.2,"lines"))}
ggsave(file.path(out,"hybrid_observed_vs_classification_corrected_male_ridges.png"),hybrid_ridge_plot("Male"),width=14,height=15,dpi=250)
ggsave(file.path(out,"hybrid_observed_vs_classification_corrected_female_ridges.png"),hybrid_ridge_plot("Female"),width=18,height=13,dpi=250)

# Cross-correlation between hybrid recruitment and survey-space abundance.
# Positive lag means survey abundance occurs after recruitment by that many
# years. Both the fitted observed-label stream and the true-hybrid expectation
# before classification are shown.
hseries<-rbind(
 data.frame(year=d$hybrid_male$data$years,sex="Male",recruitment=r$hm_rec,
   raw_small=ifelse(d$hybrid_male$data$survey_observed,rowSums(d$hybrid_male$data$survey_at_size[,seq_len(min(3,length(d$hybrid_male$data$sizes))),drop=FALSE]),NA_real_),
   fitted_observed=rowSums(r$obsHM[d$hybrid_male$year_map,,drop=FALSE]),
   corrected_survey=rowSums(r$HM[d$hybrid_male$year_map,,drop=FALSE])),
 data.frame(year=d$hybrid_female$data$years,sex="Female",recruitment=r$hf_rec,
   raw_small=ifelse(d$hybrid_female$data$survey_observed,apply(d$hybrid_female$data$survey_at_size[,seq_len(min(3,length(d$hybrid_female$data$sizes))),,drop=FALSE],1,sum),NA_real_),
   fitted_observed=apply(r$obsHF[d$hybrid_female$year_map,,,drop=FALSE],1,sum),
   corrected_survey=apply(r$HF[d$hybrid_female$year_map,,,drop=FALSE],1,sum)))
lag_cor<-function(rec,survey,lag){n<-length(rec)
 if(lag>=0){i<-seq_len(n-lag);j<-i+lag}else{j<-seq_len(n+lag);i<-j-lag}
 ok<-is.finite(rec[i])&is.finite(survey[j]);c(n=sum(ok),correlation=cor(rec[i][ok],survey[j][ok]))}
cc<-list();ii<-0
for(sx in c("Female","Male")){z<-hseries[hseries$sex==sx,]
 for(src in c("raw_small","fitted_observed","corrected_survey"))for(lag in -10:10){ii<-ii+1;q<-lag_cor(z$recruitment,z[[src]],lag);cc[[ii]]<-data.frame(sex=sx,survey_series=src,lag=lag,n=q["n"],correlation=q["correlation"])}}
cc<-do.call(rbind,cc);cc$survey_series<-factor(cc$survey_series,levels=c("raw_small","fitted_observed","corrected_survey"),labels=c("Raw observed recruitment-bin abundance","Fitted observed-hybrid label","Classification-corrected survey abundance"))
write.csv(cc,file.path(out,"hybrid_survey_recruitment_crosscorrelation.csv"),row.names=FALSE)
p<-ggplot(cc,aes(lag,correlation,colour=survey_series))+geom_hline(yintercept=0,colour="grey70")+geom_vline(xintercept=0,colour="grey70")+geom_line(linewidth=.9)+geom_point(size=1.5)+facet_wrap(~sex,ncol=1)+scale_colour_manual(values=c("Raw observed recruitment-bin abundance"="black","Fitted observed-hybrid label"="grey55","Classification-corrected survey abundance"=unname(stockcols["Hybrid"])))+scale_x_continuous(breaks=-10:10)+coord_cartesian(ylim=c(-1,1))+labs(title="Cross-correlation of hybrid recruitment and survey abundance",subtitle="Raw observations use the 27.5-37.5 mm bins; positive lag means abundance follows recruitment",x="Lag (years)",y="Pearson correlation",colour=NULL)+theme_bw()+theme(legend.position="bottom")
ggsave(file.path(out,"hybrid_survey_recruitment_crosscorrelation.png"),p,width=10,height=8,dpi=220)

# Compare estimated recruitment with abundance in the three modeled recruitment
# bins. All series are normalized within sex and series to emphasize timing;
# absolute values are retained in the CSV.
small_series<-list();ii<-0
add_small<-function(year,sex,series,value,observed=FALSE){data.frame(year=year,sex=sex,series=series,value=as.numeric(value),observed=observed)}
bm<-d$hybrid_male;dm<-bm$data;um<-which(dm$survey_observed);sm<-seq_len(min(3,length(dm$sizes)))
small_series[[1]]<-add_small(dm$years[um],"Male","Observed hybrid label",rowSums(dm$survey_at_size[um,sm,drop=FALSE]),TRUE)
small_series[[2]]<-add_small(dm$years,"Male","Fitted observed-hybrid label",rowSums(r$obsHM[bm$year_map,bm$size_map[sm],drop=FALSE]))
small_series[[3]]<-add_small(dm$years,"Male","Classification-corrected survey abundance",rowSums(r$HM[bm$year_map,bm$size_map[sm],drop=FALSE]))
small_series[[4]]<-add_small(dm$years,"Male","Estimated recruitment",r$hm_rec)
bf<-d$hybrid_female;df<-bf$data;uf<-which(df$survey_observed);sf<-seq_len(min(3,length(df$sizes)))
small_series[[5]]<-add_small(df$years[uf],"Female","Observed hybrid label",apply(df$survey_at_size[uf,sf,,drop=FALSE],1,sum),TRUE)
small_series[[6]]<-add_small(df$years,"Female","Fitted observed-hybrid label",apply(r$obsHF[bf$year_map,bf$size_map[sf],,drop=FALSE],1,sum))
small_series[[7]]<-add_small(df$years,"Female","Classification-corrected survey abundance",apply(r$HF[bf$year_map,bf$size_map[sf],,drop=FALSE],1,sum))
small_series[[8]]<-add_small(df$years,"Female","Estimated recruitment",r$hf_rec)
small<-do.call(rbind,small_series);small$key<-interaction(small$sex,small$series,drop=TRUE);small$normalized<-small$value/ave(small$value,small$key,FUN=function(z)max(z,na.rm=TRUE));small$key<-NULL
small$line_group<-year_run(small,c("sex","series"))
write.csv(small,file.path(out,"hybrid_small_abundance_vs_recruitment.csv"),row.names=FALSE)
seriescols<-c("Observed hybrid label"="black","Fitted observed-hybrid label"="grey55","Classification-corrected survey abundance"=unname(stockcols["Hybrid"]),"Estimated recruitment"="#009E73")
p<-ggplot(small,aes(year,normalized,colour=series,group=line_group))+geom_line(data=subset(small,!observed),linewidth=.9)+geom_point(data=subset(small,observed),size=1.8)+facet_wrap(~sex,ncol=1)+scale_colour_manual(values=seriescols)+coord_cartesian(ylim=c(0,1.05))+labs(title="Small hybrid survey abundance and estimated recruitment",subtitle="Survey abundance sums 27.5-37.5 mm bins; each sex/series is divided by its own maximum",x="Year",y="Normalized value",colour=NULL)+theme_bw()+theme(legend.position="bottom")
ggsave(file.path(out,"hybrid_small_abundance_vs_recruitment.png"),p,width=12,height=8,dpi=220)

if(!is.null(r$common_rec)) {
  drivers<-rbind(data.frame(year=d$years[-length(d$years)],process="Recruitment",value=r$common_rec),
                 data.frame(year=d$years[-length(d$years)],process="Natural mortality",value=r$common_M))
  write.csv(drivers,file.path(out,"shared_annual_process_drivers.csv"),row.names=FALSE)
  p<-ggplot(drivers,aes(year,value,colour=process))+geom_hline(yintercept=0,colour="grey70")+
    geom_line(linewidth=.85)+facet_wrap(~process,ncol=1,scales="free_y")+
    scale_colour_manual(values=c("Recruitment"="#0072B2","Natural mortality"="#D55E00"))+
    labs(title="Estimated shared annual process drivers",x="Year",y="Standardized annual effect",colour=NULL)+
    theme_bw()+theme(legend.position="bottom")
  ggsave(file.path(out,"shared_annual_process_drivers.png"),p,width=10,height=7,dpi=220)
  rec_labels<-c("Snow / Male","Tanner / Male","Hybrid / Male","Snow / Female","Tanner / Female","Hybrid / Female")
  M_labels<-c("Snow / Male / Immature","Snow / Male / Mature","Tanner / Male / Immature","Tanner / Male / Mature","Hybrid / Male / Total","Snow / Female / Immature","Snow / Female / Mature","Tanner / Female / Immature","Tanner / Female / Mature","Hybrid / Female / Immature","Hybrid / Female / Mature")
  variance<-rbind(data.frame(process="Recruitment",series=rec_labels,shared_fraction=r$rec_variance_explained,
    residual_sd=exp(r$log_rec_sigma),shared_loading=exp(r$log_rec_loading)),
    data.frame(process="Natural mortality",series=M_labels,shared_fraction=r$M_variance_explained,
    residual_sd=exp(r$log_M_sigma),shared_loading=exp(r$log_M_loading)))
  write.csv(variance,file.path(out,"shared_process_variance.csv"),row.names=FALSE)
  p<-ggplot(variance,aes(shared_fraction,reorder(series,shared_fraction),fill=process))+geom_col()+
    facet_wrap(~process,ncol=1,scales="free_y")+scale_x_continuous(limits=c(0,1))+
    scale_fill_manual(values=c("Recruitment"="#0072B2","Natural mortality"="#D55E00"))+
    labs(title="Fraction of annual-process variance attributed to the shared driver",x="Shared fraction",y=NULL,fill=NULL)+
    theme_bw()+theme(legend.position="bottom")
  ggsave(file.path(out,"shared_process_variance_explained.png"),p,width=11,height=9,dpi=220)
}
hd<-x$hessian_diagnostics
write.csv(data.frame(metric=c("optimizer_code","objective","max_gradient","hessian_positive_definite","minimum_hessian_eigenvalue","hessian_condition_number","fit_type","uncertainty_type"),value=c(x$fit$convergence,x$fit$objective,x$max_gradient,if(is.null(hd))NA else hd$positive_definite,if(is.null(hd))NA else hd$minimum_eigenvalue,if(is.null(hd))NA else hd$condition_number,"joint penalized mode / MAP",if(has_uncertainty)x$uncertainty_type else "Point estimates; Hessian uncertainty not available")),file.path(out,"fit_diagnostics.csv"),row.names=FALSE)
cat("JMR diagnostics written to",out,"\n")


