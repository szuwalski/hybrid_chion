maturation_curves<-function(r,d){
 list(sm=d$snow_male$data$prop_term_molt,tm=d$tanner_male$data$prop_term_molt,
 hm=matrix(r$hm_ptm,length(d$hybrid_male$data$years),length(r$hm_ptm),byrow=TRUE),
 sf=r$sf_pm,tf=r$tf_pm,hf=r$hf_pm)
}

select_maturation_sizes<-function(r,d,spec){
 curves<-maturation_curves(r,d);audit<-list()
 for(nm in names(spec)){
  s<-spec[[nm]];b<-d[[paste0(tolower(s$stock),"_",tolower(s$sex))]]$data
  z<-curves[[nm]];eligible<-matrix(s$years!=2020 & seq_along(s$years)>1,nrow(z),ncol(z))
  support<-"Model-input years; specimen counts unavailable"
  if(s$sex=="Female"){
   eligible<-eligible & is.finite(b$raw_maturity_probability) & is.finite(b$maturity_specimen_n) & b$maturity_specimen_n>10
   support<-"Observation years with >10 specimens at this size"
  }
  rows<-lapply(seq_along(s$sizes),function(j){
   v<-z[eligible[,j]&is.finite(z[,j]),j]
   data.frame(key=nm,stock=s$stock,sex=s$sex,size_mm=s$sizes[j],
    supported_years=length(v),mean_probability=if(length(v))mean(v) else NA_real_,
    annual_sd=if(length(v)>1)sd(v) else NA_real_,
    minimum=if(length(v))min(v) else NA_real_,maximum=if(length(v))max(v) else NA_real_,
    support_rule=support)
  })
  tab<-do.call(rbind,rows)
  ok<-with(tab,supported_years>=10 & is.finite(annual_sd) & mean_probability>=.2 & mean_probability<=.8)
  tab$eligible_intermediate_size<-ok
  if(!any(ok))stop("No supported intermediate maturation size for ",nm)
  fixed<-nm %in% c("hm","hf")
  cand<-which(ok)
  pick<-if(fixed)cand[which.min(abs(tab$mean_probability[cand]-.5))] else cand[which.max(tab$annual_sd[cand])]
  if(nm=="sm")pick<-match(77.5,tab$size_mm)
  tab$selected<-seq_len(nrow(tab))==pick
  tab$selection_rule<-if(fixed)"Fixed annual curve: closest to probability 0.5" else "Maximum annual probability SD; mean 0.2-0.8; >=10 supported years"
  if(nm=="sm")tab$selection_rule<-"User-selected 75-mm area: 77.5-mm model bin"
  tab$annual_model_variation<-if(fixed)"Fixed across years" else "Annual input ogive"
  audit[[nm]]<-tab
 }
 do.call(rbind,audit)
}

plot_projection_maturation<-function(process_out,maturation_selection,out,accepted,maturation_audit){
 z<-subset(process_out,process=="Terminal molt probability")
 z$panel<-paste(z$stock,z$sex,z$component,sep=" / ")
 sel<-maturation_selection
 panels<-paste(sel$stock,sel$sex,paste0(sel$size_mm," mm"),sep=" / ")
 z$panel<-factor(z$panel,levels=panels)
 write.csv(z,file.path(out,"terminal_molt_probability_trajectories.csv"),row.names=FALSE)
 p<-ggplot(z,aes(year,estimate,group=interaction(panel,period)))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#E69F00",alpha=.22)+
  geom_line(aes(linetype=period),linewidth=.7)+facet_wrap(~panel,ncol=3)+
  coord_cartesian(ylim=c(0,1))+scale_y_continuous(breaks=seq(0,1,.25))+
  geom_vline(xintercept=min(z$year[z$period=="Projected"]),linetype=3,colour="grey50")+
  theme_bw()+theme(legend.position="bottom")+
  labs(title="Terminal-molt / maturation probability at representative sizes",
   subtitle=paste(accepted,"trial draws; historical ogive resampling; hybrid curves fixed across years"),
   caption="Post-growth size; conditional maturation probability, not proportion mature in the population.\nRibbons describe projection variability; observation uncertainty in fixed ogive inputs is not estimated.",
   x="Process year (observation-year ogive minus one)",y="Probability",linetype=NULL)
 ggsave(file.path(out,"terminal_molt_probability_trajectories.png"),p,width=14,height=8,dpi=220)
 a<-maturation_audit;a$panel<-paste(a$stock,a$sex,sep=" / ")
 a$panel<-factor(a$panel,levels=unique(a$panel))
 psize<-ggplot(a,aes(size_mm,annual_sd))+
  geom_line(colour="grey60",na.rm=TRUE)+
  geom_point(aes(colour=eligible_intermediate_size),na.rm=TRUE)+
  geom_point(data=subset(a,selected),shape=21,size=3.5,fill="#E69F00",colour="black")+
  facet_wrap(~panel,ncol=3)+theme_bw()+theme(legend.position="bottom")+
  labs(title="Selecting a representative maturation size",
   subtitle="Orange circles: selected size; hybrid selections use the fixed curve's midpoint",
   x="Post-growth carapace width (mm)",y="Across-year SD of probability",colour="Eligible intermediate size")
 ggsave(file.path(out,"terminal_molt_size_selection.png"),psize,width=14,height=8,dpi=220)
}
