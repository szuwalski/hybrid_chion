# Replay every fitted transition through the projection functions. This checks
# agreement with the saved fitted states, including the final transition.
check_projection_transitions<-function(r,d){
 blocks<-c(sm="snow_male",tm="tanner_male",hm="hybrid_male",
           sf="snow_female",tf="tanner_female",hf="hybrid_female")
 checks<-list()
 for(pre in names(blocks)){
  b<-d[[blocks[[pre]]]];n<-length(b$data$years)-1L
  getr<-function(suffix)r[[paste0(pre,"_",suffix)]]
  hi<-getr(if(pre=="hm")"U" else "Ni")
  hm<-getr(if(pre=="hm")"T" else "Nm")
  if(pre %in% c("sm","tm")){
   if(!is.null(getr("F_directed"))){
    pr<-project_parent_male_two_fishery(hi[1,],hm[1,],getr("G"),getr("Mi"),getr("Mm"),
      getr("F_directed"),getr("F_incidental"),getr("rec"),getr("rec_prop"),
      getr("fish_total_directed"),getr("fish_retain_directed"),
      getr("fish_total_incidental"),getr("fish_retain_incidental"),
      b$data$prop_term_molt[-1,,drop=FALSE],n)
   }else pr<-project_parent_male(b,hi[1,],hm[1,],getr("G"),getr("Mi"),getr("Mm"),
      getr("F"),getr("rec"),getr("rec_prop"),getr("fish_total"),
      getr("fish_retain"),b$data$prop_term_molt[-1,,drop=FALSE],n)
  }else if(pre=="hm"){
   pr<-project_hybrid_male(hi[1,],hm[1,],getr("G"),getr("M"),
     getr("rec"),getr("rec_prop"),getr("ptm"),n,spring=isTRUE(b$data$spring_recruitment))
  }else{
   pr<-project_female(hi[1,],hm[1,],getr("G"),getr("Mi"),getr("Mm"),
     getr("rec"),getr("rec_prop"),getr("pm")[-1,,drop=FALSE],n,spring=isTRUE(b$data$spring_recruitment))
  }
  err<-max(abs(pr$i-hi)/pmax(1,abs(hi)),abs(pr$m-hm)/pmax(1,abs(hm)))
  growth_error<-max(abs(colSums(getr("G"))-1))
  stopifnot(is.finite(err),err<1e-8,growth_error<1e-8,
    all(is.finite(pr$i)),all(pr$i>=0),all(pr$m>=0))
  checks[[pre]]<-data.frame(component=blocks[[pre]],transitions=n,
    size_bins=ncol(hi),max_relative_state_error=err,
    max_growth_column_sum_error=growth_error)
 }
 do.call(rbind,checks)
}

# Recruitment-only pulse: growth is deliberately nonidentity and F is large.
# New recruits must be neither grown nor fished, and must survive only 0.41 M.
check_spring_recruitment_pulse<-function(){
 G<-matrix(c(0,1,0,0,0,1,1,0,0),3,3)
 rp<-c(.75,.20,.05);pm<-c(.1,.4,.8);rec<-100;Mi<-.3;Mm<-.6
 expected_i<-rec*rp*(1-pm)*exp(-.41*Mi)
 expected_m<-rec*rp*pm*exp(-.41*Mm)
 b<-list(data=list(catch=list(discard_survival=.2)))
 parent<-project_parent_male(b,rep(0,3),rep(0,3),G,Mi,Mm,10,rec,rp,
   matrix(1,1,3),matrix(1,1,3),matrix(pm,1,3),1)
 female<-project_female(rep(0,3),rep(0,3),G,Mi,Mm,rec,rp,matrix(pm,1,3),1,spring=TRUE)
 hybrid<-project_hybrid_male(rep(0,3),rep(0,3),G,Mi,rec,rp,pm,1,spring=TRUE)
 errors<-c(parent=max(abs(parent$i[2,]-expected_i),abs(parent$m[2,]-expected_m)),
   female=max(abs(female$i[2,]-expected_i),abs(female$m[2,]-expected_m)),
   hybrid=max(abs(hybrid$i[2,]-expected_i),
     abs(hybrid$m[2,]-rec*rp*pm*exp(-.41*Mi))))
 stopifnot(all(errors<1e-12))
 data.frame(block=names(errors),recruitment_pulse_error=as.numeric(errors))
}
