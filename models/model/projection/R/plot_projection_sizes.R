# size_out contains intervals aggregated within each draw, preserving covariance.
sz<-size_out
sz$stock<-factor(sz$stock,c("Snow","Tanner","Hybrid"))
sz$size_group<-factor(sz$size_group,c("<76 mm","76-101 mm",">101 mm"))
for(sx in c("Female","Male")){
 z<-subset(sz,sex==sx&component=="Total")
 p<-ggplot(z,aes(year,estimate/1e6,group=interaction(stock,size_group,period)))+
  geom_ribbon(aes(ymin=lower/1e6,ymax=upper/1e6),fill="#7B3294",alpha=.18)+
  geom_line(aes(linetype=period),linewidth=.65)+
  facet_wrap(vars(stock,size_group),ncol=3,scales="free_y")+
  geom_vline(xintercept=max(d$years),linetype=3,colour="grey50")+
  theme_bw()+theme(legend.position="bottom")+
  labs(title=paste(sx,"latent abundance by size group"),
   subtitle=paste(accepted,"trial draws; medians and 95% intervals; all maturity states combined"),
   x="Survey year",y="Abundance (millions)",linetype=NULL)
 ggsave(file.path(out,paste0("latent_abundance_by_size_group_",tolower(sx),".png")),p,width=14,height=10,dpi=220)
}
large<-subset(sz,sex=="Male"&size_group==">101 mm")
write.csv(large,file.path(out,"male_abundance_gt101mm.csv"),row.names=FALSE)
large$state<-ifelse(large$component %in% c("Immature","Nonterminal"),"Immature / nonterminal",
 ifelse(large$component=="Total","Total","Mature / terminal"))
large$state<-factor(large$state,c("Total","Immature / nonterminal","Mature / terminal"))
p<-ggplot(large,aes(year,estimate/1e6,group=interaction(stock,state,period)))+
 geom_ribbon(aes(ymin=lower/1e6,ymax=upper/1e6),fill="#2166AC",alpha=.18)+
 geom_line(aes(linetype=period),linewidth=.7)+
 facet_wrap(vars(stock,state),ncol=3,scales="free_y")+
 geom_vline(xintercept=max(d$years),linetype=3,colour="grey50")+
 theme_bw()+theme(legend.position="bottom")+
 labs(title="Large male latent abundance: >101 mm carapace width",
  subtitle=paste(accepted,"trial draws; medians and full 95% intervals; common size cutoff across stocks"),
  x="Survey year",y="Abundance (millions)",linetype=NULL)
ggsave(file.path(out,"male_abundance_gt101mm.png"),p,width=14,height=10,dpi=220)
