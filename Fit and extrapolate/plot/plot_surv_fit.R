dfHazEst<-read.csv("dfhazest_nivpfs_3y.csv")
dfHazEst<-dfHazEst[,2:79]
ltHaz<-read.csv("lthaz_nivpfs_3y.csv")
ltHaz<-ltHaz[,2:31]
ltHaz_out<-read.csv("lthaz_nivpfs_6.5y.csv")
ltHaz_out<-ltHaz_out[,2:31]
md<-c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","genf","FP1-1","FP1-2",
      "FP2-1","FP2-2","RCS1","RCS2","RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")

#calculate survival over time
###################################################################
lthaz_plot<-data.frame(ltHaz$Time,ltHaz[,11:30])
dfsurv_lthaz<-as.data.frame(array(dim=c(36,21)))
colnames(dfsurv_lthaz)<-c("Time",md)
dfsurv_lthaz$Time<-lthaz_plot$ltHaz.Time
for (i in 1:20) {
  temp_dfhaz<-data.frame(lthaz_plot$ltHaz.Time,lthaz_plot[,i+1])
  colnames(temp_dfhaz)<-c("time","haz")
  temp_dfsurv<-temp_dfhaz %>% 
    dplyr::arrange(time) %>% 
    dplyr::mutate(cumhaz = cumsum(haz)) %>% 
    dplyr::mutate(survProp = exp(-1*cumhaz))
  dfsurv_lthaz[,i+1]<-temp_dfsurv$survProp
}

dfFigSurv_lthaz = dfsurv_lthaz %>%
  gather(key = "Model", value = "survProp", -Time) %>% mutate(Model = factor(Model))
###################################################################

dfhar_out<-ltHaz_out[,11:30]
colnames(dfhar_out) <- md
dfhar_out <- cbind(data.frame(ltHaz_out$Time), dfhar_out)
dfhar_out <- rename(dfhar_out,"Time"="ltHaz_out.Time")

dfsurv_out<-as.data.frame(array(dim=c(78,21)))
colnames(dfsurv_out)<-c("Time",md)
dfsurv_out$Time<-dfhar_out$Time
###  20 model hazard  ###
for (i in 1:20) {
  temp_dfhaz<-data.frame(dfhar_out$Time,dfhar_out[,i+1])
  colnames(temp_dfhaz)<-c("time","haz")
  temp_dfsurv_out<-temp_dfhaz %>% 
    dplyr::arrange(time) %>% 
    dplyr::mutate(cumhaz = cumsum(haz)) %>% 
    dplyr::mutate(survProp = exp(-1*cumhaz))
  dfsurv_out[,i+1]<-temp_dfsurv_out$survProp
}
dfFig_outSurv = dfsurv_out %>%
  gather(key = "Model", value = "survProp", -Time) %>% mutate(Model = factor(Model))
# ###################################################################
# ##plot##
#1
# f_surv_nivos3y= ggplot() +
#   geom_line(data=dfFigSurv_lthaz, aes(x=Time, y=survProp, group=Model, colour=Model), size=0.6,show.legend = FALSE) +
#   geom_line(data = ltHaz,aes(x=Time,y=surv),size=1)+
#   scale_color_discrete(name="Model")+
#   expand_limits(y=c(0,1),x=c(0,3)) +
#   scale_x_continuous(breaks = c(seq(from=0, to=3,by = 0.5))) +
#   ylab("Overall survival") +
#   xlab("Time(months)") +
#   guides(color = guide_legend(ncol = 1))  +
#   theme(legend.position = "bottom") +
#   theme_bw() +
#   labs(title="Survival plot of NI group",
#        subtitle = "OS 3Ys")
# f_surv_nivos3y
# ##plot##
# #2
# f_surv_nivos6.5y= ggplot() +
#   geom_line(data=dfFig_outSurv, aes(x=Time, y=survProp, group=Model, colour=Model), size=0.5) +
#   geom_line(data = ltHaz_out,aes(x=Time,y=surv),size=1)+
#   scale_color_discrete(name="Model")+
#   expand_limits(y=c(0,1),x=c(0,6.5)) +
#   scale_x_continuous(breaks = c(seq(from=0, to=6.5,by = 0.5))) +
#   ylab("Overall survival") +
#   xlab("Time(months)") +
#   guides(color = guide_legend(ncol = 1))  +
#   theme(legend.key.width = unit(0.4,"cm"),legend.key.height = unit(0.4,"cm")) +
#   labs(title="Survival plot of NI group",
#        subtitle = "OS 6.5Ys")
# f_surv_nivos6.5y
###################################################################
# ##plot##
# #3
# f_surv_ipios3y= ggplot() +
#   geom_line(data=dfFigSurv_lthaz, aes(x=Time, y=survProp, group=Model, colour=Model), size=0.6,show.legend = FALSE) +
#   geom_line(data = ltHaz,aes(x=Time,y=surv),size=1)+
#   scale_color_discrete(name="Model")+
#   expand_limits(y=c(0,1),x=c(0,3)) +
#   scale_x_continuous(breaks = c(seq(from=0, to=3,by = 0.5))) +
#   ylab("Overall survival") +
#   xlab("Time(months)") +
#   guides(color = guide_legend(ncol = 1))  +
#   theme(legend.position = "bottom") +
#   theme_bw() +
#   labs(title="Survival plot of I group",
#        subtitle = "OS 3Ys")
# f_surv_ipios3y
# ##plot##
# #4
# f_surv_ipios6.5y= ggplot() +
#   geom_line(data=dfFig_outSurv, aes(x=Time, y=survProp, group=Model, colour=Model), size=0.5) +
#   geom_line(data = ltHaz_out,aes(x=Time,y=surv),size=1)+
#   scale_color_discrete(name="Model")+
#   expand_limits(y=c(0,1),x=c(0,6.5)) +
#   scale_x_continuous(breaks = c(seq(from=0, to=6.5,by = 0.5))) +
#   ylab("Overall survival") +
#   xlab("Time(months)") +
#   guides(color = guide_legend(ncol = 1))  +
#   theme(legend.key.width = unit(0.4,"cm"),legend.key.height = unit(0.4,"cm")) +
#   labs(title="Survival plot of I group",
#        subtitle = "OS 6.5Ys")
# f_surv_ipios6.5y
###################################################################
##plot##
## 5
# f_surv_nivpfs3y= ggplot() +
#   geom_line(data=dfFigSurv_lthaz, aes(x=Time, y=survProp, group=Model, colour=Model), size=0.6,show.legend = FALSE) +
#   geom_line(data = ltHaz,aes(x=Time,y=surv),size=1)+
#   scale_color_discrete(name="Model")+
#   expand_limits(y=c(0,1),x=c(0,3)) +
#   scale_x_continuous(breaks = c(seq(from=0, to=3,by = 0.5))) +
#   ylab("Progression-free survival") +
#   xlab("Time(months)") +
#   guides(color = guide_legend(ncol = 1))  +
#   theme(legend.position = "bottom") +
#   theme_bw() +
#   labs(title="Survival plot of NI group",
#        subtitle = "PFS 3Ys")
# f_surv_nivpfs3y
# ##plot##
# #6
# f_surv_nivpfs6.5y= ggplot() +
#   geom_line(data=dfFig_outSurv, aes(x=Time, y=survProp, group=Model, colour=Model), size=0.5) +
#   geom_line(data = ltHaz_out,aes(x=Time,y=surv),size=1)+
#   scale_color_discrete(name="Model")+
#   expand_limits(y=c(0,1),x=c(0,6.5)) +
#   scale_x_continuous(breaks = c(seq(from=0, to=6.5,by = 0.5))) +
#   ylab("Progression-free survival") +
#   xlab("Time(months)") +
#   guides(color = guide_legend(ncol = 1))  +
#   theme(legend.key.width = unit(0.4,"cm"),legend.key.height = unit(0.4,"cm")) +
#   labs(title="Survival plot of NI group",
#        subtitle = "PFS 6.5Ys")
# f_surv_nivpfs6.5y
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
dfHazEst<-read.csv("dfhazest_ipipfs_3y.csv")
dfHazEst<-dfHazEst[,2:79]
ltHaz<-read.csv("lthaz_ipipfs_3y.csv")
ltHaz<-ltHaz[,2:30]
ltHaz_out<-read.csv("lthaz_ipipfs_6.5y.csv")
ltHaz_out<-ltHaz_out[,2:30]
md<-c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","FP1-1","FP1-2",
      "FP2-1","FP2-2","RCS1","RCS2","RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")

#calculate survival over time
###################################################################
lthaz_plot<-data.frame(ltHaz$Time,ltHaz[,11:29])
dfsurv_lthaz<-as.data.frame(array(dim=c(36,20)))
colnames(dfsurv_lthaz)<-c("Time",md)
dfsurv_lthaz$Time<-lthaz_plot$ltHaz.Time
for (i in 1:19) {
  temp_dfhaz<-data.frame(lthaz_plot$ltHaz.Time,lthaz_plot[,i+1])
  colnames(temp_dfhaz)<-c("time","haz")
  temp_dfsurv<-temp_dfhaz %>% 
    dplyr::arrange(time) %>% 
    dplyr::mutate(cumhaz = cumsum(haz)) %>% 
    dplyr::mutate(survProp = exp(-1*cumhaz))
  dfsurv_lthaz[,i+1]<-temp_dfsurv$survProp
}

dfFigSurv_lthaz = dfsurv_lthaz %>%
  gather(key = "Model", value = "survProp", -Time) %>% mutate(Model = factor(Model))
###################################################################

dfhar_out<-ltHaz_out[,11:29]
colnames(dfhar_out) <- md
dfhar_out <- cbind(data.frame(ltHaz_out$Time), dfhar_out)
dfhar_out <- rename(dfhar_out,"Time"="ltHaz_out.Time")

dfsurv_out<-as.data.frame(array(dim=c(78,20)))
colnames(dfsurv_out)<-c("Time",md)
dfsurv_out$Time<-dfhar_out$Time
###  20 model hazard  ###
for (i in 1:19) {
  temp_dfhaz<-data.frame(dfhar_out$Time,dfhar_out[,i+1])
  colnames(temp_dfhaz)<-c("time","haz")
  temp_dfsurv_out<-temp_dfhaz %>% 
    dplyr::arrange(time) %>% 
    dplyr::mutate(cumhaz = cumsum(haz)) %>% 
    dplyr::mutate(survProp = exp(-1*cumhaz))
  dfsurv_out[,i+1]<-temp_dfsurv_out$survProp
}
dfFig_outSurv = dfsurv_out %>%
  gather(key = "Model", value = "survProp", -Time) %>% mutate(Model = factor(Model))
###################################################################
##plot##
#7
f_surv_ipipfs3y= ggplot() +
  geom_line(data=dfFigSurv_lthaz, aes(x=Time, y=survProp, group=Model, colour=Model), size=0.6,show.legend = FALSE) +
  geom_line(data = ltHaz,aes(x=Time,y=surv),size=1)+
  scale_color_discrete(name="Model")+
  expand_limits(y=c(0,1),x=c(0,3)) +
  scale_x_continuous(breaks = c(seq(from=0, to=3,by = 0.5))) +
  ylab("Progression-free survival") +
  xlab("Time(months)") +
  guides(color = guide_legend(ncol = 1))  +
  theme(legend.key.width = unit(0.4,"cm"),legend.key.height = unit(0.4,"cm")) +
  labs(title="Survival plot of I group",
       subtitle = "PFS 3Ys")
f_surv_ipipfs3y
##plot##
#8
f_surv_ipipfs6.5y= ggplot() +
  geom_line(data=dfFig_outSurv, aes(x=Time, y=survProp, group=Model, colour=Model), size=0.5) +
  geom_line(data = ltHaz_out,aes(x=Time,y=surv),size=1)+
  scale_color_discrete(name="Model")+
  expand_limits(y=c(0,1),x=c(0,6.5)) +
  scale_x_continuous(breaks = c(seq(from=0, to=6.5,by = 0.5))) +
  ylab("Progression-free survival") +
  xlab("Time(months)") +
  guides(color = guide_legend(ncol = 1))  +
  theme(legend.key.width = unit(0.4,"cm"),legend.key.height = unit(0.4,"cm")) +
  labs(title="Survival plot of I group",
       subtitle = "PFS 6.5Ys")
f_surv_ipipfs6.5y


surv_fit<-grid.arrange(f_surv_nivos3y,f_surv_nivos6.5y,
                       f_surv_ipios3y,f_surv_ipios6.5y,
                       f_surv_nivpfs3y,f_surv_nivpfs6.5y,
                       f_surv_ipipfs3y,f_surv_ipipfs6.5y,
                      ncol=2, widths=c(18, 18))
ggsave("surv_fit.png",plot=surv_fit,width = 18,height = 18,dpi = 1200)
