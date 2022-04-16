dfHazEst<-read.csv("dfhazest_nivpfs_3y.csv")
dfHazEst<-dfHazEst[,2:79]
ltHaz<-read.csv("lthaz_nivpfs_3y.csv")
ltHaz<-ltHaz[,2:31]
ltHaz_out<-read.csv("lthaz_nivpfs_6.5y.csv")
ltHaz_out<-ltHaz_out[,2:31]
####------------------------------------------------------------####
follow_up <- 36
Newtime <- data.frame(Time = ltHaz_out$Time, AtRisk = 1)
Newtime$MyId <- 1:dim(Newtime)[1]
Newtime$MyId <- ifelse(Newtime$MyId > follow_up, follow_up, Newtime$MyId)  # Random effects: Using last observed ID for extrapolation
Newtime$EventsL <- 0
Newtime$EventsL[1:follow_up] <- lag(ltHaz$Events)
Newtime$EventsL[1] <- 0
Newtime$EventsL <- ifelse(Newtime$MyId > follow_up, 0, Newtime$EventsL) # AR: Using last observed event count for extrapolation
Newtime$timedelta<-Newtime$Time[2]-Newtime$Time[1]
Newtime$lnTime<-log(Newtime$Time)

md<-c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","genf","FP1-1","FP1-2",
      "FP2-1","FP2-2","RCS1","RCS2","RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")

#calculate survival over time
dfhar<-t(dfHazEst)
colnames(dfhar) <- md
dfhar <- cbind(data.frame(Newtime$Time), dfhar)
dfhar <- rename(dfhar,"Time"="Newtime.Time")

dfsurv<-as.data.frame(array(dim=c(78,21)))
colnames(dfsurv)<-c("Time",md)
dfsurv$Time<-dfhar$Time
###  20 model hazard  ###
for (i in 1:20) {
  temp_dfhaz<-data.frame(dfhar$Time,dfhar[,i+1])
  colnames(temp_dfhaz)<-c("time","haz")
  temp_dfsurv<-temp_dfhaz %>% 
    dplyr::arrange(time) %>% 
    dplyr::mutate(cumhaz = cumsum(haz)) %>% 
    dplyr::mutate(survProp = exp(-1*cumhaz))
  dfsurv[,i+1]<-temp_dfsurv$survProp
}
###combine surv data ####
dfFigSurv = dfsurv %>%
  gather(key = "Model", value = "survProp", -Time) %>% mutate(Model = factor(Model))
## 1
# f_surv_nivos= ggplot() +
#   geom_line(data=dfFigSurv, aes(x=Time, y=survProp, group=Model, colour=Model), size=1) +
#   geom_line(data = ltHaz_out,aes(x=Time,y=surv),size=1)+
#   geom_vline(aes(xintercept = 3),colour="#990000", linetype="dashed",size=1)+
#   geom_text(aes(x= 3.3, y = 1, label = "Start to extrapolate"),fontface='bold',size=4)+
#   scale_color_discrete(name="Model")+
#   expand_limits(y=c(0,1),x=c(0,6.5)) +
#   scale_x_continuous(breaks = c(seq(from=0, to=6.5,by = 0.5))) +
#   ylab("Overall survival") +
#   xlab("Time(months)") +
#   guides(color = guide_legend(ncol = 1))  +
#   theme(legend.key.width = unit(0.4,"cm"),legend.key.height = unit(0.4,"cm")) +
#   labs(title="Survival plot of NI group",
#        subtitle = "3Ys extrapolate to 6.5Ys")
# f_surv_nivos
# #2
# f_surv_ipios= ggplot() +
#   geom_line(data=dfFigSurv, aes(x=Time, y=survProp, group=Model, colour=Model), size=1) +
#   geom_line(data = ltHaz_out,aes(x=Time,y=surv),size=1)+
#   geom_vline(aes(xintercept = 3),colour="#990000", linetype="dashed",size=1)+
#   geom_text(aes(x= 3.3, y = 1, label = "Start to extrapolate"),fontface='bold',size=4)+
#   scale_color_discrete(name="Model")+
#   expand_limits(y=c(0,1),x=c(0,6.5)) +
#   scale_x_continuous(breaks = c(seq(from=0, to=6.5,by = 0.5))) +
#   ylab("Overall survival") +
#   xlab("Time(months)") +
#   guides(color = guide_legend(ncol = 1))  +
#   theme(legend.key.width = unit(0.4,"cm"),legend.key.height = unit(0.4,"cm")) +
#   labs(title="Survival plot of I group",
#        subtitle = "3Ys extrapolate to 6.5Ys")
# f_surv_ipios
#3
f_surv_nivpfs= ggplot() +
  geom_line(data=dfFigSurv, aes(x=Time, y=survProp, group=Model, colour=Model), size=1) +
  geom_line(data = ltHaz_out,aes(x=Time,y=surv),size=1)+
  geom_vline(aes(xintercept = 3),colour="#990000", linetype="dashed",size=1)+
  geom_text(aes(x= 3.3, y = 1, label = "Start to extrapolate"),fontface='bold',size=4)+
  scale_color_discrete(name="Model")+
  expand_limits(y=c(0,1),x=c(0,6.5)) + 
  scale_x_continuous(breaks = c(seq(from=0, to=6.5,by = 0.5))) +
  ylab("Progression-free survival") +
  xlab("Time(months)") +
  guides(color = guide_legend(ncol = 1))  +
  theme(legend.key.width = unit(0.4,"cm"),legend.key.height = unit(0.4,"cm")) +
  labs(title="Survival plot of NI group",
       subtitle = "3Ys extrapolate to 6.5Ys")
f_surv_nivpfs

####-------------------------------------------------------------------------------------------------------------------------###
####-------------------------------------------------------------------------------------------------------------------------###
####-------------------------------------------------------------------------------------------------------------------------###

dfHazEst<-read.csv("dfhazest_ipipfs_3y.csv")
dfHazEst<-dfHazEst[,2:79]
ltHaz<-read.csv("lthaz_ipipfs_3y.csv")
ltHaz<-ltHaz[,2:30]
ltHaz_out<-read.csv("lthaz_ipipfs_6.5y.csv")
ltHaz_out<-ltHaz_out[,2:30]
####------------------------------------------------------------####
follow_up <- 36
Newtime <- data.frame(Time = ltHaz_out$Time, AtRisk = 1)
Newtime$MyId <- 1:dim(Newtime)[1]
Newtime$MyId <- ifelse(Newtime$MyId > follow_up, follow_up, Newtime$MyId)  # Random effects: Using last observed ID for extrapolation
Newtime$EventsL <- 0
Newtime$EventsL[1:follow_up] <- lag(ltHaz$Events)
Newtime$EventsL[1] <- 0
Newtime$EventsL <- ifelse(Newtime$MyId > follow_up, 0, Newtime$EventsL) # AR: Using last observed event count for extrapolation
Newtime$timedelta<-Newtime$Time[2]-Newtime$Time[1]
Newtime$lnTime<-log(Newtime$Time)

md<-c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","FP1-1","FP1-2",
      "FP2-1","FP2-2","RCS1","RCS2","RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")

#calculate survival over time
dfhar<-t(dfHazEst)
colnames(dfhar) <- md
dfhar <- cbind(data.frame(Newtime$Time), dfhar)
dfhar <- rename(dfhar,"Time"="Newtime.Time")

dfsurv<-as.data.frame(array(dim=c(78,20)))
colnames(dfsurv)<-c("Time",md)
dfsurv$Time<-dfhar$Time
###  20 model hazard  ###
for (i in 1:19) {
  temp_dfhaz<-data.frame(dfhar$Time,dfhar[,i+1])
  colnames(temp_dfhaz)<-c("time","haz")
  temp_dfsurv<-temp_dfhaz %>% 
    dplyr::arrange(time) %>% 
    dplyr::mutate(cumhaz = cumsum(haz)) %>% 
    dplyr::mutate(survProp = exp(-1*cumhaz))
  dfsurv[,i+1]<-temp_dfsurv$survProp
}
###combine surv data ####
dfFigSurv = dfsurv %>%
  gather(key = "Model", value = "survProp", -Time) %>% mutate(Model = factor(Model))

f_surv_ipipfs= ggplot() +
  geom_line(data=dfFigSurv, aes(x=Time, y=survProp, group=Model, colour=Model), size=1) +
  geom_line(data = ltHaz_out,aes(x=Time,y=surv),size=1)+
  geom_vline(aes(xintercept = 3),colour="#990000", linetype="dashed",size=1)+
  geom_text(aes(x= 3.3, y = 1, label = "Start to extrapolate"),fontface='bold',size=4)+
  scale_color_discrete(name="Model")+
  expand_limits(y=c(0,1),x=c(0,6.5)) + 
  scale_x_continuous(breaks = c(seq(from=0, to=6.5,by = 0.5))) +
  ylab("Progression-free survival") +
  xlab("Time(months)") +
  guides(color = guide_legend(ncol = 1))  +
  theme(legend.key.width = unit(0.4,"cm"),legend.key.height = unit(0.4,"cm")) +
  labs(title="Survival plot of NI group",
       subtitle = "3Ys extrapolate to 6.5Ys")
f_surv_ipipfs

surv_extraploate<-grid.arrange(f_surv_nivos,f_surv_ipios,
                               f_surv_nivpfs,f_surv_ipipfs,
                       ncol=2, widths=c(18, 18))
ggsave("surv_extraploate.png",plot=surv_extraploate,width = 18,height = 12,dpi = 1200)
