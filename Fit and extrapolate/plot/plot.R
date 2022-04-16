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
####-----------------process data and plot 1--------------------####
# 3y
dfFig <- ltHaz[,11:30]
colnames(dfFig) <- md
dfFig <- cbind(data.frame(ltHaz$Time), dfFig)
dfFig = dfFig %>% mutate(Time =ltHaz.Time) %>% select(-ltHaz.Time) %>%
  gather(key = "Model", value = "Haz", -Time) %>% mutate(Model = factor(Model))
dfFig$Model = fct_recode(dfFig$Model, "lognormal"="lnorm", "loglogistic"="llogis")
# 6.5y
dfFig_out <- ltHaz_out[,11:30]
colnames(dfFig_out) <- md
dfFig_out <- cbind(data.frame(ltHaz_out$Time), dfFig_out)
dfFig_out = dfFig_out %>% mutate(Time =ltHaz_out.Time) %>% select(-ltHaz_out.Time) %>%
  gather(key = "Model", value = "Haz", -Time) %>% mutate(Model = factor(Model))
dfFig_out$Model = fct_recode(dfFig_out$Model, "lognormal"="lnorm", "loglogistic"="llogis")

###  hazard ratio   ###
# # Standard models
# fig1a <- ggplot(data=filter(dfFig_out, Model %in% c("exp","weibull","gamma","lognormal","gompertz","loglogistic","gengamma","genf")), aes(x=Time, y=Haz)) +
#   geom_point(data=ltHaz_out, aes(x=Time, y=hazLT, size=AtRisk), shape = 1) + geom_line(aes(color=Model)) + 
#   labs(x="Time", y="Hazard")  + guides(size="none") +
#   theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.07), xlim=c(0,7)) 

# ##  hazard ratio   ###
# #1
# fig_niv3yos <- ggplot(data=filter(dfFig, Model %in% c("exp","weibull","gamma","lognormal","gompertz","loglogistic","gengamma","genf",
#                                                       "FP1-1","FP1-2","FP2-1","FP2-2","RCS1","RCS2",
#                                                       "RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")), aes(x=Time, y=Haz)) +
#   geom_point(data=ltHaz, aes(x=Time, y=hazLT, size=AtRisk,), shape = 1) + geom_line(aes(color=Model),size=0.8,show.legend = FALSE) +
#   labs(x="Time", y="Hazard")  + guides(size="none") +
#   theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.06), xlim=c(0,3))+
#   labs(title="Hazard plot of NI group",
#        subtitle = "OS 3Ys")
# fig_niv3yos
# #2
# fig_niv6.5yos <- ggplot(data=filter(dfFig_out, Model %in% c("exp","weibull","gamma","lognormal","gompertz","loglogistic","gengamma","genf",
#                                                       "FP1-1","FP1-2","FP2-1","FP2-2","RCS1","RCS2",
#                                                       "RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")), aes(x=Time, y=Haz)) +
#   geom_point(data=ltHaz_out, aes(x=Time, y=hazLT, size=AtRisk,), shape = 1) + geom_line(aes(color=Model),size=0.8,show.legend = FALSE) +
#   labs(x="Time", y="Hazard")  + guides(size="none") +
#   theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.06), xlim=c(0,6.5))+
#   labs(title="Hazard plot of NI group",
#        subtitle = "OS 6.5Ys")
# fig_niv6.5yos
# #3
# fig_ipi3yos <- ggplot(data=filter(dfFig, Model %in% c("exp","weibull","gamma","lognormal","gompertz","loglogistic","gengamma","genf",
#                                                       "FP1-1","FP1-2","FP2-1","FP2-2","RCS1","RCS2",
#                                                       "RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")), aes(x=Time, y=Haz)) +
#   geom_point(data=ltHaz, aes(x=Time, y=hazLT, size=AtRisk,), shape = 1) + geom_line(aes(color=Model),size=0.8,show.legend = FALSE) +
#   labs(x="Time", y="Hazard")  + guides(size="none") +
#   theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.07), xlim=c(0,3))+
#   labs(title="Hazard plot of I group",
#        subtitle = "OS 3Ys")
# fig_ipi3yos
# #4
# fig_ipi6.5yos <- ggplot(data=filter(dfFig_out, Model %in% c("exp","weibull","gamma","lognormal","gompertz","loglogistic","gengamma","genf",
#                                                             "FP1-1","FP1-2","FP2-1","FP2-2","RCS1","RCS2",
#                                                             "RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")), aes(x=Time, y=Haz)) +
#   geom_point(data=ltHaz_out, aes(x=Time, y=hazLT, size=AtRisk,), shape = 1) + geom_line(aes(color=Model),size=0.8,show.legend = FALSE) +
#   labs(x="Time", y="Hazard")  + guides(size="none") +
#   theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.07), xlim=c(0,6.5))+
#   labs(title="Hazard plot of I group",
#        subtitle = "OS 6.5Ys")
# fig_ipi6.5yos
# #5
# fig_niv3ypfs <- ggplot(data=filter(dfFig, Model %in% c("exp","weibull","gamma","lognormal","gompertz","loglogistic","gengamma","genf",
#                                                       "FP1-1","FP1-2","FP2-1","FP2-2","RCS1","RCS2",
#                                                       "RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")), aes(x=Time, y=Haz)) +
#   geom_point(data=ltHaz, aes(x=Time, y=hazLT, size=AtRisk,), shape = 1) + geom_line(aes(color=Model),size=0.8) +
#   labs(x="Time", y="Hazard")  + guides(size="none") +
#   theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.25), xlim=c(0,3))+
#   labs(title="Hazard plot of NI group",
#        subtitle = "PFS 3Ys")
# fig_niv3ypfs
# #6
# fig_niv6.5ypfs <- ggplot(data=filter(dfFig_out, Model %in% c("exp","weibull","gamma","lognormal","gompertz","loglogistic","gengamma","genf",
#                                                             "FP1-1","FP1-2","FP2-1","FP2-2","RCS1","RCS2",
#                                                             "RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")), aes(x=Time, y=Haz)) +
#   geom_point(data=ltHaz_out, aes(x=Time, y=hazLT, size=AtRisk,), shape = 1) + geom_line(aes(color=Model),size=0.8) +
#   labs(x="Time", y="Hazard")  + guides(size="none") +
#   theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.2), xlim=c(0,6.5))+
#   labs(title="Hazard plot of NI group",
#        subtitle = "PFS 6.5Ys")
# fig_niv6.5ypfs
####-----------------------------------------------------------------------------------------------------------------------------------------####
####-----------------------------------------------------------------------------------------------------------------------------------------####
####-----------------------------------------------------------------------------------------------------------------------------------------####
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
####-----------------process data and plot 1--------------------####
# 3y
dfFig <- ltHaz[,11:29]
colnames(dfFig) <- md
dfFig <- cbind(data.frame(ltHaz$Time), dfFig)
dfFig = dfFig %>% mutate(Time =ltHaz.Time) %>% select(-ltHaz.Time) %>%
  gather(key = "Model", value = "Haz", -Time) %>% mutate(Model = factor(Model))
dfFig$Model = fct_recode(dfFig$Model, "lognormal"="lnorm", "loglogistic"="llogis")
# 6.5y
dfFig_out <- ltHaz_out[,11:29]
colnames(dfFig_out) <- md
dfFig_out <- cbind(data.frame(ltHaz_out$Time), dfFig_out)
dfFig_out = dfFig_out %>% mutate(Time =ltHaz_out.Time) %>% select(-ltHaz_out.Time) %>%
  gather(key = "Model", value = "Haz", -Time) %>% mutate(Model = factor(Model))
dfFig_out$Model = fct_recode(dfFig_out$Model, "lognormal"="lnorm", "loglogistic"="llogis")
#7
fig_ipi3ypfs <- ggplot(data=filter(dfFig, Model %in% c("exp","weibull","gamma","lognormal","gompertz","loglogistic","gengamma","genf",
                                                      "FP1-1","FP1-2","FP2-1","FP2-2","RCS1","RCS2",
                                                      "RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")), aes(x=Time, y=Haz)) +
  geom_point(data=ltHaz, aes(x=Time, y=hazLT, size=AtRisk,), shape = 1) + geom_line(aes(color=Model),size=0.8) +
  labs(x="Time", y="Hazard")  + guides(size="none") +
  theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.7), xlim=c(0,3))+
  labs(title="Hazard plot of I group",
       subtitle = "PFS 3Ys")
fig_ipi3ypfs
#8
fig_ipi6.5ypfs <- ggplot(data=filter(dfFig_out, Model %in% c("exp","weibull","gamma","lognormal","gompertz","loglogistic","gengamma","genf",
                                                            "FP1-1","FP1-2","FP2-1","FP2-2","RCS1","RCS2",
                                                            "RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")), aes(x=Time, y=Haz)) +
  geom_point(data=ltHaz_out, aes(x=Time, y=hazLT, size=AtRisk,), shape = 1) + geom_line(aes(color=Model),size=0.8) +
  labs(x="Time", y="Hazard")  + guides(size="none") +
  theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.7), xlim=c(0,6.5))+
  labs(title="Hazard plot of I group",
       subtitle = "PFS 6.5Ys")
fig_ipi6.5ypfs
####-----------------------------------------------------------------------------------------------------------------------------------------####

haz_fit<-grid.arrange(fig_niv3yos,fig_niv6.5yos,
                      fig_ipi3yos,fig_ipi6.5yos,
                      fig_niv3ypfs,fig_niv6.5ypfs,
                      fig_ipi3ypfs,fig_ipi6.5ypfs,
                      ncol=2, widths=c(15, 15))
ggsave("haz_fit.png",plot=haz_fit,width = 12,height = 18,dpi = 1200)



########==============================================================================================###
