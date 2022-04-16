dfHazEst<-read.csv("dfhazest_nivpfs_3y.csv")
dfHazEst<-dfHazEst[,2:79]
ltHaz<-read.csv("lthaz_nivpfs_3y.csv")
ltHaz<-ltHaz[,2:31]
ltHaz_out<-read.csv("lthaz_nivpfs_6.5y.csv")
ltHaz_out<-ltHaz_out[,2:31]
####------------------------------------------------------------####

# Get estimates of lifetime mean survival ####
fun_LE = function(haz, time){ # Function to estimate LE
  df = tibble(x = time, y = haz) %>%
    mutate(tau = x - lag(x, default=0),
           cum_y = cumsum(y * tau),
           surv =  exp(-cum_y),
           AUC = surv * tau)
  LE = sum(df$AUC)
  return(LE)
}

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

dfFig <- t(dfHazEst)
colnames(dfFig) <- md
dfFig <- cbind(data.frame(Newtime$Time), dfFig)
dfFig = dfFig %>% mutate(Time =Newtime.Time) %>% select(-Newtime.Time) %>%
  gather(key = "Model", value = "Haz", -Time) %>% mutate(Model = factor(Model))
dfFig$Model = fct_recode(dfFig$Model, "lognormal"="lnorm", "loglogistic"="llogis")
dfLE_dfHazEst<-dfFig
lthaz_plot<-data.frame(ltHaz$Time,ltHaz[,11:30])
lthaz_plot<-rename(lthaz_plot,"Time"="ltHaz.Time")
dfLE_lthaz<-lthaz_plot
dfLE_lthaz = dfLE_lthaz %>% 
  gather(key = "Model", value = "Haz", -Time) %>% mutate(Model = factor(Model))
dfLE_lthaz$Model = fct_recode(dfLE_lthaz$Model, "lognormal"="lnorm", "loglogistic"="llogis")

####------------------------------------------------------------####
####------------------------------------------------------------####
####------------------------------------------------------------####

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

dfFig <- t(dfHazEst)
colnames(dfFig) <- md
dfFig <- cbind(data.frame(Newtime$Time), dfFig)
dfFig = dfFig %>% mutate(Time =Newtime.Time) %>% select(-Newtime.Time) %>%
  gather(key = "Model", value = "Haz", -Time) %>% mutate(Model = factor(Model))
dfFig$Model = fct_recode(dfFig$Model, "lognormal"="lnorm", "loglogistic"="llogis")
dfLE_dfHazEst<-dfFig
lthaz_plot<-data.frame(ltHaz$Time,ltHaz[,11:29])
lthaz_plot<-rename(lthaz_plot,"Time"="ltHaz.Time")
dfLE_lthaz<-lthaz_plot
dfLE_lthaz = dfLE_lthaz %>% 
  gather(key = "Model", value = "Haz", -Time) %>% mutate(Model = factor(Model))
dfLE_lthaz$Model = fct_recode(dfLE_lthaz$Model, "lognormal"="lnorm", "loglogistic"="llogis")

###
#Origin
Ly_3y_nivos<- fun_LE(ltHaz$hazKM, ltHaz$Time)
Ly_6.5y_nivos<-fun_LE(ltHaz_out$hazKM, ltHaz_out$Time)

Ly_3y_ipios<- fun_LE(ltHaz$hazKM, ltHaz$Time)
Ly_6.5y_ipios<-fun_LE(ltHaz_out$hazKM, ltHaz_out$Time)

Ly_3y_nivpfs<- fun_LE(ltHaz$hazKM, ltHaz$Time)
Ly_6.5y_nivpfs<-fun_LE(ltHaz_out$hazKM, ltHaz_out$Time)

Ly_3y_ipipfs<- fun_LE(ltHaz$hazKM, ltHaz$Time)
Ly_6.5y_ipipfs<-fun_LE(ltHaz_out$hazKM, ltHaz_out$Time)

###

# #1
# LE_6.5y_nivos = dfLE_dfHazEst %>% group_by(Model) %>% summarise(LE = fun_LE(Haz, Time))
# LE_3y_nivos = dfLE_lthaz %>% group_by(Model) %>% summarise(LE = fun_LE(Haz, Time))
# #2
# LE_6.5y_ipios = dfLE_dfHazEst %>% group_by(Model) %>% summarise(LE = fun_LE(Haz, Time))
# LE_3y_ipios = dfLE_lthaz %>% group_by(Model) %>% summarise(LE = fun_LE(Haz, Time))

# #3
# LE_6.5y_nivpfs = dfLE_dfHazEst %>% group_by(Model) %>% summarise(LE = fun_LE(Haz, Time))
# LE_3y_nivpfs = dfLE_lthaz %>% group_by(Model) %>% summarise(LE = fun_LE(Haz, Time))
# #4
# LE_6.5y_ipipfs = dfLE_dfHazEst %>% group_by(Model) %>% summarise(LE = fun_LE(Haz, Time))
# LE_3y_ipipfs = dfLE_lthaz %>% group_by(Model) %>% summarise(LE = fun_LE(Haz, Time))
###
LE<-as.data.frame(array(dim = c(20,5)))
LE$V1<-LE_6.5y_nivos$Model
LE$V2<-LE_3y_nivos$LE
LE$V3<-LE_3y_ipios$LE
LE$V4<-LE_6.5y_nivos$LE
LE$V5<-LE_6.5y_ipios$LE
colnames(LE)<-c("Model","3ys_NI_OS","3ys_I_OS","6.5ys_NI_OS","6.5ys_I_OS")

write.csv(LE,"Lifetime_estimate.csv")
###
LE<-as.data.frame(array(dim = c(20,5)))
LE$V1<-LE_6.5y_nivpfs$Model
LE$V2<-LE_3y_nivpfs$LE
LE[1:6,3]<-LE_3y_ipipfs[1:6,2]
LE[8:20,3]<-LE_3y_ipipfs[7:19,2]
LE$V4<-LE_6.5y_nivpfs$LE
LE[1:6,5]<-LE_6.5y_ipipfs[1:6,2]
LE[8:20,5]<-LE_6.5y_ipipfs[7:19,2]

colnames(LE)<-c("Model","3ys_NI_PFS","3ys_I_PFS","6.5ys_NI_PFS","6.5ys_I_PFS")

write.csv(LE,"Lifetime_estimate_pfs.csv")
