library("ggplot2")   # Mainly for plotting
library("tidyverse")   # Mainly for plotting
library("flexsurv")    # RPs (also loads survival) and BC case-study
library("gridExtra")   # Plotting
library("survHE")   # Plotting
library("discSurv")    # Create life tables
library("mgcv")        # GAM/RCS
library("survminer")   #KM curve


theme_set(theme_light())  # GGplot theme

## Load data base

niv5y <- read.delim("IPD_os_niv_6.5y.txt")
niv5y<-data.frame(niv5y$event,niv5y$time)
niv5y<-rename(niv5y,"censrec"="niv5y.event","recyrs"="niv5y.time")
niv5y$recyrs<-niv5y$recyrs/12
niv5y$rectime<-as.integer(niv5y$recyrs*365.24) 

niv5y = niv5y %>% mutate(censrec = case_when(niv5y$recyrs > 6.5 ~ integer(1), TRUE ~ niv5y$censrec),
                         recyrs = case_when(niv5y$recyrs > 6.5 ~ 6.5, TRUE ~ niv5y$recyrs),
                         rectime = recyrs * 365)

niv3y <- read.delim("IPD_os_niv_3y.txt")
niv3y<-data.frame(niv3y$event,niv3y$time)
niv3y<-rename(niv3y,"censrec"="niv3y.event","recyrs"="niv3y.time")
niv3y$recyrs<-niv3y$recyrs/12
niv3y$rectime<-as.integer(niv3y$recyrs*365.24) 

niv3y = niv3y %>% mutate(censrec = case_when(niv3y$recyrs > 3 ~ integer(1), TRUE ~ niv3y$censrec),
                         recyrs = case_when(niv3y$recyrs > 3 ~ 3, TRUE ~ niv3y$recyrs),
                         rectime = recyrs * 365)

ipi5y <- read.delim("IPD_os_ipi_6.5y.txt")
ipi5y<-data.frame(ipi5y$event,ipi5y$time)
ipi5y<-rename(ipi5y,"censrec"="ipi5y.event","recyrs"="ipi5y.time")
ipi5y$recyrs<-ipi5y$recyrs/12
ipi5y$rectime<-as.integer(ipi5y$recyrs*365.24) 

ipi5y = ipi5y %>% mutate(censrec = case_when(ipi5y$recyrs > 6.5 ~ integer(1), TRUE ~ ipi5y$censrec),
                         recyrs = case_when(ipi5y$recyrs > 6.5 ~ 6.5, TRUE ~ ipi5y$recyrs),
                         rectime = recyrs * 365)

ipi3y <- read.delim("IPD_os_ipi_3y.txt")
ipi3y<-data.frame(ipi3y$event,ipi3y$time)
ipi3y<-rename(ipi3y,"censrec"="ipi3y.event","recyrs"="ipi3y.time")
ipi3y$recyrs<-ipi3y$recyrs/12
ipi3y$rectime<-as.integer(ipi3y$recyrs*365.24) 

ipi3y = ipi3y %>% mutate(censrec = case_when(ipi3y$recyrs > 3 ~ integer(1), TRUE ~ ipi3y$censrec),
                         recyrs = case_when(ipi3y$recyrs > 3 ~ 3, TRUE ~ ipi3y$recyrs),
                         rectime = recyrs * 365)

###----------------------------------------------------------- Process NIV data ------------------------------------------------------------------------------###
###   Process NIV data   ####
## Learn about niv

table(niv5y$censrec)
niv5y %>%
  summarise(Min_surv = min(recyrs, na.rm = TRUE),
            Max_surv = max(recyrs, na.rm = TRUE),
            Mean_surv = mean(recyrs, na.rm = TRUE),
            SD_surv = sd(recyrs, na.rm = TRUE))

table(niv3y$censrec)
niv3y %>%
  summarise(Min_surv = min(recyrs, na.rm = TRUE),
            Max_surv = max(recyrs, na.rm = TRUE),
            Mean_surv = mean(recyrs, na.rm = TRUE),
            SD_surv = sd(recyrs, na.rm = TRUE))

##-----niv data 3Y Monthly life table estimates of hazard-----###
niv3y$rectime2 <- as.integer(niv3y$rectime/(365.24/12)) + 1
niv5y$rectime2 <- as.integer(niv5y$rectime/(365.24/12)) + 1
#+1 above as integer rounds down, we want to round up
ltBC <- lifeTable(niv3y, timeColumn = "rectime2", censColumn = "censrec")
ltHaz0 <- data.frame(hazKM = ltBC$Output$hazard, Time = (seq(1:length(ltBC$Output[,1]))-0.5)/12,
                    AtRisk = ltBC$Output$atRisk, Events = ltBC$Output$events)
# The above hazard is the product-limit (KM) estimate. Also calculate the life-table (acturial) estimate
ltHaz0$hazLT = ltHaz0$Events / (ltHaz0$AtRisk - ltHaz0$Events/2)
# Generate log-time
ltHaz0$lnTime <- log(ltHaz0$Time)
# For random effects add an ID for each time period
ltHaz0$MyId <- 1:dim(ltHaz0)[1] # Generate id variable 
# For AR(1) model get outcomes lagged by one.
ltHaz0$EventsL <- lag(ltHaz0$Events)
# Set first lagged value = 0 (usually would discard, but retain so IC are comparable. Can be justified as a prior value)
ltHaz0$EventsL[1] <- 0
#Set surv data
ltHaz0$surv <- ltBC$Output$S
#timedelta
ltHaz0$timedelta<-ltHaz0$Time[2]-ltHaz0$Time[1]

##-----niv data 5Y Monthly life table estimates of hazard-----### 
ltBC_out <- lifeTable(niv5y, timeColumn = "rectime2", censColumn = "censrec")
ltHaz0_out <- data.frame(hazKM = ltBC_out$Output$hazard, Time = (seq(1:length(ltBC_out$Output[,1]))-0.5)/12,
                        AtRisk = ltBC_out$Output$atRisk, Events = ltBC_out$Output$events)
# The above hazard is the product-limit (KM) estimate. Also calculate the life-table (acturial) estimate
ltHaz0_out$hazLT = ltHaz0_out$Events / (ltHaz0_out$AtRisk - ltHaz0_out$Events/2)
# Generate log-time
ltHaz0_out$lnTime <- log(ltHaz0_out$Time)
# For random effects add an ID for each time period
ltHaz0_out$MyId <- 1:dim(ltHaz0_out)[1] # Generate id variable 
# For AR(1) model get outcomes lagged by one.
ltHaz0_out$EventsL <- lag(ltHaz0_out$Events)
# Set first lagged value = 0 (usually would discard, but retain so IC are comparable. Can be justified as a prior value)
ltHaz0_out$EventsL[1] <- 0
#Set surv data
ltHaz0_out$surv <- ltBC_out$Output$S
#timedelta
ltHaz0_out$timedelta<-ltHaz0_out$Time[2]-ltHaz0_out$Time[1]

###------------------------------------------------------------ Process ipi data -----------------------------------------------------------------------------####
###  Process ipi data  ####
## Learn about ipi

table(ipi5y$censrec)
ipi5y %>%
  summarise(Min_surv = min(recyrs, na.rm = TRUE),
            Max_surv = max(recyrs, na.rm = TRUE),
            Mean_surv = mean(recyrs, na.rm = TRUE),
            SD_surv = sd(recyrs, na.rm = TRUE))

table(ipi3y$censrec)
ipi3y %>%
  summarise(Min_surv = min(recyrs, na.rm = TRUE),
            Max_surv = max(recyrs, na.rm = TRUE),
            Mean_surv = mean(recyrs, na.rm = TRUE),
            SD_surv = sd(recyrs, na.rm = TRUE))


##-----ipi data 3Y Monthly life table estimates of hazard-----##
ipi3y$rectime2 <- as.integer(ipi3y$rectime/(365.24/12)) + 1
ipi5y$rectime2 <- as.integer(ipi5y$rectime/(365.24/12)) + 1
#+1 above as integer rounds down, we want to round up
ltBC1 <- lifeTable(ipi3y, timeColumn = "rectime2", censColumn = "censrec")
ltHaz1 <- data.frame(hazKM = ltBC1$Output$hazard, Time = (seq(1:length(ltBC1$Output[,1]))-0.5)/12,
                    AtRisk = ltBC1$Output$atRisk, Events = ltBC1$Output$events)
# The above hazard is the product-limit (KM) estimate. Also calculate the life-table (acturial) estimate
ltHaz1$hazLT = ltHaz1$Events / (ltHaz1$AtRisk - ltHaz1$Events/2)
# Generate log-time
ltHaz1$lnTime <- log(ltHaz1$Time)
# For random effects add an ID for each time period
ltHaz1$MyId <- 1:dim(ltHaz1)[1] # Generate id variable 
# For AR(1) model get outcomes lagged by one.
ltHaz1$EventsL <- lag(ltHaz1$Events)
# Set first lagged value = 0 (usually would discard, but retain so IC are comparable. Can be justified as a prior value)
ltHaz1$EventsL[1] <- 0
#Set surv data
ltHaz1$surv <- ltBC1$Output$S
#timedelta
ltHaz1$timedelta<-ltHaz1$Time[2]-ltHaz1$Time[1]


##-----ipi data 5Y Monthly life table estimates of hazard----- ##
ltBC1_out <- lifeTable(ipi5y, timeColumn = "rectime2", censColumn = "censrec")
ltHaz1_out <- data.frame(hazKM = ltBC1_out$Output$hazard, Time = (seq(1:length(ltBC1_out$Output[,1]))-0.5)/12,
                        AtRisk = ltBC1_out$Output$atRisk, Events = ltBC1_out$Output$events)
# The above hazard is the product-limit (KM) estimate. Also calculate the life-table (acturial) estimate
ltHaz1_out$hazLT = ltHaz1_out$Events / (ltHaz1_out$AtRisk - ltHaz1_out$Events/2)
# Generate log-time
ltHaz1_out$lnTime <- log(ltHaz1_out$Time)
# For random effects add an ID for each time period
ltHaz1_out$MyId <- 1:dim(ltHaz1_out)[1] # Generate id variable 
# For AR(1) model get outcomes lagged by one.
ltHaz1_out$EventsL <- lag(ltHaz1_out$Events)
# Set first lagged value = 0 (usually would discard, but retain so IC are comparable. Can be justified as a prior value)
ltHaz1_out$EventsL[1] <- 0
#Set surv data
ltHaz1_out$surv <- ltBC1_out$Output$S
#timedelta
ltHaz1_out$timedelta<-ltHaz1_out$Time[2]-ltHaz1_out$Time[1]

###--------------------------------------------------------------------   log cum haz  ---------------------------------------------------------------------####
###   log cum haz  ####
#3y
lch_niv_3y<-ltHaz0 %>%
  dplyr::mutate(logTime = log(Time)) %>%
  dplyr::mutate(cumhaz = cumsum(hazKM)) %>% 
  dplyr::mutate(logcumhaz = log(cumhaz)) %>%
  dplyr::mutate(survProp = exp(-1*cumhaz))
lch_ipi_3y<-ltHaz1 %>%
  dplyr::mutate(logTime = log(Time)) %>%
  dplyr::mutate(cumhaz = cumsum(hazKM)) %>% 
  dplyr::mutate(logcumhaz = log(cumhaz)) %>%
  dplyr::mutate(survProp = exp(-1*cumhaz))

f_surv1= ggplot() +
  geom_line(data=lch_niv_3y, aes(x=logTime, y=logcumhaz), size=1) +
  expand_limits(y=c(-4,1),x=c(-4,1)) + 
  scale_x_continuous(breaks = c(seq(from=-4, to=1,by = 1))) +
  ylab("logcumhaz") +
  xlab("logTime") +
  guides(color = guide_legend(ncol = 1))  +
  theme(legend.position = "bottom") + 
  theme_bw() 
f_surv1

f_surv2= ggplot() +
  geom_line(data=lch_ipi_3y, aes(x=logTime, y=logcumhaz), size=1) +
  expand_limits(y=c(-4,1),x=c(-4,1)) + 
  scale_x_continuous(breaks = c(seq(from=-4, to=1,by = 1))) +
  ylab("logcumhaz") +
  xlab("logTime") +
  guides(color = guide_legend(ncol = 1))  +
  theme(legend.position = "bottom") + 
  theme_bw() 
f_surv2

f_surv_c1= ggplot() +
  geom_line(data=lch_niv_3y, aes(x=logTime, y=logcumhaz), size=1) +
  geom_line(data=lch_ipi_3y, aes(x=logTime, y=logcumhaz), size=1) +
  expand_limits(y=c(-4,1),x=c(-4,1)) + 
  scale_x_continuous(breaks = c(seq(from=-4, to=1,by = 1))) +
  ylab("logcumhaz") +
  xlab("logTime") +
  guides(color = guide_legend(ncol = 1))  +
  theme(legend.position = "bottom") + 
  theme_bw() 
f_surv_c1
#5y
lch_niv_5y<-ltHaz0_out %>%
  dplyr::mutate(logTime = log(Time)) %>%
  dplyr::mutate(cumhaz = cumsum(hazKM)) %>% 
  dplyr::mutate(logcumhaz = log(cumhaz)) %>%
  dplyr::mutate(survProp = exp(-1*cumhaz))
lch_ipi_5y<-ltHaz1_out %>%
  dplyr::mutate(logTime = log(Time)) %>%
  dplyr::mutate(cumhaz = cumsum(hazKM)) %>% 
  dplyr::mutate(logcumhaz = log(cumhaz)) %>%
  dplyr::mutate(survProp = exp(-1*cumhaz))

f_surv3= ggplot() +
  geom_line(data=lch_niv_5y, aes(x=logTime, y=logcumhaz), size=1) +
  expand_limits(y=c(-4,1),x=c(-4,1)) + 
  scale_x_continuous(breaks = c(seq(from=-4, to=1,by = 1))) +
  ylab("logcumhaz") +
  xlab("logTime") +
  guides(color = guide_legend(ncol = 1))  +
  theme(legend.position = "bottom") + 
  theme_bw() 
f_surv3

f_surv4= ggplot() +
  geom_line(data=lch_ipi_5y, aes(x=logTime, y=logcumhaz), size=1) +
  expand_limits(y=c(-4,1),x=c(-4,1)) + 
  scale_x_continuous(breaks = c(seq(from=-4, to=1,by = 1))) +
  ylab("logcumhaz") +
  xlab("logTime") +
  guides(color = guide_legend(ncol = 1))  +
  theme(legend.position = "bottom") + 
  theme_bw() 
f_surv4

f_surv_c2= ggplot() +
  geom_line(data=lch_niv_5y, aes(x=logTime, y=logcumhaz), size=1) +
  geom_line(data=lch_ipi_5y, aes(x=logTime, y=logcumhaz), size=1) +
  expand_limits(y=c(-4,1),x=c(-4,1)) + 
  scale_x_continuous(breaks = c(seq(from=-4, to=1,by = 1))) +
  ylab("logcumhaz") +
  xlab("logTime") +
  guides(color = guide_legend(ncol = 1))  +
  theme(legend.position = "bottom") + 
  theme_bw() 
f_surv_c2
###--------------------------------------------------------------------------------------------------------------------------------------------------------####
##check survival data
fit_niv3y<-survfit(Surv(niv3y$recyrs,niv3y$censrec)~1,data=niv3y)
ggsurvplot(fit_niv3y,data=niv3y)
fit_niv5y<-survfit(Surv(niv5y$recyrs,niv5y$censrec)~1,data=niv5y)
ggsurvplot(fit_niv5y,data=niv5y)
fit_ipi3y<-survfit(Surv(ipi3y$recyrs,ipi3y$censrec)~1,data=ipi3y)
ggsurvplot(fit_ipi3y,data=ipi3y)
fit_ipi5y<-survfit(Surv(ipi5y$recyrs,ipi5y$censrec)~1,data=ipi5y)
ggsurvplot(fit_ipi5y,data=ipi5y)

###--------------------------------------------------------------------------------------------------------------------------------------------------------####
##  data for process  models ####
##data for sd
## niv ##
ltHaz<-ltHaz0
ltHaz_out<-ltHaz0_out
sd_bc<-data.frame("recyrs"=niv3y$recyrs,"censrec"=niv3y$censrec)
sd_bc_out<-data.frame("recyrs"=niv5y$recyrs,"censrec"=niv5y$censrec)
## ipi ##
# ltHaz<-ltHaz1
# ltHaz_out<-ltHaz1_out
# sd_bc<-data.frame("recyrs"=ipi3y$recyrs,"censrec"=ipi3y$censrec)
# sd_bc_out<-data.frame("recyrs"=ipi5y$recyrs,"censrec"=ipi5y$censrec)


###--------------------------------------------------------------------------------------------------------------------------------------------------------####
####----New Data----####
follow_up <- 36
numMod <- 20 # Models considered
MyTH <- 6.5 # Time Horizon (years)
MyStep <- 12 # Number of obs. per year
MyN <- MyTH*MyStep # Total time points (observed & extrapolated)
dfHazEst <- array(dim=c(numMod, MyN))
Newtime <- data.frame(Time = ltHaz_out$Time, AtRisk = 1)
Newtime$MyId <- 1:dim(Newtime)[1]
Newtime$MyId <- ifelse(Newtime$MyId > follow_up, follow_up, Newtime$MyId)  # Random effects: Using last observed ID for extrapolation
Newtime$EventsL <- 0
Newtime$EventsL[1:follow_up] <- lag(ltHaz$Events)
Newtime$EventsL[1] <- 0
Newtime$EventsL <- ifelse(Newtime$MyId > follow_up, 0, Newtime$EventsL) # AR: Using last observed event count for extrapolation
Newtime$timedelta<-Newtime$Time[2]-Newtime$Time[1]
Newtime$lnTime<-log(Newtime$Time)
# Also have 1x GOF matrix. Rows = Methods, Columns = Method, LL, AIC
dfGOF <- data.frame(matrix(, nrow=20, ncol=4))
dfGOF_out <- data.frame(matrix(, nrow=20, ncol=4))
colnames(dfGOF) <- c("Model","LnL","Params","AIC")
colnames(dfGOF_out) <- c("Model","LnL","Params","AIC")
# Below is constant for when have to derive log-likelihood
llCons <- sum(ltHaz$Events*log(ltHaz$AtRisk) - log(factorial(ltHaz$Events)))
llCons_out<-sum(ltHaz_out$Events*log(ltHaz_out$AtRisk) - log(factorial(ltHaz_out$Events)))
# Names of models to consider
modnames <- list("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","genf","FP1-1","FP1-2",
                 "FP2-1","FP2-2","RCS1","RCS2","RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")
md<-c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","genf","FP1-1","FP1-2",
     "FP2-1","FP2-2","RCS1","RCS2","RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")
###--------------------------------------------------------------------------------------------------------------------------------------------------------####
#=============================================================================================================#
#=============================================================================================================#

#                                       2y res    +    2y to 5y Extrapolation

#=============================================================================================================#
#=============================================================================================================#

##########################
#######    RCS     #######
##########################
MODi <- 1
# First need knot locations for up to 5 internal knots.
# Basing these on equally-spaced percentiles of the observed (uncensored) death times.
bc2 <- subset(niv3y, censrec==1)
myLnL <- array(dim=5)
myAIC <- array(dim=5)
for (i in 1:5){
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz)
  myLnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[i] <- extractAIC(glmTemp)[2]
  dfHazEst[MODi,]<-predict(glmTemp, newdata=Newtime, type="response")
  dfGOF[MODi,2] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF[MODi,3] <- extractAIC(glmTemp)[1]
  MODi<-MODi+1
}
RCSres <- data.frame(c("One","Two","Three","Four","Five"))
RCSres <- cbind(RCSres,myLnL,myAIC)
colnames(RCSres) <- c("Int.Knots","LnL","AIC")
RCSres<-arrange(RCSres)


bc2 <- subset(niv3y, censrec==1)
myLnL1 <- array(dim=5)
myAIC1 <- array(dim=5)
for (i in 1:5){
  glmTemp <- gam(Events ~ s(lnTime, bs="cr", k=i+2, fx=TRUE) + offset(log(AtRisk)), knots=list(lnTime=quantile(log(bc2$recyrs), seq(from=0, to=1, by=1/(1+i))), length=i+2), family=poisson, data=ltHaz)
  
  myLnL1[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC1[i] <- extractAIC(glmTemp)[2]
  dfHazEst[MODi,]<-predict(glmTemp, newdata=Newtime, type="response")
  dfGOF[MODi,2] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF[MODi,3] <- extractAIC(glmTemp)[1]
  MODi<-MODi+1
}
RCSres1 <- data.frame(c("One","Two","Three","Four","Five"))
RCSres1 <- cbind(RCSres,myLnL1,myAIC1)
colnames(RCSres1) <- c("Int.Knots","LnL","AIC")
RCSres<-arrange(RCSres1)


# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# RCSres[1,]
# RCSres[2,]
# i<-1
# glmTemp1 <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2),family=binomial(link=cloglog), data=ltHaz)
# i<-2
# glmTemp2 <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2),family=binomial(link=cloglog), data=ltHaz)
# 
# dfGOF[MODi,1] <- modnames[[MODi]]
# dfGOF[MODi,2] <- (extractAIC(glmTemp1)[2] - 2*extractAIC(glmTemp1)[1])*(-0.5)
# dfGOF[MODi,3] <- extractAIC(glmTemp1)[1]
# # Hazard estimates
# dfHazEst[MODi,] <- predict(glmTemp1, newdata=Newtime, type="response")
# ltHaz[modnames[[MODi]]] <- predict(glmTemp1, newdata=ltHaz, type="response")
# 
# MODi<-14
# dfGOF[MODi,1] <- modnames[[MODi]]
# dfGOF[MODi,2] <- (extractAIC(glmTemp2)[2] - 2*extractAIC(glmTemp2)[1])*(-0.5)
# dfGOF[MODi,3] <- extractAIC(glmTemp2)[1]
# # Hazard estimates
# dfHazEst[MODi,] <- predict(glmTemp2, newdata=Newtime, type="response")
# ltHaz[modnames[[MODi]]] <- predict(glmTemp2, newdata=ltHaz, type="response")
# 
# 
# ##################------------surv rate calculate  ###  dfhazest   ###---------------##############
# #calculate survival over time
# dfhar<-t(dfHazEst)
# colnames(dfhar) <- dfGOF[,1]
# dfhar <- cbind(data.frame(Newtime$Time), dfhar)
# dfhar <- rename(dfhar,"Time"="Newtime.Time")
# 
# dfsurv<-as.data.frame(array(dim=c(78,21)))
# colnames(dfsurv)<-c("Time",modnames)
# dfsurv$Time<-dfhar$Time
# ###  20 model hazard  ###
# for (i in 1:20) {
#   temp_dfhaz<-data.frame(dfhar$Time,dfhar[,i+1])
#   colnames(temp_dfhaz)<-c("time","haz")
#   temp_dfsurv<-temp_dfhaz %>% 
#     dplyr::arrange(time) %>% 
#     dplyr::mutate(cumhaz = cumsum(haz)) %>% 
#     dplyr::mutate(survProp = exp(-1*cumhaz))
#   dfsurv[,i+1]<-temp_dfsurv$survProp
# }
# 
# 
# ###combine surv data ####
# dfFigSurv = dfsurv %>%
#   gather(key = "Model", value = "survProp", -Time) %>% mutate(Model = factor(Model))
# ##plot##
# f_surv= ggplot() +
#   geom_line(data=dfFigSurv, aes(x=Time, y=survProp, group=Model, colour=Model), size=1) +
#   geom_point(data = ltHaz_out,aes(x=Time,y=surv),size=1,shape=1)+
#   scale_color_discrete(name="Model")+
#   expand_limits(y=c(0,1),x=c(0,7)) + 
#   facet_wrap(~Model,nrow=4)+
#   scale_x_continuous(breaks = c(seq(from=0, to=7,by = 1))) +
#   ylab("Overall survival") +
#   xlab("Time(months)") +
#   guides(color = guide_legend(ncol = 1))  +
#   theme(legend.position = "bottom") + 
#   theme_bw() 
# f_surv
# 
# ggsave("fig_surv_niv_3y_os.png",f_surv,width = 18,height = 12,dpi = 800)
# 
# ##################------------surv rate calculate  ###  lthaz   ###---------------##############
# lthaz_plot<-data.frame(ltHaz$Time,ltHaz[,11:30])
# dfsurv_lthaz<-as.data.frame(array(dim=c(36,21)))
# colnames(dfsurv_lthaz)<-c("Time",modnames)
# dfsurv_lthaz$Time<-lthaz_plot$ltHaz.Time
# for (i in 1:20) {
#   temp_dfhaz<-data.frame(lthaz_plot$ltHaz.Time,lthaz_plot[,i+1])
#   colnames(temp_dfhaz)<-c("time","haz")
#   temp_dfsurv<-temp_dfhaz %>% 
#     dplyr::arrange(time) %>% 
#     dplyr::mutate(cumhaz = cumsum(haz)) %>% 
#     dplyr::mutate(survProp = exp(-1*cumhaz))
#   dfsurv_lthaz[,i+1]<-temp_dfsurv$survProp
# }
# 
# dfFigSurv_lthaz = dfsurv_lthaz %>%
#   gather(key = "Model", value = "survProp", -Time) %>% mutate(Model = factor(Model))
# ##plot##
# f_surv2= ggplot() +
#   geom_line(data=dfFigSurv_lthaz, aes(x=Time, y=survProp, group=Model, colour=Model), size=1) +
#   geom_line(data = ltHaz,aes(x=Time,y=surv),size=0.5)+
#   scale_color_discrete(name="Model")+
#   expand_limits(y=c(0,1),x=c(0,3)) + 
#   facet_wrap(~Model,nrow=4)+
#   scale_x_continuous(breaks = c(seq(from=0, to=3,by = 0.5))) +
#   ylab("Overall survival") +
#   xlab("Time(months)") +
#   guides(color = guide_legend(ncol = 1))  +
#   theme(legend.position = "bottom") + 
#   theme_bw() 
# f_surv2
# 
# ggsave("fig_surv_niv_3y_os_3ys.png",f_surv2,width = 18,height = 12,dpi = 800)
# 
# ########==============================================================================================###
# 
# 
# 
# 
# 
# 








##########################
#######    RCS     #######
##########################

# First need knot locations for up to 5 internal knots.
# Basing these on equally-spaced percentiles of the observed (uncensored) death times.
bc2_out <- subset(niv5y, censrec==1)
myLnL <- array(dim=5)
myAIC <- array(dim=5)
for (i in 1:5){
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2_out$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_out)
  myLnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[i] <- extractAIC(glmTemp)[2]
  dfHazEst[MODi,]<-predict(glmTemp, newdata=Newtime, type="response")
  dfGOF[MODi,2] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF[MODi,3] <- extractAIC(glmTemp)[1]
  MODi<-MODi+1
}
RCSres <- data.frame(c("One","Two","Three","Four","Five"))
RCSres <- cbind(RCSres,myLnL,myAIC)
colnames(RCSres) <- c("Int.Knots","LnL","AIC")
RCSres<-arrange(RCSres)



bc2_out <- subset(niv5y, censrec==1)
myLnL1 <- array(dim=5)
myAIC1 <- array(dim=5)
for (i in 1:5){

  glmTemp <- gam(Events ~ s(lnTime, bs="cr", k=i+2, fx=TRUE) + offset(log(AtRisk)), knots=list(lnTime=quantile(log(bc2$recyrs), seq(from=0, to=1, by=1/(1+i))), length=i+2), family=poisson, data=ltHaz)
  
  myLnL1[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC1[i] <- extractAIC(glmTemp)[2]
  dfHazEst[MODi,]<-predict(glmTemp, newdata=Newtime, type="response")
  dfGOF[MODi,2] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF[MODi,3] <- extractAIC(glmTemp)[1]
  MODi<-MODi+1
}
RCSres1 <- data.frame(c("One","Two","Three","Four","Five"))
RCSres1 <- cbind(RCSres1,myLnL,myAIC)
colnames(RCSres1) <- c("Int.Knots","LnL","AIC")
RCSres<-arrange(RCSres1)




# 
# 
# 
# 
# 
# 
# 
# 
# RCSres[1,]
# RCSres[2,]
# i<-1
# glmTemp1 <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2_out$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2),family=binomial(link=cloglog), data=ltHaz_out)
# i<-2
# glmTemp2 <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2_out$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2),family=binomial(link=cloglog), data=ltHaz_out)
# 
# dfGOF_out[MODi,1] <- modnames[[MODi]]
# dfGOF_out[MODi,2] <- (extractAIC(glmTemp1)[2] - 2*extractAIC(glmTemp1)[1])*(-0.5)
# dfGOF_out[MODi,3] <- extractAIC(glmTemp1)[1]
# # Hazard estimates
# ltHaz_out[modnames[[MODi]]] <- predict(glmTemp1, newdata=ltHaz_out, type="response")
# 
# MODi<-14
# dfGOF_out[MODi,1] <- modnames[[MODi]]
# dfGOF_out[MODi,2] <- (extractAIC(glmTemp2)[2] - 2*extractAIC(glmTemp2)[1])*(-0.5)
# dfGOF_out[MODi,3] <- extractAIC(glmTemp2)[1]
# # Hazard estimates
# ltHaz_out[modnames[[MODi]]] <- predict(glmTemp2, newdata=ltHaz_out, type="response")
# 
# 
# ########==============================================================================================###
# 
# ####-----------------process data and plot 1--------------------####
# 
# dfFig_out <- ltHaz_out[,11:30]
# colnames(dfFig_out) <- dfGOF_out[,1]
# dfFig_out <- cbind(data.frame(ltHaz_out$Time), dfFig_out)
# dfFig_out = dfFig_out %>% mutate(Time =ltHaz_out.Time) %>% select(-ltHaz_out.Time) %>%
#   gather(key = "Model", value = "Haz", -Time) %>% mutate(Model = factor(Model))
# dfFig_out$Model = fct_recode(dfFig_out$Model, "lognormal"="lnorm", "loglogistic"="llogis")
# 
# ###  hazard ratio   ###
# # Standard models
# fig1a <- ggplot(data=filter(dfFig_out, Model %in% c("exp","weibull","gamma","lognormal","gompertz","loglogistic","gengamma","genf")), aes(x=Time, y=Haz)) +
#   geom_point(data=ltHaz_out, aes(x=Time, y=hazLT, size=AtRisk), shape = 1) + geom_line(aes(color=Model)) + 
#   labs(x="Time", y="Hazard")  + guides(size="none") +
#   theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.07), xlim=c(0,7)) 
# #+ theme(axis.text.y = element_blank())
# #fig1a
# 
# # FP models
# fig1b <- ggplot(data=filter(dfFig_out, Model %in% c("FP1-1","FP1-2","FP2-1","FP2-2")), aes(x=Time, y=Haz)) +
#   geom_point(data=ltHaz_out, aes(x=Time, y=hazLT, size=AtRisk), shape = 1) + geom_line(aes(color=Model)) + 
#   labs(x="Time", y="Hazard") + guides(size="none") +
#   theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.07), xlim=c(0,7)) 
# #fig1b
# 
# # RCS models
# fig1c <- ggplot(data=filter(dfFig_out, Model %in% c("RCS1","RCS2")), aes(x=Time, y=Haz)) +
#   geom_point(data=ltHaz_out, aes(x=Time, y=hazLT, size=AtRisk), shape = 1) + geom_line(aes(color=Model)) + 
#   labs(x="", y="Hazard")+ guides(size="none") +
#   theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.07), xlim=c(0,7)) 
# #fig1c
# 
# # RP models
# fig1d <- ggplot(data=filter(dfFig_out, Model %in% c("RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")), aes(x=Time, y=Haz)) +
#   geom_point(data=ltHaz_out, aes(x=Time, y=hazLT, size=AtRisk), shape = 1) + geom_line(aes(color=Model)) + 
#   labs(x="", y="Hazard")  + guides(size="none") +
#   theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.07), xlim=c(0,7)) 
# #fig1d
# 
# fig_haz_niv_5y<-grid.arrange(fig1a, fig1b, fig1c, fig1d, ncol=2, widths=c(15, 15))
# ggsave("fig_haz_niv_6.5y_os.png",plot = fig_haz_niv_5y,width = 18,height = 12,dpi = 800)



##################------------surv rate calculate  ###  dfhazest   ###---------------##############
#calculate survival over time

dfhar_out<-t(dfHazEst)
rcsname<-c("3y_t-1","3y_t-2","3y_t-3","3y_t-4","3y_t-5",
           "3y_lnt-1","3y_lnt-2","3y_lnt-3","3y_lnt-4","3y_lnt-5",
           "6.5y_t-1","6.5y_t-2","6.5y_t-3","6.5y_t-4","6.5y_t-5",
           "6.5y_lnt-1","6.5y_lnt-2","6.5y_lnt-3","6.5y_lnt-4","6.5y_lnt-5")
colnames(dfhar_out) <- rcsname
dfhar_out <- cbind(data.frame(ltHaz_out$Time), dfhar_out)
dfhar_out <- rename(dfhar_out,"Time"="ltHaz_out.Time")

dfsurv_out<-as.data.frame(array(dim=c(78,21)))
colnames(dfsurv_out)<-c("Time",rcsname)
dfsurv_out$Time<-dfhar_out$Time

dfGOF[,1]<-rcsname
dfGOF$AIC <- -2*dfGOF$LnL + 2*dfGOF$Params

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

###combine surv data ####
dfFig_outSurv = dfsurv_out %>%
  gather(key = "Model", value = "survProp", -Time) %>% mutate(Model = factor(Model))
##plot##
f_surv_out= ggplot() +
  geom_line(data=dfFig_outSurv, aes(x=Time, y=survProp, group=Model, colour=Model), size=1) +
  geom_line(data = ltHaz_out,aes(x=Time,y=surv))+
  scale_color_discrete(name="Model")+
  expand_limits(y=c(0,1),x=c(0,7)) + 
  facet_wrap(~Model,nrow=4)+
  scale_x_continuous(breaks = c(seq(from=0, to=7,by = 1))) +
  ylab("Overall survival") +
  xlab("Time(months)") +
  guides(color = guide_legend(ncol = 1))  +
  theme(legend.position = "bottom") + 
  theme_bw() 
f_surv_out

# ggsave("fig_surv_niv_6.5y_os.png",f_surv_out,width = 18,height = 12,dpi = 800)
# 
# 
# #########################################
# #             Results export            # 
# #########################################
# FP1res[1,1]
# FP1res[2,1]
# FP2res[1,1]
# FP2res[2,1]
# RCSres[1,1]
# RCSres[2,1]
# rp_input
# 
# write.csv(dfsurv,"surv_niv_3_explo.csv") 
# write.csv(dfsurv_out,"surv_niv_6.5.csv") 


