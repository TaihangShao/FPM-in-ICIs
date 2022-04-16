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

niv5y <- read.delim("IPD_pfs_niv_6.5y.txt")
niv5y<-data.frame(niv5y$event,niv5y$time)
niv5y<-rename(niv5y,"censrec"="niv5y.event","recyrs"="niv5y.time")
niv5y$recyrs<-niv5y$recyrs/12
niv5y$rectime<-as.integer(niv5y$recyrs*365.24) 

niv5y = niv5y %>% mutate(censrec = case_when(niv5y$recyrs > 6.5 ~ integer(1), TRUE ~ niv5y$censrec),
                         recyrs = case_when(niv5y$recyrs > 6.5 ~ 6.5, TRUE ~ niv5y$recyrs),
                         rectime = recyrs * 365)

niv3y <- read.delim("IPD_pfs_niv_3y.txt")
niv3y<-data.frame(niv3y$event,niv3y$time)
niv3y<-rename(niv3y,"censrec"="niv3y.event","recyrs"="niv3y.time")
niv3y$recyrs<-niv3y$recyrs/12
niv3y$rectime<-as.integer(niv3y$recyrs*365.24) 

niv3y = niv3y %>% mutate(censrec = case_when(niv3y$recyrs > 3 ~ integer(1), TRUE ~ niv3y$censrec),
                         recyrs = case_when(niv3y$recyrs > 3 ~ 3, TRUE ~ niv3y$recyrs),
                         rectime = recyrs * 365)

ipi5y <- read.delim("IPD_pfs_ipi_6.5y.txt")
ipi5y<-data.frame(ipi5y$event,ipi5y$time)
ipi5y<-rename(ipi5y,"censrec"="ipi5y.event","recyrs"="ipi5y.time")
ipi5y$recyrs<-ipi5y$recyrs/12
ipi5y$rectime<-as.integer(ipi5y$recyrs*365.24) 

ipi5y = ipi5y %>% mutate(censrec = case_when(ipi5y$recyrs > 6.5 ~ integer(1), TRUE ~ ipi5y$censrec),
                         recyrs = case_when(ipi5y$recyrs > 6.5 ~ 6.5, TRUE ~ ipi5y$recyrs),
                         rectime = recyrs * 365)

ipi3y <- read.delim("IPD_pfs_ipi_3y.txt")
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
# ltHaz<-ltHaz0
# ltHaz_out<-ltHaz0_out
# sd_bc<-data.frame("recyrs"=niv3y$recyrs,"censrec"=niv3y$censrec)
# sd_bc_out<-data.frame("recyrs"=niv5y$recyrs,"censrec"=niv5y$censrec)
## ipi ##
ltHaz<-ltHaz1
ltHaz_out<-ltHaz1_out
sd_bc<-data.frame("recyrs"=ipi3y$recyrs,"censrec"=ipi3y$censrec)
sd_bc_out<-data.frame("recyrs"=ipi5y$recyrs,"censrec"=ipi5y$censrec)


###--------------------------------------------------------------------------------------------------------------------------------------------------------####
####----New Data----####
follow_up <- 36
numMod <- 19 # Models considered
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

MyTH2 <- 40 # Time Horizon (years)
MyStep2 <- 12 # Number of obs. per year
MyN2 <- MyTH2*MyStep2 # Total time points (observed & extrapolated)
dfHazEst2 <- array(dim=c(numMod, MyN2))
dfHazEst3 <- array(dim=c(numMod, MyN2))
Newtime2 <- data.frame(Time = seq(from=1/MyStep2, to=MyTH2, by=1/MyStep2), AtRisk = 1)
Newtime2$MyId <- 1:dim(Newtime2)[1]
Newtime2$MyId <- ifelse(Newtime2$MyId > follow_up, follow_up, Newtime2$MyId)  # Random effects: Using last observed ID for extrapolation
Newtime2$EventsL <- 0
Newtime2$EventsL[1:follow_up] <- lag(ltHaz$Events)
Newtime2$EventsL[1] <- 0
Newtime2$EventsL <- ifelse(Newtime2$MyId > follow_up, 0, Newtime2$EventsL) # AR: Using last observed event count for extrapolation
Newtime2$timedelta<-Newtime2$Time[2]-Newtime2$Time[1]
Newtime2$lnTime<-log(Newtime2$Time)

# Also have 1x GOF matrix. Rows = Methods, Columns = Method, LL, AIC
dfGOF <- data.frame(matrix(, nrow=19, ncol=4))
dfGOF_out <- data.frame(matrix(, nrow=19, ncol=4))
colnames(dfGOF) <- c("Model","LnL","Params","AIC")
colnames(dfGOF_out) <- c("Model","LnL","Params","AIC")
# Below is constant for when have to derive log-likelihood
llCons <- sum(ltHaz$Events*log(ltHaz$AtRisk) - log(factorial(ltHaz$Events)))
llCons_out<-sum(ltHaz_out$Events*log(ltHaz_out$AtRisk) - log(factorial(ltHaz_out$Events)))
# Names of models to consider
modnames <- list("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","FP1-1","FP1-2",
                 "FP2-1","FP2-2","RCS1","RCS2","RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")
md<-c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","FP1-1","FP1-2",
     "FP2-1","FP2-2","RCS1","RCS2","RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")
###--------------------------------------------------------------------------------------------------------------------------------------------------------####
#=============================================================================================================#
#=============================================================================================================#

#                                       2y res    +    2y to 5y Extrapolation

#=============================================================================================================#
#=============================================================================================================#

####################################
#######    Standard Dist     #######
####################################
MODi <- 1 # Model index
MyDists <- list("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma")
for (i in 1:7){
  glmTemp <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data = sd_bc, dist = MyDists[[i]])
  dfHazEst[MODi,] <- summary(glmTemp, t=Newtime$Time, type="hazard")[[1]]$est/12
  dfHazEst2[MODi,] <- summary(glmTemp, t=Newtime2$Time, type="hazard")[[1]]$est/12
  ltHaz[[MyDists[[i]]]] <- summary(glmTemp, t=ltHaz$Time, type="hazard")[[1]]$est/12
  dfGOF[MODi,1] <- MyDists[[i]]
  dfGOF[MODi,2] <- sum(ltHaz$Events*log(ltHaz[[MyDists[[i]]]]) - ltHaz[[MyDists[[i]]]]*ltHaz$AtRisk) + llCons
  dfGOF[MODi,3] <- glmTemp$npars
  MODi<-MODi+1
}

## genf dist does not finite(non-finite finite-difference value)

#########################
#######    FP     #######
#########################

#-----FP1 -----
myLnL <- array(dim=8)
myAIC <- array(dim=8)
MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
for (i in 1:7){
  glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)
  myLnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[i] <- extractAIC(glmTemp)[2]
}
### run for 0
glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)
myLnL[8] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
myAIC[8] <- extractAIC(glmTemp)[2]

FP1res <- data.frame(c("-2","-1","-0.5","0","0.5","1","2","3"))
FP1res <- cbind(FP1res,myLnL,myAIC)
colnames(FP1res) <- c("Powers","LnL","AIC")
FP1res <-arrange(FP1res,AIC)
FP1res[1,]
FP1res[2,]

#-----FP2 -----
myLnL <- array(dim=36)
myAIC <- array(dim=36)
MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
index <- 1
for (i in 1:7){
  for (j in 1:7){
    if (j > i) {
      glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][j])+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)# 
      myLnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
      myAIC[index] <- extractAIC(glmTemp)[2]
      index <- index + 1
    }
  }
}
for (i in 1:7) {
  glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][i]*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)# 
  myLnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[index] <- extractAIC(glmTemp)[2]
  index <- index + 1
}

for (i in 1:7) {
  glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)# 
  myLnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[index] <- extractAIC(glmTemp)[2]
  index <- index + 1
}

glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)# 
myLnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
myAIC[index] <- extractAIC(glmTemp)[2]

FP2res <- data.frame(c("-2,-1","-2,-0.5","-2.0.5","-2,1","-2,2","-2,3","-1,-0.5","-1,0.5","-1,1","-1,2","-1,3","-0.5,0.5","-0.5,1","-0.5,2","-0.5,3",
                       "0.5,1","0.5,2","0.5,3","1,2","1,3","2,3","-2,-2","-1,-1","-0.5,-0.5","0.5,0.5","1,1","2,2","3,3",
                       "-2,0","-1,0","-0.5,0","0.5,0","1,0","2,0","3,0","0,0"))
FP2res <- cbind(FP2res,myLnL,myAIC)
colnames(FP2res) <- c("Powers","LnL","AIC")
FP2res <-arrange(FP2res,AIC)
FP2res[1,]
FP2res[2,]


##FP  choose model####
# # Now for the closed-test procedure for function selection
# # Note currently identifying best models manually
# # modify power
# modFP2 <- glm(Events ~ I(lnTime^(-1)) + I(lnTime^1) + offset(log(AtRisk)), family=poisson, data=ltHaz)#
# modFP1 <- glm(Events ~ I(lnTime^(-1))+ offset(log(AtRisk)) , family=poisson, data=ltHaz)
# modLin <- glm(Events ~ lnTime + offset(log(AtRisk)), family=poisson, data=ltHaz)
# modNULL <- glm(Events ~ 1 + offset(log(AtRisk)), family=poisson, data=ltHaz)  
# # Overall association of the outcome with time (Sig result = include time)
# anova(modNULL, modFP2, test="LRT") # p-value 6.798e-15
# # Evidence for non-linearity (Sig result = non-linear model)
# anova(modLin, modFP2, test="LRT") # p-value 2.265e-11
# # Simpler or more complex non-linear model?  (Sig result = FP2, else FP1)
# anova(modFP1, modFP2, test="LRT") # p-value 2.937e-08  
# # Hence use FP2 model

###combine FP results####
#(FP2-best)

modFP1_1 <- glm(cbind(Events,AtRisk-Events) ~ I(Time^(0.5))+ offset(log(timedelta)) , family=binomial(link=cloglog), data=ltHaz)
modFP1_2 <- glm(cbind(Events,AtRisk-Events) ~ I(Time^1)+ offset(log(timedelta)) , family=binomial(link=cloglog), data=ltHaz)
modFP2_1 <- glm(cbind(Events,AtRisk-Events) ~ I(Time^(-0.5)) + I(Time^(-0.5)*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)
modFP2_2 <- glm(cbind(Events,AtRisk-Events) ~ I(Time^(-0.5)) + I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz)

MODi <- 8
dfGOF[MODi,1] <- modnames[[MODi]]
dfGOF[MODi,2] <- (extractAIC(modFP1_1)[2] - 2*extractAIC(modFP1_1)[1])*(-0.5)
dfGOF[MODi,3] <- extractAIC(modFP1_1)[1]
# Hazard estimates
dfHazEst[MODi,] <- predict(modFP1_1, newdata=Newtime, type="response") # Extrapolated
dfHazEst2[MODi,] <- predict(modFP1_1, newdata=Newtime2, type="response") # Extrapolated
ltHaz[modnames[[MODi]]] <- predict(modFP1_1, newdata=ltHaz, type="response")  # Within-sample



MODi<-9
dfGOF[MODi,1] <- modnames[[MODi]]
dfGOF[MODi,2] <- (extractAIC(modFP1_2)[2] - 2*extractAIC(modFP1_2)[1])*(-0.5)
dfGOF[MODi,3] <- extractAIC(modFP1_2)[1]
# Hazard estimates
dfHazEst[MODi,] <- predict(modFP1_2, newdata=Newtime, type="response") # Extrapolated
dfHazEst2[MODi,] <- predict(modFP1_2, newdata=Newtime2, type="response") # Extrapolated
ltHaz[modnames[[MODi]]] <- predict(modFP1_2, newdata=ltHaz, type="response")  # Within-sample

MODi<-10
dfGOF[MODi,1] <- modnames[[MODi]]
dfGOF[MODi,2] <- (extractAIC(modFP2_1)[2] - 2*extractAIC(modFP2_1)[1])*(-0.5)
dfGOF[MODi,3] <- extractAIC(modFP2_1)[1]
# Hazard estimates
dfHazEst[MODi,] <- predict(modFP2_1, newdata=Newtime, type="response") # Extrapolated
dfHazEst2[MODi,] <- predict(modFP2_1, newdata=Newtime2, type="response") # Extrapolated
ltHaz[modnames[[MODi]]] <- predict(modFP2_1, newdata=ltHaz, type="response")  # Within-sample

MODi<-11
dfGOF[MODi,1] <- modnames[[MODi]]
dfGOF[MODi,2] <- (extractAIC(modFP2_2)[2] - 2*extractAIC(modFP2_2)[1])*(-0.5)
dfGOF[MODi,3] <- extractAIC(modFP2_2)[1]
# Hazard estimates
dfHazEst[MODi,] <- predict(modFP2_2, newdata=Newtime, type="response") # Extrapolated
dfHazEst2[MODi,] <- predict(modFP2_2, newdata=Newtime2, type="response") # Extrapolated
ltHaz[modnames[[MODi]]] <- predict(modFP2_2, newdata=ltHaz, type="response")  # Within-sample


##########################
#######    RCS     #######
##########################
MODi <- 12
# First need knot locations for up to 5 internal knots.
# Basing these on equally-spaced percentiles of the observed (uncensored) death times.
bc2 <- subset(ipi3y, censrec==1)
myLnL <- array(dim=5)
myAIC <- array(dim=5)
for (i in 1:5){
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz)
  myLnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[i] <- extractAIC(glmTemp)[2]
}
RCSres <- data.frame(c("One","Two","Three","Four","Five"))
RCSres <- cbind(RCSres,myLnL,myAIC)
colnames(RCSres) <- c("Int.Knots","LnL","AIC")
RCSres<-arrange(RCSres)
RCSres[1,]
RCSres[2,]
i<-1
glmTemp1 <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2),family=binomial(link=cloglog), data=ltHaz)
i<-2
glmTemp2 <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2),family=binomial(link=cloglog), data=ltHaz)

dfGOF[MODi,1] <- modnames[[MODi]]
dfGOF[MODi,2] <- (extractAIC(glmTemp1)[2] - 2*extractAIC(glmTemp1)[1])*(-0.5)
dfGOF[MODi,3] <- extractAIC(glmTemp1)[1]
# Hazard estimates
dfHazEst[MODi,] <- predict(glmTemp1, newdata=Newtime, type="response")
dfHazEst2[MODi,] <- predict(glmTemp1, newdata=Newtime2, type="response")
ltHaz[modnames[[MODi]]] <- predict(glmTemp1, newdata=ltHaz, type="response")

MODi<-13
dfGOF[MODi,1] <- modnames[[MODi]]
dfGOF[MODi,2] <- (extractAIC(glmTemp2)[2] - 2*extractAIC(glmTemp2)[1])*(-0.5)
dfGOF[MODi,3] <- extractAIC(glmTemp2)[1]
# Hazard estimates
dfHazEst[MODi,] <- predict(glmTemp2, newdata=Newtime, type="response")
dfHazEst2[MODi,] <- predict(glmTemp2, newdata=Newtime2, type="response")
ltHaz[modnames[[MODi]]] <- predict(glmTemp2, newdata=ltHaz, type="response")

#########################
#######    RP     #######
#########################
MODi <- 14
MyAIC <- array(dim=c(6,3))
MyScale <- list("hazard","odds","normal")
for (i in 1:3){
  for (j in 0:5){
    fit<-try(MyTemp <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data = sd_bc_out, k = j, scale = MyScale[[i]]))
    if("try-error" %in% class(fit)) {
      MyAIC[[(i-1)*6+1+j]] <- "error"
    }
    else{
      flexsurvspline(Surv(recyrs, censrec) ~ 1, data = sd_bc_out, k = j, scale = MyScale[[i]])
      MyAIC[[(i-1)*6+1+j]] <- (-2*MyTemp$loglik+2*MyTemp$npars)
    }
  }
}
  
MyAICResults <- as.data.frame(cbind(seq(1:6)-1,MyAIC))
colnames(MyAICResults) <- c("Int.Knots","Hazard","Odds","Normal")
best_rp_hazard<-data.frame(MyAICResults$Int.Knots,MyAICResults$Hazard)
colnames(best_rp_hazard)<-c("knots","AIC")
best_rp_odds<-data.frame(MyAICResults$Int.Knots,MyAICResults$Odds)
colnames(best_rp_odds)<-c("knots","AIC")
best_rp_normal<-data.frame(MyAICResults$Int.Knots,MyAICResults$Normal)
colnames(best_rp_normal)<-c("knots","AIC")
best_rp_hazard<-arrange(best_rp_hazard,AIC)
best_rp_odds<-arrange(best_rp_odds,AIC)
best_rp_normal<-arrange(best_rp_normal,AIC)

rp_input<-as.data.frame(array(dim=c(6,2)))
rp_input$V1<- c("hazard","hazard","odds","odds","normal","normal")
rp_input$V2<-c(best_rp_hazard[1,1],
best_rp_hazard[2,1],
best_rp_odds[1,1],
best_rp_odds[2,1],
best_rp_normal[1,1],
best_rp_normal[2,1])
colnames(rp_input)<-c("scale","knots")
rp_input$knots<-as.numeric(rp_input$knots)

for (i in 1:6) {
  rp_scale<-rp_input[i,1]
  rp_k<-rp_input[i,2]
  rpTemp <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data = sd_bc, k = rp_k, scale = rp_scale)
  rpAIC <- (-2*rpTemp$loglik+2*rpTemp$npars)
  dfHazEst[MODi,] <- summary(rpTemp, t=Newtime$Time, type="hazard")[[1]]$est/12
  dfHazEst2[MODi,] <- summary(rpTemp, t=Newtime2$Time, type="hazard")[[1]]$est/12
  ltHaz[modnames[[MODi]]] <- summary(rpTemp, t=ltHaz$Time, type="hazard")[[1]]$est/12
  
  dfGOF[MODi,1] <- modnames[[MODi]]
  dfGOF[MODi,2] <- sum(ltHaz$Events*log(ltHaz[modnames[[MODi]]]) - ltHaz[[modnames[[MODi]]]]*ltHaz$AtRisk) + llCons
  dfGOF[MODi,3] <- rpTemp$npars  
  
  MODi<-MODi+1
}  


#######################################
#             AIC & export            # 
#######################################
dfGOF$AIC <- -2*dfGOF$LnL + 2*dfGOF$Params
write.csv(dfGOF,"niv-3y-pfs-AIC.csv")

FP1res[1,1]
FP1res[2,1]
FP2res[1,1]
FP2res[2,1]
RCSres[1,1]
RCSres[2,1]
rp_input

write.csv(dfHazEst,"dfhazest_ipipfs_3y.csv")
write.csv(ltHaz,"lthaz_ipipfs_3y.csv")
########==============================================================================================###

####-----------------process data and plot 1--------------------####

dfFig <- t(dfHazEst)
colnames(dfFig) <- dfGOF[,1]
dfFig <- cbind(data.frame(Newtime$Time), dfFig)
dfFig = dfFig %>% mutate(Time =Newtime.Time) %>% select(-Newtime.Time) %>%
  gather(key = "Model", value = "Haz", -Time) %>% mutate(Model = factor(Model))
dfFig$Model = fct_recode(dfFig$Model, "lognormal"="lnorm", "loglogistic"="llogis")

###  hazard ratio   ###
# Standard models
fig1a <- ggplot(data=filter(dfFig, Model %in% c("exp","weibull","gamma","lognormal","gompertz","loglogistic","gengamma","genf")), aes(x=Time, y=Haz)) +
  geom_point(data=ltHaz, aes(x=Time, y=hazLT, size=AtRisk), shape = 1) + geom_line(aes(color=Model)) + 
  labs(x="Time", y="Hazard")  + guides(size="none") +
  theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.7), xlim=c(0,3)) 
#+ theme(axis.text.y = element_blank())
#fig1a

# FP models
fig1b <- ggplot(data=filter(dfFig, Model %in% c("FP1-1","FP1-2","FP2-1","FP2-2")), aes(x=Time, y=Haz)) +
  geom_point(data=ltHaz, aes(x=Time, y=hazLT, size=AtRisk), shape = 1) + geom_line(aes(color=Model)) + 
  labs(x="Time", y="Hazard") + guides(size="none") +
  theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.7), xlim=c(0,3)) 
#fig1b

# RCS models
fig1c <- ggplot(data=filter(dfFig, Model %in% c("RCS1","RCS2")), aes(x=Time, y=Haz)) +
  geom_point(data=ltHaz, aes(x=Time, y=hazLT, size=AtRisk), shape = 1) + geom_line(aes(color=Model)) + 
  labs(x="", y="Hazard")+ guides(size="none") +
  theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.7), xlim=c(0,3)) 
#fig1c

# RP models
fig1d <- ggplot(data=filter(dfFig, Model %in% c("RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")), aes(x=Time, y=Haz)) +
  geom_point(data=ltHaz, aes(x=Time, y=hazLT, size=AtRisk), shape = 1) + geom_line(aes(color=Model)) + 
  labs(x="", y="Hazard")  + guides(size="none") +
  theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.7), xlim=c(0,3)) 
#fig1d

fig_haz_ipi_2y<-grid.arrange(fig1a, fig1b, fig1c, fig1d, ncol=2, widths=c(15, 15))
ggsave("fig_haz_ipi_3y_pfs.png",plot = fig_haz_ipi_2y,width = 18,height = 12,dpi = 800)


#############-----------other figs--------------#############
# f1<-ggplot() +
#   geom_line(data=filter(dfFig, Model %in% c("exp","weibull","gamma","lognormal","gompertz","loglogistic","gengamma","genf")),aes(x=Time, y=Haz, group=Model,colour=Model), size=1) +
#   #geom_hline(yintercept=1, lty=2) + 
#   #geom_vline(xintercept=14, lty=2) + 
#   facet_wrap(~Model,nrow=3)+
#   scale_color_discrete(name="Model")+
#   #scale_y_log10(limits = c(0.1, 10), breaks = c(0.1,0.25, 0.5,  1,  2, 4, 10)) +
#   #scale_x_continuous(breaks = c(seq(from=0, to=60,by = 12))) +
#   ylab("Hazard Ratio") +
#   xlab("Time(months)") +
#   theme(legend.position = "bottom") + 
#   theme_bw() 
# f1
# 
# f2<-ggplot() +
#   geom_line(data=filter(dfFig, Model %in% c("FP","RCS","RP")),aes(x=Time, y=Haz, group=Model,colour=Model), size=1) +
#   #geom_hline(yintercept=1, lty=2) + 
#   #geom_vline(xintercept=14, lty=2) + 
#   facet_wrap(~Model,nrow=3)+
#   scale_color_discrete(name="Model")+
#   #scale_y_log10(limits = c(0.1, 10), breaks = c(0.1,0.25, 0.5,  1,  2, 4, 10)) +
#   #scale_x_continuous(breaks = c(seq(from=0, to=60,by = 12))) +
#   ylab("Hazard Ratio") +
#   xlab("Time(months)") +
#   theme(legend.position = "bottom") + 
#   theme_bw() 
# f2


##################------------surv rate calculate  ###  dfhazest   ###---------------##############
#calculate survival over time
dfhar<-t(dfHazEst)
colnames(dfhar) <- dfGOF[,1]
dfhar <- cbind(data.frame(Newtime$Time), dfhar)
dfhar <- rename(dfhar,"Time"="Newtime.Time")

dfsurv<-as.data.frame(array(dim=c(78,20)))
colnames(dfsurv)<-c("Time",modnames)
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
##plot##
f_surv= ggplot() +
  geom_line(data=dfFigSurv, aes(x=Time, y=survProp, group=Model, colour=Model), size=1) +
  geom_point(data = ltHaz_out,aes(x=Time,y=surv),size=1,shape=1)+
  scale_color_discrete(name="Model")+
  expand_limits(y=c(0,1),x=c(0,7)) + 
  facet_wrap(~Model,nrow=4)+
  scale_x_continuous(breaks = c(seq(from=0, to=7,by = 1))) +
  ylab("Overall survival") +
  xlab("Time(months)") +
  guides(color = guide_legend(ncol = 1))  +
  theme(legend.position = "bottom") + 
  theme_bw() 
f_surv

ggsave("fig_surv_ipi_3y_os.png",f_surv,width = 18,height = 12,dpi = 800)

# f_surv1= ggplot() +
#   geom_line(data=dfFigSurv, aes(x=Time, y=survProp, group=Model, colour=Model), size=1) +
#   geom_point(data = ltHaz_out,aes(x=Time,y=surv),size=3,shape=1)+
#   scale_color_discrete(name="Model")+
#   expand_limits(y=c(0,1),x=c(0,6)) + 
#   scale_x_continuous(breaks = c(seq(from=0, to=6,by = 1))) +
#   ylab("Proportion surviving") +
#   xlab("Time(months)") +
#   guides(color = guide_legend(ncol = 1))  +
#   theme(legend.position = "bottom") + 
#   theme_bw() 
# f_surv1

##################------------surv rate calculate  ###  lthaz   ###---------------##############
lthaz_plot<-data.frame(ltHaz$Time,ltHaz[,11:29])
dfsurv_lthaz<-as.data.frame(array(dim=c(36,20)))
colnames(dfsurv_lthaz)<-c("Time",modnames)
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
##plot##
f_surv2= ggplot() +
  geom_line(data=dfFigSurv_lthaz, aes(x=Time, y=survProp, group=Model, colour=Model), size=1) +
  geom_line(data = ltHaz,aes(x=Time,y=surv),size=0.5)+
  scale_color_discrete(name="Model")+
  expand_limits(y=c(0,1),x=c(0,3)) + 
  facet_wrap(~Model,nrow=4)+
  scale_x_continuous(breaks = c(seq(from=0, to=3,by = 0.5))) +
  ylab("Overall survival") +
  xlab("Time(months)") +
  guides(color = guide_legend(ncol = 1))  +
  theme(legend.position = "bottom") + 
  theme_bw() 
f_surv2

ggsave("fig_surv_ipi_3y_os_3ys.png",f_surv2,width = 18,height = 12,dpi = 800)

# f_surv3= ggplot() +
#   geom_line(data=dfFigSurv_lthaz, aes(x=Time, y=survProp, group=Model, colour=Model), size=1) +
#   geom_line(data = ltHaz,aes(x=Time,y=surv),size=1,colour="black")+
#   scale_color_discrete(name="Model")+
#   expand_limits(y=c(0,1),x=c(0,2)) + 
#   scale_x_continuous(breaks = c(seq(from=0, to=2,by = 0.5))) +
#   scale_y_continuous(breaks = c(seq(from=0.5, to=1,by = 0.1)),limits = c(0.5,1)) +
#   ylab("Proportion surviving") +
#   xlab("Time(months)") +
#   guides(color = guide_legend(ncol = 1))  +
#   theme(legend.position = "bottom") + 
#   theme_bw() 
# f_surv3


########==============================================================================================###

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

dfLE_dfHazEst<-dfFig
lthaz_plot<-rename(lthaz_plot,"Time"="ltHaz.Time")
dfLE_lthaz<-lthaz_plot
dfLE_lthaz = dfLE_lthaz %>% 
  gather(key = "Model", value = "Haz", -Time) %>% mutate(Model = factor(Model))
dfLE_lthaz$Model = fct_recode(dfLE_lthaz$Model, "lognormal"="lnorm", "loglogistic"="llogis")


LE_5y_predict = dfLE_dfHazEst %>% group_by(Model) %>% summarise(LE = fun_LE(Haz, Time))
LE_2y = dfLE_lthaz %>% group_by(Model) %>% summarise(LE = fun_LE(Haz, Time)) 

ggplot(data=LE_5y_predict,mapping=aes(x=Model,y=LE,fill=LE,group=factor(1)))+
  geom_bar(stat="identity",width=0.8)+
  geom_text(aes(label = round(LE,digits = 3), vjust = -0.8, hjust = 0.5), show.legend = TRUE) +
  expand_limits(y=c(4,6))+
  scale_y_continuous(breaks = c(seq(from=4, to=6,by = 0.2)))+
  coord_cartesian(ylim = c(4,6))


ggplot(data=LE_2y,mapping=aes(x=Model,y=LE,fill=LE,group=factor(1)))+
  geom_bar(stat="identity",width=0.8)+
  geom_text(aes(label = round(LE,digits = 3), vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  expand_limits(y=c(2.3,2.6))+
  scale_y_continuous(breaks = c(seq(from=2.3, to=2.6,by = 0.02)))+
  coord_cartesian(ylim = c(2.3,2.6))


### SSE calculate haz #####
#KM surv data combined with predicted data
haz_ref<-ltHaz_out$hazLT
haz_ref<-haz_ref[37:78]

SSE_haz<-as.data.frame(t(dfHazEst))
colnames(SSE_haz)<-modnames
SSE_haz<-SSE_haz[37:78,]
#SD1
SSE<-data.frame(matrix(nrow=19,ncol=2))
colnames(SSE)<-c("modnames","SSE")
SSE$modnames<-md


for (j in 1:19){
  SSE_temp<-0
  for (i in 1:42){
    SSE_temp <- SSE_temp + (haz_ref[i]-SSE_haz[i,j])^2
  }
  SSE[j,2]<-SSE_temp
}
SSE
Res_SSE<-SSE
Res_SSE$SSE<-SSE$SSE*10000
Res_SSE

Res_SSE<-arrange(Res_SSE,SSE)

sse_plot<-ggplot(data=Res_SSE,mapping=aes(x=modnames,y=SSE,fill=modnames,group=factor(1)))+
  geom_bar(stat="identity",width=0.8)+
  geom_text(aes(label = round(SSE,digits = 2), vjust = -0.8, hjust = 0.5), show.legend = FALSE)

ggsave("fig_sse_ipi_3y_pfs.png",plot = sse_plot,width = 18,height = 12,dpi = 800)


### SAE calculate surv #####
surv_ref<-ltHaz_out$surv
surv_ref<-surv_ref[37:78]

SAE_surv<-dfsurv[37:78,2:20]
colnames(SAE_surv)<-modnames
#SD1
SAE<-data.frame(matrix(nrow=19,ncol=2))
colnames(SAE)<-c("modnames","SAE")
SAE$modnames<-md


for (j in 1:20){
  SAE_temp<-0
  for (i in 1:42){
    SAE_temp <- SAE_temp + abs(surv_ref[i]-SAE_surv[i,j])
  }
  SAE[j,2]<-SAE_temp
}
SAE
Res_SAE<-SAE
Res_SAE

Res_SAE<-arrange(Res_SAE,SAE)

sae_plot<-ggplot(data=Res_SAE,mapping=aes(x=modnames,y=SAE,fill=modnames,group=factor(1)))+
  geom_bar(stat="identity",width=0.8)+
  geom_text(aes(label = round(SAE,digits = 2), vjust = -0.8, hjust = 0.5), show.legend = FALSE)
ggsave("fig_sae_niv_3y_os.png",plot = sae_plot,width = 18,height = 12,dpi = 800)


#=============================================================================================================#

#                                               3y lifetime data

#=============================================================================================================#

dfhar2<-t(dfHazEst2)
colnames(dfhar2) <- dfGOF[,1]
dfhar2 <- cbind(data.frame(Newtime2$Time), dfhar2)
dfhar2 <- rename(dfhar2,"Time"="Newtime2.Time")

dfsurvLE<-as.data.frame(array(dim=c(480,20)))
colnames(dfsurvLE)<-c("Time",modnames)
dfsurvLE$Time<-dfhar2$Time
###  20 model hazard  ###
for (i in 1:19) {
  temp_dfhaz<-data.frame(dfhar2$Time,dfhar2[,i+1])
  colnames(temp_dfhaz)<-c("time","haz")
  temp_dfsurvLE<-temp_dfhaz %>% 
    dplyr::arrange(time) %>% 
    dplyr::mutate(cumhaz = cumsum(haz)) %>% 
    dplyr::mutate(survProp = exp(-1*cumhaz))
  dfsurvLE[,i+1]<-temp_dfsurvLE$survProp
}

write.csv(dfsurvLE,"ipi-pfs-3y-LE.csv")


#=============================================================================================================#
#=============================================================================================================#

#                                               5y Res

#=============================================================================================================#
#=============================================================================================================#





####################################
#######    Standard Dist     #######
####################################
MODi <- 1 # Model index
MyDists <- list("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma")
for (i in 1:7){
  glmTemp <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data = sd_bc_out, dist = MyDists[[i]])
  ltHaz_out[[MyDists[[i]]]] <- summary(glmTemp, t=ltHaz_out$Time, type="hazard")[[1]]$est/12
  dfHazEst3[MODi,] <- summary(glmTemp, t=Newtime2$Time, type="hazard")[[1]]$est/12
  dfGOF_out[MODi,1] <- MyDists[[i]]
  dfGOF_out[MODi,2] <- sum(ltHaz_out$Events*log(ltHaz_out[[MyDists[[i]]]]) - ltHaz_out[[MyDists[[i]]]]*ltHaz_out$AtRisk) + llCons_out
  dfGOF_out[MODi,3] <- glmTemp$npars
  MODi<-MODi+1
}

#########################
#######    FP     #######
#########################

#-----FP1 -----
myLnL <- array(dim=8)
myAIC <- array(dim=8)
MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
for (i in 1:7){
  glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_out)
  myLnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[i] <- extractAIC(glmTemp)[2]
}
### run for 0
glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_out)
myLnL[8] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
myAIC[8] <- extractAIC(glmTemp)[2]

FP1res <- data.frame(c("-2","-1","-0.5","0","0.5","1","2","3"))
FP1res <- cbind(FP1res,myLnL,myAIC)
colnames(FP1res) <- c("Powers","LnL","AIC")
FP1res <-arrange(FP1res,AIC)
FP1res[1,]
FP1res[2,]

#-----FP2 -----
myLnL <- array(dim=36)
myAIC <- array(dim=36)
MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
index <- 1
for (i in 1:7){
  for (j in 1:7){
    if (j > i) {
      glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][j])+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_out)# 
      myLnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
      myAIC[index] <- extractAIC(glmTemp)[2]
      index <- index + 1
    }
  }
}
for (i in 1:7) {
  glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][i]*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_out)# 
  myLnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[index] <- extractAIC(glmTemp)[2]
  index <- index + 1
}

for (i in 1:7) {
  glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_out)# 
  myLnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[index] <- extractAIC(glmTemp)[2]
  index <- index + 1
}

glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_out)# 
myLnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
myAIC[index] <- extractAIC(glmTemp)[2]

FP2res <- data.frame(c("-2,-1","-2,-0.5","-2.0.5","-2,1","-2,2","-2,3","-1,-0.5","-1,0.5","-1,1","-1,2","-1,3","-0.5,0.5","-0.5,1","-0.5,2","-0.5,3",
                       "0.5,1","0.5,2","0.5,3","1,2","1,3","2,3","-2,-2","-1,-1","-0.5,-0.5","0.5,0.5","1,1","2,2","3,3",
                       "-2,0","-1,0","-0.5,0","0.5,0","1,0","2,0","3,0","0,0"))
FP2res <- cbind(FP2res,myLnL,myAIC)
colnames(FP2res) <- c("Powers","LnL","AIC")
FP2res <-arrange(FP2res,AIC)
FP2res[1,]
FP2res[2,]


##FP  choose model####
# # Now for the closed-test procedure for function selection
# # Note currently identifying best models manually
# # modify power
# modFP2 <- glm(Events ~ I(lnTime^(-1)) + I(lnTime^1) + offset(log(AtRisk)), family=poisson, data=ltHaz_out)#
# modFP1 <- glm(Events ~ I(lnTime^(-1))+ offset(log(AtRisk)) , family=poisson, data=ltHaz_out)
# modLin <- glm(Events ~ lnTime + offset(log(AtRisk)), family=poisson, data=ltHaz_out)
# modNULL <- glm(Events ~ 1 + offset(log(AtRisk)), family=poisson, data=ltHaz_out)  
# # Overall association of the outcome with time (Sig result = include time)
# anova(modNULL, modFP2, test="LRT") # p-value 6.798e-15
# # Evidence for non-linearity (Sig result = non-linear model)
# anova(modLin, modFP2, test="LRT") # p-value 2.265e-11
# # Simpler or more complex non-linear model?  (Sig result = FP2, else FP1)
# anova(modFP1, modFP2, test="LRT") # p-value 2.937e-08  
# # Hence use FP2 model

###combine FP results####
#(FP2-best)

modFP1_1 <- glm(cbind(Events,AtRisk-Events) ~ I(Time^(0.5))+ offset(log(timedelta)) , family=binomial(link=cloglog), data=ltHaz_out)
modFP1_2 <- glm(cbind(Events,AtRisk-Events) ~ I(Time^1)+ offset(log(timedelta)) , family=binomial(link=cloglog), data=ltHaz_out)
modFP2_1 <- glm(cbind(Events,AtRisk-Events) ~ I(Time^(-0.5)) + I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_out)
modFP2_2 <- glm(cbind(Events,AtRisk-Events) ~ I(Time^(-0.5)) + I(Time^(-0.5)*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_out)

MODi <- 8
dfGOF_out[MODi,1] <- modnames[[MODi]]
dfGOF_out[MODi,2] <- (extractAIC(modFP1_1)[2] - 2*extractAIC(modFP1_1)[1])*(-0.5)
dfGOF_out[MODi,3] <- extractAIC(modFP1_1)[1]
# Hazard estimates
dfHazEst3[MODi,] <- predict(modFP1_1, newdata=Newtime2, type="response") # Extrapolated
ltHaz_out[modnames[[MODi]]] <- predict(modFP1_1, newdata=ltHaz_out, type="response")  # Within-sample

predict(modFP1_1, newdata=ltHaz_out_newdata, type="response")
predict(modFP1_1, newdata=ltHaz_out_olddata, type="response")

MODi<-9
dfGOF_out[MODi,1] <- modnames[[MODi]]
dfGOF_out[MODi,2] <- (extractAIC(modFP1_2)[2] - 2*extractAIC(modFP1_2)[1])*(-0.5)
dfGOF_out[MODi,3] <- extractAIC(modFP1_2)[1]
# Hazard estimates
dfHazEst3[MODi,] <- predict(modFP1_2, newdata=Newtime2, type="response") # Extrapolated
ltHaz_out[modnames[[MODi]]] <- predict(modFP1_2, newdata=ltHaz_out, type="response")  # Within-sample

MODi<-10
dfGOF_out[MODi,1] <- modnames[[MODi]]
dfGOF_out[MODi,2] <- (extractAIC(modFP2_1)[2] - 2*extractAIC(modFP2_1)[1])*(-0.5)
dfGOF_out[MODi,3] <- extractAIC(modFP2_1)[1]
# Hazard estimates
dfHazEst3[MODi,] <- predict(modFP2_1, newdata=Newtime2, type="response") # Extrapolated
ltHaz_out[modnames[[MODi]]] <- predict(modFP2_1, newdata=ltHaz_out, type="response")  # Within-sample

MODi<-11
dfGOF_out[MODi,1] <- modnames[[MODi]]
dfGOF_out[MODi,2] <- (extractAIC(modFP2_2)[2] - 2*extractAIC(modFP2_2)[1])*(-0.5)
dfGOF_out[MODi,3] <- extractAIC(modFP2_2)[1]
# Hazard estimates
dfHazEst3[MODi,] <- predict(modFP2_2, newdata=Newtime2, type="response") # Extrapolated
ltHaz_out[modnames[[MODi]]] <- predict(modFP2_2, newdata=ltHaz_out, type="response")  # Within-sample


##########################
#######    RCS     #######
##########################
MODi <- 12
# First need knot locations for up to 5 internal knots.
# Basing these on equally-spaced percentiles of the observed (uncensored) death times.
bc2_out <- subset(ipi5y, censrec==1)
myLnL <- array(dim=5)
myAIC <- array(dim=5)
for (i in 1:5){
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2_out$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_out)
  myLnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  myAIC[i] <- extractAIC(glmTemp)[2]
}
RCSres <- data.frame(c("One","Two","Three","Four","Five"))
RCSres <- cbind(RCSres,myLnL,myAIC)
colnames(RCSres) <- c("Int.Knots","LnL","AIC")
RCSres<-arrange(RCSres)
RCSres[1,]
RCSres[2,]
i<-1
glmTemp1 <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2_out$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2),family=binomial(link=cloglog), data=ltHaz_out)
i<-2
glmTemp2 <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2_out$recyrs, seq(from=0, to=1, by=1/(1+i))), length=i+2),family=binomial(link=cloglog), data=ltHaz_out)

dfGOF_out[MODi,1] <- modnames[[MODi]]
dfGOF_out[MODi,2] <- (extractAIC(glmTemp1)[2] - 2*extractAIC(glmTemp1)[1])*(-0.5)
dfGOF_out[MODi,3] <- extractAIC(glmTemp1)[1]
# Hazard estimates
dfHazEst3[MODi,] <- predict(glmTemp1, newdata=Newtime2, type="response")
ltHaz_out[modnames[[MODi]]] <- predict(glmTemp1, newdata=ltHaz_out, type="response")

MODi<-13
dfGOF_out[MODi,1] <- modnames[[MODi]]
dfGOF_out[MODi,2] <- (extractAIC(glmTemp2)[2] - 2*extractAIC(glmTemp2)[1])*(-0.5)
dfGOF_out[MODi,3] <- extractAIC(glmTemp2)[1]
# Hazard estimates
dfHazEst3[MODi,] <- predict(glmTemp2, newdata=Newtime2, type="response")
ltHaz_out[modnames[[MODi]]] <- predict(glmTemp2, newdata=ltHaz_out, type="response")

#########################
#######    RP     #######
#########################
MODi <- 14
MyAIC <- array(dim=c(6,3))
MyBIC <- array(dim=c(6,3))
MyScale <- list("hazard","odds","normal")
for (i in 1:3){
  for (j in 0:5){
    fit<-try(MyTemp <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data = sd_bc_out, k = j, scale = MyScale[[i]]))
    if("try-error" %in% class(fit)) {
      MyAIC[[(i-1)*6+1+j]] <- "error"
    }
    else{
      flexsurvspline(Surv(recyrs, censrec) ~ 1, data = sd_bc_out, k = j, scale = MyScale[[i]])
      MyAIC[[(i-1)*6+1+j]] <- (-2*MyTemp$loglik+2*MyTemp$npars)
    }
  }
}
MyAICResults <- as.data.frame(cbind(seq(1:6)-1,MyAIC))
colnames(MyAICResults) <- c("Int.Knots","Hazard","Odds","Normal")
best_rp_hazard<-data.frame(MyAICResults$Int.Knots,MyAICResults$Hazard)
colnames(best_rp_hazard)<-c("knots","AIC")
best_rp_odds<-data.frame(MyAICResults$Int.Knots,MyAICResults$Odds)
colnames(best_rp_odds)<-c("knots","AIC")
best_rp_normal<-data.frame(MyAICResults$Int.Knots,MyAICResults$Normal)
colnames(best_rp_normal)<-c("knots","AIC")
best_rp_hazard<-arrange(best_rp_hazard,AIC)
best_rp_odds<-arrange(best_rp_odds,AIC)
best_rp_normal<-arrange(best_rp_normal,AIC)

rp_input<-as.data.frame(array(dim=c(6,2)))
rp_input$V1<- c("hazard","hazard","odds","odds","normal","normal")
rp_input$V2<-c(best_rp_hazard[1,1],
               best_rp_hazard[2,1],
               best_rp_odds[1,1],
               best_rp_odds[2,1],
               best_rp_normal[1,1],
               best_rp_normal[2,1])
colnames(rp_input)<-c("scale","knots")
rp_input$knots<-as.numeric(rp_input$knots)

for (i in 1:6) {
  rp_scale<-rp_input[i,1]
  rp_k<-rp_input[i,2]
  rpTemp <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data = sd_bc_out, k = rp_k, scale = rp_scale)
  rpAIC <- (-2*rpTemp$loglik+2*rpTemp$npars)
  ltHaz_out[modnames[[MODi]]] <- summary(rpTemp, t=ltHaz_out$Time, type="hazard")[[1]]$est/12
  dfHazEst3[MODi,] <- summary(rpTemp, t=Newtime2$Time, type="hazard")[[1]]$est/12
  
  dfGOF_out[MODi,1] <- modnames[[MODi]]
  dfGOF_out[MODi,2] <- sum(ltHaz_out$Events*log(ltHaz_out[modnames[[MODi]]]) - ltHaz_out[[modnames[[MODi]]]]*ltHaz_out$AtRisk) + llCons_out
  dfGOF_out[MODi,3] <- rpTemp$npars  
  
  MODi<-MODi+1
}  


#######################################
#             AIC & export            # 
#######################################
dfGOF_out$AIC <- -2*dfGOF_out$LnL + 2*dfGOF_out$Params
write.csv(dfGOF_out,"ipi-6.5y-pfs-AIC.csv")

write.csv(ltHaz_out,"lthaz_ipipfs_6.5y.csv")

########==============================================================================================###

####-----------------process data and plot 1--------------------####

dfFig_out <- ltHaz_out[,11:29]
colnames(dfFig_out) <- dfGOF_out[,1]
dfFig_out <- cbind(data.frame(ltHaz_out$Time), dfFig_out)
dfFig_out = dfFig_out %>% mutate(Time =ltHaz_out.Time) %>% select(-ltHaz_out.Time) %>%
  gather(key = "Model", value = "Haz", -Time) %>% mutate(Model = factor(Model))
dfFig_out$Model = fct_recode(dfFig_out$Model, "lognormal"="lnorm", "loglogistic"="llogis")

###  hazard ratio   ###
# Standard models
fig1a <- ggplot(data=filter(dfFig_out, Model %in% c("exp","weibull","gamma","lognormal","gompertz","loglogistic","gengamma","genf")), aes(x=Time, y=Haz)) +
  geom_point(data=ltHaz_out, aes(x=Time, y=hazLT, size=AtRisk), shape = 1) + geom_line(aes(color=Model)) + 
  labs(x="Time", y="Hazard")  + guides(size="none") +
  theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.6), xlim=c(0,7)) 
#+ theme(axis.text.y = element_blank())
#fig1a

# FP models
fig1b <- ggplot(data=filter(dfFig_out, Model %in% c("FP1-1","FP1-2","FP2-1","FP2-2")), aes(x=Time, y=Haz)) +
  geom_point(data=ltHaz_out, aes(x=Time, y=hazLT, size=AtRisk), shape = 1) + geom_line(aes(color=Model)) + 
  labs(x="Time", y="Hazard") + guides(size="none") +
  theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.6), xlim=c(0,7)) 
#fig1b

# RCS models
fig1c <- ggplot(data=filter(dfFig_out, Model %in% c("RCS1","RCS2")), aes(x=Time, y=Haz)) +
  geom_point(data=ltHaz_out, aes(x=Time, y=hazLT, size=AtRisk), shape = 1) + geom_line(aes(color=Model)) + 
  labs(x="", y="Hazard")+ guides(size="none") +
  theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.6), xlim=c(0,7)) 
#fig1c

# RP models
fig1d <- ggplot(data=filter(dfFig_out, Model %in% c("RP-hazard-1","RP-hazard-2","RP-odds-1","RP-odds-2","RP-normal-1","RP-normal-2")), aes(x=Time, y=Haz)) +
  geom_point(data=ltHaz_out, aes(x=Time, y=hazLT, size=AtRisk), shape = 1) + geom_line(aes(color=Model)) + 
  labs(x="", y="Hazard")  + guides(size="none") +
  theme(legend.position = "bottom") + coord_cartesian(ylim=c(0,0.6), xlim=c(0,7)) 
#fig1d

fig_haz_niv_5y<-grid.arrange(fig1a, fig1b, fig1c, fig1d, ncol=2, widths=c(15, 15))
ggsave("fig_haz_ipi_6.5y_pfs.png",plot = fig_haz_niv_5y,width = 18,height = 12,dpi = 800)



##################------------surv rate calculate  ###  dfhazest   ###---------------##############
#calculate survival over time
dfhar_out<-ltHaz_out[,11:29]
colnames(dfhar_out) <- dfGOF_out[,1]
dfhar_out <- cbind(data.frame(ltHaz_out$Time), dfhar_out)
dfhar_out <- rename(dfhar_out,"Time"="ltHaz_out.Time")

dfsurv_out<-as.data.frame(array(dim=c(78,20)))
colnames(dfsurv_out)<-c("Time",modnames)
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

###combine surv data ####
dfFig_outSurv = dfsurv_out %>%
  gather(key = "Model", value = "survProp", -Time) %>% mutate(Model = factor(Model))
##plot##
f_surv_out= ggplot() +
  geom_line(data=dfFig_outSurv, aes(x=Time, y=survProp, group=Model, colour=Model), size=1) +
  geom_point(data = ltHaz_out,aes(x=Time,y=surv),size=1,shape=1)+
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

ggsave("fig_surv_ipi_6.5y_pfs.png",f_surv_out,width = 18,height = 12,dpi = 800)

#########################################
#             Results export            # 
#########################################
FP1res[1,1]
FP1res[2,1]
FP2res[1,1]
FP2res[2,1]
RCSres[1,1]
RCSres[2,1]
rp_input

write.csv(dfsurv,"surv_ipi_3_explo.csv") 
write.csv(dfsurv_out,"surv_ipi_6.5.csv") 

#=============================================================================================================#

#                                               6.5y lifetime data

#=============================================================================================================#

dfhar3<-t(dfHazEst3)
colnames(dfhar3) <- dfGOF[,1]
dfhar3 <- cbind(data.frame(Newtime2$Time), dfhar3)
dfhar3 <- rename(dfhar3,"Time"="Newtime2.Time")

dfsurvLE2<-as.data.frame(array(dim=c(480,20)))
colnames(dfsurvLE2)<-c("Time",modnames)
dfsurvLE2$Time<-dfhar3$Time
###  20 model hazard  ###
for (i in 1:19) {
  temp_dfhaz<-data.frame(dfhar3$Time,dfhar3[,i+1])
  colnames(temp_dfhaz)<-c("time","haz")
  temp_dfsurvLE<-temp_dfhaz %>% 
    dplyr::arrange(time) %>% 
    dplyr::mutate(cumhaz = cumsum(haz)) %>% 
    dplyr::mutate(survProp = exp(-1*cumhaz))
  dfsurvLE2[,i+1]<-temp_dfsurvLE$survProp
}

write.csv(dfsurvLE2,"ipi-pfs-6.5y-LE.csv")
