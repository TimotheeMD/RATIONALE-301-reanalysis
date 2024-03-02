## How to reconstruct data from published KM curves

### everything explained here: https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-021-01308-8

### install package
install.packages("IPDfromKM")
library(IPDfromKM)

## First step, go on digitilzest and for each arm (important, use 100 as the maximum in y axis), you export an csv or excel files on your Desktop (EXP for the experimental arm, and CON for the control arm)

## first the experimental arm
##read the data from csv file : 
E <- read.csv(file="/Users/Desktop/EXP.csv", sep = ",", header=T)
## or just use the Import Dataset function in R Studio to import excel files. 
## then put your EXP data in E : 
E <- EXP
## attribute new column names to E
colnames(E) <-c("time", "survival probability")

## the same with the control arm (called C)
C <- read.csv(file="/Users/Desktop/CON.csv", sep = ",", header=T)
C <- CON
colnames(C) <-c("time", "survival probability")

##time risk, separate by comma (no space allowed)
trisk <- c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54)

## number at risk experimental arm (no space allowed)
nrisk.E <- c(342,307,259,228,191,170,155,137,126,111,101,98,77,53,33,18,4,0,0)

## number at risk control arm (no space allowed)
nrisk.C <- c(332,291,247,208,179,147,136,113,96,84,77,66,52,39,29,13,4,1,0)

## preprocess the data
pre_E <- preprocess(dat=E, trisk=trisk, nrisk=nrisk.E, maxy=100)
pre_C <- preprocess(dat=C, trisk=trisk, nrisk=nrisk.C, maxy=100)

## then individual patient data = IPD, (experimental treatment is 1, control group treatment is 0)
est_E_OS <- getIPD(prep=pre_E, armID=1, tot.events = NULL)
est_C_OS <- getIPD(prep=pre_C, armID=0, tot.events = NULL)

## you can isolate the IPD part (time, status, treat) from est_E_OS and est_C_OS
est_E_OS$IPD
est_C_OS$IPD

## the summary function allow to have the data of events, censoring, etc. 
summary(est_E_OS)
summary (est_C_OS)

#### SURVIVAL ANALYIS
install.packages("rms")
library(rms)

## here you can create a KM curve based on the novel IPD data generated based on the data that have been digitilzed
## you can also recapitulate the Cox analysis. 
surv_data<- rbind(est_E_OS$IPD, est_C_OS$IPD)
names(surv_data) <- c("time", "status", "arm")
x <- surv_data$time
y <- surv_data$status
z <- surv_data$arm

## Kaplan-Meier Curve
par(mar = c(1, 1, 1, 1), xaxs = "i", yaxs = "i")
plot(survfit(Surv(x,y)~z), xlim = c(0, 30), ylim = c(0, 1), 
     col=c("#66CC66", "#6699CC"), lwd=3, xaxt='n', bty = "l", 
     las = 1, cex.axis = 1.5, tcl  = -0.5)
axis(side=1, at=seq(0, 30, by=3), cex.axis = 1.5)

## Cox, HR and p-value
summary(coxph(Surv(x,y)~z))