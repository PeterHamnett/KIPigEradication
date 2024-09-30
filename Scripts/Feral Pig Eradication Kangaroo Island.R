# Peter Hamnett, Frédérik Saltré and Corey Bradshaw
# Global Ecology, Flinders University — globalecologyflinders.com
# feral pig (sus scrofa) eradication on Kangaroo Island
# https://github.com/PeterHamnett/KIPigEradication
### update 16/9/2024
# base stochastic model modified from https://github.com/KathrynVenning/FeralCatEradication


## remove everything
rm(list = ls())
library(readr)
library(dplyr)
library(plotly)
library(tidyr)
library(ggpubr)


# functions
# Set functions
AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# if you've downloaded this script from https://github.com/PeterHamnett/KIPigEradication/settings...
## change file paths to your own file locations
setwd("C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/Scripts")
source("matrixOperators.r")



#### Functional Response Estimates ####

pigs_effort_no_LU <- read_csv("pigeffort3.csv")
effort <- data.frame(pigs_effort_no_LU)
head(effort)
effort.orig <- effort 

# tidy effort data
effort_dates <- effort

# convert Date column from chr to date format
effort_dates$Date <- as.Date(effort$Date, format = "%d/%m/%y")
effort <- effort_dates

# sort by date so oldest records are at the top of the data frame
effort <- effort[order(effort$Date),] 

# remove records where effort is na
effort <- na.omit(effort)

# create a new column with the inverse of pigs/hr i.e., hrs per pig
effort$hrsPig <- effort$efforthrs/effort$numkilled
head(effort)

#create continuous variables for model fitting
lpropRemaining.cont <- log(seq(0.01, 1, 0.01))
propRemaining.cont <- exp(lpropRemaining.cont)


# create data subsets for each control technique
shot_efrt <- effort %>% filter (controlType == "shot")
shot_efrt <- shot_efrt %>% filter (operator != "PJ") # excludes shooting by PJ = 8 out of 82 shooting events (all events recorded as 100hrs - imprecise and invalid for fitting functional response model)
TAAC_efrt <-effort %>% filter (controlType == "TAAC")
TAAC_efrt <- TAAC_efrt[which(is.infinite(TAAC_efrt$hrsPig)==F),] # remove entries where TAC effort was infinite i.e., flights where no pigs were killed (0 kills/hours = infinite value)
trap_efrt <-effort %>% filter (controlType == "trapped")
trap_efrt <-trap_efrt %>% filter (org != "DEW") # excludes DEW outliers 2 out of 10 trapping events
poison_efrt <- effort %>% filter (controlType == "poisoned")
poison_efrt <- poison_efrt %>% filter (org != "LB") # excludes outliers - poisoning done by the landscape board 5 out of 15 poisoning events

# plot logarithmic and exponential effort~prop remaining relationships for each control type and fit models
par(mfrow=c(2,2))
par(mar = c(3,3,3,3))
# Thesis figure 4 - results - Functional response estimates

#shot
plot(hrsPig~propRemaining, data = shot_efrt, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,110))
mtext("a", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
  # logarithmic fit
    fitShot <- lm(hrsPig ~ log(propRemaining), data=shot_efrt, na.action =) # logarithmic fit
    summary(fitShot)
    linreg.ER(log(shot_efrt$propRemaining),shot_efrt$hrsPig)
    hrsPigShot.pred <- coef(fitShot)[1] + coef(fitShot)[2]*lpropRemaining.cont
    lines(exp(lpropRemaining.cont), hrsPigShot.pred, lty=2, col="red")
  # exponential: y = a * exp(-b*x)
    s.param.init <- c(200, 5)
    fitShotExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                      data = shot_efrt,
                      #algorithm = "port",
                      start = c(a = s.param.init[1], b = s.param.init[2]),
                      trace = TRUE,      
                      nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
    fitShotExp.summ <- summary(fitShotExp)
    hrsPigShotExp.pred <- coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * propRemaining.cont)
    lines(propRemaining.cont,hrsPigShotExp.pred,lty=3,lwd=2,col="blue")

## TAAC
plot(hrsPig~propRemaining, data = TAAC_efrt, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,5))
mtext("b", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
  # logarithmic fit
    fitTAAC <- lm(hrsPig ~ log(propRemaining), data=TAAC_efrt) # logarithmic fit
    summary(fitTAAC)
    linreg.ER(log(TAAC_efrt$propRemaining),TAAC_efrt$hrsPig)
    hrsPigTAAC.pred <- coef(fitTAAC)[1] + coef(fitTAAC)[2]*lpropRemaining.cont
    lines(exp(lpropRemaining.cont), hrsPigTAAC.pred, lty=2, col="red")
  # exponential: y = a * exp(-b*x)
    s.param.init <- c(5, 1)
    fitTAACExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                      data = TAAC_efrt,
                      algorithm = "port",
                      start = c(a = s.param.init[1], b = s.param.init[2]),
                      trace = TRUE,      
                      nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
    fitTAACExp.summ <- summary(fitTAACExp)
    hrsPigTAACExp.pred <- coef(fitTAACExp)[1] * exp(-coef(fitTAACExp)[2] * propRemaining.cont)
    lines(propRemaining.cont,hrsPigTAACExp.pred,lty=3,lwd=2,col="blue")

#trapped
plot(hrsPig~propRemaining, data = trap_efrt, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,40))
mtext("c", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
  # logarithmic fit
  fitTrap <- lm(hrsPig ~ log(propRemaining), data=trap_efrt, na.action = ) # logarithmic fit
  summary(fitTrap)
  linreg.ER(log(trap_efrt$propRemaining),trap_efrt$hrsPig)
  hrsPigTrap.pred <- coef(fitTrap)[1] + coef(fitTrap)[2]*lpropRemaining.cont
  lines(exp(lpropRemaining.cont), hrsPigTrap.pred, lty=2, col="red")
  # exponential: y = a * exp(-b*x)
  s.param.init <- c(12, 2)
  fitTrapExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                    data = trap_efrt,
                    algorithm = "port",
                    start = c(a = s.param.init[1], b = s.param.init[2]),
                    trace = TRUE,      
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
  fitTrapExp.summ <- summary(fitTrapExp)
  hrsPigTrapExp.pred <- coef(fitTrapExp)[1] * exp(-coef(fitTrapExp)[2] * propRemaining.cont)
  lines(propRemaining.cont,hrsPigTrapExp.pred,lty=3,lwd=2,col="blue")

#poisoned
plot(hrsPig~propRemaining, data = poison_efrt, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,75))
mtext("d", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
  # logarithmic fit: y = a + b*log(x)
    fitPoison <- lm(hrsPig ~ log(propRemaining), data=poison_efrt) # logarithmic fit
    summary(fitPoison)
    linreg.ER(log(poison_efrt$propRemaining),poison_efrt$hrsPig)
    lpropRemaining.cont <- log(seq(0.01, 1, 0.01))
    hrsPigPoison.pred <- coef(fitPoison)[1] + coef(fitPoison)[2]*lpropRemaining.cont
    lines(exp(lpropRemaining.cont), hrsPigPoison.pred, lty=2, col="red")
  # exponential: y = a * exp(-b*x)
    s.param.init <- c(200, 5)
    fitPoisonExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                        data = poison_efrt,
                        algorithm = "port",
                        start = c(a = s.param.init[1], b = s.param.init[2]),
                        trace = TRUE,      
                        nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
    fitPoisonExp.summ <- summary(fitPoisonExp)
    propRemaining.cont <- exp(lpropRemaining.cont)
    hrsPigPoisonExp.pred <- coef(fitPoisonExp)[1] * exp(-coef(fitPoisonExp)[2] * propRemaining.cont)
    lines(propRemaining.cont,hrsPigPoisonExp.pred,lty=3,lwd=2,col="blue")

# review EXP model fit parameters for each control type
fitShotExp.summ
fitTAACExp.summ    
fitTrapExp.summ
fitPoisonExp.summ 


# Appendix S1: Figure S5   
  # make data frame to compare efficiencies of each method and plot FRs without plotting the individual events (as per figure above)
control.effs <- data.frame(exp(lpropRemaining.cont),hrsPigShot.pred,hrsPigTAAC.pred,hrsPigTrap.pred,hrsPigPoison.pred)
  colnames(control.effs) <-c("prop.rem","shot", "TAAC","trap", "poison")
  
  # Plot control.effs
  plot(control.effs$prop.rem, control.effs$shot, type = "l", col = "red", xlab = "Proportion of pigs remaining", ylab = "Hours per pig", ylim = range(control.effs), lty = 1, lwd = 2)
  
  # Add lines for the other control types
  lines(control.effs$prop.rem, control.effs$TAAC, col = "blue", lty = 2, lwd = 2)
  lines(control.effs$prop.rem, control.effs$trap, col = "green", lty = 3,lwd = 2)
  lines(control.effs$prop.rem, control.effs$poison, col = "black", lty = 4,lwd = 2)
  legend("topright", legend = c("shot", "TAAC", "trap", "poison"), col = c("red", "blue", "green", "black"), lty = c(1,2,3,4),lwd = 2)

 #export for plotting in other programs 
  write.csv(control.effs, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/CombinedExpFRs.csv")    
# linear and intercept only model fits for comparison
# linear fit
  
  fitShotLin <- lm(hrsPig ~ propRemaining, data=shot_efrt) # fits the Linear model
  summary(fitShotLin) # provides summary statistics of the LM for AIC comparison, not plotted
  
  fitTAACLin <- lm(hrsPig ~ propRemaining, data=TAAC_efrt) # fits the Linear model
  summary(fitTAACLin) # provides summary statistics of the LM for AIC comparison, not plotted
  
  fitTrapLin <- lm(hrsPig ~ propRemaining, data=trap_efrt) # fits the Linear model
  summary(fitTrapLin) # provides summary statistics of the LM for AIC comparison, not plotted
  
  fitPoisonLin <- lm(hrsPig ~ propRemaining, data=poison_efrt) # fits the Linear model
  summary(fitPoisonLin) # provides summary statistics of the LM for AIC comparison, not plotted

# intercept only
  fitShotInt <- lm(hrsPig ~ 1, data=shot_efrt) 
  summary(fitShotInt)
  
  fitTAACInt <- lm(hrsPig ~ 1, data=TAAC_efrt) 
  summary(fitTAACInt)
  
  fitTrapInt <- lm(hrsPig ~ 1, data=trap_efrt) 
  summary(fitTrapInt)
  
  fitPoisonInt <- lm(hrsPig ~ 1, data=poison_efrt) 
  summary(fitPoisonInt)

# goodness of fit comparison for different models
# wAIC values used for comparison and for calculating evidence ratios
shot.AIC.vec <- c(AIC(fitShot), AIC(fitShotExp), AIC(fitShotLin))
shot.d.vec <- delta.AIC(shot.AIC.vec)
shot.w.vec <- weight.AIC(shot.d.vec)
shot.w.vec  
shot.d.vec

TAAC.AIC.vec <- c(AIC(fitTAAC), AIC(fitTAACExp), AIC(fitTAACLin))
TAAC.d.vec <- delta.AIC(TAAC.AIC.vec)
TAAC.w.vec <- weight.AIC(TAAC.d.vec)
TAAC.w.vec
TAAC.d.vec

trap.AIC.vec <- c(AIC(fitTrap), AIC(fitTrapExp), AIC(fitTrapLin))
trap.d.vec <- delta.AIC(trap.AIC.vec)
trap.w.vec <- weight.AIC(trap.d.vec)
trap.w.vec
trap.d.vec

poison.AIC.vec <- c(AIC(fitPoison), AIC(fitPoisonExp), AIC(fitPoisonLin))
poison.d.vec <- delta.AIC(poison.AIC.vec)
poison.w.vec <- weight.AIC(poison.d.vec)
poison.w.vec
poison.d.vec

# return to plot arrangement tp single plots rather than 2x2 plot grid
par(mfrow=c(1,1))

#### Stochastic Population Projection ####
#### create Leslie matrix ####

age.max = 14 # Snow et al. (2019) report longevity up to 12-14, so I have set longevity at 14 years

## create fertility and survival vectors 
# derived from Choquenot (1996)

#### Fertility ####
# "Average litter sizes vary from 4.9 to 6.3 piglets". Mean litter size = (4.9 + 6.3)/2 = 5.6
# " on average, females have 0.85 litters per year"
# mean fertility per female = mean litter size per female x litters per female per year = 5.6 x 0.85 = 4.76
# assuming 1:1 sex ratio...
# ...mean number females born each year = mean birth rate/2 = 4.76/2 = 2.38 females offspring born per female per year
# juvenile fertility approximated as Adult fertility/3 allowing for some juvenile females to reproduce in their first year = 2.38/3 = 0.79

f.vec <- c(0.79, 2.38, 2.38, 2.38, 2.38, 2.38, 2.38, 2.38, 2.38, 2.38, 2.38, 2.38, 2.38, 2.38)

## SD calculated as (high range value - low range value)/(2*1.96) to approximate SD under the assumptions of the Student's t distribution

# adult fertility high range value (females only) = upper limit litter size x litters per year / 2 = 6.3 x 0.85 / 2 = 3.7059
# adult fertility low range value (females only)= lower limit litter size x litters per year / 2 = 4.9 x 0.85 / 2 = 2.0825
# Adult SD = (3.7059-2.0825)/(2*1.96) = 1.6234/3.92 = 0.4141
A.f.sd <- 0.4141

# juvenile fertility  high range = (upper limit litter size/3) x litters per year /2 = 2.1 * 0.85 / 2 = 0.8925
# juvenile fertility  low range = (upper limit litter size/3) x litters per year /2 = 1.63 * 0.85 / 2 = 0.6942
# juvenile SD = (0.8925-0.6942)/(2*1.96) = 0.1983/3.92 = 0.0506
J.f.sd <- 0.0506 

# vector of fertility SDs
f.sd.vec <- c(J.f.sd, A.f.sd, A.f.sd, A.f.sd, A.f.sd, A.f.sd, A.f.sd, A.f.sd, A.f.sd, A.f.sd, A.f.sd, A.f.sd, A.f.sd, A.f.sd) #mean and standard deviations vector, juvenile and adult fertility 

## visualise change in mean fertility as a function of age
plot(0:13,f.vec, pch=19, type="b")

#### Survival ####
# "piglet mortality ranging from 10 –15% (0.9 survival) when food supplies and weather are favourable, up to 100% (0.0 survival) when conditions are poor"
# mean juvenile survival = upper range + lower range / 2 = 0.9 + 0.0 / 2 = 0.9 / 2 = 0.45
# "Adult mortality can vary from 15 (0.85 survival) to 50% (0.5 survival) with few pigs living beyond five years of age"
# mean Adult survival = upper range + lower range / 2 = 0.85 + 0.5 / 2 = 1.35 / 2 = 0.675
s.vec <- c(0.45, 0.675, 0.675, 0.675, 0.675, 0.675, 0.675, 0.675, 0.675, 0.675, 0.675, 0.675, 0.675, 0.675) ##feral pig survival values from Choquenot (1996)

# survival SD
# Adult SD = (upper - lower range) / (2*1.96) (0.85 - 0.5)/(2*1.96) = 0.35/3.92 = 0.0893
A.s.sd <- 0.0893

# Adult SD = (upper - lower range) / (2*1.96) (0.9- 0.0)/(2*1.96) = 0.9/3.92 = 0.2296
J.s.sd <- 0.2296 #mean and standard deviations, juvenile survival
 
# vector of survival SDs
s.sd.vec <- c(J.s.sd, A.s.sd, A.s.sd, A.s.sd, A.s.sd, A.s.sd, A.s.sd, A.s.sd, A.s.sd, A.s.sd, A.s.sd, A.s.sd, A.s.sd, A.s.sd) #mean and standard deviations vector, juvenile and adult survival

## visualise change in mean survival as a function of age
plot(0:13,s.vec,pch=19,type="b")

# create matrix
popmat <- matrix(data = 0, nrow=age.max, ncol=age.max) ##creates a matrix where data = 0 and dimensions age.max x age.max; 14 x 14
diag(popmat[2:age.max,]) <- s.vec[1:13] ## diagonally in popmat from row 2, col 1 to row 6, col 5 populated with s.vec
popmat[age.max,age.max] <- s.vec[14] ## position [14,14] = 0
popmat[1,] <- f.vec ## row 1 of popmat populated with f.vec
popmat.orig <- popmat ## save original matrix

## matrix properties. Functions from "matrixOperators.R"
max.lambda(popmat) ## 1-yr lambda. 
max.r(popmat) # rate of population change, 1-yr
ssd <- stable.stage.dist(popmat) ## stable stage distribution
plot(0:13,ssd,pch=19,type="b")
R.val(popmat, age.max) # reproductive value
gen.l <- G.val(popmat, age.max) # mean generation length

## initial population vector
pop.found <- 500/2 # +/- 50 founding population size personal communication, B. Page, PIRSA,  2021 . 
## Divided by 2 for females only
ssd <- stable.stage.dist(popmat) ## sum of all values in ssd = 1
init.vec <- ssd * pop.found #initial population vector. Sum of all values here = 500 as suggested by PIRSA
## (can change to whatever we want our founding pop to be)


#### Deterministic population projection ####
# i.e., project population growth from the leslie matrix without incorporating stochasticity
## set time limit for projection in 1-yr increments
yr.now <- 2020 
#************************

yr.end <- 2030 #end year for our projection time frame
#************************
t <- (yr.end - yr.now) #timeframe
## time frame can be adjusted to explore different scenarios, but set to 3 here in line with PIRSA target time frame for pig eradication

tot.F <- sum(popmat.orig[1,]) ## this is the sum of fertility values in row 1 of the deterministic matrix. 
## Initially this is the same as f.vec, but popmat will change for each iteration as values are stochastically resampled 
popmat <- popmat.orig #resets matrix 
yr.vec <- seq(yr.now,yr.end) #year vector, 2020, 2021, 2022.... Increases in increments of one by default so increment size doesn't need to be specified

## set population storage matrices
n.mat <- matrix(0, nrow=age.max,ncol=(t+1)) #matrix with values = 0, number of columns = to max.age (=6) and number of columns = 
n.mat[,1] <- init.vec #fill first matrix column with initial population vector

## set up projection loop
for (i in 1:t) { #  for every year of the specified timeframe...
  n.mat[,i+1] <- popmat %*% n.mat[,i] #  the corresponding column of n.mat...
  #is populated with the product of current year column multiplied by the deterministic matrix popmat
}

n.pigs <- colSums(n.mat) #number of pigs in any year is the sum of all values in the corresponding column 
## no density reduction treatment or carrying capacity at this stage
yrs <- seq(yr.now, yr.end, 1)
plot(yrs, (n.pigs),type="l",lty=2,pch=19,xlab="year",ylab="N")

#### compensatory density feedback: SURVIVAL ####
# K = carry capacity
K.max <- 2500
## Southgate, R. (2018). Feral Pig Sus Scrofa and Feral Cat Felis catus monitoring on Kangaroo Island: 2014 to 2017 - found pig density to be stable during study period
## Masters (2011) pig numbers estimated from densities observed elsewhere in Australia (0.4 - 4 km2 )
## pre-fire pig population estimated to be ~5000 (PIRSA personal communication)
## Assume pop remained stable between 2017 study and 2019 bushfire, so set k.max as 5000/2 (0.5 to account for females only) 

K.s.vec <- c(1,K.max*0.2, K.max*0.45, K.max*0.65, K.max*0.75, K.max*.85, K.max*0.95) ##describes the x axis of the reduction curve
red.s.vec <- c(1,0.995,0.98,0.93,0.8,0.58,0.28) ## describes the y axis of the reduction curve

plot(K.s.vec,red.s.vec,pch=19)
Kred.s.dat <- data.frame(K.s.vec,red.s.vec)

# logistic power function a/(1+(x/b)^c) #fits logistic power function to population relative to carry capacity, K
s.param.init <- c(1, K.max, 3) ## These parameters are arbitrary and arrived at after testing various combinations until expected feedback response achieved.
## parameter adjusts the initial value of the modifier and should be close to 1.
## parameter b changes how quickly the modifier approaches K
## parameter c adjusts the shape of the relationship between the modifier and K.
fit.s.lp <- nls(red.s.vec ~ a/(1+(K.s.vec/b)^c), 
                data = Kred.s.dat,
                algorithm = "port",
                start = c(a = s.param.init[1], b = s.param.init[2], c = s.param.init[3]),
                trace = TRUE,      
                nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fit.s.lp.summ <- summary(fit.s.lp)
plot(K.s.vec,red.s.vec,pch=19,xlab="N",ylab="s reduction factor")
K.s.vec.cont <- seq(1,12*pop.found,1)
pig.s.lp.fx <- coef(fit.s.lp)[1]/(1+(K.s.vec.cont/coef(fit.s.lp)[2])^coef(fit.s.lp)[3])
lines(K.s.vec.cont,pig.s.lp.fx,lty=2,lwd=3,col="red")

s.a.lp <- coef(fit.s.lp)[1]
s.b.lp <- coef(fit.s.lp)[2]
s.c.lp <- coef(fit.s.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
n.mat[,1] <- init.vec
popmat <- popmat.orig

## set up projection: deterministic pop with S density feedback
for (i in 1:t) {
  
  totN.i <- sum(n.mat[,i])
  pigs.s.red <- s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp)
  diag(popmat[2:age.max,]) <- s.vec[1:13]*pigs.s.red
  popmat[age.max,age.max] <- s.vec[14]*pigs.s.red
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pigs <- colSums(n.mat)
plot(yrs, n.pigs,type="l",main = "deterministic population projection", sub = " S density feedback", lty=1,pch=19,xlab="year",ylab="N",ylim=c(0,1.8*K.max)) #untreated population increases, rate of increase relative to K, no stochastic sampling
abline(h=K.max, lty=2, col="red") #carry capacity

#### iterations ####
iter <- 1000 #final model run at 10 000
itdiv <- iter/10 #final model rate at iter/1000


#### untreated population projection ####

## stochastic projection with survival density feedback
## set storage matrices & vectors
# set t to 10 years so projection shows response relative to K

yr.now <- 2020
#************************

yr.end <- 2120 #end year for our projection timeframe
yrs <- seq(yr.now, yr.end, 1)
#************************
t <- (yr.end - yr.now)

n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1)) #storage matrix with 1000 rows and t+1 columns based on t year timeframe

for (e in 1:iter) {
  popmat <- popmat.orig
  
  n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
  n.mat[,1] <- init.vec
  
  for (i in 1:t) {
    # stochastic survival values
    s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
    s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
    s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
    
    # stochastic fertility sampler (gaussian)
    
    fert.stoch <- rtruncnorm(length(popmat.orig[,1]), a=0, b=Inf, mean=popmat.orig[1,], sd=f.sd.vec)
    
    totN.i <- sum(n.mat[,i])
    pigs.s.red <- as.numeric(s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp))
    #a/(1+(K.vec/b)^c)
    
    popmat[1,] <- fert.stoch
    diag(popmat[2:age.max,]) <- s.stoch[1:13]*pigs.s.red
    popmat[age.max,age.max] <- s.stoch[14]*pigs.s.red
    #popmat[age.max,age.max] <- 0
    
    
    n.mat[,i+1] <- popmat %*% n.mat[,i]
    
  } # end i loop
  
  n.sums.mat[e,] <- ((as.vector(colSums(n.mat))))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # median pop size over all iterations
n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # 97.5th percentile over all iterations
n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # 2.5th percentile over all iterations

# produce dataframe recording median and CIs for pop size in each year
datN <- data.frame(yrs, n.md, n.lo, n.up)
tail(datN)

## plot pop change in ggplot2
ggplot(data = datN, mapping = aes(x=yrs)) +
  geom_line(aes(y=n.md), color = "black") + 
  ylim(0,5000) +
  scale_x_continuous(breaks = seq(2020, 2120, by = 10)) +
  geom_line(aes(y=n.lo),color = "red", linetype = 2) + 
  geom_line(aes(y=n.up),color = "red", linetype = 2) + 
  theme_bw() +
  geom_hline(yintercept = 2500, linetype = 3) +
  labs(x = "Years", y = "Population (females)")

write.csv(datN, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/Choq_stochastic.pop.10yrS.csv")

#### CONSTANT PROPORTIONAL YEARLY HARVEST - NO COSTS####
# This section calculates minimum proportion remaining after t years of constant proportional cull.

#  3 years
yr.now <- 2020 
#************************

yr.end <- 2023 #end year for our projection timeframe
yrs <- seq(yr.now, yr.end, 1)
#************************
t <- (yr.end - yr.now)

# # harvest rate

harv.prop.consist <- seq(0.2,0.99,0.05) #sequence harvest/culling quotas, e.g remove 0.2-.99 proportion of founding pop in increasing increments of 0.05

# define our quasi-extinction probability storage vector
min.med.n3 <- min.lo.n3 <- min.up.n3 <- rep(0,length(harv.prop.consist))

for (s in 1:length(harv.prop.consist)) {

  # set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1))
  #storage matrix with 1000 rows and 4 columns (years) based on timeframe

  for (e in 1:iter) {
    popmat <- popmat.orig

    n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
    n.mat[,1] <- init.vec

    for (i in 1:t) {
      # stochastic survival values
      s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
      s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
      s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)

      # stochastic fertility sampler (gaussian)
      fert.stoch <- rtruncnorm(length(popmat.orig[,1]), a=0, b=Inf, mean=popmat.orig[1,], sd=f.sd.vec)

      totN.i <- sum(n.mat[,i])
      pigs.s.red <- as.numeric(s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp))
      #a/(1+(K.vec/b)^c)

      popmat[1,] <- fert.stoch
      diag(popmat[2:age.max,]) <- s.stoch[1:13]*pigs.s.red
      popmat[age.max,age.max] <- s.stoch[14]*pigs.s.red
      #popmat[age.max,age.max] <- 0

      n.mat[,i+1] <- popmat %*% n.mat[,i]

      ## harvest things here
      n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.consist[s], 0), 0)


      if (length(which(n.mat[,i+1] < 0)) > 0) {
        n.mat[which(n.mat[,i+1] < 0), i+1] <- 0

      }

    }# end i loop

    n.sums.mat[e,] <- as.vector((colSums(n.mat))/pop.found) # / pop.mat for min proportion remaining population

    if (e %% itdiv==0) print(e)

  } # end e loop

  # calculate minimum population size

  min.pop.vec <- apply(n.sums.mat, MARGIN=1, min)
  min.med.n3[s] <- median(min.pop.vec, na.rm=T)
  min.lo.n3[s] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
  min.up.n3[s] <- quantile(min.pop.vec, probs=0.975, na.rm=T)

  #
  # # print("##############")
  print(paste("harvest proportion = ", harv.prop.consist[s], sep=""))
  # print("##############")

} # ends S loop

minn.prop.pop3 <- data.frame(harv.prop.consist, min.med.n3, min.lo.n3, min.up.n3)
minn.prop.pop3$med.pop <- minn.prop.pop3$min.med.n3*250

## plot in ggplot (thesis figure 3)
harv3plot <- ggplot(data = minn.prop.pop3, mapping = aes(x=harv.prop.consist)) +
  geom_line(aes(y=min.med.n3), color = "black") + 
  geom_line(aes(y=min.lo.n3), color = "red", linetype = 2) + 
  geom_line(aes(y=min.up.n3), color = "red", linetype = 2) + 
  geom_vline(xintercept = 0.8, linetype = 3, color = "black") +
  theme_bw() +
  labs(x = "Harvest Proportion", y = "Min N")

write.csv(minn.prop.pop3, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/Choq_constant.prop.cull3.csv")

### repeat for 10 year projection interval
#  10 years
yr.now <- 2020 
#************************

yr.end <- 2030 #end year for our projection timeframe
yrs <- seq(yr.now, yr.end, 1)
#************************
t <- (yr.end - yr.now)

# # harvest rate

harv.prop.consist <- seq(0.2,0.99,0.05) #sequence harvest/culling quotas, e.g remove 0.2-.99 proportion of founding pop in increasing increments of 0.05

# define our quasi-extinction probability storage vector
min.med.n <- min.lo.n <- min.up.n <- rep(0,length(harv.prop.consist))

for (s in 1:length(harv.prop.consist)) {
  
  # set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1))
  #storage matrix with 1000 rows and 4 columns (years) based on timeframe
  
  for (e in 1:iter) {
    popmat <- popmat.orig
    
    n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
    n.mat[,1] <- init.vec
    
    for (i in 1:t) {
      # stochastic survival values
      s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
      s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
      s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
      
      # stochastic fertility sampler (gaussian)
      fert.stoch <- rtruncnorm(length(popmat.orig[,1]), a=0, b=Inf, mean=popmat.orig[1,], sd=f.sd.vec)
      
      totN.i <- sum(n.mat[,i])
      pigs.s.red <- as.numeric(s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp))
      #a/(1+(K.vec/b)^c)
      
      popmat[1,] <- fert.stoch
      diag(popmat[2:age.max,]) <- s.stoch[1:13]*pigs.s.red
      popmat[age.max,age.max] <- s.stoch[14]*pigs.s.red
      #popmat[age.max,age.max] <- 0
      
      n.mat[,i+1] <- popmat %*% n.mat[,i]
      
      ## harvest things here
      n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.consist[s], 0), 0)
      
      
      if (length(which(n.mat[,i+1] < 0)) > 0) {
        n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
        
      }
      
    }# end i loop
    
    n.sums.mat[e,] <- as.vector((colSums(n.mat))/pop.found) # / pop.mat for min proportion remaining population
    
    if (e %% itdiv==0) print(e)
    
  } # end e loop
  
  # calculate minimum population size
  
  min.pop.vec <- apply(n.sums.mat, MARGIN=1, min)
  min.med.n[s] <- median(min.pop.vec, na.rm=T)
  min.lo.n[s] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
  min.up.n[s] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
  
  # # print("##############")
  print(paste("harvest proportion = ", harv.prop.consist[s], sep=""))
  # print("##############")
  
} # ends S loop

plot(harv.prop.consist, min.med.n, type="l", pch=19, xlab="harvest proportion", ylab="min N", ylim=c(min(min.lo.n),max(min.up.n)))
lines(harv.prop.consist, min.lo.n, col="red", lty=2)
lines(harv.prop.consist, min.up.n, col="red", lty=2)

minn.prop.pop <- data.frame(harv.prop.consist, min.med.n, min.lo.n, min.up.n)
minn.prop.pop$med.pop <- minn.prop.pop$min.med.n*250

#plot the lines individually
harv10plot <- ggplot(data = minn.prop.pop, mapping = aes(x=harv.prop.consist)) +
  geom_line(aes(y=min.med.n), color = "black") + 
  geom_line(aes(y=min.lo.n), color = "red", linetype = 2) + 
  geom_line(aes(y=min.up.n), color = "red", linetype = 2) + 
  geom_vline(xintercept = 0.6, linetype = 3, color = "black") +
  theme_bw() +
  labs(x = "Harvest Proportion", y = "Min N")

write.csv(minn.prop.pop, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/Choq_constant.prop.cull10.csv")

# multipanel of constant prop harvest over 3 and 10 years (Thesis figure 3 - results)
harvplot <- ggarrange(harv3plot, harv10plot, common.legend = TRUE, labels = c("a", "b"), legend = "bottom", nrow = 1, ncol = 2)
harvplot

#### CONSTANT PROPORTIONAL CULL COSTS####
# fixed cost components (not recalculated each loop)
labour.ph <- 36.87 ## $36.87 per hour, based on SA Public Sector Award OPS4 classification. 
Vehicle.ph<- 10 ## arbitrary per hour vehicle cost to reflect vehicle use for shooting, trapping and baiting
labour.vehicle <- labour.ph+Vehicle.ph # combined labour and vehicle cost for shooting, trapping and baiting
bullets.pp <- 4   ## 2 bullets per pig @ $2 per bullet
freefeed.ph <- 14 ## based on daily deployment of 10kg grain at $1.4 per kg. Deployment of grain takes 1 hr per day 

# TAAC specific costs
TAAC.initial <- 12246 # total initialisation cost $12,246.00
TAAC.crew.pd <- 2482.66 # daily cost for flight time and crew wages. includes daily fuel consumption. Contractor charges daily, not by the hour
TACC.crew.accomm <- 420.00 # daily crew food and accommodation
TAAC.mrksmn.accomm <- 199.75 # marksman meals and accommodation (marksman supplied by SA Gov - costs based on daily meals and accommodation allowances under Commissioners Determination - 3.2)
TAAC.daily.ft <- 3.7 #average total daily flight time (hrs/day) = 3.7 therefore TAAC.total.pd accrued for every 3.7hrs effort  
TAAC.resupply <- 1246.00 # cost for additional fuel every 30 days = 30 x 3.7hrs =  111hrs   
TAAC.pd <- TAAC.crew.pd + TACC.crew.accomm + (labour.ph*7.5) + TAAC.mrksmn.accomm # Total daily TAAC costs

#shot specific costs
freefeed.ph.shot <- freefeed.ph * 0.588 # assumes 0.588 of shooting done with free feeding

# trap specific costs
trap.initial <- labour.vehicle*20 # additional unreported effort of ~20 hrs per trap for trap construction
trap.unit <- 9500.00 # cost per unit for the Jager Pro trap
trap.ff<-freefeed.ph*mean(trap_efrt$efforthrs) # cost of grain per event
#trap.ff.labour <-labour.ph*mean(trap_efrt$efforthrs) # feeding labour cost per event
# mean(trap_efrt$efforthrs) = 9.785 
trap.ff.labour <-labour.vehicle*10 
#rounded up to 10 as 9.875 is very precise compared to remaining estimated 20hrs setup per trap. 
trap.ff.cost<-trap.ff.labour + trap.ff # combined cost of grain and labour per event

# poison specific costs
poison.placebo <- 224 # 6 placebo baits x 4 days
poison.toxic <- 156 # 6 toxic baits x 1 day
poison.dispenser.unit <- 485
avg.poison.effort <- mean(poison_efrt$efforthrs)
poison.labour <- labour.vehicle*avg.poison.effort
poison.ff <- freefeed.ph*(avg.poison.effort-5)
max.annual.poison.events <-365/avg.poison.effort

# set t back to 3. Proportional cull loops repeated below for t=10  
##Will set up as loop if time allows
t<-3
# set storage matrices & vectors

n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1)) 
cost.mat <- matrix(data=0, nrow=iter, ncol=t)
time.mat <- matrix(data=0, nrow=iter, ncol=t)
trap.number <- rep(NA, t)
trap.purchase.number <- rep(NA, t)
dispenser.no <- rep(NA, t)
#relative prop cost and total cost vectors
tot.TAAC.rcost <- tot.trap.rcost <- tot.poison.rcost <- tot.shot.rcost <-  
  tot.TAAC.cost <-tot.trap.cost <- tot.poison.cost <- tot.shot.cost <- tot.cost <-  
  tot.TAAC.hrs <-tot.trap.hrs <- tot.poison.hrs <- tot.shot.hrs <- tot.time <- rep(NA, t)
prop.rcost.vec <-  inv.prop.rcost.vec <- weight.rcost.vec <- rep(NA, 4)

# create matrix for calling different cull scenarios. Each row provides different set of proportions
cull.props.mat <- matrix(0, 6, 4)
cull.props.mat[1,] <- weight.rcost.vec # relative cost proportional cull
cull.props.mat[2,] <- rep(0.25,4) # even prop cull
diag(cull.props.mat[3:6,]) <- c(1,1,1,1)
# cull.props.mat[3,] <- TAAC only
# cull.props.mat[4,] <- Shot only
# cull.props.mat[5,] <- Trap only
# cull.props.mat[6,] <- Poison only 

harv.prop.consist <- seq(0.35,0.95,0.05) 
#sequence harvest/culling quotas, e.g remove 0.35-.95 proportion of founding pop in increasing increments of 0.05

# cost, time and minimum n storage vector
tot.med <- tot.lo <- tot.up <- min.med.n <- min.lo.n <- min.up.n <- tot.time.med <- tot.time.lo <- tot.time.up <- rep(0,length(harv.prop.consist))

# storage for outcomes of each scenario at each harvest proportion
scenarios.mat <- matrix(NA, t, 10)
harv.prop.cost.mat <- matrix(NA, t, 10)
scenarios.tot <- rep(NA,11)

# names for each sculling scenario to be used later when plotting outcomes
scenario.names <- c("Relative cost proportional allocation", 
                    " 25pc harvest per method", "Thermal-assisted aerial cull", "Shoot", "Trap", "Poison")
par(mfrow = c(1,1))


for (f in 1:nrow(cull.props.mat)){
  ## f loop runs projection and cost calculations for each cull scenario

for (s in 1:length(harv.prop.consist)) {
  
    for (e in 1:iter) {
      popmat <- popmat.orig
      
      harvested.tot <- tot.cost <- tot.time <- rep(NA,t)
      n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
      n.mat[,1] <- init.vec
      
      for (i in 1:t) {
        # stochastic survival values
        s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
        s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
        s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
        
        # stochastic fertility sampler (Gaussian)
        fert.stoch <- rtruncnorm(length(popmat.orig[,1]), a=0, b=Inf, mean=popmat.orig[1,], sd=f.sd.vec)
        
        totN.i <- sum(n.mat[,i])
        pigs.s.red <- as.numeric(s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp))
        #a/(1+(K.vec/b)^c)
        
        popmat[1,] <- fert.stoch
        diag(popmat[2:age.max,]) <- s.stoch[1:13]*pigs.s.red
        popmat[age.max,age.max] <- s.stoch[14]*pigs.s.red
        #popmat[age.max,age.max] <- 0
        
        n.mat[,i+1] <- popmat %*% n.mat[,i]
        
        ## harvest things here
        # create harvested component matrix
        harvested.tot[i] <- round(sum(n.mat[,i+1])*harv.prop.consist[s], 0)
        
        # do the harvest
        n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * harvested.tot[i], 0)
        
        # prevent n.mat from containing negative values. Cannot have less than 0 pigs
        if (length(which(n.mat[,i+1] < 0)) > 0) {
          n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
        }
        
        # calculate proportion remaining after each iteration. 
        ## proportions are relative to number used to calculate functional response curves (males and females), 
            ## not initial pop in projection(n=250: females only)
        remaining.preharv <- sum(n.mat[,i])
        prop.remain <- remaining.preharv/957*2
        
        # effort per control type relative to proportion remaining
        hrsPigTAAC <- ifelse(coef(fitTAACExp)[1] * exp(-coef(fitTAACExp)[2] * prop.remain)>=min(TAAC_efrt$hrsPig),
                             coef(fitTAACExp)[1] * exp(-coef(fitTAACExp)[2] * prop.remain),
                             min(TAAC_efrt$hrsPig))
        hrsPigShot <- ifelse(coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * prop.remain)>=min(shot_efrt$hrsPig), 
                             coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * prop.remain),
                             min(shot_efrt$hrsPig))
        hrsPigTrap <- ifelse(coef(fitTrapExp)[1] * exp(-coef(fitTrapExp)[2] * prop.remain)>=min(trap_efrt$hrsPig),
                             coef(fitTrapExp)[1] * exp(-coef(fitTrapExp)[2] * prop.remain),
                             min(trap_efrt$hrsPig))
        hrsPoison <- ifelse(coef(fitPoisonExp)[1] * exp(-coef(fitPoisonExp)[2] * prop.remain)>=min(poison_efrt$hrsPig), 
                            coef(fitPoisonExp)[1] * exp(-coef(fitPoisonExp)[2] * prop.remain),
                            min(poison_efrt$hrsPig))
        
        
        #### calculating and assigning proportions #####
        
        
        if (f == 1){ ##  calculate relative proportional cost allocation ####
                # TAAC
                tot.TAAC.hrs <- hrsPigTAAC*harvested.tot[i]
                ## cost accrued in increments of 111 hrs TAAC work completed.
                tot.TAAC.rcost[i] <- TAAC.initial + ((TAAC.pd*tot.TAAC.hrs)/TAAC.daily.ft) + bullets.pp*harvested.tot[i] + ((TAAC.resupply*tot.TAAC.hrs)/111)
                
                # shot
                tot.shot.hrs <- hrsPigShot*harvested.tot[i]
                labour.cost <- (labour.vehicle  + freefeed.ph.shot) * tot.shot.hrs
                bullet.cost <- bullets.pp * harvested.tot[i]
                tot.shot.rcost[i] <- labour.cost + bullet.cost
                
                # Trap
                tot.trap.hrs <- hrsPigTrap*harvested.tot[i]
                trap.number[i] <-ceiling((tot.trap.hrs/(10*(20+mean(trap_efrt$efforthrs)))))
                trap.purchase.number[i] <- trap.number[i]# ceiling used so number of traps is an integer
                # need to deduct number of traps purchased in previous years
                if (i > 1) {
                  trap.purchase.number[i] <- ifelse(trap.purchase.number[i] < trap.purchase.number[1], 0, trap.purchase.number[i] - trap.purchase.number[1])
                }
                tot.trap.rcost[i] <-(trap.purchase.number[i]*trap.unit) + ## cost to purchase no. traps required in year [i]
                  (trap.number[i] * 10 * (trap.initial + trap.ff.cost)) + # cost to set up and feed no. traps required in year [i]. Assuming 10 events per trap per year
                  bullets.pp*harvested.tot[i]
                
                #poison
                tot.poison.hrs <- hrsPoison*harvested.tot[i]
                total.poison.events.pa <- tot.poison.hrs/avg.poison.effort # no of events
                dispenser.no[i] <- ceiling((total.poison.events.pa/max.annual.poison.events) * 1.45)
                if (i > 1) {
                  dispenser.no[i] <- ifelse(dispenser.no[i] < dispenser.no[1], 0, dispenser.no[i] - dispenser.no[1])
                }
                tot.poison.rcost[i] <- (total.poison.events.pa * (poison.labour + poison.ff)) + 
                (total.poison.events.pa * (1.45 * (poison.placebo + poison.toxic))) + 
                  (dispenser.no[i]*poison.dispenser.unit)
                
                tot.rcost <- tot.TAAC.rcost[i] + tot.shot.rcost[i] + tot.trap.rcost[i] + tot.poison.rcost[i]
                tot.rcost.vec <- c(tot.TAAC.rcost[i],tot.shot.rcost[i],tot.trap.rcost[i],tot.poison.rcost[i])
                cull.props.mat[1,] <- (1/(tot.rcost.vec/tot.rcost))/sum(1/(tot.rcost.vec/tot.rcost))
                   # populates cull.props.mat with relative cost weighting
                }
        
        harvest.TAAC <- harvested.tot[i]*cull.props.mat[f,1]
        harvest.shot <- harvested.tot[i]*cull.props.mat[f,2]
        harvest.trap <- harvested.tot[i]*cull.props.mat[f,3]
        harvest.poison <- harvested.tot[i]*cull.props.mat[f,4]
        
        #### actual cost calculation #### 
        # relative cost prop weighting
        # TAAC
        tot.TAAC.hrs[i] <- hrsPigTAAC*harvest.TAAC
        tot.TAAC.cost[i] <- TAAC.initial + ((TAAC.pd*tot.TAAC.hrs[i])/TAAC.daily.ft) + bullets.pp*harvest.TAAC + ((TAAC.resupply*tot.TAAC.hrs[i])/111)
        
        # shot
        tot.shot.hrs[i] <- hrsPigShot*harvest.shot
        labour.cost <- (labour.vehicle  + freefeed.ph.shot) * tot.shot.hrs[i]
        bullet.cost <- bullets.pp * harvest.shot
        tot.shot.cost[i] <- labour.cost + bullet.cost
        
        # Trap
        tot.trap.hrs[i] <- hrsPigTrap*harvest.trap
        trap.number[i] <- ceiling((tot.trap.hrs[i]/(10*(20+mean(trap_efrt$efforthrs))))) # ceiling used so number of traps is an integer
        trap.purchase.number[i] <- trap.number[i]
        # need to deduct number of traps purchased in previous years
        if (i > 1) {
          trap.purchase.number[i] <- ifelse(trap.purchase.number[i] < trap.purchase.number[1], 0, trap.purchase.number[i] - trap.purchase.number[1])
        }
        tot.trap.cost[i] <-(trap.purchase.number[i]*trap.unit) + ## cost to purchase no. traps required in year [i]
          (trap.number[i] * 10 * (trap.initial + trap.ff.cost)) + # cost to set up and feed no. traps required in year [i]. Assuming 10 events per trap per year
          bullets.pp*harvest.trap
        
        #poison
        tot.poison.hrs[i] <- hrsPoison*harvest.poison
        total.poison.events.pa <- tot.poison.hrs[i]/avg.poison.effort # no of events
        dispenser.no[i] <- ceiling((total.poison.events.pa/max.annual.poison.events) * 1.45)
        if (i > 1) {
          dispenser.no[i] <- ifelse(dispenser.no[i] < dispenser.no[1], 0, dispenser.no[i] - dispenser.no[1])
        }
        tot.poison.cost[i] <- (total.poison.events.pa * (poison.labour + poison.ff)) + 
          total.poison.events.pa * (1.45 * (poison.placebo + poison.toxic)) + 
          (dispenser.no[i]*poison.dispenser.unit)
        
        tot.cost[i] <- sum(c(tot.TAAC.cost[i],tot.shot.cost[i],tot.trap.cost[i],tot.poison.cost[i]))
        tot.time[i] <- sum(c(tot.TAAC.hrs[i], tot.shot.hrs[i], tot.trap.hrs[i], tot.poison.hrs[i]))
        
      } # end i loop
      
      cost.mat[e,] <- tot.cost
      time.mat[e,] <- tot.time
      
      n.sums.mat[e,] <- as.vector((colSums(n.mat))/pop.found) # / pop.mat for min proportion remaining population 
      
      if (e %% itdiv==0) print(e) 
      
    } # end e loop
    
    # Calculate total median cost + upper and lower 95% confidence intervals
    
    median.cost <- apply(cost.mat, MARGIN=2, median, na.rm=T)
    up.cost <-  apply(cost.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
    lo.cost <-  apply(cost.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
    
    tot.med[s] <- sum(median.cost)
    tot.up[s] <- sum(up.cost)
    tot.lo[s] <- sum(lo.cost)
    print(c(tot.lo[s], tot.med[s], tot.up[s]))
    
    #calculate total median time + upper and lower 95% confidence intervals
    median.time <- apply(time.mat, MARGIN=2, median, na.rm=T)
    up.time <-  apply(time.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
    lo.time <-  apply(time.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
    
    tot.time.med[s] <- sum(median.time)
    tot.time.up[s] <- sum(up.time)
    tot.time.lo[s] <- sum(lo.time)
    print(c(tot.time.lo[s], tot.time.med[s], tot.time.up[s]))
   
    
    # plot costper year for each scenario and cull proportion combination
    plot(yrs[1:(t)], median.cost, type="l", lwd=2, xlab="year", ylab="cost ($)", 
         main = paste("Annual costs:", scenario.names[f]), 
         sub = paste("cull proportion", harv.prop.consist[s]), 
         ylim=(c(min(lo.cost),max(up.cost))))
    lines(yrs[1:(t)], up.cost, lty=2, col="red")  
    lines(yrs[1:(t)], lo.cost, lty=2, col="red")
    
    # calculate minimum population size for each 
    min.pop.vec <- apply(n.sums.mat, MARGIN=1, min)
    min.med.n[s] <- median(min.pop.vec, na.rm=T)
    min.lo.n[s] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
    min.up.n[s] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
    
    # calculate median prop remaining after each year
    med.N.Yr.vec <- apply(n.sums.mat, MARGIN = 2, median)
    
    print(paste("harvest proportion = ", harv.prop.consist[s], sep=""))
    
    # matrix of lower, median and upper cost per year for each scenario and proportion
    harv.prop.cost.mat[,1] <- scenario.names[f]
    harv.prop.cost.mat[,2] <- harv.prop.consist[s]
    harv.prop.cost.mat[,3] <- yrs[1:(t)]
    harv.prop.cost.mat[,4] <- lo.cost
    harv.prop.cost.mat[,5] <- median.cost
    harv.prop.cost.mat[,6] <- up.cost
    harv.prop.cost.mat[,7] <- lo.time
    harv.prop.cost.mat[,8] <- median.time
    harv.prop.cost.mat[,9] <- up.time
    harv.prop.cost.mat[1:3,10] <- med.N.Yr.vec[2:4]
    
    scenarios.mat <-rbind(scenarios.mat, harv.prop.cost.mat)
    
    # matrix of lower, median and upper total cost and minimum N for each scenario and proportion
    
    harv.prop.annual <- c(scenario.names[f], harv.prop.consist[s],
                          tot.lo[s], tot.med[s], tot.up[s], tot.time.lo[s], tot.time.med[s], tot.time.up[s], min.lo.n[s], min.med.n[s], min.up.n[s])

    scenarios.tot <-rbind(scenarios.tot, harv.prop.annual)
  
    } # end s loop
} # end f loop

scenarios.df <- as.data.frame(scenarios.mat[-(1:3),])
names(scenarios.df) <- c("scenario", "harvest.prop", "year","low.cost", "median.cost", "upper.cost", 
                         "low.time", "median.time", "upper.time", "median.N") 
rownames(scenarios.df)<-NULL
scenarios.df[ ,-1] <-lapply(scenarios.df[ ,-1], as.numeric)
View(scenarios.df)

scenarios.tot.df<-as.data.frame(scenarios.tot[-1,])
names(scenarios.tot.df) <- c("scenario", "harvest.prop", "tot.cost.lo", "tot.cost.med", "tot.cost.up", 
                             "tot.time.low", "tot.time.med", "tot.time.up", "min.n.lower", "min.n.median", "min.n.upper")
rownames(scenarios.tot.df)<-NULL
scenarios.tot.df[ ,-1] <-lapply(scenarios.tot.df[ ,-1], as.numeric)
View(scenarios.tot.df)

##### repeat cull scenarios projection for 10 years instead of 3 to compare costs and effort####

## set time limit for projection in 1-yr increments
yr.now <- 2020 
yr.end <- 2030 #end year for our projection timeframe

t <- (yr.end - yr.now) #timeframe

yrs <- seq(yr.now, yr.end, 1)

harv.prop.consist <- seq(0.35,0.95,0.05) # range of effective cull props for 10 year constant proportional cull

n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1)) 
cost.mat <- time.mat <- matrix(data=0, nrow=iter, ncol=t)
trap.number <- rep(NA, t)
trap.purchase.number <- rep(NA, t)
dispenser.no <- rep(NA, t)
#relative prop cost and total cost vectors
tot.TAAC.rcost <- tot.trap.rcost <- tot.poison.rcost <- tot.shot.rcost <-  
  tot.TAAC.cost <-tot.trap.cost <- tot.poison.cost <- tot.shot.cost <- tot.cost <-  
  tot.TAAC.hrs <-tot.trap.hrs <- tot.poison.hrs <- tot.shot.hrs <- tot.time <-  rep(NA, t)
prop.rcost.vec <-  inv.prop.rcost.vec <- weight.rcost.vec <- rep(NA, 4)

# cost, time and minimum n storage vector

tot.med <- tot.lo <- tot.up <- min.med.n <- min.lo.n <- min.up.n <- tot.time.med <- tot.time.lo <- tot.time.up <- rep(0,length(harv.prop.consist))


# storage matrices and vector for outcomes of each scenario at each harvest proportion
scenarios.mat <- matrix(NA, t, 10)
harv.prop.cost.mat <- matrix(NA, t, 10)

scenarios.tot <- rep(NA,11)

scenario.names <- c("Relative cost proportional allocation", 
                    " 25pc harvest per method", "Thermal-assisted aerial cull", "Shoot", "Trap", "Poison")

par(mfrow = c(1,1))

for (f in 1:nrow(cull.props.mat)){
  ## f loop runs projection and cost calculations for each cull scenario
  
  for (s in 1:length(harv.prop.consist)) {
    ## s loop runs
    for (e in 1:iter) {
      popmat <- popmat.orig
      
      harvested.tot <- tot.cost <- rep(NA,t)
      n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
      n.mat[,1] <- init.vec
      
      for (i in 1:t) {
        # stochastic survival values
        s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
        s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
        s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
        
        # stochastic fertility sampler (Gaussian)
        fert.stoch <- rtruncnorm(length(popmat.orig[,1]), a=0, b=Inf, mean=popmat.orig[1,], sd=f.sd.vec)
        
        totN.i <- sum(n.mat[,i])
        pigs.s.red <- as.numeric(s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp))
        #a/(1+(K.vec/b)^c)
        
        popmat[1,] <- fert.stoch
        diag(popmat[2:age.max,]) <- s.stoch[1:13]*pigs.s.red
        popmat[age.max,age.max] <- s.stoch[14]*pigs.s.red
        #popmat[age.max,age.max] <- 0
        
        n.mat[,i+1] <- popmat %*% n.mat[,i]
        
        ## harvest things here
        # create harvested component matrix
        harvested.tot[i] <- round(sum(n.mat[,i+1])*harv.prop.consist[s], 0)
        
        # do the harvest
        n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * harvested.tot[i], 0)
        
        # prevent n.mat from containing negative values. Cannot have less than 0 pigs
        if (length(which(n.mat[,i+1] < 0)) > 0) {
          n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
        }
        
        # calculate proportion remaining after each iteration
        remaining.preharv <- sum(n.mat[,i])
        prop.remain <- remaining.preharv/957*2
        
        # effort per control type relative to proportion remaining
        hrsPigTAAC <- ifelse(coef(fitTAACExp)[1] * exp(-coef(fitTAACExp)[2] * prop.remain)>=min(TAAC_efrt$hrsPig),
                             coef(fitTAACExp)[1] * exp(-coef(fitTAACExp)[2] * prop.remain),
                             min(TAAC_efrt$hrsPig))
        hrsPigShot <- ifelse(coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * prop.remain)>=min(shot_efrt$hrsPig), 
                             coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * prop.remain),
                             min(shot_efrt$hrsPig))
        hrsPigTrap <- ifelse(coef(fitTrapExp)[1] * exp(-coef(fitTrapExp)[2] * prop.remain)>=min(trap_efrt$hrsPig),
                             coef(fitTrapExp)[1] * exp(-coef(fitTrapExp)[2] * prop.remain),
                             min(trap_efrt$hrsPig))
        hrsPoison <- ifelse(coef(fitPoisonExp)[1] * exp(-coef(fitPoisonExp)[2] * prop.remain)>=min(poison_efrt$hrsPig), 
                            coef(fitPoisonExp)[1] * exp(-coef(fitPoisonExp)[2] * prop.remain),
                            min(poison_efrt$hrsPig))
        
        #### calculating and assigning proportions for relative cost method #####
        
        
        if (f == 1){ ##  calculate relative proportional cost allocation ####
          # TAAC
          tot.TAAC.hrs <- hrsPigTAAC*harvested.tot[i]
          ## cost accrued in increments of 111 hrs TAAC work completed.
          tot.TAAC.rcost[i] <- TAAC.initial + ((TAAC.pd*tot.TAAC.hrs)/TAAC.daily.ft) + bullets.pp*harvested.tot[i] + ((TAAC.resupply*tot.TAAC.hrs)/111)
          
          # shot
          tot.shot.hrs <- hrsPigShot*harvested.tot[i]
          labour.cost <- (labour.vehicle  + freefeed.ph.shot) * tot.shot.hrs
          bullet.cost <- bullets.pp * harvested.tot[i]
          tot.shot.rcost[i] <- labour.cost + bullet.cost
          
          # Trap
          tot.trap.hrs <- hrsPigTrap*harvested.tot[i]
          trap.number[i] <-ceiling((tot.trap.hrs/(10*(20+mean(trap_efrt$efforthrs)))))
          trap.purchase.number[i] <- trap.number[i]# ceiling used so number of traps is an integer
          # need to deduct number of traps purchased in previous years
          if (i > 1) {
            trap.purchase.number[i] <- ifelse(trap.purchase.number[i] < trap.purchase.number[1], 0, trap.purchase.number[i] - trap.purchase.number[1])
          }
          tot.trap.rcost[i] <-(trap.purchase.number[i]*trap.unit) + ## cost to purchase no. traps required in year [i]
            (trap.number[i] * 10 * (trap.initial + trap.ff.cost)) + # cost to set up and feed no. traps required in year [i]. Assuming 10 events per trap per year
            bullets.pp*harvested.tot[i]
          
          #poison
          tot.poison.hrs <- hrsPoison*harvested.tot[i]
          total.poison.events.pa <- tot.poison.hrs/avg.poison.effort # no of events
          dispenser.no[i] <- ceiling((total.poison.events.pa/max.annual.poison.events) * 1.45)
          if (i > 1) {
            dispenser.no[i] <- ifelse(dispenser.no[i] < dispenser.no[1], 0, dispenser.no[i] - dispenser.no[1])
          }
          tot.poison.rcost[i] <- (total.poison.events.pa * (poison.labour + poison.ff)) + 
            (total.poison.events.pa * (1.45 * (poison.placebo + poison.toxic))) + 
            (dispenser.no[i]*poison.dispenser.unit)
          
          tot.rcost <- tot.TAAC.rcost[i] + tot.shot.rcost[i] + tot.trap.rcost[i] + tot.poison.rcost[i]
          tot.rcost.vec <- c(tot.TAAC.rcost[i],tot.shot.rcost[i],tot.trap.rcost[i],tot.poison.rcost[i])
          cull.props.mat[1,] <- (1/(tot.rcost.vec/tot.rcost))/sum(1/(tot.rcost.vec/tot.rcost))
          # populates cull.props.mat with relative cost weighting
        }
        
        harvest.TAAC <- harvested.tot[i]*cull.props.mat[f,1]
        harvest.shot <- harvested.tot[i]*cull.props.mat[f,2]
        harvest.trap <- harvested.tot[i]*cull.props.mat[f,3]
        harvest.poison <- harvested.tot[i]*cull.props.mat[f,4]
        
        #### actual cost calculation #### 
        # relative cost prop weighting
        # TAAC
        tot.TAAC.hrs[i] <- hrsPigTAAC*harvest.TAAC
        tot.TAAC.cost[i] <- TAAC.initial + ((TAAC.pd*tot.TAAC.hrs[i])/TAAC.daily.ft) + bullets.pp*harvest.TAAC + ((TAAC.resupply*tot.TAAC.hrs[i])/111)
        
        # shot
        tot.shot.hrs[i] <- hrsPigShot*harvest.shot
        labour.cost <- (labour.vehicle  + freefeed.ph.shot) * tot.shot.hrs[i]
        bullet.cost <- bullets.pp * harvest.shot
        tot.shot.cost[i] <- labour.cost + bullet.cost
        
        # Trap
        tot.trap.hrs[i] <- hrsPigTrap*harvest.trap
        trap.number[i] <- ceiling((tot.trap.hrs[i]/(10*(20+mean(trap_efrt$efforthrs))))) # ceiling used so number of traps is an integer
        trap.purchase.number[i] <- trap.number[i]
        # if loop to deduct number of traps purchased in previous years as these can be reused in subsequent years
        if (i > 1) {
          trap.purchase.number[i] <- ifelse(trap.purchase.number[i] < trap.purchase.number[1], 0, trap.purchase.number[i] - trap.purchase.number[1])
        }
        tot.trap.cost[i] <-(trap.purchase.number[i]*trap.unit) + ## cost to purchase no. traps required in year [i]
          (trap.number[i] * 10 * (trap.initial + trap.ff.cost)) + # cost to set up and feed no. traps required in year [i]. Assuming 10 events per trap per year
          bullets.pp*harvest.trap
        
        #poison
        tot.poison.hrs[i] <- hrsPoison*harvest.poison
        total.poison.events.pa <- tot.poison.hrs[i]/avg.poison.effort # no of events
        dispenser.no[i] <- ceiling((total.poison.events.pa/max.annual.poison.events) * 1.45)
        if (i > 1) {
          dispenser.no[i] <- ifelse(dispenser.no[i] < dispenser.no[1], 0, dispenser.no[i] - dispenser.no[1])
        }
        tot.poison.cost[i] <- (total.poison.events.pa * (poison.labour + poison.ff)) + 
          total.poison.events.pa * (1.45 * (poison.placebo + poison.toxic)) + 
          (dispenser.no[i]*poison.dispenser.unit)
        
        tot.cost[i] <- sum(c(tot.TAAC.cost[i],tot.shot.cost[i],tot.trap.cost[i],tot.poison.cost[i]))
        tot.time[i] <- sum(c(tot.TAAC.hrs[i], tot.shot.hrs[i], tot.trap.hrs[i], tot.poison.hrs[i]))
        
      } # end i loop
      
      cost.mat[e,] <- tot.cost
      time.mat[e,] <- tot.time
      
      n.sums.mat[e,] <- as.vector((colSums(n.mat))/pop.found) # / pop.mat for min proportion remaining population 
      
      if (e %% itdiv==0) print(e) 
      
    } # end e loop
    
    # Calculate total median cost, upper and lower quantile costs
    
    median.cost <- apply(cost.mat, MARGIN=2, median, na.rm=T)
    up.cost <-  apply(cost.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
    lo.cost <-  apply(cost.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
    
    tot.med[s] <- sum(median.cost)
    tot.up[s] <- sum(up.cost)
    tot.lo[s] <- sum(lo.cost)
    print(c(tot.lo[s], tot.med[s], tot.up[s]))
    
    #calculate total median time + upper and lower 95% confidence intervals
    median.time <- apply(time.mat, MARGIN=2, median, na.rm=T)
    up.time <-  apply(time.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
    lo.time <-  apply(time.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
    
    tot.time.med[s] <- sum(median.time)
    tot.time.up[s] <- sum(up.time)
    tot.time.lo[s] <- sum(lo.time)
    print(c(tot.time.lo[s], tot.time.med[s], tot.time.up[s]))
    
    # plot costper year for each scenario and cull proportion combination
    plot(yrs[1:(t)], median.cost, type="l", lwd=2, xlab="year", ylab="cost ($)", 
         main = paste("Annual costs:", scenario.names[f]), 
         sub = paste("cull proportion", harv.prop.consist[s]), 
         ylim=(c(min(lo.cost),max(up.cost))))
    lines(yrs[1:(t)], up.cost, lty=2, col="red")  
    lines(yrs[1:(t)], lo.cost, lty=2, col="red")
    
    # calculate minimum population size for each 
    min.pop.vec <- apply(n.sums.mat, MARGIN=1, min)
    min.med.n[s] <- median(min.pop.vec, na.rm=T)
    min.lo.n[s] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
    min.up.n[s] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
    
    # calculate median prop remaining after each year
    med.N.Yr.vec <- apply(n.sums.mat, MARGIN = 2, median)
    
    print(paste("harvest proportion = ", harv.prop.consist[s], sep=""))
    
    # matrix of lower, median and upper cost per year for each scenario and proportion
    harv.prop.cost.mat[,1] <- scenario.names[f]
    harv.prop.cost.mat[,2] <- harv.prop.consist[s]
    harv.prop.cost.mat[,3] <- yrs[1:(t)]
    harv.prop.cost.mat[,4] <- lo.cost
    harv.prop.cost.mat[,5] <- median.cost
    harv.prop.cost.mat[,6] <- up.cost
    harv.prop.cost.mat[,7] <- lo.time
    harv.prop.cost.mat[,8] <- median.time
    harv.prop.cost.mat[,9] <- up.time
    harv.prop.cost.mat[1:10,10] <- med.N.Yr.vec[2:11]
    
    scenarios.mat <-rbind(scenarios.mat, harv.prop.cost.mat)
    
    # matrix of lower, median and upper total cost and minimum N for each scenario and proportion
    
    harv.prop.annual <- c(scenario.names[f], harv.prop.consist[s],
                          tot.lo[s], tot.med[s], tot.up[s], 
                          tot.time.lo[s], tot.time.med[s], tot.time.up[s], 
                          min.lo.n[s], min.med.n[s], min.up.n[s])
    
    scenarios.tot <-rbind(scenarios.tot, harv.prop.annual)
    
  } # end s loop
} # end f loop

scenarios10.df <- as.data.frame(scenarios.mat[-(1:10),])
names(scenarios10.df) <- c("scenario", "harvest.prop", "year","low.cost", "median.cost", "upper.cost", 
                           "low.time", "median.time", "upper.time", "median.N") 
rownames(scenarios10.df)<-NULL
scenarios10.df[ ,-1] <-lapply(scenarios10.df[ ,-1], as.numeric)
View(scenarios10.df)

scenarios10.tot.df<-as.data.frame(scenarios.tot[-1,])
names(scenarios10.tot.df) <- c("scenario", "harvest.prop", "tot.cost.lo", "tot.cost.med", "tot.cost.up", 
                               "tot.time.low", "tot.time.med", "tot.time.up", "min.n.lower", "min.n.median", "min.n.upper")
rownames(scenarios10.tot.df)<-NULL
scenarios10.tot.df[ ,-1] <-lapply(scenarios10.tot.df[ ,-1], as.numeric)
View(scenarios10.tot.df)

############### Section not used ##################
# add columns to scenarios.df and scenarios10.df calculating ongoing costs from pig damages as a function of N remaining @ $185.19/pig
# assumes costs are linear
# scenarios.df <- scenarios.df %>%
#   mutate(
#     `female pigs remaining` = median.N * 250,
#     `ongoing annual pig damages` = median.N * 500 * 185.19
#   )
# print(scenarios.df)
# 
# scenarios10.df <- scenarios10.df %>%
#   mutate(
#     `female pigs remaining` = median.N * 250,
#     `ongoing annual pig damages` = median.N * 500 * 185.19
#   )
# print(scenarios10.df)

################ 

write_csv(scenarios.df, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/ChoScen.df.csv")
write_csv(scenarios.tot.df, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/ChoScen.tot.df.csv")
write_csv(scenarios10.df, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/ChoScen10.df.csv")
write_csv(scenarios10.tot.df, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/ChoScen10.tot.df.csv")

# use scenarios10.df and scenarios10.tot.df to produce combined plots fof cost per year for all scenarios per cull proportion
par(mfrow=c(1,1))
# set a colourblind friendly palette
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#D55E00", "#CC79A7", "#0072B2")
library(scales)

# plot of total cost per scenario after 3 years
p3C<-ggplot(scenarios.tot.df, aes(x = harvest.prop, y = tot.cost.med, group=scenario, colour=scenario)) +  
  geom_line() +  theme_bw() +scale_colour_manual(values=cbp2) +
  labs(title = "", subtitle =  "", y = "", x ="") + geom_vline(xintercept = 0.8, linetype = 3, color = "black") +
  scale_y_continuous(labels = comma)

# plot of total time per scenario after 3 years
p3T<- ggplot(scenarios.tot.df, aes(x = harvest.prop, y = tot.time.med, group=scenario, colour=scenario)) +  
  geom_line() + theme_bw() + scale_colour_manual(values=cbp2) +
  labs(title = "", subtitle =  "", y = "", x ="") + geom_vline(xintercept = 0.8, linetype = 3, color = "black")

# plot of total cost per scenario after 10 years
p10C<-ggplot(scenarios10.tot.df, aes(x = harvest.prop, y = tot.cost.med, group=scenario, colour=scenario)) +  
  geom_line() + theme_bw() + scale_colour_manual(values=cbp2) +
        labs(title = "", subtitle =  "", y = "", x ="") + geom_vline(xintercept = 0.6, linetype = 3, color = "black") +
        scale_y_continuous(labels = comma)

# plot of total time per scenario after 10 years

p10T <-ggplot(scenarios10.tot.df, aes(x = harvest.prop, y = tot.time.med, group=scenario, colour=scenario)) +  
  geom_line() + theme_bw() + scale_colour_manual(values=cbp2) +
        labs(title = "", subtitle =  "", y = "", x ="") + geom_vline(xintercept = 0.5, linetype = 3, color = "black")

# Thesis figures 5 and 7 - results
totcostfig <- ggarrange(p3C, p10C, common.legend = TRUE, labels = c("a", "b"), legend = "bottom", nrow = 1, ncol = 2)
toteffortfig <- ggarrange(p3T, p10T, common.legend = TRUE, labels = c("a", "b"), legend = "bottom", nrow = 1, ncol = 2)


print(totcostfig)
print(toteffortfig)

# identify least cost scenario above the minimum cull threshold
# 3-year projection interval
threshold3<- 0.7
filtered_df3 <- scenarios.tot.df %>%
  filter(harvest.prop >= threshold3) %>%  # Filter rows based on the harvest.prop column
  group_by(scenario) %>%
  filter(tot.cost.med == min(tot.cost.med)) %>%
  ungroup()

print(filtered_df3)

# 10 year projection interval
threshold10<- 0.5
filtered_df10 <- scenarios10.tot.df %>%
  filter(harvest.prop >= threshold10) %>%  # Filter rows based on the harvest.prop column
  group_by(scenario) %>%
  filter(tot.cost.med == min(tot.cost.med)) %>%
  ungroup()

print(filtered_df10)

write_csv(filtered_df3, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/Cho_least_cost3.csv")
write_csv(filtered_df10, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/Cho_least_cost10.df.csv")


########################################################################################################
# Appendix S1: plots of effort per pig killed as a function of pig population density (proportion of pig population remaining). 

# First, we fit the models to raw data for each control methods,
# but found that functional response curves for shooting and trapping didn't conform to the expected shape
# poisoning revealed only a weak response
# highlight data contributed by individuals or organisations


# 1.1: create subsets for raw data for each control technique
shot_efrt_raw <- effort %>% filter (controlType == "shot")
TAAC_efrt_raw <-effort %>% filter (controlType == "TAAC")
TAAC_efrt_raw <- TAAC_efrt_raw[which(is.infinite(TAAC_efrt_raw$hrsPig)==F),] # remove entries where TAC effort was infinite i.e., flights where no pigs were killed (0 kills/hours = infinite value)
trap_efrt_raw <-effort %>% filter (controlType == "trapped")
poison_efrt_raw <- effort %>% filter (controlType == "poisoned")


# fit logarithmic and exponential models and plot 
# Appendix S1: fig. S1

par(mfrow=c(1,3))
par(mar = c(3,2,3,2))

#shot
plot(hrsPig~propRemaining, data = shot_efrt_raw, xlab = '', ylab = '', pch=19, cex=0.7, col = ifelse(operator == "PJ","red","black") , xlim=c(0,1), ylim=c(0,110))
mtext("a", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit
fitShot <- lm(hrsPig ~ log(propRemaining), data=shot_efrt_raw, na.action =) # logarithmic fit
summary(fitShot)
linreg.ER(log(shot_efrt_raw$propRemaining),shot_efrt_raw$hrsPig)
hrsPigShot.pred <- coef(fitShot)[1] + coef(fitShot)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigShot.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)j
s.param.init <- c(200, 5)
fitShotExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                  data = shot_efrt_raw,
                  #algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,      
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fitShotExp.summ <- summary(fitShotExp)
fitShotExp.summ
hrsPigShotExp.pred <- coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigShotExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
shot_raw_mods <-data.frame(propRemaining.cont, hrsPigShot.pred, hrsPigShotExp.pred)
colnames(shot_raw_mods) <- c("pRemain", "logfit", "expfit")
write.csv(shot_raw_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_raw_mods.csv")
shot_raw <-data.frame(shot_efrt_raw$propRemaining,shot_efrt_raw$hrsPig)
colnames(shot_raw)<- c("pRemain", "effort")
write.csv(shot_raw, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_raw.csv")

#trapped
plot(hrsPig~propRemaining, data = trap_efrt_raw, xlab = '', ylab = '', pch=19, cex=0.7, col = ifelse(org == "DEW", "red","black"),xlim=c(0,1), ylim=c(0,40))
mtext("b", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit
fitTrap <- lm(hrsPig ~ log(propRemaining), data=trap_efrt_raw, na.action = ) # logarithmic fit
summary(fitTrap)
linreg.ER(log(trap_efrt_raw$propRemaining),trap_efrt_raw$hrsPig)
hrsPigTrap.pred <- coef(fitTrap)[1] + coef(fitTrap)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigTrap.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(12, 2)
fitTrapExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                  data = trap_efrt_raw,
                  algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,      
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fitTrapExp.summ <- summary(fitTrapExp)
hrsPigTrapExp.pred <- coef(fitTrapExp)[1] * exp(-coef(fitTrapExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigTrapExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
trap_raw_mods <-data.frame(propRemaining.cont, hrsPigTrap.pred, hrsPigTrapExp.pred)
colnames(trap_raw_mods) <- c("pRemain", "logfit", "expfit")
write.csv(trap_raw_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/trap_raw_mods.csv")
trap_raw <-data.frame(trap_efrt_raw$propRemaining,trap_efrt_raw$hrsPig)
colnames(trap_raw)<- c("pRemain", "effort")
write.csv(trap_raw, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/trap_raw.csv")

#poisoned
plot(hrsPig~propRemaining, data = poison_efrt_raw, xlab = '', ylab = '', pch=19, cex=0.7, col = ifelse(org == "LB","red", "black"),xlim=c(0,1), ylim=c(0,75))
mtext("c", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit: y = a + b*log(x)
fitPoison <- lm(hrsPig ~ log(propRemaining), data=poison_efrt_raw) # logarithmic fit
summary(fitPoison)
linreg.ER(log(poison_efrt_raw$propRemaining),poison_efrt_raw$hrsPig)
lpropRemaining.cont <- log(seq(0.01, 1, 0.01))
hrsPigPoison.pred <- coef(fitPoison)[1] + coef(fitPoison)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigPoison.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(200, 5)
fitPoisonExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                    data = poison_efrt_raw,
                    algorithm = "port",
                    start = c(a = s.param.init[1], b = s.param.init[2]),
                    trace = TRUE,      
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fitPoisonExp.summ <- summary(fitPoisonExp)
propRemaining.cont <- exp(lpropRemaining.cont)
hrsPigPoisonExp.pred <- coef(fitPoisonExp)[1] * exp(-coef(fitPoisonExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigPoisonExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
poison_raw_mods <-data.frame(propRemaining.cont, hrsPigPoison.pred, hrsPigPoisonExp.pred)
colnames(poison_raw_mods) <- c("pRemain", "logfit", "expfit")
write.csv(poison_raw_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/poison_raw_mods.csv")
poison_raw <-data.frame(poison_efrt_raw$propRemaining,poison_efrt_raw$hrsPig)
colnames(poison_raw)<- c("pRemain", "effort")
write.csv(poison_raw, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/poison_raw.csv")

#2: stratify control types by operator or organisation
# note: only 1 operator for TAAC so no subset created

#2.1: shooting
shot_npws <-  shot_efrt_raw %>% filter (org == "NPWS")
shot_other <-  shot_efrt_raw %>% filter (org == "other")
shot_farmer <-  shot_efrt_raw %>% filter (org == "farmer")
shot_op_other <-  shot_efrt_raw %>% filter (operator == "other")
shot_BF <- shot_efrt_raw %>% filter (operator == "BF")
shot_DJ <- shot_efrt_raw %>% filter (operator == "DJ")
shot_PJ <- shot_efrt_raw %>% filter (operator == "PJ")

# Fig. S2: shooting by operator and organisation/group
par(mfrow=c(3,3))
par(mar = c(3,3,3,3))
# Figures for appendix

#shot_DJ
plot(hrsPig~propRemaining, data = shot_DJ, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,110))
mtext("DJ", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit
fitShot <- lm(hrsPig ~ log(propRemaining), data=shot_DJ, na.action =) # logarithmic fit
summary(fitShot)
linreg.ER(log(shot_DJ$propRemaining),shot_DJ$hrsPig)
hrsPigShot.pred <- coef(fitShot)[1] + coef(fitShot)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigShot.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(200, 5)
fitShotExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                  data = shot_PJ,
                  #algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,      
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/2048))
fitShotExp.summ <- summary(fitShotExp)
fitShotExp.summ
hrsPigShotExp.pred <- coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigShotExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
shot_DJ_mods <-data.frame(propRemaining.cont, hrsPigShot.pred, hrsPigShotExp.pred)
colnames(shot_DJ_mods) <- c("pRemain", "logfit", "expfit")
write.csv(shot_DJ_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_DJ_mods.csv")
shot_DJ_ef <-data.frame(shot_DJ$propRemaining,shot_DJ$hrsPig)
colnames(shot_DJ_ef)<- c("pRemain", "effort")
write.csv(shot_DJ_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_DJ_ef.csv")

#shot_BF
plot(hrsPig~propRemaining, data = shot_BF, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,110))
mtext("BF", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit
fitShot <- lm(hrsPig ~ log(propRemaining), data=shot_BF, na.action =) # logarithmic fit
summary(fitShot)
linreg.ER(log(shot_BF$propRemaining),shot_BF$hrsPig)
hrsPigShot.pred <- coef(fitShot)[1] + coef(fitShot)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigShot.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(200, 5)
fitShotExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                  data = shot_BF,
                  #algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,      
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/2048))
fitShotExp.summ <- summary(fitShotExp)
fitShotExp.summ
hrsPigShotExp.pred <- coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigShotExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
shot_BF_mods <-data.frame(propRemaining.cont, hrsPigShot.pred, hrsPigShotExp.pred)
colnames(shot_BF_mods) <- c("pRemain", "logfit", "expfit")
write.csv(shot_BF_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_BF_mods.csv")
shot_BF_ef <-data.frame(shot_BF$propRemaining,shot_BF$hrsPig)
colnames(shot_BF_ef)<- c("pRemain", "effort")
write.csv(shot_BF_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_BF_ef.csv")

#shot_PJ
plot(hrsPig~propRemaining, data = shot_PJ, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,110))
mtext("PJ", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit
fitShot <- lm(hrsPig ~ log(propRemaining), data=shot_DJ, na.action =) # logarithmic fit
summary(fitShot)
linreg.ER(log(shot_DJ$propRemaining),shot_DJ$hrsPig)
hrsPigShot.pred <- coef(fitShot)[1] + coef(fitShot)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigShot.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(200, 5)
fitShotExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                  data = shot_DJ,
                  #algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,      
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/2048))
fitShotExp.summ <- summary(fitShotExp)
fitShotExp.summ
hrsPigShotExp.pred <- coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigShotExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
shot_PJ_mods <-data.frame(propRemaining.cont, hrsPigShot.pred, hrsPigShotExp.pred)
colnames(shot_PJ_mods) <- c("pRemain", "logfit", "expfit")
write.csv(shot_PJ_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_PJ_mods.csv")
shot_PJ_ef <-data.frame(shot_PJ$propRemaining,shot_PJ$hrsPig)
colnames(shot_PJ_ef)<- c("pRemain", "effort")
write.csv(shot_PJ_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_PJ_ef.csv")

#shot_op_other
plot(hrsPig~propRemaining, data = shot_op_other, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,110))
mtext("other_individuals", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit
fitShot <- lm(hrsPig ~ log(propRemaining), data=shot_op_other, na.action =) # logarithmic fit
summary(fitShot)
linreg.ER(log(shot_op_other$propRemaining),shot_op_other$hrsPig)
hrsPigShot.pred <- coef(fitShot)[1] + coef(fitShot)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigShot.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(200, 5)
fitShotExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                  data = shot_op_other,
                  #algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,      
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/2048))
fitShotExp.summ <- summary(fitShotExp)
fitShotExp.summ
hrsPigShotExp.pred <- coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigShotExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
shot_op_other_mods <-data.frame(propRemaining.cont, hrsPigShot.pred, hrsPigShotExp.pred)
colnames(shot_op_other_mods) <- c("pRemain", "logfit", "expfit")
write.csv(shot_op_other_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_op_other_mods.csv")
shot_op_other_ef <-data.frame(shot_op_other$propRemaining,shot_op_other$hrsPig)
colnames(shot_op_other_ef)<- c("pRemain", "effort")
write.csv(shot_op_other_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_op_other_ef.csv")

#shot_npws
plot(hrsPig~propRemaining, data = shot_npws, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,110))
mtext("NPWS", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit
fitShot <- lm(hrsPig ~ log(propRemaining), data=shot_npws, na.action =) # logarithmic fit
summary(fitShot)
linreg.ER(log(shot_npws$propRemaining),shot_npws$hrsPig)
hrsPigShot.pred <- coef(fitShot)[1] + coef(fitShot)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigShot.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(200, 5)
fitShotExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                  data = shot_npws,
                  #algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,      
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/2048))
fitShotExp.summ <- summary(fitShotExp)
fitShotExp.summ
hrsPigShotExp.pred <- coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigShotExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
shot_npws_mods <-data.frame(propRemaining.cont, hrsPigShot.pred, hrsPigShotExp.pred)
colnames(shot_npws_mods) <- c("pRemain", "logfit", "expfit")
write.csv(shot_npws_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_npws_mods.csv")
shot_npws_ef <-data.frame(shot_npws$propRemaining,shot_npws$hrsPig)
colnames(shot_npws_ef)<- c("pRemain", "effort")
write.csv(shot_npws_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_npws_ef.csv")

#shot_farmer
plot(hrsPig~propRemaining, data = shot_farmer, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,110))
mtext("farmers", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit
fitShot <- lm(hrsPig ~ log(propRemaining), data=shot_farmer, na.action =) # logarithmic fit
summary(fitShot)
linreg.ER(log(shot_farmer$propRemaining),shot_farmer$hrsPig)
hrsPigShot.pred <- coef(fitShot)[1] + coef(fitShot)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigShot.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(200, 5)
fitShotExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                  data = shot_farmer,
                  #algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,      
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/2048))
fitShotExp.summ <- summary(fitShotExp)
fitShotExp.summ
hrsPigShotExp.pred <- coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigShotExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
shot_farmer_mods <-data.frame(propRemaining.cont, hrsPigShot.pred, hrsPigShotExp.pred)
colnames(shot_farmer_mods) <- c("pRemain", "logfit", "expfit")
write.csv(shot_farmer_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_farmer_mods.csv")
shot_farmer_ef <-data.frame(shot_farmer$propRemaining,shot_farmer$hrsPig)
colnames(shot_farmer_ef)<- c("pRemain", "effort")
write.csv(shot_farmer_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_farmer_ef.csv")

#shot_other
plot(hrsPig~propRemaining, data = shot_other, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,110))
mtext("other organisations/operator types", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit
fitShot <- lm(hrsPig ~ log(propRemaining), data=shot_other, na.action =) # logarithmic fit
summary(fitShot)
linreg.ER(log(shot_other$propRemaining),shot_other$hrsPig)
hrsPigShot.pred <- coef(fitShot)[1] + coef(fitShot)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigShot.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(200, 5)
fitShotExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                  data = shot_other,
                  #algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,      
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/2048))
fitShotExp.summ <- summary(fitShotExp)
fitShotExp.summ
hrsPigShotExp.pred <- coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigShotExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
shot_other_mods <-data.frame(propRemaining.cont, hrsPigShot.pred, hrsPigShotExp.pred)
colnames(shot_other_mods) <- c("pRemain", "logfit", "expfit")
write.csv(shot_other_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_other_mods.csv")
shot_other_ef <-data.frame(shot_other$propRemaining,shot_other$hrsPig)
colnames(shot_other_ef)<- c("pRemain", "effort")
write.csv(shot_other_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/shot_other_ef.csv")

#2.2: trapping
trap_DEW <-trap_efrt_raw %>% filter (org == "DEW") 
trap_NPWS <-trap_efrt_raw %>% filter (org == "NPWS")
trap_PIRSA <-trap_efrt_raw %>% filter (org == "PIRSA")

# fig. S3: trapping by operator and organisation/group
par(mfrow=c(1,3))
par(mar = c(3,3,3,3))

#trapped NPWS
plot(hrsPig~propRemaining, data = trap_NPWS, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,40))
mtext("NPWS", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit
fitTrap <- lm(hrsPig ~ log(propRemaining), data=trap_NPWS, na.action = ) # logarithmic fit
summary(fitTrap)
linreg.ER(log(trap_NPWS$propRemaining),trap_NPWS$hrsPig)
hrsPigTrap.pred <- coef(fitTrap)[1] + coef(fitTrap)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigTrap.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(12, 2)
fitTrapExp <- nls(hrsPig ~ a * exp(-b*propRemaining),
                  data = trap_NPWS,
                  algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fitTrapExp.summ <- summary(fitTrapExp)
hrsPigTrapExp.pred <- coef(fitTrapExp)[1] * exp(-coef(fitTrapExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigTrapExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
trap_NPWS_mods <-data.frame(propRemaining.cont, hrsPigTrap.pred, hrsPigTrapExp.pred)
colnames(trap_NPWS_mods) <- c("pRemain", "logfit", "expfit")
write.csv(trap_NPWS_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/trap_NPWS_mods.csv")
trap_NPWS_ef <-data.frame(trap_NPWS$propRemaining,trap_NPWS$hrsPig)
colnames(trap_NPWS_ef)<- c("pRemain", "effort")
write.csv(trap_NPWS_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/trap_NPWS_ef.csv")

#trapped PIRSA
plot(hrsPig~propRemaining, data = trap_PIRSA, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,40))
mtext("PIRSA", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit
fitTrap <- lm(hrsPig ~ log(propRemaining), data=trap_PIRSA, na.action = ) # logarithmic fit
summary(fitTrap)
linreg.ER(log(trap_PIRSA$propRemaining),trap_PIRSA$hrsPig)
hrsPigTrap.pred <- coef(fitTrap)[1] + coef(fitTrap)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigTrap.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(12, 2)
fitTrapExp <- nls(hrsPig ~ a * exp(-b*propRemaining),
                  data = trap_PIRSA,
                  algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fitTrapExp.summ <- summary(fitTrapExp)
hrsPigTrapExp.pred <- coef(fitTrapExp)[1] * exp(-coef(fitTrapExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigTrapExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
trap_PIRSA_mods <-data.frame(propRemaining.cont, hrsPigTrap.pred, hrsPigTrapExp.pred)
colnames(trap_PIRSA_mods) <- c("pRemain", "logfit", "expfit")
write.csv(trap_PIRSA_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/trap_PIRSA_mods.csv")
trap_PIRSA_ef <-data.frame(trap_PIRSA$propRemaining,trap_PIRSA$hrsPig)
colnames(trap_PIRSA_ef)<- c("pRemain", "effort")
write.csv(trap_PIRSA_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/trap_PIRSA_ef.csv")

#trapped DEW
plot(hrsPig~propRemaining, data = trap_DEW, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,40))
mtext("DEW", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit
fitTrap <- lm(hrsPig ~ log(propRemaining), data=trap_DEW, na.action = ) # logarithmic fit
summary(fitTrap)
linreg.ER(log(trap_efrt$propRemaining),trap_DEW$hrsPig)
hrsPigTrap.pred <- coef(fitTrap)[1] + coef(fitTrap)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigTrap.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(12, 2)
fitTrapExp <- nls(hrsPig ~ a * exp(-b*propRemaining),
                  data = trap_DEW,
                  algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fitTrapExp.summ <- summary(fitTrapExp)
hrsPigTrapExp.pred <- coef(fitTrapExp)[1] * exp(-coef(fitTrapExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigTrapExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
trap_DEW_mods <-data.frame(propRemaining.cont, hrsPigTrap.pred, hrsPigTrapExp.pred)
colnames(trap_DEW_mods) <- c("pRemain", "logfit", "expfit")
write.csv(trap_DEW_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/trap_DEW_mods.csv")
trap_DEW_ef <-data.frame(trap_DEW$propRemaining,trap_DEW$hrsPig)
colnames(trap_DEW_ef)<- c("pRemain", "effort")
write.csv(trap_DEW_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/trap_DEW_ef.csv")

#2.3: poisoning
poison_LB <- poison_efrt_raw %>% filter (org == "LB") 
poison_NPWS <- poison_efrt_raw %>% filter (org == "NPWS")
poison_DEW <- poison_efrt_raw %>% filter (org == "DEW")
poison_other <- poison_efrt_raw %>% filter (org == "other")

#fig. S4: poisoning by operator and organisation.
par(mfrow=c(2,2))
par(mar = c(3,3,3,3))

#Poisoned NPWS
plot(hrsPig~propRemaining, data = poison_NPWS, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,75))
mtext("NPWS", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit: y = a + b*log(x)
fitPoison <- lm(hrsPig ~ log(propRemaining), data=poison_NPWS) # logarithmic fit
summary(fitPoison)
linreg.ER(log(poison_NPWS$propRemaining),poison_NPWS$hrsPig)
lpropRemaining.cont <- log(seq(0.01, 1, 0.01))
hrsPigPoison.pred <- coef(fitPoison)[1] + coef(fitPoison)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigPoison.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(200, 5)
fitPoisonExp <- nls(hrsPig ~ a * exp(-b*propRemaining),
                    data = poison_NPWS,
                    algorithm = "port",
                    start = c(a = s.param.init[1], b = s.param.init[2]),
                    trace = TRUE,
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fitPoisonExp.summ <- summary(fitPoisonExp)
propRemaining.cont <- exp(lpropRemaining.cont)
hrsPigPoisonExp.pred <- coef(fitPoisonExp)[1] * exp(-coef(fitPoisonExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigPoisonExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
poison_NPWS_mods <-data.frame(propRemaining.cont, hrsPigPoison.pred, hrsPigPoisonExp.pred)
colnames(poison_NPWS_mods) <- c("pRemain", "logfit", "expfit")
write.csv(poison_NPWS_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/poison_NPWS_mods.csv")
poison_NPWS_ef <-data.frame(poison_NPWS$propRemaining,poison_NPWS$hrsPig)
colnames(poison_NPWS_ef)<- c("pRemain", "effort")
write.csv(poison_NPWS_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/poison_NPWS_ef.csv")

#Poisoned DEW
plot(hrsPig~propRemaining, data = poison_DEW, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,75))
mtext("DEW", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit: y = a + b*log(x)
fitPoison <- lm(hrsPig ~ log(propRemaining), data=poison_DEW) # logarithmic fit
summary(fitPoison)
linreg.ER(log(poison_DEW$propRemaining),poison_DEW$hrsPig)
lpropRemaining.cont <- log(seq(0.01, 1, 0.01))
hrsPigPoison.pred <- coef(fitPoison)[1] + coef(fitPoison)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigPoison.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(200, 5)
fitPoisonExp <- nls(hrsPig ~ a * exp(-b*propRemaining),
                    data = poison_DEW,
                    algorithm = "port",
                    start = c(a = s.param.init[1], b = s.param.init[2]),
                    trace = TRUE,
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fitPoisonExp.summ <- summary(fitPoisonExp)
propRemaining.cont <- exp(lpropRemaining.cont)
hrsPigPoisonExp.pred <- coef(fitPoisonExp)[1] * exp(-coef(fitPoisonExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigPoisonExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
poison_DEW_mods <-data.frame(propRemaining.cont, hrsPigPoison.pred, hrsPigPoisonExp.pred)
colnames(poison_DEW_mods) <- c("pRemain", "logfit", "expfit")
write.csv(poison_DEW_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/poison_DEW_mods.csv")
poison_DEW_ef <-data.frame(poison_DEW$propRemaining,poison_DEW$hrsPig)
colnames(poison_DEW_ef)<- c("pRemain", "effort")
write.csv(poison_DEW_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/poison_DEW_ef.csv")

#Poisoned LB
plot(hrsPig~propRemaining, data = poison_LB, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,75))
mtext("KILB", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit: y = a + b*log(x)
fitPoison <- lm(hrsPig ~ log(propRemaining), data=poison_LB) # logarithmic fit
summary(fitPoison)
linreg.ER(log(poison_LB$propRemaining),poison_LB$hrsPig)
lpropRemaining.cont <- log(seq(0.01, 1, 0.01))
hrsPigPoison.pred <- coef(fitPoison)[1] + coef(fitPoison)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigPoison.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(200, 5)
fitPoisonExp <- nls(hrsPig ~ a * exp(-b*propRemaining),
                    data = poison_LB,
                    algorithm = "port",
                    start = c(a = s.param.init[1], b = s.param.init[2]),
                    trace = TRUE,
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fitPoisonExp.summ <- summary(fitPoisonExp)
propRemaining.cont <- exp(lpropRemaining.cont)
hrsPigPoisonExp.pred <- coef(fitPoisonExp)[1] * exp(-coef(fitPoisonExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigPoisonExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
poison_LB_mods <-data.frame(propRemaining.cont, hrsPigPoison.pred, hrsPigPoisonExp.pred)
colnames(poison_LB_mods) <- c("pRemain", "logfit", "expfit")
write.csv(poison_LB_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/poison_LB_mods.csv")
poison_LB_ef <-data.frame(poison_LB$propRemaining,poison_LB$hrsPig)
colnames(poison_LB_ef)<- c("pRemain", "effort")
write.csv(poison_LB_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/poison_LB_ef.csv")

#Poisoned Other
plot(hrsPig~propRemaining, data = poison_other, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,75))
mtext("other", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
# logarithmic fit: y = a + b*log(x)
fitPoison <- lm(hrsPig ~ log(propRemaining), data=poison_other) # logarithmic fit
summary(fitPoison)
linreg.ER(log(poison_other$propRemaining),poison_other$hrsPig)
lpropRemaining.cont <- log(seq(0.01, 1, 0.01))
hrsPigPoison.pred <- coef(fitPoison)[1] + coef(fitPoison)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigPoison.pred, lty=2, col="red")
# exponential: y = a * exp(-b*x)
s.param.init <- c(200, 5)
fitPoisonExp <- nls(hrsPig ~ a * exp(-b*propRemaining),
                    data = poison_other,
                    algorithm = "port",
                    start = c(a = s.param.init[1], b = s.param.init[2]),
                    trace = TRUE,
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fitPoisonExp.summ <- summary(fitPoisonExp)
propRemaining.cont <- exp(lpropRemaining.cont)
hrsPigPoisonExp.pred <- coef(fitPoisonExp)[1] * exp(-coef(fitPoisonExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigPoisonExp.pred,lty=3,lwd=2,col="blue")

#export models and points for plotting 
poison_other_mods <-data.frame(propRemaining.cont, hrsPigPoison.pred, hrsPigPoisonExp.pred)
colnames(poison_other_mods) <- c("pRemain", "logfit", "expfit")
write.csv(poison_other_mods, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/poison_other_mods.csv")
poison_other_ef <-data.frame(poison_other$propRemaining,poison_other$hrsPig)
colnames(poison_other_ef)<- c("pRemain", "effort")
write.csv(poison_other_ef, "C:/Users/hamn0008/OneDrive - Flinders/Documents/Honours/Pigs_test/data/poison_other_ef.csv")

