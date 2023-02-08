# Peter Hamnett, Kathryn Venning, Frédérik Saltré and Corey Bradshaw
# Global Ecology, Flinders University — globalecologyflinders.com
# feral pig (sus scrofa) eradication on Kangaroo Island
# https://github.com/PeterHamnett/KIPigEradication
### update 8/2/2023
# base stochastic model modified from https://github.com/KathrynVenning/FeralCatEradication

## this script presents outputs presented in Appendix 4 of the manuscript ##

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
setwd("/Users/peterhamnett/Documents - Peter’s MacBook Pro/R Resources/Projects /Pigs_test")
source("~/Documents - Peter’s MacBook Pro/R Resources/Projects /Pigs_test/Scipts/matrixOperators.r")



#### Functional Response Estimates ####

# as above, change file path to local file location
pigs_effort_no_LU <- read_csv("~/Documents - Peter’s MacBook Pro/R Resources/Projects /Pigs_test/data/final/pigeffort3.csv")
effort <- data.frame(pigs_effort_no_LU)
head(effort)
effort.orig <- effort # back up version of effort to revert to just in case

# create a copy of effort in case this goes wrong
effort_dates <- effort

# convert Date column from chr to date format
effort_dates$Date <- as.Date(effort$Date, format = "%d/%m/%y")

# it worked so overwrite effort with the new data
effort <- effort_dates

# sort by date so oldest records are at the top of the data frame
effort <- effort[order(effort$Date),] 

# remove records where effort is na
effort <- na.omit(effort)

# create a new column with the inverse of pigs/hr i.e., hrs per pig
## this will allow us to quantify changes in cost relative to proportion of population remaining, once we have calculated cost per hour for each control type.
effort$hrsPig <- effort$efforthrs/effort$numkilled

#Check the new column was added
head(effort)

#create continuous variables for model fitting
lpropRemaining.cont <- log(seq(0.01, 1, 0.01))
propRemaining.cont <- exp(lpropRemaining.cont)


# create data subsets for each control technique
shot_efrt <- effort %>% filter (controlType == "shot")
shot_efrt <- shot_efrt %>% filter (operator != "PJ") # excludes shooting by PJ 8 out of 82 shooting events (all events recorded as 100hrs - imprecise and invalid for fitting functional response model)
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
plot(hrsPig~KProp, data = shot_efrt, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,110))
mtext("a", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
  # logarithmic fit
    fitShot <- lm(hrsPig ~ log(KProp), data=shot_efrt, na.action =) # logarithmic fit
    summary(fitShot)
    linreg.ER(log(shot_efrt$KProp),shot_efrt$hrsPig)
    hrsPigShot.pred <- coef(fitShot)[1] + coef(fitShot)[2]*lpropRemaining.cont
    lines(exp(lpropRemaining.cont), hrsPigShot.pred, lty=2, col="red")
  # exponential: y = a * exp(-b*x)
    s.param.init <- c(200, 5)
    fitShotExp <- nls(hrsPig ~ a * exp(-b*KProp), 
                      data = shot_efrt,
                      #algorithm = "port",
                      start = c(a = s.param.init[1], b = s.param.init[2]),
                      trace = TRUE,      
                      nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
    fitShotExp.summ <- summary(fitShotExp)
    fitShotExp.summ
    hrsPigShotExp.pred <- coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * propRemaining.cont)
    lines(propRemaining.cont,hrsPigShotExp.pred,lty=3,lwd=2,col="blue")

## TAAC
plot(hrsPig~KProp, data = TAAC_efrt, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,5))
mtext("b", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
  # logarithmic fit
    fitTAAC <- lm(hrsPig ~ log(KProp), data=TAAC_efrt) # logarithmic fit
    summary(fitTAAC)
    linreg.ER(log(TAAC_efrt$KProp),TAAC_efrt$hrsPig)
    hrsPigTAAC.pred <- coef(fitTAAC)[1] + coef(fitTAAC)[2]*lpropRemaining.cont
    lines(exp(lpropRemaining.cont), hrsPigTAAC.pred, lty=2, col="red")
  # exponential: y = a * exp(-b*x)
    s.param.init <- c(5, 1)
    fitTAACExp <- nls(hrsPig ~ a * exp(-b*KProp), 
                      data = TAAC_efrt,
                      algorithm = "port",
                      start = c(a = s.param.init[1], b = s.param.init[2]),
                      trace = TRUE,      
                      nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
    fitTAACExp.summ <- summary(fitTAACExp)
    hrsPigTAACExp.pred <- coef(fitTAACExp)[1] * exp(-coef(fitTAACExp)[2] * propRemaining.cont)
    lines(propRemaining.cont,hrsPigTAACExp.pred,lty=3,lwd=2,col="blue")

#trapped
plot(hrsPig~KProp, data = trap_efrt, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,40))
mtext("c", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
  # logarithmic fit
  fitTrap <- lm(hrsPig ~ log(KProp), data=trap_efrt, na.action = ) # logarithmic fit
  summary(fitTrap)
  linreg.ER(log(trap_efrt$KProp),trap_efrt$hrsPig)
  hrsPigTrap.pred <- coef(fitTrap)[1] + coef(fitTrap)[2]*lpropRemaining.cont
  lines(exp(lpropRemaining.cont), hrsPigTrap.pred, lty=2, col="red")
  # exponential: y = a * exp(-b*x)
  s.param.init <- c(12, 2)
  fitTrapExp <- nls(hrsPig ~ a * exp(-b*KProp), 
                    data = trap_efrt,
                    algorithm = "port",
                    start = c(a = s.param.init[1], b = s.param.init[2]),
                    trace = TRUE,      
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
  fitTrapExp.summ <- summary(fitTrapExp)
  hrsPigTrapExp.pred <- coef(fitTrapExp)[1] * exp(-coef(fitTrapExp)[2] * propRemaining.cont)
  lines(propRemaining.cont,hrsPigTrapExp.pred,lty=3,lwd=2,col="blue")

#poisoned
plot(hrsPig~KProp, data = poison_efrt, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,75))
mtext("d", side = 3, line =-1.5, adj = 0.05, cex = 1)
abline(h=0, lty=3, lwd=0.5)
  # logarithmic fit: y = a + b*log(x)
    fitPoison <- lm(hrsPig ~ log(KProp), data=poison_efrt) # logarithmic fit
    summary(fitPoison)
    linreg.ER(log(poison_efrt$KProp),poison_efrt$hrsPig)
    lpropRemaining.cont <- log(seq(0.01, 1, 0.01))
    hrsPigPoison.pred <- coef(fitPoison)[1] + coef(fitPoison)[2]*lpropRemaining.cont
    lines(exp(lpropRemaining.cont), hrsPigPoison.pred, lty=2, col="red")
  # exponential: y = a * exp(-b*x)
    s.param.init <- c(200, 5)
    fitPoisonExp <- nls(hrsPig ~ a * exp(-b*KProp), 
                        data = poison_efrt,
                        algorithm = "port",
                        start = c(a = s.param.init[1], b = s.param.init[2]),
                        trace = TRUE,      
                        nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
    fitPoisonExp.summ <- summary(fitPoisonExp)
    propRemaining.cont <- exp(lpropRemaining.cont)
    hrsPigPoisonExp.pred <- coef(fitPoisonExp)[1] * exp(-coef(fitPoisonExp)[2] * propRemaining.cont)
    lines(propRemaining.cont,hrsPigPoisonExp.pred,lty=3,lwd=2,col="blue")

# used in final presentation, but not in thesis. Will resurrect for manuscript
  # all control types combined for overlaying on a single plot
  # comb_efrt <-rbind(shot_efrt,TAAC_efrt,trap_efrt,poison_efrt)
  # plot with all 4 log models presented together in same frame
  # plot(hrsPig~propRemaining, data = comb_efrt, xlab = '', ylab = '', pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,35))
  # abline(h=0, lty=3, lwd=0.5)
  # lines(exp(lpropRemaining.cont), hrsPigShot.pred, lty=1, col="red")
  # lines(exp(lpropRemaining.cont), hrsPigTAAC.pred, lty=2, col="red")
  # lines(exp(lpropRemaining.cont), hrsPigTrap.pred, lty=3, col="red")
  # lines(exp(lpropRemaining.cont), hrsPigPoison.pred, lty=4, col="red")

# linear and intercept only model fits for comparison
# linear fit
  fitShotLin <- lm(hrsPig ~ KProp, data=shot_efrt) # fits the Linear model
  summary(fitShotLin) # provides summary statistics of the LM for AIC comparison, not plotted
  
  fitTAACLin <- lm(hrsPig ~ KProp, data=TAAC_efrt) # fits the Linear model
  summary(fitTAACLin) # provides summary statistics of the LM for AIC comparison, not plotted
  
  fitTrapLin <- lm(hrsPig ~ KProp, data=trap_efrt) # fits the Linear model
  summary(fitTrapLin) # provides summary statistics of the LM for AIC comparison, not plotted
  
  fitPoisonLin <- lm(hrsPig ~ KProp, data=poison_efrt) # fits the Linear model
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

age.max = 6 # Choquenot, D (1996) says few pigs live longer than 5 in Australia, but much greater longevity reported both in Australia and elsewhere

## create fertility and survival vectors 
#### Fertility ####
## fertility values for each age class derived from Bieber and Ruf (2005)
f.vec <- c(0.8, 2.3375, 2.925, 2.925, 2.925, 2.835)
## alternative f.vec <-c(0.79, 2.38, 2.38, 2.38, 2.38, 2.38) values calculated from Choquenot (1996), but frequency distribution not known for calculating SD

## visualise change in mean fertility as a function of age
plot(0:5,f.vec, pch=19, type="b")

# fertility errors based on Bieber and Ruf (2005). 
## SD calculated as the mean(c(((intermediate value - poor  value)/2),((good value - intermediate value)/2)))
J.f.sd <- mean(c(((0.8 - 0.525)/2),((1.125 - 0.8)/2))) #mean and standard deviations, juvenile fertility
Y.f.sd <- mean(c(((1.8 - 1.625)/2) ,((2.3375 - 1.8)/2))) #mean and standard deviations, yearling fertility
A.f.sd <- mean(c(((2.835 - 1.7)/2),((2.925 - 1.7)/2))) #mean and standard deviations, adult fertility
f.sd.vec <- c(J.f.sd, Y.f.sd, A.f.sd, A.f.sd, A.f.sd, A.f.sd) #mean and standard deviations vector, juvenile and adult fertility 

#### Survival ####
s.vec <- c(0.33, 0.40, 0.66, 0.66, 0.66, 0.58) ##feral pig survival # intermediate values from Bieber and Ruf (2005), except 6 which is the poor value.
## s.vec <- c(0.45, 0.675, 0.675, 0.675, 0.675, 0.675) ## alternative values from Choquenot (1996) 

## SD calculated as the mean(c(((mid range value - low range value)/2),((high range value - mid range value)/2)))
# survival errors based on Bieber and Ruf (2005). 
J.s.sd <- mean(c(((0.33 - 0.25)/2),((0.52 - 0.33)/2))) #mean and standard deviations, juvenile survival
Y.s.sd <- mean(c(((0.40 - 0.31)/2),((0.60 - 0.40)/2))) #mean and standard deviations, yearling survival
A.s.sd <- mean(c(((0.66 -0.58)/2),((0.71 -0.66)/2))) #mean and standard deviations, adult survival
s.sd.vec <- c(J.s.sd, Y.s.sd, A.s.sd, A.s.sd, A.s.sd, A.s.sd) #mean and standard deviations vector, juvenile and adult survival

## visualise change in mean survival as a function of age
plot(0:5,s.vec,pch=19,type="b")

# create matrix
popmat <- matrix(data = 0, nrow=age.max, ncol=age.max) ##creates a matrix where data = 0 and dimensions age.max x age.max; 6 x 6
diag(popmat[2:age.max,]) <- s.vec[1:5] ## diagonally in popmat from row 2, col 1 to row 6, col 5 populated with s.vec
popmat[age.max,age.max] <- s.vec[6] ## position [6,6] = 0
popmat[1,] <- f.vec ## row 1 of popmat populated with f.vec
popmat.orig <- popmat ## save original matrix

## matrix properties. Functions from "matrixOperators.R"
max.lambda(popmat) ## 1-yr lambda. 
max.r(popmat) # rate of population change, 1-yr
ssd <- stable.stage.dist(popmat) ## stable stage distribution
plot(0:5,ssd,pch=19,type="b")
R.val(popmat, age.max) # reproductive value
gen.l <- G.val(popmat, age.max) # mean generation length

## initial population vector
pop.found <- 5000/2 # +/- 50 founding population size personal communication, B. Page, PIRSA,  2021 . 
## Divided by 2 for females only
ssd <- stable.stage.dist(popmat) ## sum of all values in ssd = 1
init.vec <- ssd * pop.found #initial population vector. Sum of all values here = 500 as suggested by PIRSA
## (can change to whatever we want our founding pop to be)


#### Deterministic population projection ####
# i.e., project population growth from the leslie matrix without incorporating stochasticity
## set time limit for projection in 1-yr increments
yr.now <- 2020 
#************************

yr.end <- 2023 #end year for our projection time frame
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

K.s.vec <- c(1,K.max*0.4,K.max*0.5,0.7*K.max,0.85*K.max,0.95*K.max) ##describes the x axis of the reduction curve
red.s.vec <- c(1,0.995,0.98,0.96,0.92,0.85) ## describes the y axis of the reduction curve
plot(K.s.vec,red.s.vec,pch=19,type="b")
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

#### compensatory density feedback: FERTILITY ####
K.f.vec <- c(1,K.max*0.6,0.85*K.max,0.92*K.max,0.95*K.max) ## describes the x axis of the reduction curve
red.f.vec <- c(1,0.993,0.987,0.985,0.984)

# K.f.vec <- c(1,K.max*0.6,0.85*K.max,0.92*K.max,0.95*K.max) ## describes the x axis of the reduction curve
# red.f.vec <- c(1,0.995,0.988,0.985,0.9835)
plot(K.f.vec,red.f.vec,pch=19,type="b")
Kred.f.dat <- data.frame(K.f.vec,red.f.vec)

# logistic power function a/(1+(x/b)^c) #fits logistic power function to population relative to carry capacity, K
f.param.init <- c(1, K.max, 3) ## These values are arbitrary 
fit.f.lp <- nls(red.f.vec ~ a/(1+(K.f.vec/b)^c), 
                data = Kred.f.dat,
                algorithm = "port",
                start = c(a = f.param.init[1], b = f.param.init[2], c = f.param.init[3]),
                trace = TRUE,      
                nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fit.f.lp.summ <- summary(fit.f.lp)
plot(K.f.vec,red.f.vec,pch=19,xlab="N",ylab="f reduction factor")
K.f.vec.cont <- seq(1,12*pop.found,1)
pig.f.lp.fx <- coef(fit.f.lp)[1]/(1+(K.f.vec.cont/coef(fit.f.lp)[2])^coef(fit.f.lp)[3])
lines(K.f.vec.cont,pig.f.lp.fx,lty=2,lwd=3,col="red")

f.a.lp <- coef(fit.f.lp)[1]
f.b.lp <- coef(fit.f.lp)[2]
f.c.lp <- coef(fit.f.lp)[3]
## compensatory density-feedback deterministic model
## set population storage matrices
n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
n.mat[,1] <- init.vec
popmat <- popmat.orig


## set up projection: deterministic pop with both S and F density feedback
for (i in 1:t) {
  
  totN.i <- sum(n.mat[,i])
  pigs.s.red <- s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp)
  pigs.f.red <- f.a.lp/(1+(totN.i/f.b.lp)^f.c.lp)
  diag(popmat[2:age.max,]) <- s.vec[1:5]*pigs.s.red
  popmat[age.max,age.max] <- s.vec[6]*pigs.s.red
  popmat[1,1:6]<- f.vec*pigs.f.red
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pigs <- colSums(n.mat)
plot(yrs, n.pigs,type="l",main = "deterministic population projection", sub = " S and F density feedback", lty=1,pch=19,xlab="year",ylab="N",ylim=c(0,1.8*K.max)) #untreated population increases, rate of increase relative to K, no stochastic sampling
abline(h=K.max, lty=2, col="red") #carry capacity

#### iterations ####
iter <- 1000 #final model run at 10 000
itdiv <- iter/100 #final model rate at iter/1000


#### untreated population projection ####

## stochastic projection with survival and fertility density feedback
## set storage matrices & vectors
# set t to 10 years so projection shows response relative to K

yr.now <- 2020 
#************************

yr.end <- 2030 #end year for our projection timeframe
yrs <- seq(yr.now, yr.end, 1)
#************************
t <- (yr.end - yr.now)
n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1)) #storage matrix with 1000 rows and 11 columns based on 10 year timeframe

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
    fert.stch <- rnorm(length(popmat[,1]), popmat[1,], f.sd.vec)
    fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
    
    totN.i <- sum(n.mat[,i])
    pigs.s.red <- as.numeric(s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp))
    pigs.f.red <- as.numeric(f.a.lp/(1+(totN.i/f.b.lp)^f.c.lp))
    #a/(1+(K.vec/b)^c)
    
    popmat[1,] <- fert.stoch*pigs.f.red
    diag(popmat[2:age.max,]) <- s.stoch[1:5]*pigs.s.red
    popmat[age.max,age.max] <- s.stoch[6]*pigs.s.red
    #popmat[age.max,age.max] <- 0
    
    n.mat[,i+1] <- popmat %*% n.mat[,i]
    
  } # end i loop
  
  n.sums.mat[e,] <- ((as.vector(colSums(n.mat))))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# produce dataframe recording median and CIs for pop size in each year
datN <- data.frame(yrs, n.md, n.lo, n.up)
tail(datN)

## plot pop change in ggplot2 (thesis figure 2 - results)
ggplot(data = datN, mapping = aes(x=yrs)) +
  geom_line(aes(y=n.md), color = "black") + 
  geom_line(aes(y=n.lo),color = "red", linetype = 2) + 
  geom_line(aes(y=n.up),color = "red", linetype = 2) + 
  theme_bw() +
  geom_hline(yintercept = 2500, linetype = 3) +
  labs(x = "Years", y = "Population (females)")


## untreated pop. projection with S density reduction only (no F) for inclusion in supplementary material 
## demonstrates unrestrained at upper confidence interval if density dependent feedback not applied to fertility

n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1)) #storage matrix with 1000 rows and 11 columns based on 10 year timeframe

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
    fert.stch <- rnorm(length(popmat[,1]), popmat[1,], f.sd.vec)
    fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
    
    totN.i <- sum(n.mat[,i])
    pigs.s.red <- as.numeric(s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp))
    #a/(1+(K.vec/b)^c)
    
    popmat[1,] <- fert.stoch
    diag(popmat[2:age.max,]) <- s.stoch[1:5]*pigs.s.red
    popmat[age.max,age.max] <- s.stoch[6]*pigs.s.red
    #popmat[age.max,age.max] <- 0
    
    n.mat[,i+1] <- popmat %*% n.mat[,i]
    
  } # end i loop
  
  n.sums.mat[e,] <- ((as.vector(colSums(n.mat))))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# produce dataframe recording median and CIs for pop size in each year
datN <- data.frame(yrs, n.md, n.lo, n.up)
tail(datN)

## plot pop change in ggplot2 (thesis Figure S1 - Appendix 2)
ggplot(data = datN, mapping = aes(x=yrs)) +
  geom_line(aes(y=n.md), color = "black") + 
  geom_line(aes(y=n.lo),color = "red", linetype = 2) + 
  geom_line(aes(y=n.up),color = "red", linetype = 2) + 
  theme_bw() +
  geom_hline(yintercept = 2500, linetype = 3) +
  labs(x = "Years", y = "Population (females)")

### Additional analysis investigating cull rates and costs for eradication in a population at carrying capacity

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
      fert.stch <- rnorm(length(popmat[,1]), popmat[1,], f.sd.vec/2)
      fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
      
      totN.i <- sum(n.mat[,i])
      pigs.s.red <- as.numeric(s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp))
      pigs.f.red <- as.numeric(f.a.lp/(1+(totN.i/f.b.lp)^f.c.lp))
      #a/(1+(K.vec/b)^c)
      
      popmat[1,] <- fert.stoch*pigs.f.red
      diag(popmat[2:age.max,]) <- s.stoch[1:5]*pigs.s.red
      popmat[age.max,age.max] <- s.stoch[6]*pigs.s.red
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
  
  
  
  # plot(yrs,n.md,type="l",xlab="year", ylab="minimum N", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
  # lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
  # lines(yrs,n.up,lty=2,col="red",lwd=1.5)
  #
  # # print("##############")
  print(paste("harvest proportion = ", harv.prop.consist[s], sep=""))
  # print("##############")
  
} # ends S loop

minn.prop.pop3 <- data.frame(harv.prop.consist, min.med.n3, min.lo.n3, min.up.n3)
minn.prop.pop3$med.pop <- minn.prop.pop3$min.med.n3*pop.found

## plot in ggplot (thesis figure 3)
harv3plot <- ggplot(data = minn.prop.pop3, mapping = aes(x=harv.prop.consist)) +
  geom_line(aes(y=min.med.n3), color = "black") + 
  geom_line(aes(y=min.lo.n3), color = "red", linetype = 2) + 
  geom_line(aes(y=min.up.n3), color = "red", linetype = 2) + 
  geom_vline(xintercept = 0.9, linetype = 3, color = "black") +
  theme_bw() +
  labs(x = "Harvest Proportion", y = "Min N")

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
      fert.stch <- rnorm(length(popmat[,1]), popmat[1,], f.sd.vec/2)
      fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
      
      totN.i <- sum(n.mat[,i])
      pigs.s.red <- as.numeric(s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp))
      pigs.f.red <- as.numeric(f.a.lp/(1+(totN.i/f.b.lp)^f.c.lp))
      #a/(1+(K.vec/b)^c)
      
      popmat[1,] <- fert.stoch*pigs.f.red
      diag(popmat[2:age.max,]) <- s.stoch[1:5]*pigs.s.red
      popmat[age.max,age.max] <- s.stoch[6]*pigs.s.red
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
minn.prop.pop$med.pop <- minn.prop.pop$min.med.n*pop.found

#plot the lines individually
harv10plot <- ggplot(data = minn.prop.pop, mapping = aes(x=harv.prop.consist)) +
  geom_line(aes(y=min.med.n), color = "black") + 
  geom_line(aes(y=min.lo.n), color = "red", linetype = 2) + 
  geom_line(aes(y=min.up.n), color = "red", linetype = 2) + 
  geom_vline(xintercept = 0.6, linetype = 3, color = "black") +
  theme_bw() +
  labs(x = "Harvest Proportion", y = "Min N")

# multipanel of constant prop harvest over 3 and 10 years (Thesis figure 3 - results)
harvplot <- ggarrange(harv3plot, harv10plot, common.legend = TRUE, labels = c("a", "b"), legend = "bottom", nrow = 1, ncol = 2)
harvplot


#### CONSTANT PROPORTIONAL CULL COSTS####
# fixed cost components (not recalculated each loop)
labour.ph <- 36.87 ## $36.87 per hour, based on SA Public Sector Award OPS4 classification. 
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
trap.initial <- labour.ph*20 # additional unreported effort of ~20 hrs per trap for trap construction
trap.unit <- 9500.00 # cost per unit for the Jager Pro trap
trap.ff<-freefeed.ph*mean(trap_efrt$efforthrs) # cost of grain per event
trap.ff.labour <-labour.ph*mean(trap_efrt$efforthrs) # feeding labour cost per event
trap.ff.cost<-trap.ff.labour + trap.ff # combined cost of grain and labour per event

# poison specific costs
poison.placebo <- 224 # 6 placebo baits x 4 days
poison.toxic <- 156 # 6 toxic baits x 1 day
poison.dispenser.unit <- 485
avg.poison.effort <- mean(poison_efrt$efforthrs)
poison.labour <- labour.ph*avg.poison.effort
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
        fert.stch <- rnorm(length(popmat[,1]), popmat[1,], f.sd.vec/2)
        fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
        
        totN.i <- sum(n.mat[,i])
        pigs.s.red <- as.numeric(s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp))
        pigs.f.red <- as.numeric(f.a.lp/(1+(totN.i/f.b.lp)^f.c.lp))
        #a/(1+(K.vec/b)^c)
        
        popmat[1,] <- fert.stoch*pigs.f.red
        diag(popmat[2:age.max,]) <- s.stoch[1:5]*pigs.s.red
        popmat[age.max,age.max] <- s.stoch[6]*pigs.s.red
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
            ## not initial pop in projection(250 - females only)
        remaining.preharv <- sum(n.mat[,i])
        prop.remain <- remaining.preharv/K.max
        
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
                labour.cost <- (labour.ph  + freefeed.ph.shot) * tot.shot.hrs
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
        labour.cost <- (labour.ph  + freefeed.ph.shot) * tot.shot.hrs[i]
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

scenariosNK.df <- as.data.frame(scenarios.mat[-(1:3),])
names(scenariosNK.df) <- c("scenario", "harvest.prop", "year","low.cost", "median.cost", "upper.cost", 
                         "low.time", "median.time", "upper.time", "median.N") 
rownames(scenariosNK.df)<-NULL
scenariosNK.df[ ,-1] <-lapply(scenariosNK.df[ ,-1], as.numeric)
View(scenariosNK.df)

scenariosNK.tot.df<-as.data.frame(scenarios.tot[-1,])
names(scenariosNK.tot.df) <- c("scenario", "harvest.prop", "tot.cost.lo", "tot.cost.med", "tot.cost.up", 
                             "tot.time.low", "tot.time.med", "tot.time.up", "min.n.lower", "min.n.median", "min.n.upper")
rownames(scenariosNK.tot.df)<-NULL
scenariosNK.tot.df[ ,-1] <-lapply(scenariosNK.tot.df[ ,-1], as.numeric)
View(scenariosNK.tot.df)

##### repeat cull scenarios projection for 10 years instead of 3 to compare costs and effort####
### Ideally this should be changed so that time frames t = 3 and t = 10 are run in loops

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

# names for each culling scenario to be used later when plotting outcomes
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
        fert.stch <- rnorm(length(popmat[,1]), popmat[1,], f.sd.vec/2)
        fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
        
        totN.i <- sum(n.mat[,i])
        pigs.s.red <- as.numeric(s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp))
        pigs.f.red <- as.numeric(f.a.lp/(1+(totN.i/f.b.lp)^f.c.lp))
        #a/(1+(K.vec/b)^c)
        
        popmat[1,] <- fert.stoch*pigs.f.red
        diag(popmat[2:age.max,]) <- s.stoch[1:5]*pigs.s.red
        popmat[age.max,age.max] <- s.stoch[6]*pigs.s.red
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
        prop.remain <- remaining.preharv/K.max
        
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
          labour.cost <- (labour.ph  + freefeed.ph.shot) * tot.shot.hrs
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
        labour.cost <- (labour.ph  + freefeed.ph.shot) * tot.shot.hrs[i]
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

scenarios10NK.df <- as.data.frame(scenarios.mat[-(1:10),])
names(scenarios10NK.df) <- c("scenario", "harvest.prop", "year","low.cost", "median.cost", "upper.cost", 
                           "low.time", "median.time", "upper.time", "median.N") 
rownames(scenarios10NK.df)<-NULL
scenarios10NK.df[ ,-1] <-lapply(scenarios10NK.df[ ,-1], as.numeric)
View(scenarios10NK.df)

scenarios10NK.tot.df<-as.data.frame(scenarios.tot[-1,])
names(scenarios10NK.tot.df) <- c("scenario", "harvest.prop", "tot.cost.lo", "tot.cost.med", "tot.cost.up", 
                               "tot.time.low", "tot.time.med", "tot.time.up", "min.n.lower", "min.n.median", "min.n.upper")
rownames(scenarios10NK.tot.df)<-NULL
scenarios10NK.tot.df[ ,-1] <-lapply(scenarios10NK.tot.df[ ,-1], as.numeric)
View(scenarios10NK.tot.df)

write_csv(scenariosNK.df, "/Users/peterhamnett/Documents - Peter’s MacBook Pro/R Resources/Projects /Pigs_test\\scenariosNK.df.csv")
write_csv(scenariosNK.tot.df, "/Users/peterhamnett/Documents - Peter’s MacBook Pro/R Resources/Projects /Pigs_test\\scenariosNK.tot.df.csv")
write_csv(scenarios10NK.df, "/Users/peterhamnett/Documents - Peter’s MacBook Pro/R Resources/Projects /Pigs_test\\scenarios10NK.df.csv")
write_csv(scenarios10NK.tot.df, "/Users/peterhamnett/Documents - Peter’s MacBook Pro/R Resources/Projects /Pigs_test\\scenarios10NK.tot.df.csv")

# use scenarios10NK.df and scenarios10NK.tot.df to produce combined plots fof cost per year for all scenarios per cull proportion
par(mfrow=c(1,1))
# set a colourblind friendly palette
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#D55E00", "#CC79A7", "#0072B2")
library(scales)

# plot of total cost per scenario after 3 years
p3C<-ggplot(scenariosNK.tot.df, aes(x = harvest.prop, y = tot.cost.med, group=scenario, colour=scenario)) +  
  geom_line() +  theme_bw() +scale_colour_manual(values=cbp2) +
  labs(title = "", subtitle =  "", y = "", x ="") + geom_vline(xintercept = 0.9, linetype = 3, color = "black") +
  scale_y_continuous(labels = comma)

# plot of total time per scenario after 3 years
p3T<- ggplot(scenariosNK.tot.df, aes(x = harvest.prop, y = tot.time.med, group=scenario, colour=scenario)) +  
  geom_line() + theme_bw() + scale_colour_manual(values=cbp2) +
  labs(title = "", subtitle =  "", y = "", x ="") + geom_vline(xintercept = 0.9, linetype = 3, color = "black")

# plot of total cost per scenario after 10 years
p10C<-ggplot(scenarios10NK.tot.df, aes(x = harvest.prop, y = tot.cost.med, group=scenario, colour=scenario)) +  
  geom_line() + theme_bw() + scale_colour_manual(values=cbp2) +
        labs(title = "", subtitle =  "", y = "", x ="") + geom_vline(xintercept = 0.6, linetype = 3, color = "black") +
        scale_y_continuous(labels = comma)

# plot of total time per scenario after 10 years

p10T <-ggplot(scenarios10NK.tot.df, aes(x = harvest.prop, y = tot.time.med, group=scenario, colour=scenario)) +  
  geom_line() + theme_bw() + scale_colour_manual(values=cbp2) +
        labs(title = "", subtitle =  "", y = "", x ="") + geom_vline(xintercept = 0.6, linetype = 3, color = "black")

# Thesis figures 5 and 7 - results
totcostfig <- ggarrange(p3C, p10C, common.legend = TRUE, labels = c("a", "b"), legend = "bottom", nrow = 1, ncol = 2)
toteffortfig <- ggarrange(p3T, p10T, common.legend = TRUE, labels = c("a", "b"), legend = "bottom", nrow = 1, ncol = 2)

print(totcostfig)
print(toteffortfig)
