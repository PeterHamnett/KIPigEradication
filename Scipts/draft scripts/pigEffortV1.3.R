rm(list = ls())

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

library(readr)
library(dplyr)
pigs_effort_no_LU <- read_csv("data/Final/pigeffort.csv")
View(pigs_effort_no_LU)
effort<-data.frame(pigs_effort_no_LU)
head(effort)
effort.orig <- effort # back up version of effort to revert to just in case

# create a copy of effort in case this goes wrong
effort_dates <- effort

# convert Date column from chr to date format
effort_dates$Date <- as.Date(effort$Date, format = "%d/%m/%y")

# it worked so overwrite effort with the new data
effort <- effort_dates

# sort by date 
effort <- effort[order(effort$Date),] 
# Oldest records are now at the top of the data frame

# Total number of pigs killed (where effort was recorded)= 
sum(effort$numkilled)
# 740

# Create column showing proportion of pigs remaining after each event in chronological order, assuming starting proportion is 1 and final proportion is 0
effort$propRemaining <- (sum(effort$numkilled) - cumsum(effort$numkilled))/sum(effort$numkilled)
# cumsum is the sum of the current value and all previous values in the column
# so here we have (total pigs killed - pigs killed to date)/total pigs killed
# give proportions in declining in order from just short of 1 (1 is before any pigs were killed) down to 0
head(effort)
tail(effort)

# some of the rows contain na for effort, we want to remove these rows for further analysis of effort
effort <- na.omit(effort)

## remove entries where efforthrs = 100 to see if this has any impact on shooting response curve
# effort100 <- effort [effort$efforthrs <100,]
# effort <- effort100 # sothis new subset can be used without changing any of the subsequent code
# Changed response curve, but prefer to keep these entries in as they appear to reflect shooting effort over a 

# create a new column with the inverse of pigs/hr i.e., hrs per pig
## this will allow us to quantify changes in cost relative to proportion of population remaining, once we have calculated cost per hour for each control type.

effort$hrsPig <- effort$efforthrs/effort$numkilled
#Check the new column was added
head(effort)

# Create subsets of effort dataframe for each control type

poison_efrt <-effort %>% filter (controlType == "poisoned")
poison_efrt <-poison_efrt %>% filter (org != "LB") # excludes outliers - poisoning done by the landscape board
shot_efrt <-effort %>% filter (controlType == "shot")
shot_efrt <- shot_efrt %>% filter (operator != "PJ") # excludes shooting by BF
trap_efrt <-effort %>% filter (controlType == "trapped")
trap_efrt <-trap_efrt %>% filter (org != "DEW") # excludes DEW outliers
TAAC_efrt <-effort %>% filter (controlType == "TAAC") # only conducted by HeliSurveys so no further filtering 

# produce histograms of hrs_pig for all control types to review distribution of values
par(mfrow=c(2,2))
hist(poison_efrt$hrsPig)
hist(shot_efrt$hrsPig)
hist(trap_efrt$hrsPig)
hist(TAAC_efrt$hrsPig)

# now we can plot pigs_hr and hrs-pig against proportion remaining for each control type to see if there is a relationship between the two.
# we would expect efficiency to decrease as numbers decline
par(mfrow=c(2,2))
plot(hrsPig~propRemaining, data = poison_efrt, xlim=rev(range(propRemaining)), main = "Poison CPUE", col = factor(poison_efrt$org))
legend ("topleft", legend = levels (factor(poison_efrt$org)), pch = 19, col = factor(levels(factor(poison_efrt$org))))
plot(hrsPig~propRemaining, data = trap_efrt, xlim=rev(range(propRemaining)), main = "Trapping CPUE", col = factor(trap_efrt$org))
legend ("topright", legend = levels (factor(trap_efrt$org)), pch = 19, col = factor(levels(factor(trap_efrt$org))))
plot(hrsPig~propRemaining, data = shot_efrt, xlim=rev(range(propRemaining)), main = "Shooting CPUE", col = factor(shot_efrt$operator))
legend ("topright", legend = levels (factor(shot_efrt$operator)), pch = 19, col = factor(levels(factor(shot_efrt$operator))))
plot(hrsPig~propRemaining, data = TAAC_efrt, xlim=rev(range(propRemaining)), main = "TAAC CPUE" )

#lets try using ggplot2 and adding a smoothed trend line
#pigs_hr
# library(ggplot2)
# ggplot(poison_efrt, aes(x = propRemaining,y = pigs_hr)) + scale_x_reverse() + geom_point() + geom_smooth()  + labs(title = "Poison CPUE")
# ggplot(trap_efrt, aes(x = propRemaining,y = pigs_hr)) + scale_x_reverse() + geom_point() + geom_smooth() + labs(title = "Trapping CPUE")
# ggplot(shot_efrt, aes(x = propRemaining,y = pigs_hr)) + scale_x_reverse() + geom_point() + geom_smooth() + labs(title = "Shooting CPUE")
# ggplot(TAAC_efrt, aes(x = propRemaining,y = pigs_hr)) + scale_x_reverse() + geom_point() + geom_smooth() + labs(title = "TAAC CPUE")
# 
# # and  hrs_pig
# ggplot(poison_efrt, aes(x = propRemaining,y = hrsPig)) + scale_x_reverse() + geom_point() + geom_smooth()  + labs(title = "Poison: Effort per pig")
# ggplot(trap_efrt, aes(x = propRemaining,y = hrsPig)) + scale_x_reverse() + geom_point() + geom_smooth() + labs(title = "Trapping: Effort per pig")
# ggplot(shot_efrt, aes(x = propRemaining,y = hrsPig)) + scale_x_reverse() + geom_point() + geom_smooth() + labs(title = "Shooting: Effort per pig")
# ggplot(TAAC_efrt, aes(x = propRemaining,y = hrsPig)) + scale_x_reverse() + geom_point() + geom_smooth() + labs(title = "TAAC: Effort per pig")
# 
# # the above plots don't show a clear linear, logarithmic or exponential relationship

par(mfrow=c(2,2))
## Poison
plot(hrsPig~propRemaining, data = poison_efrt, pch=19, cex=0.7, main = "Poison CPUE", xlim=c(0,1), ylim=c(0,75), col = factor(poison_efrt$org))
abline(h=0, lty=3, lwd=0.5)
legend ("topleft", legend = levels (factor(poison_efrt$org)), pch = 19, col = factor(levels(factor(poison_efrt$org))))

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
lines(propRemaining.cont,hrsPigPoisonExp.pred,lty=3,lwd=2,col="green")



## Trap
plot(hrsPig~propRemaining, data = trap_efrt, pch=19, cex=0.7, main = "Trapping CPUE", xlim=c(0,1), ylim=c(0,40), col = factor(trap_efrt$org))
abline(h=0, lty=3, lwd=0.5)
legend ("topleft", legend = levels (factor(trap_efrt$org)), pch = 19, col = factor(levels(factor(trap_efrt$org))))
# logarithmic fit
# remove zeros
trap_efrtNZ <- trap_efrt[-(which(trap_efrt$propRemaining==0)), ]
fitTrap <- lm(hrsPig ~ log(propRemaining), data=trap_efrtNZ, na.action = ) # logarithmic fit
summary(fitTrap)
linreg.ER(log(trap_efrtNZ$propRemaining),trap_efrtNZ$hrsPig)
hrsPigTrap.pred <- coef(fitTrap)[1] + coef(fitTrap)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigTrap.pred, lty=2, col="red")

# exponential: y = a * exp(-b*x)
s.param.init <- c(10, 2)
fitTrapExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                    data = trap_efrt,
                    algorithm = "port",
                    start = c(a = s.param.init[1], b = s.param.init[2]),
                    trace = TRUE,      
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fitTrapExp.summ <- summary(fitTrapExp)
hrsPigTrapExp.pred <- coef(fitTrapExp)[1] * exp(-coef(fitTrapExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigTrapExp.pred,lty=3,lwd=2,col="green")


## Shot
plot(hrsPig~propRemaining, data = shot_efrt, pch=19, cex=0.7, main = "Shot CPUE", xlim=c(0,1), ylim=c(0,110), col = factor(shot_efrt$operator))
abline(h=0, lty=3, lwd=0.5)
legend ("topleft", legend = levels (factor(shot_efrt$operator)), pch = 19, col = factor(levels(factor(shot_efrt$operator))))
# logarithmic fit
fitShot <- lm(hrsPig ~ log(propRemaining), data=shot_efrt) # logarithmic fit
summary(fitShot)
linreg.ER(log(shot_efrt$propRemaining),shot_efrt$hrsPig)
hrsPigShot.pred <- coef(fitShot)[1] + coef(fitShot)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigShot.pred, lty=2, col="red")

# exponential: y = a * exp(-b*x)
s.param.init <- c(200, 5)
fitShotExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                  data = shot_efrt,
                  algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,      
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fitShotExp.summ <- summary(fitShotExp)
hrsPigShotExp.pred <- coef(fitShotExp)[1] * exp(-coef(fitShotExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigShotExp.pred,lty=3,lwd=2,col="green")


## TAAC
plot(hrsPig~propRemaining, data = TAAC_efrt, pch=19, cex=0.7, main = "TAAC CPUE", xlim=c(0,1), ylim=c(0,5))
abline(h=0, lty=3, lwd=0.5)

# remove infinite
TAAC_efrtNI <- TAAC_efrt[which(is.infinite(TAAC_efrt$hrsPig)==F),]
# logarithmic fit
fitTAAC <- lm(hrsPig ~ log(propRemaining), data=TAAC_efrtNI) # logarithmic fit
summary(fitTAAC)
linreg.ER(log(TAAC_efrtNI$propRemaining),TAAC_efrtNI$hrsPig)
hrsPigTAAC.pred <- coef(fitTAAC)[1] + coef(fitTAAC)[2]*lpropRemaining.cont
lines(exp(lpropRemaining.cont), hrsPigTAAC.pred, lty=2, col="red")

# exponential: y = a * exp(-b*x)
s.param.init <- c(5, 1)
fitTAACExp <- nls(hrsPig ~ a * exp(-b*propRemaining), 
                  data = TAAC_efrtNI,
                  algorithm = "port",
                  start = c(a = s.param.init[1], b = s.param.init[2]),
                  trace = TRUE,      
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fitTAACExp.summ <- summary(fitTAACExp)
hrsPigTAACExp.pred <- coef(fitTAACExp)[1] * exp(-coef(fitTAACExp)[2] * propRemaining.cont)
lines(propRemaining.cont,hrsPigTAACExp.pred,lty=3,lwd=2,col="green")
par(mfrow=c(1,1))

