# Peter Hamnett, Kathryn Venning, Frédérik Saltré and Corey Bradshaw
# Global Ecology, Flinders University — globalecologyflinders.com
# feral pig (sus scrofa) eradication on Kangaroo Island
# https://github.com/PeterHamnett/Pigs_test
# requires library - Plotly
### update 2/12/2021
## update includes: modification of code from https://github.com/KathrynVenning/FeralCatEradication for application to sus scrofa


## remove everything
rm(list = ls())

# libraries
library(plotly)
options(scipen = 1000) ##

## functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

## source/matrix operators
source("~/Documents - Peter’s MacBook Pro/R Resources/Projects /Pigs_test/Scipts/matrixOperators.r")

# create Leslie matrix
age.max = 6 # Choquenot, D (1996)

#fertility and survival values after Bieber and Ruf (2005), Table 2.
## all intermediate values (table 2 lists good, intermediate and poor values reflecting variation in S and F according to environmental conditions)
## with the exception of last age class which is assigned the poor value to reflect senescence 
## F =  mean litter size/2*proportion of individuals reproducing rather than stated fertility values. Assumes 1:1 sex ratio

## Although Bieber and Ruf (2005) is a european study, fertility and survival values are comparable to those listed in Choquenot (1996):
    ##Adult mortality 15-50%. Probability of survival range = 0.85 -0.5, avg = 0.675
    ##Juvenile mortality as 10-15% up to 100% in bad years. probability of survival range = 0.9-0, avg = 0.45
    ##Adult fertility 0.85 litters per year. Range 4.9-5.6 piglets per litter. Avg female offspring per year  = 2.38
    ##juvenile fertility: divide adult fertility by 3 to account for proportion of females breeding in first year of life

## create vectors 

f.vec <- c(0.8, 2.3375, 2.925, 2.925, 2.925, 2.835) 
## values from Choquenot (1996) ## f.vec <-c(0.79, 2.38, 2.38, 2.38, 2.38, 2.38)

## KI feral pig birth rates matrix, data for female offspring produced each year.

plot(0:5,f.vec, pch=19, type="b")

# fertility errors based on Bieber and Ruf (2005). Unsure how to (or if I can) produce SD from Choquenot (1996)
## SD calculated as the mean(c(((intermediate value - poor  value)/2),((good value - intermediate value)/2)))
J.f.sd <- mean(c(((0.8 - 0.525)/2),((1.125 - 0.8)/2))) #mean and standard deviations, juvenile fertility
Y.f.sd <- mean(c(((1.8 - 1.625)/2) ,((2.3375 - 1.8)/2))) #mean and standard deviations, yearling fertility
A.f.sd <- mean(c(((2.835 - 1.7)/2),((2.925 - 1.7)/2))) #mean and standard deviations, adult fertility
f.sd.vec <- c(J.f.sd, Y.f.sd, A.f.sd, A.f.sd, A.f.sd, A.f.sd) #mean and standard deviations vector, juvenile and adult fertility 

#survival
s.vec <- c(0.33, 0.40, 0.66, 0.66, 0.66, 0.58) ##feral pig survival # intermediate values from Bieber and Ruf (2005), except 6 which is the poor value.
## alternative values from Choquenot (1996) ## s.vec <- c(0.45, 0.675, 0.675, 0.675, 0.675, 0.675)

## SD calculated as the mean(c(((mid range value - low range value)/2),((high range value - mid range value)/2)))
# survival errors based on Bieber and Ruf (2005). Unsure how to (or if I can) produce SD from Choquenot (1996)
J.s.sd <- mean(c(((0.33 - 0.25)/2),((0.52 - 0.33)/2))) #mean and standard deviations, juvenile survival
Y.s.sd <- mean(c(((0.40 - 0.31)/2),((0.60 - 0.40)/2))) #mean and standard deviations, yearling survival
A.s.sd <- mean(c(((0.66 -0.58)/2),((0.71 -0.66)/2))) #mean and standard deviations, adult survival
s.sd.vec <- c(J.s.sd, Y.s.sd, A.s.sd, A.s.sd, A.s.sd, A.s.sd) #mean and standard deviations vector, juvenile and adult survival

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
pop.found <- 500/2 # +/- 50 founding population size personal communication, B. Page, PIRSA,  2021 . 
##Divided by 2 for females only
ssd <- stable.stage.dist(popmat) ## sum of all values in ssd = 1
init.vec <- ssd * pop.found #initial population vector. Sum of all values here = 500 as suggested by PIRSA
## (can change to whatever we want our founding pop to be)

#################
## project
## set time limit for projection in 1-yr increments
yr.now <- 2020 # Change to start year to 2021?
#************************

yr.end <- 2220 #end year for our projection timeframe
#************************
t <- (yr.end - yr.now) #timeframe
## we can adjust this later to explore different scenarios 
## e.g., probability of achieving eradication by various control techniques or proportional rate of offtake within 3 years

tot.F <- sum(popmat.orig[1,]) ## this is the sum of fertility values in row 1 of the deterministic matrix. 
## Initially this is the same as f.vec, but popmat will change for each iteration as values are stochastically resampled 
popmat <- popmat.orig #resets matrix 
yr.vec <- seq(yr.now,yr.end) #year vector, 2020, 2021, 2022.... Increases in increments of one by default so increment size doesn't need to be specified

## set population storage matrices
n.mat <- matrix(0, nrow=age.max,ncol=(t+1)) #matrix with values = 0, number of columns = to max.age (=6) and number of columns = 
n.mat[,1] <- init.vec #fill first matrix column with initial population vector

## set up projection loop
for (i in 1:t) { #  for every year of the specified timeframe...
  n.mat[,i+1] <- popmat %*% n.mat[,i] #  the corresponding column of n.mat 
  #is populated with the product of current year column multiplied by the deterministic matrix popmat
}

n.pigs <- colSums(n.mat) #number of pigs in any year is the sum of all values in the corresponding column 
## no density reduction treatment or carrying capacity at this stage
yrs <- seq(yr.now, yr.end, 1)
plot(yrs, (n.pigs),type="l",lty=2,pch=19,xlab="year",ylab="N")

# compensatory density feedback: SURVIVAL
# K = carry capacity
# population rate of increase relative to carry capacity. Larger distance between population and K = faster population growth
K.max <- 2500
## Southgate, R. (2018). Feral Pig Sus Scrofa and Feral Cat Felis catus monitoring on Kangaroo Island: 2014 to 2017 - found pig density to be stable during study period
## pre-fire pig population estimated to be ~5000 (PIRSA personal communication)
## Assume pop remained stable between 2017 study and 2019 bushfire, so set k.max as 5000/2 (0.5 to account for females only) 

K.s.vec <- c(1,K.max*0.4,K.max*0.5,0.7*K.max,0.85*K.max,0.95*K.max) ##describes the x axis of the reduction curve
red.s.vec <- c(1,0.995,0.98,0.96,0.92,0.85) ## describes the y axis of the reduction curve
plot(K.s.vec,red.s.vec,pch=19,type="b")
Kred.s.dat <- data.frame(K.s.vec,red.s.vec)

# logistic power function a/(1+(x/b)^c) #fits logistic power function to population relative to carry capacity, K
s.param.init <- c(1, K.max, 3) ## These values are arbitrary 
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

# compensatory density feedback: FERTILITY
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


## set up projection loop for deterministic population with S density feedback only
## this is to help me visualise the effect of S feedback - will delete from final code
for (i in 1:t) {
  totN.i <- sum(n.mat[,i])
  pigs.s.red <- s.a.lp/(1+(totN.i/s.b.lp)^s.c.lp)
  diag(popmat[2:age.max,]) <- s.vec[1:5]*pigs.s.red
  popmat[age.max,age.max] <- s.vec[6]*pigs.s.red
  popmat[1,1:6]<- f.vec
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pigs <- colSums(n.mat)
plot(yrs, n.pigs,type="l",main= "deterministic pop projection",sub = "S density feedback", lty=1,pch=19,xlab="year",ylab="N",ylim=c(0,1.8*K.max)) #untreated population increases, rate of increase relative to K, no stochastic sampling
abline(h=K.max, lty=2, col="red") #dashed red line indicating carry capacity

## set up projection loop for deterministic population with F density feedback
## this is to help me visualise the effect of F feedback - will delete from final code
for (i in 1:t) {
  totN.i <- sum(n.mat[,i])
  pigs.f.red <- f.a.lp/(1+(totN.i/f.b.lp)^f.c.lp)
  diag(popmat[2:age.max,]) <- s.vec[1:5]
  popmat[age.max,age.max] <- s.vec[6]
  popmat[1,1:6]<- f.vec*pigs.f.red
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pigs <- colSums(n.mat)
plot(yrs, n.pigs,type="l",main= "deterministic pop projection", sub = "F density feedback", lty=1,pch=19,xlab="year",ylab="N",ylim=c(0,1.8*K.max)) #untreated population increases, rate of increase relative to K, no stochastic sampling
abline(h=K.max, lty=2, col="red") #carry capacity

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


#################################################### 
## iterations and quasi ext for each following model
####################################################
iter <- 1000 #final model run at 10 000
itdiv <- iter/100 #final model rate at iter/1000

################################################################################################################
## untreated population
###############################################################################################################
## stochastic projection with survival and fertility density feedback
## set storage matrices & vectors

n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1)) #storage matrix with 1000 rows and 16 columns based on timeframe

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

plot(yrs,n.md,type="l", main= "untreated stochastic population projection", sub = "S and F density feedback", xlab="year", ylab="N", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,n.up,lty=2,col="red",lwd=1.5)

datN <- data.frame(yrs, n.md, n.lo, n.up)
tail(datN)
## not entirely happy with this plot. 
## Short 21 year time frame hides the fact that population is slowly declining rather than reaching equilibrium near K

###############################################################################################################################
## constant proportional yearly harvest
###############################################################################################################################

# harvest rate
harv.prop.consist <- seq(0.2,0.99,0.05) #sequence harvest/culling quotas, e.g remove 0.2-.99 proportion of founding pop in increasing increments of 0.05

# define our quasi-extinction probability storage vector
min.med.n <- min.lo.n <- min.up.n <- rep(0,length(harv.prop.consist))

for (s in 1:length(harv.prop.consist)) {
   
  # set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1)) 
  #storage matrix with 1000 rows and 16 columns (years) based on timeframe
 
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
  
  n.md <- apply((n.sums.mat), MARGIN=2, mean, na.rm=T) # minimum over all iterations
  n.up <- apply((n.sums.mat), MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply((n.sums.mat), MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  plot(yrs,n.md,type="l",xlab="year", ylab="minimum N", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
  lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
  lines(yrs,n.up,lty=2,col="red",lwd=1.5)
  
  print("##############")
  print(paste("harvest proportion = ", harv.prop.consist[s], sep=""))
  print("##############")
  
} # ends S loop

plot(harv.prop.consist, min.med.n, type="l", pch=19, xlab="harvest proportion", ylab="min N", ylim=c(min(min.lo.n),max(min.up.n)))
lines(harv.prop.consist, min.lo.n, col="red", lty=2)
lines(harv.prop.consist, min.up.n, col="red", lty=2)

minn.prop.pop <- data.frame(harv.prop.consist, min.med.n, min.lo.n, min.up.n)


##################################################################################################################################################
## high harvest for first 2 years, constant proportional harvest in remaining years
#####################################################################################################################################################

# harvest rate
harv.prop.init <- seq(0.5,0.9,0.05)
harv.prop.maint <- seq(0.1,0.5,0.05)

# storage
minn.med.mat <- minn.lo.mat <- minn.up.mat <- pmin.med.mat <- pmin.lo.mat <- pmin.up.mat <- matrix(data=NA, ncol=length(harv.prop.maint), nrow=length(harv.prop.init)) #storage matrices

for (m in 1:length(harv.prop.maint)) {
  
  for (n in 1:length(harv.prop.init)) {
    
    # storage
    n.sums.mat <- p.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
    
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
        
        ## harvest 
        if (i < 3) {
          n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.init[n], 0), 0)
        } else {
          n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.maint[m], 0), 0)
        }
        
        if (length(which(n.mat[,i+1] < 0)) > 0) {
          n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
        }
        
      } # end i loop
      
      n.sums.mat[e,] <- as.vector(colSums(n.mat))
      p.sums.mat[e,] <- n.sums.mat[e,] / pop.found
      
      if (e %% itdiv==0) print(e) 
    } # end e loop (stochastic iterations)
    
    min.pop.vec <- apply(n.sums.mat, MARGIN=1, min, na.rm=T)
    min.ppop.vec <- apply(p.sums.mat, MARGIN=1, min, na.rm=T)
    
    # median, lower & upper minimum population sizes
    minn.med.mat[n, m] <- median(min.pop.vec, na.rm=T) 
    minn.lo.mat[n, m] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
    minn.up.mat[n, m] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
    
    # median, lower & upper minimum proportional population sizes
    pmin.med.mat[n, m] <- median(min.ppop.vec, na.rm=T)
    pmin.lo.mat[n, m] <- quantile(min.ppop.vec, probs=0.025, na.rm=T) 
    pmin.up.mat[n, m] <- quantile(min.ppop.vec, probs=0.975, na.rm=T)
    
    
    print("##############################")
    print(paste("init harvest proportion = ", harv.prop.init[n], sep=""))
    print("##############################")
    
  } # end n loop (initial harvest rate)
  
  print("##############################")
  print(paste("maint harvest proportion = ", harv.prop.maint[m], sep=""))
  print("##############################")
  
} # end m loop (maintenance harvest rate)

## plot 3D surfaces
f1 <- list(
  family = "Avenir Light",
  size = 26,
  color = "black"
)
f2 <- list(
  family = "Avenir Light",
  size = 18,
  color = "black"
)
f3 <- list(
  family = "Avenir Light",
  size = 16,
  color = "black"
)

# minimum proportional population size (median)
par(mar=c(5,5,2,8))
pminmed3d <- plot_ly(z = ~pmin.med.mat, autocontour=F, type="contour", line = list(smoothing = 0.90), contours = list(start=0.01, end=0.32, size=0.025, showlabels = TRUE, labelfont=list(
  size=18, family="Avenir Light", face="bold", color="white"))) %>%
  colorbar(title = "med min pN1", titlefont=f2, tickfont=f2) %>%
  layout(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)))
pminmed3d

twophase.med <- data.frame(pmin.med.mat)
colnames(twophase.med) <- harv.prop.maint
rownames(twophase.med) <- harv.prop.init

twophase.lo <- data.frame(pmin.lo.mat)
colnames(twophase.lo) <- harv.prop.maint
rownames(twophase.lo) <- harv.prop.init

twophase.up <- data.frame(pmin.up.mat)
colnames(twophase.up) <- harv.prop.maint
rownames(twophase.up) <- harv.prop.init

pmin3d <- plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~pmin.med.mat) %>%
  add_surface(z = ~pmin.lo.mat, opacity = 0.55) %>%
  add_surface(z = ~pmin.up.mat, opacity = 0.55) %>%
  layout(scene = list(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)),
    zaxis = list(title="min pN1", tickfont=f3, titlefont=f1)))
pmin3d

# quasi ext (median)

      ## error below - object qext.mat not found, therefore object minmed3d cannot be created
      ## should qext.mat be created along with other storage matrices in row 395?
par(mar=c(5,5,2,8))
minmed3d <- plot_ly(z = ~qext.mat, autocontour=T, type="contour", line = list(smoothing = 0.90), contours = list(showlabels = TRUE, labelfont=list(
  size=18, family="Avenir Light", face="bold", color="white"))) %>%
  colorbar(title = "qE", titlefont=f2, tickfont=f2) %>%
  layout(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)))
minmed3d

min3d <- plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~minn.med.mat) %>%
  add_surface(z = ~minn.lo.mat, opacity = 0.55) %>%
  add_surface(z = ~minn.up.mat, opacity = 0.55) %>%
  layout(scene = list(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)),
    zaxis = list(title="min N1", tickfont=f3, titlefont=f1)))
min3d


#######################################################################################################################################################################################################
######################################################################################### COSTS #################################################################################################
########################################################################################################################################################################################################
## high harvest for initial 2 years, consistent harvest remaining years
############################################################################################################################################################################################################
## contributed by CJA Bradshaw
###########################################################################################

Dudley.area <- 375*100 #ha
KI.area <- 4405*100  #ha

## cost parameters
felixer.unit <- 13000 # AU$ #cost from felix vs felixer report 
trap.unit <- c(157,297) # AU$ # cost per trap from traps.com.au
shoot.ph <- 518.54/20 # ammo & labour (total AU$ over 20 hours) # Holmes et al 2015
only.bait.unit <- (2.07 + 0.2) # From Curiosity correspondence. $2.07 per bait + $0.20  + $250 administration fee per order, + freight fee. 500 baits per pack

# Felixers data from Moseby et al (2020)
# 20 felixers for "felixer paddock", Arid Recovery (n1 = 48), killed 31 cats (n2 = 17), over 41 days 
num.felixer <- round((20/48) * pop.found, 0)
pfelixer.killr <- (31/48 * (1/(41/365)))


felixer.area <- 26*100 #ha; density from Arid Recovery trial, "felixer paddock" = 26 km^2
felixer.dens <- 20/felixer.area #20 felixer traps over the area, average density 0.77 felixers/km^2
KI.felixer.num <- round(KI.area * felixer.dens, 0) # number of felixers needed if same density was applied throughout Kangaroo Island
KI.felixer.num # not neccessarily reflective of the use of felixers as they are used in targeted areas and spread sporadically, as opposed to systematicaly placed like traps or baits

# traps
# 40 traps killed 21 over 148 days Hodgens 
ptrap.killr <- (21/262 * (1/(148/365)))
trap.dens <- 40/Dudley.area
KI.trap.num <- round(KI.area * trap.dens, 0)
KI.trap.num

# shooting
# 14725 person-hours killed 872 (+ 172 from wounds) cats (Marion) Parkes et al. 2014 & Bloomer & Bester 1992
# assume cats not killed by Felixers & traps shot by hunters
cats.pph <- (872+172)/14725


# baiting 
# 943 baits killed 11 cats over 18.86 km^2. Pre-baiting dens = 1.18 cats/km^2, post-baiting = 0.58 cats/km^2. Ref, 'Dudley peninsula feral cat eradication operations plan: summary may 2020 - mid 2023"
# KI uses Curiosity (PAPP)
## Kangaroo Island area - 4,405 km^2  or 440 500 ha 
# can't bait built-up areas, need 500m buffer zone around towns, built up areas 362 ha, how many built up areas? 
## parndana (second largest town) approx area as circle - 2km^2 (??), + 500m buffer = area 10km^2. 5 'main towns' KI. 5*2 = 10km^2 or 1000 ha, :- approx 1000 ha can't bait urban
# can't bait beaches. KI 540 km coastline, arbitary 100m buffer around coastline = 594 km can't bait + buffer zone. 540 * 1.1 = area no bait
## Dirk Hartog Island, 15 cats collared, average density 0.701 cats/km^2 (average area = 10.515 km^2 (A = 15*0.701), 50 baits per km, baits = 50*10.515 = 525.75), 14 died following bait consumption ... 525.74/14 = 37.55 baits/cat
# Dirk Hartog Island, used eradicat (1080)
nobaitfarm <- (2303 - (2303*.94))*100 #ha; can't be baited 
nobaitcoast <- (540 * 1.1)*100 #ha; dist around costline, *1.1 for the 100m buffer around coast 
nobaittown <- 1000 #ha; can't bait town 
nobaitarea <- nobaitfarm + nobaitcoast + nobaittown # total area can't be baited
baitareaKI <- KI.area - nobaitarea #ha; area eligible for baiting 
baitdens <- 50/100 # 50 baits per km^2 converted to baits per ha
baitnum <- (baitareaKI * baitdens) #number of baits need for entire Island

baitadminfee <- 250 # administration fee, once off for baits, or twice off for two years 
baitdrop.time <- (30/60/60) * (baitareaKI/100) # 50 baits drops every 30 seconds or 1 km^2, with plane speed 240km/h - bait area/100 to convert back to km^2
trips <- baitnum/3500 #can only take 3500 baits per trip 
upandback <- seq(1,32,0.5) #ha; dist from airport to start of each bait transect
averagedist <- (sum(upandback))/(length(upandback)) #ha; average dist from airport to transect
baitreloadtime <- ((averagedist*trips)*2)/240 #52 trips needed to drop all baits, *2 for to and from airport, plan speed 240km
planehph <- 750  #cost per hour plane hire when dropping baits, inc wages of 2x pilots (1x loading and dropping baits)
planeferrycost <- (3*600)*2 # 3 hour flight William Creek to Kangaroo Island (*2 for return), $600 p/h for charter
baitreloadcost <- baitreloadtime*planehph #average extra cost for reloading baits
baitdropcost <- baitdrop.time*planehph
planecost <- planeferrycost + baitreloadcost +  baitdropcost
cost.total.bait <- planecost + baitadminfee + (only.bait.unit * baitnum) #cost of total baiting to cover entire island 

  
pbait.killr <- 14/525.74 # 14 cats killed by 525.74 baits



###########################################################################################

###########################################################################################
## Type III functional response (reduction in capture efficiency with decreasing density)
max.eff <- 1 # max efficiency 
min.eff <- 0 #min efficiency
max.pN <- 1 #max population proportion
min.pN <- 0 #min population proportion
infl.eff <- 0.5

pN.vec <- c(min.pN, 0.2, 0.4, 0.5, 0.7, 0.8, max.pN) 
eff.vec <- c(min.eff, 0.05, 0.3, infl.eff, 0.85, 0.95, max.eff)
plot(pN.vec, eff.vec, type="b", pch=19)
eff.dat <- data.frame(pN.vec, eff.vec)
colnames(eff.dat) <- c("pN", "eff")

# a/(1 + b*e^(-cx)) (logistic)
param.init <- c(1, 85, 8.9)
fit.eff <- nls(eff ~ a/(1+(b*exp(-c*pN))), 
               data = eff.dat,
               algorithm = "port",
               start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
               trace = TRUE,      
               nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fit.eff.summ <- summary(fit.eff)
plot(pN.vec,eff.vec,pch=19,xlab="pN",ylab="efficiency")
pN.vec.cont <- seq(0,1,0.01)
pred.eff.fx <- coef(fit.eff)[1]/(1+(coef(fit.eff)[2]*exp(-coef(fit.eff)[3]*pN.vec.cont)))
lines(pN.vec.cont,pred.eff.fx,lty=2,lwd=3,col="red")

a.eff <- coef(fit.eff)[1]
b.eff <- coef(fit.eff)[2]
c.eff <- coef(fit.eff)[3]

###########################################################################################


#################################################### 
## iterations and quasi ext for each following model
####################################################
iter <- 10000 #final model run at 10 000
itdiv <- iter/1000 #final model rate at iter/1000

## run choices
## make up shortfall in kill by ...
#shortfall.method <- "F" # adding Felixer units
shortfall.method <- "T" # adding traps
#shortfall.method <- "H" # increasing hunting pressure

# harvest rate
harv.prop.init <- seq(0.5,0.9,0.05)
harv.prop.maint <- seq(0.1,0.5,0.05)
q.ext <- 20

# storage
qext.mat <- minn.med.mat <- minn.lo.mat <- minn.up.mat <- pmin.med.mat <- pmin.lo.mat <- pmin.up.mat <- totcost.med <- totcost.lo <- totcost.up <- matrix(data=NA, ncol=length(harv.prop.maint), nrow=length(harv.prop.init))

for (m in 1:length(harv.prop.maint)) {
  
  for (n in 1:length(harv.prop.init)) {
    
    # storage
   init.k.sums.mat <- k.sums.mat <- n.sums.mat <- p.sums.mat <- totalcost.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
    
    for (e in 1:iter) {
      
      popmat <- popmat.orig
      
      init.k.mat <- n.mat <- k.mat <- matrix(0, nrow=age.max,ncol=(t+1))
      n.mat[,1] <- init.vec
      
      for (i in 1:t) {
        # stochastic survival values
        s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
        s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
        s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
        
        # stochastic fertility sampler (gaussian)
        fert.stch <- rnorm(length(popmat[,1]), popmat[1,], s.sd.vec)
        fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
        
        totN.i <- sum(n.mat[,i])
        pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
        
        popmat[1,] <- fert.stoch  # add new stochastically resampled fertilities
        diag(popmat[2:age.max,]) <- s.stoch*pred.red
        #popmat[age.max,age.max] <- s.stoch[age.max]*pred.red
        
        n.mat[,i+1] <- popmat %*% n.mat[,i]
        
        ## harvest things here
        if (i < 3) {
          k.mat[,i+1] <- round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.init[n], 0), 0)
          n.mat[,i+1] <- n.mat[,i+1] - k.mat[,i+1]
          init.k.mat[,i+1] <- n.mat[,i+1]
        } else {
          k.mat[,i+1] <- round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.maint[m], 0), 0)
          n.mat[,i+1] <- n.mat[,i+1] - k.mat[,i+1]
        }
        
        if (length(which(n.mat[,i+1] < 0)) > 0) {
          n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
        }
        if (length(which(k.mat[,i+1] < 0)) > 0) {
          k.mat[which(k.mat[,i+1] < 0), i+1] <- 0
        }
        
      } # end i loop
      
      init.k.sums.mat[e,] <- as.vector(colSums(init.k.mat))
      k.sums.mat[e,] <- as.vector(colSums(k.mat))
      n.sums.mat[e,] <- as.vector(colSums(n.mat))
      p.sums.mat[e,] <- n.sums.mat[e,] / pop.found
      
      # cost of cats killed here
      eff.vec.iter <- a.eff/(1+(b.eff*exp(-c.eff*p.sums.mat[e,]))) # efficiency this iteration
      
      # calculate numbers killed per year using baiting and trapping first two years
      bait.kill.base <- round(init.k.sums.mat[e,] * (eff.vec.iter*pbait.killr), 0)
      trap.kill.base <- round(k.sums.mat[e,] * (eff.vec.iter*ptrap.killr), 0)
      bt.kill.base <- trap.kill.base + bait.kill.base
      shortfall <- k.sums.mat[e,] - bt.kill.base # how many cats not being killed by these methods?
      
      #base cost
      base.cost <- (cost.total.bait*2) + (KI.trap.num*runif(1,min=trap.unit[1],max=trap.unit[2])) # at initial roll-out numbers
      
      # make up shortfall
      if (shortfall.method == "H") {
        makeup.iter <- shoot.ph*(shortfall / (cats.pph*eff.vec.iter)) # how many person-hours required to make up shortfall?
      }
      if (shortfall.method == "F") {
        makeup.iter <- felixer.unit*(shortfall / (pfelixer.killr*eff.vec.iter)) # how many person-hours required to make up shortfall?
      }
      if (shortfall.method == "T") {
        makeup.iter <- (runif(1,min=trap.unit[1],max=trap.unit[2]))*(shortfall / (ptrap.killr*eff.vec.iter)) # how many person-hours required to make up shortfall?
      }
      
      totalcost.mat[e,] <- base.cost + makeup.iter 
      
      if (e %% itdiv==0) print(e) 
    } # end e loop (stochastic iterations)
    
    min.pop.vec <- apply(n.sums.mat, MARGIN=1, min, na.rm=T)
    min.ppop.vec <- apply(p.sums.mat, MARGIN=1, min, na.rm=T)
    
    # median, lower & upper minimum population sizes
    minn.med.mat[n, m] <- median(min.pop.vec, na.rm=T) 
    minn.lo.mat[n, m] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
    minn.up.mat[n, m] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
    
    # median, lower & upper minimum proportional population sizes
    pmin.med.mat[n, m] <- median(min.ppop.vec, na.rm=T)  
    pmin.lo.mat[n, m] <- quantile(min.ppop.vec, probs=0.025, na.rm=T) 
    pmin.up.mat[n, m] <- quantile(min.ppop.vec, probs=0.975, na.rm=T)
    
    # quasi-extinction
    qext.mat[n, m] <- (sum(ifelse(round(min.pop.vec, 0) < q.ext, 1, 0)) / iter)
    
    ## costs
    totcost.vec <- apply(totalcost.mat, MARGIN=1, sum, na.rm=T)
    totcost.med[n, m] <- median(totcost.vec, na.rm=T)
    colnames(totcost.med) <- harv.prop.maint
    rownames(totcost.med) <- harv.prop.init
    
    totcost.lo[n, m] <- quantile(totcost.vec, probs=0.025, na.rm=T)
    colnames(totcost.lo) <- harv.prop.maint
    rownames(totcost.lo) <- harv.prop.init
    
    totcost.up[n, m] <- quantile(totcost.vec, probs=0.975, na.rm=T)
    colnames(totcost.up) <- harv.prop.maint
    rownames(totcost.up) <- harv.prop.init
    
    print("##############################")
    print(paste("init harvest proportion = ", harv.prop.init[n], sep=""))
    print("##############################")
    
  } # end n loop (initial harvest rate)
  
  print("##############################")
  print(paste("maint harvest proportion = ", harv.prop.maint[m], sep=""))
  print("##############################")
  
} # end m loop (maintenance harvest rate)

## plot 3D surfaces
f1 <- list(
  family = "Avenir Light",
  size = 26,
  color = "black"
)
f2 <- list(
  family = "Avenir Light",
  size = 18,
  color = "black"
)
f3 <- list(
  family = "Avenir Light",
  size = 16,
  color = "black"
)

# total cost (median)
par(mar=c(5,5,2,8))
costcontmed3d <- plot_ly(z = ~totcost.med, autocontour=T, type="contour", line = list(smoothing = 0.90), contours = list(showlabels = TRUE, labelfont=list(
  size=18, family="Avenir Light", face="bold", color="white"))) %>%
  colorbar(title = "tot $", titlefont=f2, tickfont=f2) %>%
  layout(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)))
costcontmed3d

cost3d <- plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~totcost.med) %>%
  add_surface(z = ~totcost.lo, opacity = 0.55) %>%
  add_surface(z = ~totcost.up, opacity = 0.55) %>%
  layout(scene = list(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)),
    zaxis = list(title="tot $", tickfont=f3, titlefont=f1)))
cost3d

