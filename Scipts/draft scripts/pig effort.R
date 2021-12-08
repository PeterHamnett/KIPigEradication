source("~/Documents - Peter’s MacBook Pro/R Resources/Projects /Pigs_test/data/Final/pigs_effort_no_LU.csv")
library(readr)
library(dplyr)
pigs_effort_no_LU <- read_csv("data/Final/pigs_effort_no_LU.csv")
View(pigs_effort_no_LU)
effort<-data.frame(pigs_effort_no_LU)
head(effort)

# create a copy of effort in case this goes wrong
effort_dates<-effort

# convert Date column from chr to date format
effort_dates$Date<-as.Date(effort$Date, format = "%d/%m/%y")
# it worked so overwrite effort with the new date ata
effort<-effort_dates

# sort by date 
effort <- effort[order(effort$Date),] 
# Oldest records are now at the top of the data frame

# Total number of pigs killed (where effort was recorded)= 
sum(effort$no..killed)
# 683

# Create column showing proportion of pigs remaining after each event in chronological order, assuming starting proportion is 1 and final proportion is 0
effort$Prop_remaining <- (sum(effort$no..killed) - cumsum(effort$no..killed))/sum(effort$no..killed)
# cumsum is the sum of the current value and all previous values in the column
# so here we have (total pigs killed - pigs killed to date)/total pigs killed
# give proportions in declining in order from just short of 1 (1 is before any pigs were killed) down to 0
head(effort)
tail(effort)

# create a new column with the inverse of pigs/hr i.e., hrs per pig
## this will allow us to quantify changes in cost relative to proportion of population remaining, once we have calculated cost per hour for each control type.

effort$hrs_pig<-effort$effort.hrs/effort$no..killed
#Check the new column was added
head(effort)

# Create subsets of effort dataframe for each control type
poison_efrt <-effort %>% filter (Control.Type == "poisoned")
shot_efrt <-effort %>% filter (Control.Type == "shot")
trap_efrt <-effort %>% filter (Control.Type == "trapped")
TAAC_efrt <-effort %>% filter (Control.Type == "TAAC")

# produce histograms of hrs_pig for all control types to review distribution of values
hist(poison_efrt$hrs_pig)
hist(shot_efrt$hrs_pig)
hist(trap_efrt$hrs_pig)
hist(TAAC_efrt$hrs_pig)

# now we can plot pigs_hr and hrs-pig against proportion remaining for each control type to see if there is a relationship between the two.
# we would expect efficiency to decrease as numbers decline

plot(pigs_hr~Prop_remaining, data = poison_efrt, xlim=rev(range(Prop_remaining)), main = "Poison CPUE")
plot(pigs_hr~Prop_remaining, data = trap_efrt, xlim=rev(range(Prop_remaining)), main = "Trapping CPUE")
plot(pigs_hr~Prop_remaining, data = shot_efrt, xlim=rev(range(Prop_remaining)), main = "Shooting CPUE")
plot(pigs_hr~Prop_remaining, data = TAAC_efrt, xlim=rev(range(Prop_remaining)), main = "TAAC CPUE" )

# and hrs-pig 
plot(hrs_pig~Prop_remaining, data = poison_efrt, xlim=rev(range(Prop_remaining)), main = "Poison: Effort per pig")
plot(hrs_pig~Prop_remaining, data = trap_efrt, xlim=rev(range(Prop_remaining)), main = "Trapping: Effort per pig")
plot(hrs_pig~Prop_remaining, data = shot_efrt, xlim=rev(range(Prop_remaining)), main = "Shooting: Effort per pig")
plot(hrs_pig~Prop_remaining, data = TAAC_efrt, xlim=rev(range(Prop_remaining)), main = "TAAC: Effort per pig")

#lets try using ggplot2 and adding a smoothed trend line
#pigs_hr
library(ggplot2)
ggplot(poison_efrt, aes(x = Prop_remaining,y = pigs_hr)) + scale_x_reverse() + geom_point() + geom_smooth()  + labs(title = "Poison CPUE")
ggplot(trap_efrt, aes(x = Prop_remaining,y = pigs_hr)) + scale_x_reverse() + geom_point() + geom_smooth() + labs(title = "Trapping CPUE")
ggplot(shot_efrt, aes(x = Prop_remaining,y = pigs_hr)) + scale_x_reverse() + geom_point() + geom_smooth() + labs(title = "Shooting CPUE")
ggplot(TAAC_efrt, aes(x = Prop_remaining,y = pigs_hr)) + scale_x_reverse() + geom_point() + geom_smooth() + labs(title = "TAAC CPUE")

# and  hrs_pig
ggplot(poison_efrt, aes(x = Prop_remaining,y = hrs_pig)) + scale_x_reverse() + geom_point() + geom_smooth( method = lm)  + labs(title = "Poison: Effort per pig")
ggplot(trap_efrt, aes(x = Prop_remaining,y = hrs_pig)) + scale_x_reverse() + geom_point() + geom_smooth() + labs(title = "Trapping: Effort per pig")
ggplot(shot_efrt, aes(x = Prop_remaining,y = hrs_pig)) + scale_x_reverse() + geom_point() + geom_smooth() + labs(title = "Shooting: Effort per pig")
ggplot(TAAC_efrt, aes(x = Prop_remaining,y = hrs_pig)) + scale_x_reverse() + geom_point() + geom_smooth() + labs(title = "TAAC: Effort per pig")



