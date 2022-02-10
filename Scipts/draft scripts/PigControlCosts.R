######################################################################################################################################
############################## COSTS #################################################################################################
######################################################################################################################################

# Project.area <- 192700 #ha
## not sure if area and density need to factored into calculations


# no.pigs <-      ## derived from constant proportional cull 
## proportion of pigs culled*population at start of time increment

# no.hrs <-       ## derived from no.pigs*density dependent efficiency functions for each control method.

## cost parameters
#### GENERIC COSTS####
# these apply to more than one control method
labour.ph <-36.87 ## $36.87 per hour, based on SA Public Sector Award OPS4 classification. 
## Applies to all control types
## TAAC accrues labour in 8hr blocks (daily) whereas others accrue in 1hr increments

Bullets.pp <- 4   ## 2 bullets per pig @ $2 per bullet
## applies to shooting, TAAC and trapping.
## 2 per pig allows for misses and sight zeroing.

freefeed.ph <- 14 ## based on daily deployment of 10kg grain at $1.4 per kg. Deployment of grain takes 1 hr per day 
## Same rate applies for all methods where free feeding is used: 1 hr FF per day @ $14/hr

#### METHOD SPECIFIC COSTS####

####THERMAL ASSISTED AERIAL CULL (TAAC)####
TAAC.initial <-12146 # total initialisation cost $12,246.00, comprising...
# crew and helicopter mobilisation (from Jindabyne, NSW to KI)  $7,400.00 
# Project data management $2,400.00  
# Initial fuel delivery $2,446.00
TAAC.crew.pd <-2482.66 # daily cost for flight time and crew wages. Contractor charges daily, not by the hour
TACC.crew.accomm <-420.00# daily crew food and accommodation
TAAC.mrksmn.accomm <-125 # marksman meals and accommodation (marksman supplied by SA Gov)
TAAC.pd <-TAAC.crew.pd +  TACC.crew.accomm + (labour.ph*7.5) + TAAC.mrksmn.accomm # Total daily TAAC costs
## labour*7.5 assumes marksman is allocated to TAAC in daily increments, rather than undertaking other activities when not flying
TAAC.daily.ft <-3.7#average daily flight time (hrs/day). AM + PM flights    
## Could divide daily cost by 3.7 to calculate hourly rate, but prefer to accrue costs daily as contractor charges daily rather than hourly
## therefore, TAAC.total.pd accrued for every 3.7hrs effort  
TAAC.resupply<-1246.00# cost for additional fuel every 30 days = 30 x 3.7hrs =  111hrs   
## cost accrued in increments of 111hrs TAAC work completed.

      ## HAVEN'T WORKED OUT HOW TO INCORPORATE FUEL RESUPPLY YET##


#TAAC.total.cost <- TAAC.initial + ((TAAC.pd*no.hrs)/TAAC.daily.ft) + Bullets.pp*no.pigs + TAAC.resupply
# total cost is initialisation + ((cost per day*no.hrs)/3.7) + cost per pig*no.pigs +     

## no.pigs needs to be derived from proportional cull loop
## no.hrs needs to be derived from density dependent efficiency function

####SHOOTING#####

## opportunistic and Free-fed shooting are combined so need to work out an approach for combining proportional FF costs into total shooting effort
## 83 events where O or FF could be determined
## ratio of events O:FF = 68:15 = 0.2205882 FF events for every Opportunistic event
## proportion of events: O = 15/83 = 0.1807229
## proportion of events: FF = 68/83 = 0.8192771
## 216 kills where O or FF could be determined
## ratio of kills O:FF = 193:23 = 0.119171 FF kills for every Opportunistic kill
## proportion of total kills : O = 193/216 = 0.8935185
## prop total kill: FF = 23/216 = 0.1064815
##  505 hrs total shooting effort recorded 
## proportion of effort:O = 208/505 = 0.4118812
## proportion of effort: FF = 297/505 = 0.5881188
## ratio of effort O:FF = 208:297 = 0.7003367 hrs Opportunistic effort for every 1hr FF effort
## proportion of effort seems the relevant measure as can be applied directly to the cull effort estimate. 
## FF cost will be slightly exaggerated as final hour or two of FF shooting doesn't accrue feed costs. 
## However, FF cost is negligible, only $14 per hour. 
## Effort data aren't stratified into time spent FF vs shooting, but feel this is close enough based on available data
shoot.ff.prop <- 0.5881188
shoot.ff.ph <-freefeed.ph*shoot.ff.prop
## so, total shooting estimate is labour per hour x number of hours + free feed cost per hour x no hours + cost of bullets per pig x no pigs killed.

#shoot.total.cost <- labour.ph*no.hrs + shoot.ff.ph*no.hrs + bullets.pp*no.pigs 
# replace no.hrs and no.pigs with outputs from relevant storage matrices

#### TRAPPING ####

trap.initial <- labour.ph*20
## per event
## MK (PIRSA) says reported effort does not include set up/breakdown of traps when moved between locations
## this is about 20hrs per trap site (2 staff members required for 5-10 trips to the trap site to slowly erect trap panels)

trap.ff.labour <-labour.ph*mean(trap_efrt$efforthrs)
## per event
trap.ff<-freefeed.ph*mean(trap_efrt$efforthrs)
## per event
trap.ff.cost<-trap.ff.labour + trap.ff
## per event
## Average observed effort (excluding outliers) is 9.875hrs.
## MK (PIRSA) suggested average was about 15hrs, but changed to observed average as more defensible

## Assuming trap initialisation and deployment effort are constant (based on estimated average effort for set up and FF), 
## no. pigs per trap decreases as hrs per pig increases with declining pig abundance.

trap.unit <- 9500.00 # cost per unit for the Jager Pro trap
## All effort completed using Jager Pro traps - no conventional traps used so cost per trap unit is consistent


#trap.number <-ceiling(no.hrs/10*(20+mean(trap_efrt$efforthrs))

# no. of traps required to satisfy proportional cull quota
## total effort required per year divided by the amount of effort a single trap can do per year
## traps are only deployed when needed, rather than being permanently deployed at a site
## trap site duration is about 4-6 weeks so..
## maximum trap use is about 8-12 trap sites per year. Use 10 events per trap per year
## used ceiling() to round up to nearest integer as we can only deploy whole traps
## update once functional response model has been integrated into proportional cull simulation.

# trap.total.cost.pa <-trap.number*(trap.initial + trap.ff.cost + trap.unit) + bullets.pp*no.pigs
## but! If fewer traps needed in subsequent years than in first year, first year traps can be reused without buying any more.

#####BAITING######
avg.bait.effort <- mean(poison_efrt$efforthrs)
bait.labour<-labour.ph*avg.bait.effort
#labour cost per event for deploying grain, placebo and toxic bait
## based on average observed effort (excluding outliers) ... avg.bait.effort = 12.7 hrs
## MK estimated effort of 15 days followed by 4 days placebo and 1 day toxic baits, but this is not backed up by the data 
## HogGone manufacturer recommends 2 days placebo and 2 days toxic bait, but PIRSA use a modified approach
## can change to labour.ph*20 as per MK estimate, but prefer to use real data rather than arbitrary

bait.ff<-freefeed.ph*(avg.bait.effort-5)
# cost of grain for free feeding, per event
# based on average observed effort, minus 5 hours effort for deploying placebos and toxic bait
# PIRSA do 4 days placebo and 1 day toxic baits (discussed with MK).

bait.placebo <-224
# cost is per dispenser
# $14 per placebo bait x 6 per dispenser x 4 days of placebo deployment per dispenser
## 14 X 6 X 4 = 224

bait.toxic <- 156
# cost is per dispenser
# $26 each x 6 per dispenser x 1 day of toxic bait deployment
# 26 x 6 = 156

bait.dispenser.unit<-485

## how many pigs can one dispenser kill in a single event?
# dispensers holds 6 baits
## each bait is 625g
## 100-200g of bait can kill an individual (according to manufacturer)
## bait acts quickly (1-3 hours) so consumption beyond lethal dose is not expected (according to manufacturer)
## so if bait capacity of dispenser is 6 x 625g = 3750g
## 3750g can kill between 18.75 and 37.5 pigs

## manufacturer (ACTA) deployment rate ## 
##recommends 1 dispenser per 20 pigs.

## HOWEVER 
## PIRSA deployment rate##
## PIRSA deploy roughly 1 dispenser for every 3 pigs detected.
## While baits contain enough poison to kill much higher numbers, individuals can consume much more than the lethal dose before being affected
## individual could consume several baits in a few minutes.
## Also, some individuals (e.g. boars and large sows), can prevent other individuals from accessing the bait.
## All up, project used 20 sites with 29 dispensers: average #dispensers per site is 1.45.
## this includes sites that have been omitted as outliers

## how many consecutive events or sites can a single set of dispensers do in a year?
max.bait.events <-avg.bait.effort/365
## based on 12.7 days per baiting event, at 1 hour effort per day, an individual bait dispenser can service 28.74016 events per year
## could round down to allow for staff absences/holidays etc.
## but should be able to bait year round if leave is staggered between multiple staff members.

# total.bait.events.pa <-no.hrs/avg.bait.effort
## this is the number of baiting events required per year to achieve the proportional cull quota

# dispenser.no <- Ceiling((total.bait.events.pa/max.bait.events) *1.45)
## number of bait dispensers required to satisfy annual proportional cull quota
## total.bait.events.pa/max.bait.events gives the number of distinct sites 
## multiplied by 1.45 to give number of dispensers 
## Ceiling () used to round up to nearest integer so whole number of bait dispensers used.

bait.total.cost.pa <- total.bait.events.pa * 
  (bait.labour + bait.ff + (dispenser.no * (bait.placebo + bait.toxic))) + dispenser.no(bait.dispenser.unit)

## but! If fewer dispensers needed in subsequent years, no more dispensers need to be purchased as they can be reused.
## so need to tweak this