###
# estimate age structure at equilibrium 
# assumes constant fishing at a given F 
# authors:
# Henning Winker (JRC)
# Massimiliano Cardinale (SLU) 
# Christopher Griffiths (SLU)

###
# setup
rm(list=ls())
gc()

options(scipen = 99999)

###
# packages - see README for install instructions
library(FLasher) 
library(FLBRP)
library(tidyverse)
library(viridis)

sessionInfo()
# FLasher - FLasher_0.7.0 (needs to be latest version with ffwd())
# FLBRP - FLBRP_2.5.9.9004

###
# data
load('data/FLR_stock_objects.rdata') # 81 stocks

# check for stocks with Fmsy and Blim values
stk.id = do.call(c,lapply(stks,function(x){
  if(!is.na(x@benchmark[[1]]) & !is.na(x@benchmark[[4]])){ # reference points stored in the @benchmark slot of the FLR stock object
    x@name
    }}))

# remove stocks not listed in stk.id
stks = stks[stk.id] # 72 stocks

###
# SRR for simulation to equilibrium - segmented regression with fixed ICES Blim
srs<- FLSRs(lapply(stks, function(x) { 
  return(fmle(as.FLSR(x,model=segreg),fixed=list(b=x@benchmark[["Blim"]]))) # may take a few minutes to run
  }))
# see https://flr-project.org/doc/Modelling_stock_recruitment_with_FLSR.html for how different SR functions (e.g., bevholt(), ricker() and shepherd()) can be implemented within as.FLSR

###
# estimate age structure at Fmsy and F0 for one stock
# age structure at Fmsy and F0 for all 72 stocks are provided in data/age_structure_at_Eq.rdata
# could also be calculated using the below code by looping over stocks IDs (1,...,72)

# select a stock (hke.27.3a46-8abd, her.27.5a and ple.27.21-23 used as examples in MS)
stock = 'hke.27.3a46-8abd'
stk = stks[[stock]]
sr=srs[[stock]]

# set F values (could be any F)
Fmsy= stk@benchmark[["Fmsy"]]
F0 = 0.000001 # needs to very small and not exactly 0

# compute Eq using FLBRP
bp = FLBRP(stk,sr)
fbar(bp)[,1][] = F0 
fbar(bp)[,2:101][] = Fmsy 

# extract age structure
eqstk = window(as(bp, "FLStock"),end=2) # extract age structure at equilibrium from bp
ageEq = eqstk@stock.n # age structure (numbers at age) at equilibrium under F0 (year 1) and Fmsy (year 2) as a FLQuant object
datEq = as.data.frame(FLQuants(n=stock.n(eqstk))) # coerce to a dataframe

# plot
datEq$F <- factor(datEq$year, levels = c("1", "2"), # add level for F
                  labels = c("F0", "Fmsy"))

g = ggplot(data = datEq, aes(x = as.factor(age), y = data/1000, fill = F))+
  geom_bar(stat = 'identity', col = 'black', alpha = 0.4)+
  theme_bw()+
  xlab('Age')+
  ylab('Numbers-at-age (1000s individuals)')+
  scale_fill_viridis(discrete = T, option = 'turbo')+
  ggtitle(expression(paste('Age structure at equilibrium')))+
  scale_y_continuous(limits = c(0,900), breaks = seq(0,900,100))+
  facet_wrap(~F)+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14),
        legend.position = '')
g

# save
write.csv(datEq, file = 'output/age_structure_Eq_hke.27.3a46-8abd.csv')
ggsave(g, filename = 'output/age_structure_Eq_hke.27.3a46-8abd.png')
