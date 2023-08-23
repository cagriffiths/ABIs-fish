###
# calculate Amsy and A0
# calculate Pmsy and P0 (threshold = 90th percentile of age distribution)
# calculate ABImsy and ABI0
# authors:
# Christopher Griffiths (SLU)
# Henning Winker (JRC)
# Massimiliano Cardinale (SLU) 

###
# setup
rm(list=ls())
gc()

options(scipen = 99999)

###
# packages - see README for install instructions
library(FLCore)
library(tidyverse)
library(DescTools)
library(viridis)

sessionInfo()
# FLCore - FLCore_2.6.19.9023
# DescTools - DescTools_0.99.47

###
# data
load('data/age_structure_at_Eq.rdata') # FLR stock objects (n = 72) and estimated age structures at equilibrium (code provided in age_structure_at_Eq.R)

# check FLR objects are in the same order
stks = stks[sort(stks@names)]
srs = srs[sort(srs@names)]
eqstks = eqstks[sort(eqstks@names)]

###
# set percentile threshold for Amsy and A0 calc
thres = 0.90

###
# indicator calculations
# ABImsy and ABI0 values by year and stock are also provided in data/indicator_data.rdata

num = NROW(stks@names) # vector of stock IDs for looping

for(j in 1:num){
  
  ##
  # stock age structure at equilibrium
  stock_eq = eqstks[[j]]
  dat_eq = as.data.frame(FLQuants(Numbers=stock.n(stock_eq))) # extract n at age data 
  
  ##
  # calc Amsy and Pmsy
  msy_dat = filter(dat_eq, year == 2) # extract age structure at Fmsy
  msy_dat = filter(msy_dat, !age == min(msy_dat$age)) # remove first age class to limit impact of incoming recruitment
  
  msy_thres = sum(msy_dat$data)*thres # calc n at 90 percentile
  
  msy_dat$cumsum = cumsum(msy_dat$data) # cumsum
  Amsy = filter(msy_dat, cumsum == Closest(msy_dat$cumsum, msy_thres))$age # extract age closest to 90% percentile
  
  if(Amsy == max(msy_dat$age)){ # catch to ensure that Amsy != Amax
    Amsy = Amsy-1
  }
  
  msy_old <- filter(msy_dat, age > Amsy) # filter for above Amsy
  Pmsy <- sum(msy_old$data)/sum(msy_dat$data) # calc proportion above Amsy 

  ##
  # calc A0 and P0
  NF_dat = filter(dat_eq, year == 1) # extract age structure at F0
  NF_dat = filter(NF_dat, !age == min(NF_dat$age)) # remove first age class to limit impact of incoming recruitment
  
  NF_thres = sum(NF_dat$data)*thres # calc n at 90 percentile
  
  NF_dat$cumsum = cumsum(NF_dat$data) # cumsum
  A0 = filter(NF_dat, cumsum == Closest(NF_dat$cumsum, NF_thres))$age # extract age closest to 90% percentile
  
  if(A0 == max(NF_dat$age)){ # catch to ensure that A0 != Amax
    A0 = A0-1
  }
  
  NF_old <- filter(NF_dat, age > A0) # filter for above A0
  P0 <- sum(NF_old$data)/sum(NF_dat$data) # calc proportion above Amsy 
  
  ##
  # stock assessment data
  stock = stks[[j]]
  dat_stock = as.data.frame(FLQuants(Numbers=stock.n(stock)))
  
  ##
  # build stock specific indicator data frame
  ind = data.frame(matrix(NA, ncol=10, nrow=NROW(unique(dat_stock$year))))
  colnames(ind) = c('stock', 'year', 'Amsy', 'A0', 'Pmsy', 'P0','Ptmsy','Pt0','ABImsy', 'ABI0')
  
  ##
  # add fixed variables
  ind$stock = rep(stock@name)
  ind$year = unique(dat_stock$year)
  ind$Amsy = Amsy
  ind$A0 = A0
  ind$Pmsy = Pmsy
  ind$P0 = P0
  
  ##
  # loop over years and calc ABImsy and ABI0
  for(i in 1:NROW(unique(dat_stock$year))){
    year_dat = filter(dat_stock, year == unique(dat_stock$year)[i]) # filter by year
    year_dat = filter(year_dat, !age == min(year_dat$age)) # remove first age class to limit impact of incoming recruitment
    
    year_thres = sum(year_dat$data)*thres # calc n at 90 percentile
    year_dat$cumsum = cumsum(year_dat$data) # cumsum
    
    year_msy = filter(year_dat, age > Amsy) # filter for above Amsy
    year_NF = filter(year_dat, age > A0) # filter for above A0
    
    prop_msy = sum(year_msy$data)/sum(year_dat$data) # proportion above Amsy
    prop_NF = sum(year_NF$data)/sum(year_dat$data) # proportion above A0
    
    ABImsy = prop_msy/Pmsy # calc ABImsy
    ABI0 = prop_NF/P0 # calc A0
    
    ind$Ptmsy[i] = prop_msy # store proportion above Amsy
    ind$Pt0[i] = prop_NF # store proportion above A0
    
    ind$ABImsy[i] = ABImsy # store ABImsy
    ind$ABI0[i] = ABI0 # store ABI0
  }
  
  ##
  # add SSB, F and R
  ind$SSB = as.data.frame(FLQuants(ssb=ssb(stock)))$data # store SSB
  ind$F = as.data.frame(FLQuants(f = fbar(stock)))$data # store F
  ind$R = filter(dat_stock, age == min(dat_stock$age))$data # store R (abundance at age in the smallest age bin)
  
  ##
  # export stock specific indicator data frame
  write.csv(ind, file = paste0('data/indicator data by stock/',ind$stock[1],'.csv'))
  
  ##
  # marker
  print(paste0(ind$stock[1],', no = ',j))
}

