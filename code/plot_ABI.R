###
# plot ABI0 and ABImsy through time 
# alongside other indicators of stock status 
# creates stock plots identical to Figure 3 and 4 in Griffiths et al. 2023
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
library(patchwork)
library(tidyverse)

sessionInfo()
# FLCore - FLCore_2.6.19.9023
# patchwork - patchwork_1.1.2

###
# data
load('data/indicator_data.rdata') # indicator data (code provided in calc_ABI.R)
load('data/age_structure_at_Eq.rdata') # FLR stock objects (n = 72) and estimated age structures at equilibrium (code provided in age_structure_at_Eq.R) 

###
# create plots - example for one stock 

# select a stock (her.27.5a and ple.27.21-23 used as examples in MS)
stock_id = 'her.27.5a'
out_dat = filter(indicator_data, stock == stock_id) # indicator values
stock_dat = stks[[stock_id]]  # extract stock assessment FLR object for reference points

# add loess smoothers for ABImsy and ABI0 for visualisation of trends
out_dat$loessMSY <- predict(loess(ABImsy ~ year, data=out_dat, span=0.5))
out_dat$loessNF <- predict(loess(ABI0 ~ year, data=out_dat, span=0.5))

# plots
p1 = ggplot(data = out_dat, aes(x = year, y = log(ABImsy)))+ # ABImsy
  geom_point(fill='firebrick3', size = 3, shape = 21)+
  geom_line(aes(x = year, y = log(loessMSY)),col='firebrick3', size = 1.0)+
  theme_bw()+
  xlab('Year')+
  ylab(expression(paste('log(AB',I[MSY],')')))+
  scale_x_continuous(breaks=seq(min(out_dat$year),max(out_dat$year)+1,3))+
  geom_hline(yintercept=log(1), col='black', size=1.0, linetype = 'dashed')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14))

p2 = ggplot(data = out_dat, aes(x = year, y = log(ABI0)))+ # ABI0
  geom_point(fill='royalblue1', size = 3, shape = 21)+
  geom_line(aes(x = year, y = log(loessNF)),col='royalblue1', size = 1.0)+
  theme_bw()+
  xlab('Year')+
  ylab(expression(paste('log(AB',I[0],')')))+
  scale_x_continuous(breaks=seq(min(out_dat$year),max(out_dat$year)+1,3))+
  geom_hline(yintercept=log(1), col='black', size=1.0, linetype = 'dashed')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14))

p3 = ggplot(data = out_dat, aes(x = year, y = SSB/1000))+ # SSB in 1000 tonnes
  geom_bar(stat="identity", alpha=0.3, fill = 'lightseagreen', col='black')+
  theme_bw()+
  xlab('Year')+
  ylab('Estimated SSB (1000 tonnes)')+
  geom_hline(yintercept=stock_dat@benchmark[['Blim']]/1000, col = 'red', linetype = 'dashed', size=1.0)+
  geom_hline(yintercept=stock_dat@benchmark[['Bpa']]/1000, col = 'black', linetype = 'dashed', size=1.0)+
  geom_hline(yintercept=stock_dat@benchmark[['Btrigger']]/1000, col = 'blue', linetype = 'dashed', size=1.0)+
  scale_x_continuous(breaks=seq(min(out_dat$year),max(out_dat$year)+1,3))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14))

p4 = ggplot(data = out_dat, aes(x = year, y = F))+ # F 
  geom_point(col='gray20', size = 1.5)+
  geom_line(color= "gray20", size = 1.0)+
  theme_bw()+
  xlab('Year')+
  ylab(paste0('F'))+
  geom_hline(yintercept=stock_dat@benchmark[['Flim']], col = 'red', linetype = 'dotted', size=1.0)+
  geom_hline(yintercept=stock_dat@benchmark[['Fpa']], col = 'black', linetype = 'dotted', size=1.0)+
  geom_hline(yintercept=stock_dat@benchmark[['Fmsy']], col = 'blue', linetype = 'dotted', size=1.0)+
  scale_x_continuous(breaks=seq(min(out_dat$year),max(out_dat$year)+1,3))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14))

p5 = ggplot(data = out_dat, aes(x = year, y = R/1000))+ # Recruitment in 1000 individuals
  geom_bar(stat="identity", fill= "mediumblue", alpha=0.3, col='black')+
  theme_bw()+
  xlab('Year')+
  ylab(paste0('Recruitment (age ',stock_dat@range[['min']],'; in 1000s)'))+
  scale_x_continuous(breaks=seq(min(out_dat$year),max(out_dat$year)+1,3))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14))

###
# merge plots and export
p = (p1 + p2) /(p3 + p4 + p5) + plot_layout(guides = 'keep') + plot_annotation(tag_levels = 'A')
ggsave(p, filename = 'output/stock_plot_her.27.5a.png', width = 12, height = 11)
