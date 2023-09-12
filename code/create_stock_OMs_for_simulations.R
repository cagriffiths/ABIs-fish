###
# code to create operating models for simulation testing
# considers 6 stocks: pil.27.8c9a, pra.27.3a4a, her.27.3a47d, ple.27.420, cod.27.5a and mon.27.8c9a
# authors:
# Henning Winker (JRC)
# email: henning.winker@gmail.com

###
# setup
rm(list=ls())
gc()

options(scipen = 99999)

###
# packages - see README for install instructions
library(FLCore)
library(ggplotFL)
library(mse)
library(FLRef)
library(FLBRP)
library(FLSRTMB)

sessionInfo()
# FLCore - FLCore_2.6.19.9075
# ggplotFL - ggplotFL_2.7.0 
# mse - mse_2.2.3
# FLRef - FLRef_1.0.9 
# FLBRP - FLBRP_2.5.9.9004
# FLSRTMB - FLSRTMB_1.1.4

###
# data
load('data/FLR_stock_objects.rdata') # 81 stocks

# check stock names
do.call(rbind,lapply(stks,function(x) x@name)) == stks@names 

# add life history attributes for pandalus, some of which were taken from FishLife 
# and others taken from the 2022 assessment of the stock
attr(stks$pra.27.3a4a,"fishlife") = c(linf = 2.80, k = 0.48, winf = 0.030, tmax = 8,
                                      tm = 2, m = 0.75, lm = 1.787, rho = 0.26,
                                      sigmaR = 0.43, s = 0.71, fmsy = NA,
                                      r = NA, g = NA, sd.logit.s = 0.3)

###
# Stocks for simulation testing
stksel = c( # vector of selected stocks
  "pil.27.8c9a",
  "pra.27.3a4a",
  "her.27.3a47d",
  "ple.27.420",
  "cod.27.5a",
  "mon.27.8c9a"
)

# filter stock objects by selected stocks
stks = stks[stksel]

###
# generate properties for the simulation tests

# Make sure age-0 maturity (recruits) is 0
mat(stks$mon.27.8c9a)[ac(0),] = 0

# function for generic F trajectory generation
ftrj = function(fmsy,fhi=2.5,flo=0.8,y0=125,y1=150,y2=175,steps=201){
  x = seq(1,steps,1)
  f = rep(fmsy,steps)
  f[x > y0 & x <= y1] = ((fhi*fmsy-fmsy)/(y1-y0))*(x[x > y0 & x <= y1] - y0) +fmsy
  f[x > y1 & x <= y2] = (-(fhi*fmsy-flo*fmsy)/(y2-y1))*(x[x > y1 & x <= y2]-y1) +fhi*fmsy
  f[x > y2] = fmsy*flo
  return(f)
}

# Check F trajectory
plot(ftrj(fmsy=0.2)[101:200],type="l")

# Condition a BevHolt functions with r0 prior  
srs = FLSRs(lapply(stks,function(x){
  r0 = x@eqsim[["B0"]]/mean(spr0y(x))
  # condition B-H SRR on R0 (scale) and steepness fishlife
  bh = srrTMB(as.FLSR(x,model=bevholtSV),spr0=mean(spr0y(x)),nyears=NULL,SDreport = F,r0.pr=c(r0,0.1),s=x@fishlife[["s"]],s.logitsd =x@fishlife[["sd.logit.s"]])
  # Create FLBRP 
  bh@name = x@name
  bh
}))

# Check SRfor mon.27.8c9a
plotsrs(srs$mon.27.8c9a,b0=T)

# Compute reference points
brps = FLBRPs(Map(function(x,y){
# Create FLBRP 
br = brp(FLBRP(x,y))
br@name = x@name
br
},x=stks,y=srs))

# extract refpts
rpts = FLPars(lapply(brps,function(x){
  # Get true refpts
  FLPar(
  fmsy = refpts(x)["msy","harvest"],
  bmsy = refpts(x)["msy","ssb"],
  msy = refpts(x)["msy","yield"],
  r0 = refpts(x)["msy","rec"] # changed to msy
  )
}))

###
# simulate stock OMs
iters = 1000 
stksims = Map(function(x,y){

# transform back into a stock
stk0 = stf(as(x,"FLStock"), 100) # extent by 100 yrs
nyears = length(dims(stk0)$minyear:dims(stk0)$maxyear)# Create iterations
stki = propagate(stk0,iters)

# create random rec devs
devs = FLQuants(devs04 = ar1rlnorm(rho=0, years=1:nyears, iters=iters,meanlog= 0, sdlog=0.4),
                devs07 = ar1rlnorm(rho=0, years=1:nyears, iters=iters,meanlog= 0, sdlog=0.7))


# run simulation test across the 6 selected stocks - might take 5/10 mins to run
sts = FLStocks(lapply(devs,function(z){
  bh = as(x,"FLSR")
  ft = ftrj(an(y["fmsy"]))
  sim = ffwd(stki, sr=bh,
               fbar=FLQuant(ft, dimnames=list(year=2:201)),
               deviances=z)
  sim@landings = computeLandings(sim) # manually compute landings
  sim@name = paste0(x@name)
  sim = window(sim,start=101)
  return(sim)
}))
return(sts)
},x=brps,y=rpts)

# save simulated data in data
# save(stksims,brps,srs,rpts,file="data/sim6stks.rdata") # fairly large .rdata file [~450 MB] and exceeds limits for GitHub

###
# checking and plots for a single stock
# simple/basic plot
plot(stksims$pra.27.3a4a) 

# pretty/complex plot
stksim = stksims$pra.27.3a4a$devs04
fbrp = brps$pra.27.3a4a
thinning = 10
nyears=101
probs=c(0.05,0.2,0.50,0.8,0.95)
colour  = "dodgerblue"
p <- ggplotFL::plot(stksim,metrics=list(Rec=rec,SSB=ssb,Landings=landings,F=fbar)[1:4],iter=seq(1,iters,thinning))+scale_color_manual(values=c(grey.colors(length(seq(1,iters,thinning)))))+
  theme_bw()+xlab("Year")+theme(legend.position = "none")

p <- p +ggplotFL::geom_flquantiles(fill=colour, probs=probs[c(1,3,5)], alpha=0.4) +
  ggplotFL::geom_flquantiles(fill=colour, probs=probs[c(2,3,4)], alpha=0.5)+
  facet_wrap(~qname, scales="free",ncol=2)
# add refpts
det.refs = FLPars(Rec=FLPar(R0=refpts(fbrp)["virgin","rec"]),
                  SSB=FLPar(Bmsy=refpts(fbrp)["msy","ssb"],B0=refpts(fbrp)["virgin","ssb"]),
                  Landings=FLPar(MSY = refpts(fbrp)["msy","yield"]),
                  F=FLPar(Fmsy=refpts(fbrp)["msy","harvest"]))
ggp = ggplotFL::geom_flpar(data=det.refs,
                            x=c(1.5*nyears,rep(1.5*nyears,2),1.5*nyears,1.5*nyears),
                            colour = c("blue",c("darkgreen","blue"),("darkgreen"),c("darkgreen")))
ggp[[2]]$aes_params$size=3
p + ggp+ggtitle(stksim@name,"Two-way trip")

# save plot in output
ggsave(p, filename = 'output/sim_plot_pra.27.3a4a.png')


  