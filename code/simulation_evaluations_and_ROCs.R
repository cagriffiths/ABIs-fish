###
# code for performance evaluation of ABI
# considers 6 stocks: pil.27.8c9a, pra.27.3a4a, her.27.3a47d, ple.27.420, cod.27.5a and mon.27.8c9a
# generates ROC curves and assiocated metrics (TSS and AUCs)
# also includes the code used for the threshold testing (percentile of age distribution)
# authors:
# Henning Winker (JRC)
# email: henning.winker@gmail.com

###
# setup
rm(list=ls())
gc()

options(scipen = 99999)

###
# Install packages
#----------------------------------------------
# install.packages("TMB", type="source")
# install.packages(c('FLCore','ggplotFL', 'FLBRP', 'mse', 'FLSRTMB'), repo='https://flr.r-universe.dev')
# devtools::install_github("henning-winker/FLRef")
# install.packages("pROC")
#----------------------------------------------

###
# packages - see README for install instructions
library(FLCore)
library(ggplotFL)
library(mse)
library(FLRef)
library(FLBRP)
library(FLSRTMB)
library(pROC) # For ROC curve calculations

sessionInfo()
# FLCORE - FLCore_2.6.19.9075
# ggplotFL - ggplotFL_2.7.0
# mse - mse_2.2.3
# FLRef - FLRef_1.0.9
# FLBRP - FLBRP_2.5.9.9004
# FLSRTMB - FLSRTMB_1.1.4
# pROC - pROC_1.18.0

###
# data
load('data/sim6stks.rdata') # 6 selected stocks

# few quick checks
stk = stksims$cod.27.5a$devs04
rpt = rpts$cod.27.5a
plot(iter(stk,1))

###
# compute Amsy and Pmsy (i.e., ref prop of numbers above age thresholds)
computeRP <- function(brp,thresh=0.9){
  eqstk = brp
  fbar(eqstk)[,1][] = 0.00001 # compute for eq Fmsy
  fbar(eqstk)[,1:101][] = an(refpts(brp)["msy","harvest"]) # compute for equilibrium Fmsy
  eqstk = window(as(eqstk, "FLStock"),start=2,end=2) # year1 F=0, year2=Fmsy
  eqstk@name = stk@name # name stk
  n_a = stock.n(eqstk)[-1,] # remove first age
  ages = dims(n_a)$min:dims(n_a)$max
  cums = apply(n_a,2:6,cumsum)
  n_thresh = sum(n_a*thresh)
  aref = min(ages[which((n_thresh-cums)^2==min((n_thresh-cums)^2))]+1,range(eqstk)["plusgroup"]-1)
  rp = sum(n_a[ac(aref:max(ages)),])/sum(n_a)
  return(list(rp=rp,aref=aref,bmsy=an(refpts(brp)["msy","ssb"])))
}

###
# plotting and data manipulation
# make a larger dataset for ggplot
dfmu = do.call(rbind,Map(function(x,y,z){
# compute RP
RP = computeRP(z)  
mu = rbind(
data.frame(as.data.frame(FLQuants(
F=iterMeans(fbar(x[[1]])/y[[1]]),
Yield=iterMeans(landings(x[[1]])/y[[3]]),
Recruits = iterMeans(rec(x[[1]])/y[[4]]),
SSB=iterMeans(ssb(x[[1]])/y[[2]]),
"ABI[MSY]" = iterMeans(quantSums(stock.n(x[[1]])[ac(RP$aref:range(x[[1]])[2]),])/quantSums(stock.n(x[[1]])[-1,]))/RP$rp
)),stock=x[[1]]@name, sigR="0.4"),

data.frame(as.data.frame(FLQuants(
  F=iterMeans(fbar(x[[2]])/y[[1]]),
  Yield=iterMeans(landings(x[[2]])/y[[3]]),
  Recruits = iterMeans(rec(x[[2]])/y[[4]]),
  SSB=iterMeans(ssb(x[[2]])/y[[2]]),
  "ABI[MSY]" = iterMeans(quantSums(stock.n(x[[2]])[ac(RP$aref:range(x[[2]])[2]),])/quantSums(stock.n(x[[2]])[-1,]))/RP$rp
  )),stock=x[[2]]@name,sigR="0.7")
)

lci = rbind(
  data.frame(as.data.frame(FLQuants(
    F=quantile(fbar(x[[1]])/y[[1]],0.05),
    Yield=quantile(landings(x[[1]])/y[[3]],0.05),
    Recruits = quantile(rec(x[[1]])/y[[4]],0.05),
    SSB=quantile(ssb(x[[1]])/y[[2]],0.05),
    "ABI[MSY]" = quantile(quantSums(stock.n(x[[1]])[ac(RP$aref:range(x[[1]])[2]),])/quantSums(stock.n(x[[1]])[-1,]),0.05)/RP$rp
  )),stock=x[[1]]@name, sigR="0.4"),
  
  data.frame(as.data.frame(FLQuants(
    F=quantile(fbar(x[[2]])/y[[1]],0.05),
    Yield=quantile(landings(x[[2]])/y[[3]],0.05),
    Recruits = quantile(rec(x[[2]])/y[[4]],0.05),
    SSB=quantile(ssb(x[[2]])/y[[2]],0.05),
    "ABI[MSY]"  = quantile(quantSums(stock.n(x[[2]])[ac(RP$aref:range(x[[2]])[2]),])/quantSums(stock.n(x[[2]])[-1,]),0.05)/RP$rp
  )),stock=x[[2]]@name,sigR="0.7")
)

uci = rbind(
  data.frame(as.data.frame(FLQuants(
    F=quantile(fbar(x[[1]])/y[[1]],0.95),
    Yield=quantile(landings(x[[1]])/y[[3]],0.95),
    Recruits = quantile(rec(x[[1]])/y[[4]],0.95),
    SSB=quantile(ssb(x[[1]])/y[[2]],0.95),
    "ABI[MSY]" = quantile(quantSums(stock.n(x[[1]])[ac(RP$aref:range(x[[1]])[2]),])/quantSums(stock.n(x[[1]])[-1,]),0.95)/RP$rp
  )),stock=x[[1]]@name, sigR="0.4"),
  
  data.frame(as.data.frame(FLQuants(
    F=quantile(fbar(x[[2]])/y[[1]],0.95),
    Yield=quantile(landings(x[[2]])/y[[3]],0.95),
    Recruits = quantile(rec(x[[2]])/y[[4]],0.95),
    SSB=quantile(ssb(x[[2]])/y[[2]],0.95),
    "ABI[MSY]"  = quantile(quantSums(stock.n(x[[2]])[ac(RP$aref:range(x[[2]])[2]),])/quantSums(stock.n(x[[2]])[-1,]),0.95)/RP$rp
  )),stock=x[[2]]@name,sigR="0.7")
)

data.frame(mu,lci=lci$data,uci=uci$data)

},x=stksims,y=rpts,z=brps))

# show first iteration
it = 1
dfiter = do.call(rbind,Map(function(x,y,z){
  RP = computeRP(z)  
  rbind(
    data.frame(as.data.frame(FLQuants(
      F=iter(fbar(x[[1]])/y[[1]],it),
      Yield=iter(landings(x[[1]])/y[[3]],it),
      Recruits = iter(rec(x[[1]])/y[[4]],it),
      SSB=iter(ssb(x[[1]])/y[[2]],it),
      "ABI[MSY]"= iter(quantSums(stock.n(x[[2]])[ac(RP$aref:range(x[[1]])[2]),])/quantSums(stock.n(x[[1]])[-1,]),it)/RP$rp
      )),stock=x[[1]]@name, sigR="0.4"),
    data.frame(as.data.frame(FLQuants(
      F=iter(fbar(x[[2]])/y[[1]],it),
      Yield=iter(landings(x[[2]])/y[[3]],it),
      Recruits = iter(rec(x[[2]])/y[[4]],it),
      SSB=iter(ssb(x[[2]])/y[[2]],it),
      "ABI[MSY]" = iter(quantSums(stock.n(x[[2]])[ac(RP$aref:range(x[[2]])[2]),])/quantSums(stock.n(x[[2]])[-1,]),it)/RP$rp
      )),stock=x[[1]]@name,sigR="0.7")
  )
  
},x=stksims,y=rpts,z=brps))

# adjust factor levels for plotting order
dfmu$stock = factor(dfmu$stock,levels=unique(dfmu$stock)) 
dfiter$stock = factor(dfiter$stock,levels=unique(dfiter$stock)) 
dfmu$sigR = factor(dfmu$sigR,levels=c("0.7","0.4"))
dfiter$sigR = factor(dfiter$sigR,levels=c("0.7","0.4"))

# change axis label to ABI[MSY]
flab = as_labeller(unique(dfmu$qname),
            default = label_parsed)

# plot
p_sim = ggplot(dfmu,aes(year,data,col=sigR))+facet_grid(stock~qname,labeller = flab)+theme_bw()+
  geom_ribbon(data=dfmu,aes(ymin = lci, ymax = uci,fill=sigR),alpha=0.3,color=0)+
  geom_line(data=dfiter,aes(year,data,group =iter),alpha=1,linewidth=0.4)+
  geom_hline(yintercept = 1,linetype=1)+
  geom_hline(yintercept = 0.8,linetype=2)+
  geom_line(data=dfmu[dfmu$sigR==ac(0.4),],size=0.8,col=1)+ylab("Ratio")+
  labs(color=expression(sigma[R]),fill=expression(sigma[R]))+
  scale_color_manual(values=c("brown3","blue"))+
  scale_fill_manual(values=c("brown3","blue"))+xlab("Year")+
  theme(axis.text.x=element_blank())

###
# generate ROC curves and plot
it = 1:1000
thresh = 0.8 # Bmsy threshold
roc_inp = do.call(rbind,Map(function(x,y,z){
  RP = computeRP(z)  
  
df1 = data.frame(
    as.data.frame(iter(quantSums(stock.n(x[[1]])[ac(RP$aref:range(x[[1]])[2]),])/quantSums(stock.n(x[[1]])[-1,]),it)/(thresh*RP$rp)),
    label=as.data.frame(iter(ssb(x[[1]])/(thresh*y[[2]]),it)<1)$data,
    stock=x[[1]]@name, sigR="0.4"
  )
df2=  data.frame(
    as.data.frame(iter(quantSums(stock.n(x[[2]])[ac(RP$aref:range(x[[2]])[2]),])/quantSums(stock.n(x[[2]])[-1,]),it)/(thresh*RP$rp)),
    label=as.data.frame(iter(ssb(x[[2]])/(thresh*y[[2]]),it)<1)$data,
    stock=x[[1]]@name, sigR="0.7"
  )
  
  df1$label = ifelse(df1$label==TRUE,1,0)
  r1 = roc(df1$label,df1$data)  
  pts = coords(r1)
  pts$tss = pts$sensitivity+pts$specificity-1
  pt = pts[which(abs(pts$threshold-1)==min(abs(pts$threshold-1)))[1],]
  pt$auc = an(auc(r1))
  roc1 = NULL
  roc1$curve = data.frame(fpr = 1-r1$specificities,tpr=r1$sensitivities,thresh=pts$threshold,tss=pts$sensitivity+pts$specificity-1)
  roc1$pt = pt 
  df2$label = ifelse(df2$label==TRUE,1,0)
  r2 = roc(df2$label,df2$data)  
  pts = coords(r2)
  pts$tss = pts$sensitivity+pts$specificity-1
  pt = pts[which(abs(pts$threshold-1)==min(abs(pts$threshold-1)))[1],]
  pt$auc = an(auc(r2))
  roc2 = NULL
  roc2$curve = data.frame(fpr = 1-r2$specificities,tpr=r2$sensitivities,thresh=pts$threshold,tss=pts$sensitivity+pts$specificity-1)
  roc2$pt = pt 
  
rbind(
    data.frame(stock=x[[1]]@name,sigR="0.4",roc1$curve,x=1-roc1$pt$specificity,y=roc1$pt$sensitivity,TSS=roc1$pt$tss,AUC=roc1$pt$auc),
    data.frame(stock=x[[2]]@name,sigR="0.7",roc2$curve,x=1-roc2$pt$specificity,y=roc2$pt$sensitivity,TSS=roc2$pt$tss,AUC=roc2$pt$auc)
  )
  
},x=stksims,y=rpts,z=brps))

roc_inp$stock = factor(roc_inp$stock,levels=unique(roc_inp$stock)) 
roc_inp$sigR = factor(roc_inp$sigR,levels=rev(unique(roc_inp$sigR))) 

pdat = aggregate(cbind(x,y)~stock+sigR, roc_inp,mean)
pdat$x=0
pdat$y=1

df = aggregate(cbind(x,y,AUC)~stock+sigR, roc_inp,mean)
df$TSS = round(1-df$x+df$y-1,2)
df$AUC = round(df$AUC,2)

# plot denoting TRUE SKILL SCORE (TSS)
p_roc = ggplot(data=roc_inp,aes(fpr,tpr,col=sigR))+facet_wrap(~stock,nrow=3)+theme_bw()+
  geom_line(cex=0.7)+geom_abline()+geom_point(data=df,
 aes(x,y,col=sigR),size=3)+
  geom_point(data=pdat,aes(x,y),pch=8,col=1,size=2)+
  geom_label(data=df,aes(x,y,label=TSS), size=2,show.legend =F)+
  ylab("True Positive Rate")+xlab("False Positive Rate")+
  labs(color=expression(sigma[R]))+
  scale_color_manual(values=c("brown3","blue"))+ggtitle("TSS values")

# alternative plot denoting AREA UNDER THE CURVE (AUC)
p_rocAUC = ggplot(data=roc_inp,aes(fpr,tpr,col=sigR))+facet_wrap(~stock,nrow=3)+theme_bw()+
  geom_line(cex=0.7)+geom_abline()+geom_point(data=df,
                                              aes(x,y,col=sigR),size=3)+
  geom_point(data=pdat,aes(x,y),pch=8,col=1,size=2)+
  geom_label(data=df,aes(x,y,label=AUC), size=2,show.legend =F)+
  ylab("True Positive Rate")+xlab("False Positive Rate")+
  labs(color=expression(sigma[R]))+
  scale_color_manual(values=c("brown3","blue"))+ggtitle("AUC values")

###
# save plots
# simulation plot
ggsave(p_sim, filename = 'output/simplot.png', width = 7, height = 8)

# ROCs with TSS 
ggsave(p_roc, filename = 'output/ROCs_with_TSS.png', width = 7, height = 8)

# ROCs with AUC
ggsave(p_rocAUC, filename = 'output/ROCs_with_AUCs.png', width = 7, height = 8)

################################################################################
###
# threshold experiment for optimizing cut-off age for ABImsy

# threshold range tested
threshs = c(0.5,0.6,0.7,0.8,0.85,0.9,0.95)

# get means for ggplotting across tresholds
dfthr = do.call(rbind,Map(function(x,y,z){
  
# compute reference points
RP = computeRP(z)  
do.call(rbind,lapply(threshs,function(t){
  RP = computeRP(z,t)  
    data.frame(as.data.frame(FLQuants(
      SSB=iterMeans(ssb(x[[1]])/y[[2]]),
      "ABI[MSY]" = iterMeans(quantSums(stock.n(x[[1]])[ac(RP$aref:range(x[[1]])[2]),])/quantSums(stock.n(x[[1]])[-1,]))/RP$rp
    )),stock=x[[1]]@name, sigR="0.4",threshold=t)
  }))},x=stksims,y=rpts,z=brps))

# colour generator for plots
rc4 = function (n, alpha = 1) 
{
  n = n+1
  x <- seq(0, 1, length = n)
  r <- 1/(1 + exp(20 - 35 * x))
  g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
  b <- dnorm(x, 0.25, 0.15)/max(dnorm(x, 0.25, 0.15))
  rgb.m <- matrix(c(r, g, b), ncol = 3)
  col  <- apply(rgb.m, 1, function(v) rgb(v[1], v[2], v[3], alpha = alpha))[-1]
  cols <- adjustcolor(col, alpha.f = alpha)  
}

# organise levels for plotting
dfthr$stock = factor(dfthr$stock,levels=unique(dfthr$stock))
dfthr$Threshold = factor(dfthr$threshold)

# plot comparing thresholds
p_thr = ggplot(dfthr[dfthr$qname=="ABI[MSY]",],aes(year,data,color=Threshold))+theme_bw()+
  facet_wrap(~stock,ncol=2)+
  geom_line(data=dfthr[dfthr$qname=="SSB",],aes(year,data,linetype="SSB"),color=1,size=0.9)+
  geom_line(size=0.8)+geom_hline(yintercept = 1,linetype=1)+ylab(expression(ABI[MSY]))+
  theme(axis.text.x=element_blank(),legend.title =element_blank())+
  scale_color_manual(values=rev(c(rc4(7,alpha=0.8))))

# save plot
ggsave(p_thr, filename = 'output/threshold_test.png', width = 7, height = 8)


        