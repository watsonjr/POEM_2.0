################################################################################

# Calc correlation between SAUP and POEM catch
# Historic pristine

################################################################################

rm(list=ls())

library(chron)
library(lattice)
library(ggplot2)
library(ggbiplot)
library(gridExtra)
library(corrgram)
library(Hmisc)
library(PerformanceAnalytics)

#PU laptop
source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/")

#UAF laptop
#source(file = "/Users/Colleen/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
#setwd("/Users/Colleen/Dropbox/Princeton/POEM_other/SAUP/")

cpath = "/Volumes/GFDL/NC/"
cfile = "Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05"

### load data
lmec <- read.csv(paste(cpath,cfile,"/LME_saup_catch_hist_fished_",cfile,".csv",sep=""),sep=",",header=T,stringsAsFactors=F)
lmec$logS <- log10(lmec$saup)
lmec$logP <- log10(lmec$poem)

#---------------------------------------- CORR TABLE -------------------------------------
mydata1 <- as.matrix(lmec[,2:3])
Rcorr <- rcorr(mydata1)
(Rcorr$P)
(Rcorr$r)
pdf("Hist_fished_SAUP_chart_corr.pdf")
chart.Correlation(mydata1, histogram=FALSE, pch=19)
dev.off()

ni <- which(lmec$logS!=-Inf) 
mydata2 <- as.matrix(lmec[ni,4:5])
Rcorr2 <- rcorr(mydata2)
(Rcorr2$P)
(Rcorr2$r)
pdf("Hist_fished_SAUP_chart_corr_log10.pdf")
chart.Correlation(mydata2, histogram=FALSE, pch=19)
dev.off()


### PLOTS ##################################################################
# R -------------------------------------------------------------------------
#Scatter plots
p1 <- ggplot(lmec, aes(y=poem, x=saup)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  ylab("POEM catch (MT)") + 
  xlab("SAUP catch (MT)")
p1
pdf("Hist_fished_SAUP_scatter.pdf")
p1
dev.off()

p2 <- ggplot(lmec, aes(y=logP, x=logS)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  ylab("log10 POEM catch (MT)") + 
  xlab("log10 SAUP catch (MT)")
p2
pdf("Hist_fished_SAUP_scatter_log10.pdf")
p2
dev.off()

