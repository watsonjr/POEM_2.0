# Calc recruitment variability from POEM output
# Heath's PV

rm(list=ls())

library(Hmisc)
library(pastecs)
library(fBasics)
require(combinat)

setwd("/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ/")

# load data
f <- read.csv("onelocs_hist_forage_larv.csv",sep=",", header = T,
               stringsAsFactors = F)
d <- read.csv("onelocs_hist_dem_larv.csv",sep=",", header = T,
              stringsAsFactors = F)
p <- read.csv("onelocs_hist_pisc_larv.csv",sep=",", header = T,
              stringsAsFactors = F)
nid <- ncol(f)

area <- read.csv("/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/grid_360x200_id_locs_area.csv")
A <- 0
A[1:2] <- area$Area[1]
A[3:4] <- area$Area[2]
A[5:6] <- area$Area[3]
A[7:8] <- area$Area[4]
A[9:10] <- area$Area[5]
A[11:12] <- area$Area[6]

Rstats <- as.data.frame(names(d))
for(i in 1:nid){
  #Population variability metric of Heath 2006
  #PV of R
  
  #Convert from biomass to abundance
  dnum <- round(d[,i] / 2529.822128134704 * A[i]) + 1
  pnum <- round(p[,i] / 2529.822128134704 * A[i]) + 1
  fnum <- round(f[,i] / 2.5298221281347053 * A[i]) + 1
  
  # FERNANDO'S WAY
  P <- log(na.omit(fnum));
  if (length(P)!=0){
    C <- combn( P, 2 ) # all pairwise combinations of abundances
    Z <- dim( C )[2] # number of combinations
    # operate over all combinations of pairs of observations
    Num <- abs( C[1, ] - C[2,] ) # absolute value of difference
    Denom <- apply( C, 2, max )  # maximum value
    Diff <- Num / Denom          # proportional difference
    PV <- mean(Diff, na.rm=T)
    Rstats[i,2] <- PV
    rm(list=c('P','C','Z','Num','Denom','Diff','PV'))
  }
  
  P <- log(na.omit(pnum))
  if (length(P)!=0){
    C <- combn( P, 2 ) # all pairwise combinations of abundances
    Z <- dim( C )[2] # number of combinations
    # operate over all combinations of pairs of observations
    Num <- abs( C[1, ] - C[2,] ) # absolute value of difference
    Denom <- apply( C, 2, max )  # maximum value
    Diff <- Num / Denom          # proportional difference
    PV <- mean(Diff, na.rm=T)
    Rstats[i,3] <- PV
    rm(list=c('P','C','Z','Num','Denom','Diff','PV'))
  }
  
  P <- log(na.omit(dnum))
  if (length(P)!=0){
    C <- combn( P, 2 ) # all pairwise combinations of abundances
    Z <- dim( C )[2] # number of combinations
    # operate over all combinations of pairs of observations
    Num <- abs( C[1, ] - C[2,] ) # absolute value of difference
    Denom <- apply( C, 2, max )  # maximum value
    Diff <- Num / Denom          # proportional difference
    PV <- mean(Diff, na.rm=T)
    Rstats[i,4] <- PV
    rm(list=c('P','C','Z','Num','Denom','Diff','PV'))
  }
  
}

names(Rstats) <- c("ID","F","P","D")

write.table(Rstats,"onelocs_hist_recruitment_PV.csv",row.names=F,sep=",")

