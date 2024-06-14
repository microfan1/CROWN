library(openxlsx)
library(R.matlab)

#Remember to change this path
setwd("E:/elective/RA/CUHK/Tutorial/code/methods/POET/")
source('./POET/R/zPOET.R')
Y<-t(readMat("data_temp.mat")$Return)
K<- POETKhat(Y)$K2HL
L<-POET(Y,K,0.5,'soft','vad')
Sy<-L$SigmaY
Su<-L$SigmaU
fac<-t(L$factors)
fload<-L$loadings
writeMat("data_temp.mat",Sy=Sy,Su=Su,F=fac,B=fload)
