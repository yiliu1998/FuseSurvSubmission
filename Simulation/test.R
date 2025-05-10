covar.name=c("X1","X2","X3")
site.var="site"
trt.name="A"
time.var="Y"
event="Delta"
fit.times=1:90
eval.times=c(30,60,90)
prop.SL.library=c("SL.mean", "SL.glm", "SL.glm.interaction")
event.SL.library=c("survSL.km", "survSL.coxph", "survSL.weibreg") 
cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.weibreg")
n.folds=5 
n.fed.boot=200
s=1

library(CFsurvival)
library(survSuperLearner)
library(SuperLearner)
library(dplyr)
library(glmnet)
library(caret)
source("baseFuncs.R")

load("obsdata_l2.Rdata")
data = dat.homo[[2]]
