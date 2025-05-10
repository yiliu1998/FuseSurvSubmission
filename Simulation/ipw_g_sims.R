library(survival)
library(survey)
library(dplyr)

newdat.ipw <- data.frame(A=c(0,1))
method <- c("IPW TGT", "IPW POOL", "G TGT", "G POOL")
get_surv_probs <- function(sf, times) {
  s <- summary(sf, times = times)
  mat <- matrix(s$surv, nrow = length(times))
  rowMeans(mat, na.rm = TRUE)
}

ipw_g_sim <- function(dat, M=500) {
  result <- list()
  for(i in 1:M) {
    dat.i <- dat[[i]]
    dat.i.tgt <- dat.i[dat.i$site==0,]
    
    # IPW TGT
    ps <- predict(glm(A~X1+X2+X3, data=dat.i.tgt, family=binomial), type="response")
    dat.i.tgt$w <- ifelse(dat.i.tgt$A==1, 1/ps, 1/(1-ps))
    cox_ipw.tgt <- coxph(Surv(Y, Delta) ~ A, data=dat.i.tgt, weights=w)
    surv.ipw.tgt <- summary(survfit(cox_ipw.tgt, newdata=newdat.ipw), times=90)$surv
    
    # IPW POOL
    ps <- predict(glm(A~X1+X2+X3, data=dat.i, family=binomial), type="response")
    dat.i$w <- ifelse(dat.i$A==1, 1/ps, 1/(1-ps))
    cox_ipw <- coxph(Surv(Y, Delta) ~ A, data=dat.i, weights=w)
    surv.ipw.pool <- summary(survfit(cox_ipw, newdata=newdat.ipw), times=90)$surv
    
    # G TGT
    cox.g.tgt <- coxph(Surv(Y, Delta)~A+X1+X2+X3, data=dat.i.tgt)
    newdat0 <- dat.i.tgt %>% mutate(A=0)
    newdat1 <- dat.i.tgt %>% mutate(A=1)
    sf0 <- survfit(cox.g.tgt, newdata=newdat0)
    sf1 <- survfit(cox.g.tgt, newdata=newdat1)
    surv.g.tgt <- c(get_surv_probs(sf0, times=90), 
                    get_surv_probs(sf1, times=90))
    
    # G POOL
    cox.g <- coxph(Surv(Y, Delta)~A+X1+X2+X3, data=dat.i.tgt)
    newdat0 <- dat.i %>% mutate(A=0)
    newdat1 <- dat.i %>% mutate(A=1)
    sf0 <- survfit(cox.g, newdata=newdat0)
    sf1 <- survfit(cox.g, newdata=newdat1)
    surv.g.pool <- c(get_surv_probs(sf0, times=90), 
                    get_surv_probs(sf1, times=90))
    
    result[[i]] <- data.frame(method, 
                              rbind(surv.ipw.tgt, surv.ipw.pool, surv.g.tgt, surv.g.pool))
    if(i%%50==0) print(i)
  }
  return(result=result)
}

load("obsdata_l.Rdata")
result.homo.ipwg <- ipw_g_sim(dat=dat.homo)
result.diffX.ipwg <- ipw_g_sim(dat=dat.diffX)
result.diffT.ipwg <- ipw_g_sim(dat=dat.diffT)
result.diffC.ipwg <- ipw_g_sim(dat=dat.diffC)
result.diffAll.ipwg <- ipw_g_sim(dat=dat.diffAll)

save(file="simRes_ipw_g.Rdata", 
     result.homo.ipwg, result.diffX.ipwg, result.diffT.ipwg, 
     result.diffC.ipwg, result.diffAll.ipwg)

### Result summary and comparison
load("truth.Rdata")
load("simRes_ipw_g.Rdata")
load("Res_homo_l.Rdata")
load("Res_diffX_l.Rdata")
load("Res_diffT_l.Rdata")
load("Res_diffC_l.Rdata")
load("Res_diffAll_l.Rdata")

Bias.results <- function(datlist, datlist.ipwg, M=500) {
  TGT.1 <- IVW.1 <- POOL.1 <- CCOD.1 <- FED.1 <- FED.bt.1 <- IPW.TGT.1 <- IPW.POOL.1 <- G.TGT.1 <- G.POOL.1 <- c()
  TGT.0 <- IVW.0 <- POOL.0 <- CCOD.0 <- FED.0 <- FED.bt.0 <- IPW.TGT.0 <- IPW.POOL.0 <- G.TGT.0 <- G.POOL.0 <- c()
  for(j in 1:M) {
    TGT.1 <- c(TGT.1, (datlist[[j]]$df.TGT[,"surv1"] - S1.true))
    IVW.1 <- c(IVW.1, (datlist[[j]]$df.IVW[,"surv1"] - S1.true))
    POOL.1 <- c(POOL.1, (datlist[[j]]$df.POOL[,"surv1"] - S1.true))
    CCOD.1 <- c(CCOD.1, (datlist[[j]]$df.CCOD[,"surv1"] - S1.true))
    FED.1 <- c(FED.1, (datlist[[j]]$df.FED[,"surv1"] - S1.true))
    FED.bt.1 <- c(FED.bt.1, (datlist[[j]]$df.FED.bt[,"surv1"] - S1.true))
    IPW.TGT.1 <- c(IPW.TGT.1, (datlist.ipwg[[j]][1,3] - S1.true))
    IPW.POOL.1 <- c(IPW.POOL.1, (datlist.ipwg[[j]][2,3] - S1.true))
    G.TGT.1 <- c(G.TGT.1, (datlist.ipwg[[j]][3,3] - S1.true))
    G.POOL.1 <- c(G.POOL.1, (datlist.ipwg[[j]][4,3] - S1.true))
    
    TGT.0 <- c(TGT.0, (datlist[[j]]$df.TGT[,"surv0"] - S0.true))
    IVW.0 <- c(IVW.0, (datlist[[j]]$df.IVW[,"surv0"] - S0.true))
    POOL.0 <- c(POOL.0, (datlist[[j]]$df.POOL[,"surv0"] - S0.true))
    CCOD.0 <- c(CCOD.0, (datlist[[j]]$df.CCOD[,"surv0"] - S0.true))
    FED.0 <- c(FED.0, (datlist[[j]]$df.FED[,"surv0"] - S0.true))
    FED.bt.0 <- c(FED.bt.0, (datlist[[j]]$df.FED.bt[,"surv0"] - S0.true))
    IPW.TGT.0 <- c(IPW.TGT.0, (datlist.ipwg[[j]][1,2] - S0.true))
    IPW.POOL.0 <- c(IPW.POOL.0, (datlist.ipwg[[j]][2,2] - S0.true))
    G.TGT.0 <- c(G.TGT.0, (datlist.ipwg[[j]][3,2] - S0.true))
    G.POOL.0 <- c(G.POOL.0, (datlist.ipwg[[j]][4,2] - S0.true))
  }
  rmse.TGT.1 <- mean(TGT.1) %>% abs()
  rmse.IVW.1 <- mean(IVW.1) %>% abs()
  rmse.POOL.1 <- mean(POOL.1) %>% abs()
  rmse.CCOD.1 <- mean(CCOD.1) %>% abs()
  rmse.FED.1 <- mean(FED.1) %>% abs()
  rmse.FED.bt.1 <- mean(FED.bt.1) %>% abs()
  rmse.IPW.TGT.1 <- mean(IPW.TGT.1) %>% abs()
  rmse.IPW.POOL.1 <- mean(IPW.POOL.1) %>% abs()
  rmse.G.TGT.1 <- mean(G.TGT.1) %>% abs()
  rmse.G.POOL.1 <- mean(G.POOL.1) %>% abs()
  
  rmse.TGT.0 <- mean(TGT.0) %>% abs()
  rmse.IVW.0 <- mean(IVW.0) %>% abs()
  rmse.POOL.0 <- mean(POOL.0) %>% abs()
  rmse.CCOD.0 <- mean(CCOD.0) %>% abs()
  rmse.FED.0 <- mean(FED.0) %>% abs()
  rmse.FED.bt.0 <- mean(FED.bt.0) %>% abs()
  rmse.IPW.TGT.0 <- mean(IPW.TGT.0) %>% abs()
  rmse.IPW.POOL.0 <- mean(IPW.POOL.0) %>% abs()
  rmse.G.TGT.0 <- mean(G.TGT.0) %>% abs()
  rmse.G.POOL.0 <- mean(G.POOL.0) %>% abs()
  
  return(df.RMSE = data.frame(method=c("TGT", "IVW", "POOL", "CCOD", "FED", "FED (BOOT)", 
                                       "IPW TGT", "IPW POOL", "G-form TGT", "G-form POOL"), 
                              Surv1=c(rmse.TGT.1, rmse.IVW.1, rmse.POOL.1, rmse.CCOD.1, rmse.FED.1,
                                      rmse.FED.bt.1, rmse.IPW.TGT.1, rmse.IPW.POOL.1, rmse.G.TGT.1, 
                                      rmse.G.POOL.1),
                              Surv0=c(rmse.TGT.0, rmse.IVW.0, rmse.POOL.0, rmse.CCOD.0, rmse.FED.0,
                                      rmse.FED.bt.0, rmse.IPW.TGT.0, rmse.IPW.POOL.0, rmse.G.TGT.0, 
                                      rmse.G.POOL.0)
  ))
}

library(xtable)
r1 <- Bias.results(datlist=result.homo, datlist.ipwg=result.homo.ipwg) 
r2 <- Bias.results(datlist=result.diffX, datlist.ipwg=result.diffX.ipwg)
r3 <- Bias.results(datlist=result.diffT, datlist.ipwg=result.diffT.ipwg)
r4 <- Bias.results(datlist=result.diffC, datlist.ipwg=result.diffC.ipwg)
r5 <- Bias.results(datlist=result.diffAll, datlist.ipwg=result.diffAll.ipwg)
print(xtable(cbind(r1, r2[,-1], r3[,-1], r4[,-1], r5[,-1])), 
      include.rownames=F)

RMSE.results <- function(datlist, datlist.ipwg, M=500) {
  TGT.1 <- IVW.1 <- POOL.1 <- CCOD.1 <- FED.1 <- FED.bt.1 <- IPW.TGT.1 <- IPW.POOL.1 <- G.TGT.1 <- G.POOL.1 <- c()
  TGT.0 <- IVW.0 <- POOL.0 <- CCOD.0 <- FED.0 <- FED.bt.0 <- IPW.TGT.0 <- IPW.POOL.0 <- G.TGT.0 <- G.POOL.0 <- c()
  for(j in 1:M) {
    TGT.1 <- c(TGT.1, (datlist[[j]]$df.TGT[,"surv1"] - S1.true)^2)
    IVW.1 <- c(IVW.1, (datlist[[j]]$df.IVW[,"surv1"] - S1.true)^2)
    POOL.1 <- c(POOL.1, (datlist[[j]]$df.POOL[,"surv1"] - S1.true)^2)
    CCOD.1 <- c(CCOD.1, (datlist[[j]]$df.CCOD[,"surv1"] - S1.true)^2)
    FED.1 <- c(FED.1, (datlist[[j]]$df.FED[,"surv1"] - S1.true)^2)
    FED.bt.1 <- c(FED.bt.1, (datlist[[j]]$df.FED.bt[,"surv1"] - S1.true)^2)
    IPW.TGT.1 <- c(IPW.TGT.1, (datlist.ipwg[[j]][1,3] - S1.true)^2)
    IPW.POOL.1 <- c(IPW.POOL.1, (datlist.ipwg[[j]][2,3] - S1.true)^2)
    G.TGT.1 <- c(G.TGT.1, (datlist.ipwg[[j]][3,3] - S1.true)^2)
    G.POOL.1 <- c(G.POOL.1, (datlist.ipwg[[j]][4,3] - S1.true)^2)
    
    TGT.0 <- c(TGT.0, (datlist[[j]]$df.TGT[,"surv0"] - S0.true)^2)
    IVW.0 <- c(IVW.0, (datlist[[j]]$df.IVW[,"surv0"] - S0.true)^2)
    POOL.0 <- c(POOL.0, (datlist[[j]]$df.POOL[,"surv0"] - S0.true)^2)
    CCOD.0 <- c(CCOD.0, (datlist[[j]]$df.CCOD[,"surv0"] - S0.true)^2)
    FED.0 <- c(FED.0, (datlist[[j]]$df.FED[,"surv0"] - S0.true)^2)
    FED.bt.0 <- c(FED.bt.0, (datlist[[j]]$df.FED.bt[,"surv0"] - S0.true)^2)
    IPW.TGT.0 <- c(IPW.TGT.0, (datlist.ipwg[[j]][1,2] - S0.true)^2)
    IPW.POOL.0 <- c(IPW.POOL.0, (datlist.ipwg[[j]][2,2] - S0.true)^2)
    G.TGT.0 <- c(G.TGT.0, (datlist.ipwg[[j]][3,2] - S0.true)^2)
    G.POOL.0 <- c(G.POOL.0, (datlist.ipwg[[j]][4,2] - S0.true)^2)
  }
  rmse.TGT.1 <- mean(TGT.1) %>% sqrt()
  rmse.IVW.1 <- mean(IVW.1) %>% sqrt()
  rmse.POOL.1 <- mean(POOL.1) %>% sqrt()
  rmse.CCOD.1 <- mean(CCOD.1) %>% sqrt()
  rmse.FED.1 <- mean(FED.1) %>% sqrt()
  rmse.FED.bt.1 <- mean(FED.bt.1) %>% sqrt()
  rmse.IPW.TGT.1 <- mean(IPW.TGT.1) %>% sqrt()
  rmse.IPW.POOL.1 <- mean(IPW.POOL.1) %>% sqrt()
  rmse.G.TGT.1 <- mean(G.TGT.1) %>% sqrt()
  rmse.G.POOL.1 <- mean(G.POOL.1) %>% sqrt()
  
  rmse.TGT.0 <- mean(TGT.0) %>% sqrt()
  rmse.IVW.0 <- mean(IVW.0) %>% sqrt()
  rmse.POOL.0 <- mean(POOL.0) %>% sqrt()
  rmse.CCOD.0 <- mean(CCOD.0) %>% sqrt()
  rmse.FED.0 <- mean(FED.0) %>% sqrt()
  rmse.FED.bt.0 <- mean(FED.bt.0) %>% sqrt()
  rmse.IPW.TGT.0 <- mean(IPW.TGT.0) %>% sqrt()
  rmse.IPW.POOL.0 <- mean(IPW.POOL.0) %>% sqrt()
  rmse.G.TGT.0 <- mean(G.TGT.0) %>% sqrt()
  rmse.G.POOL.0 <- mean(G.POOL.0) %>% sqrt()
  
  return(df.RMSE = data.frame(method=c("TGT", "IVW", "POOL", "CCOD", "FED", "FED (BOOT)", 
                                       "IPW TGT", "IPW POOL", "G-form TGT", "G-form POOL"), 
                              Surv1=c(rmse.TGT.1, rmse.IVW.1, rmse.POOL.1, rmse.CCOD.1, rmse.FED.1,
                                      rmse.FED.bt.1, rmse.IPW.TGT.1, rmse.IPW.POOL.1, rmse.G.TGT.1, 
                                      rmse.G.POOL.1),
                              Surv0=c(rmse.TGT.0, rmse.IVW.0, rmse.POOL.0, rmse.CCOD.0, rmse.FED.0,
                                      rmse.FED.bt.0, rmse.IPW.TGT.0, rmse.IPW.POOL.0, rmse.G.TGT.0, 
                                      rmse.G.POOL.0)
                              ))
}

library(xtable)
r1 <- RMSE.results(datlist=result.homo, datlist.ipwg=result.homo.ipwg) 
r2 <- RMSE.results(datlist=result.diffX, datlist.ipwg=result.diffX.ipwg)
r3 <- RMSE.results(datlist=result.diffT, datlist.ipwg=result.diffT.ipwg)
r4 <- RMSE.results(datlist=result.diffC, datlist.ipwg=result.diffC.ipwg)
r5 <- RMSE.results(datlist=result.diffAll, datlist.ipwg=result.diffAll.ipwg)
print(xtable(cbind(r1, r2[,-1], r3[,-1], r4[,-1], r5[,-1])), 
      include.rownames=F)






