all.bias.list <- function(datlist) {
  M <- 500
  TGT.1 <- IVW.1 <- POOL.1 <- CCOD.1 <- FedCSA.1 <- FedCSA.bt.1 <- c()
  TGT.0 <- IVW.0 <- POOL.0 <- CCOD.0 <- FedCSA.0 <- FedCSA.bt.0 <- c()
  for(j in 1:M) {
    TGT.1 <- c(TGT.1, datlist[[j]]$df.TGT[,"surv1"] - S1.true)
    IVW.1 <- c(IVW.1, datlist[[j]]$df.IVW[,"surv1"] - S1.true)
    POOL.1 <- c(POOL.1, datlist[[j]]$df.POOL[,"surv1"] - S1.true)
    CCOD.1 <- c(CCOD.1, datlist[[j]]$df.CCOD[,"surv1"] - S1.true)
    FedCSA.1 <- c(FedCSA.1, datlist[[j]]$df.FED[,"surv1"] - S1.true)
    FedCSA.bt.1 <- c(FedCSA.bt.1, datlist[[j]]$df.FED.bt[,"surv1"] - S1.true)
    
    TGT.0 <- c(TGT.0, datlist[[j]]$df.TGT[,"surv0"] - S0.true)
    IVW.0 <- c(IVW.0, datlist[[j]]$df.IVW[,"surv0"] - S0.true)
    POOL.0 <- c(POOL.0, datlist[[j]]$df.POOL[,"surv0"] - S0.true)
    CCOD.0 <- c(CCOD.0, datlist[[j]]$df.CCOD[,"surv0"] - S0.true)
    FedCSA.0 <- c(FedCSA.0, datlist[[j]]$df.FED[,"surv0"] - S0.true)
    FedCSA.bt.0 <- c(FedCSA.bt.0, datlist[[j]]$df.FED.bt[,"surv0"] - S0.true)
  }
  return(list(TGT.1, IVW.1, POOL.1, CCOD.1, FedCSA.1, FedCSA.bt.1,
              TGT.0, IVW.0, POOL.0, CCOD.0, FedCSA.0, FedCSA.bt.0))
}

library(tidyverse)
library(RColorBrewer)
load("Res_homo_s.Rdata")
load("Res_diffX_s.Rdata")
load("Res_diffT_s.Rdata")
load("Res_diffC_s.Rdata")
load("Res_diffAll_s.Rdata")
load("truth.Rdata")

homo <- all.bias.list(result.homo)
diffX <- all.bias.list(result.diffX)
diffT <- all.bias.list(result.diffT)
diffC <- all.bias.list(result.diffC)
diffAll <- all.bias.list(result.diffAll)

all.bias <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.bias)
M <- 500
t <- c(30,60,90)
time <- rep(t, n/length(t))
Case <- factor(c(rep("Homogeneous", n/5), rep("Covariate Shift", n/5), rep("Outcome Shift", n/5), 
                 rep("Censoring Shift", n/5), rep("All Shift", n/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift"))
Method <- factor(rep(c(rep("TGT", length(t)*M), rep("IVW", length(t)*M), rep("POOL", length(t)*M), 
                       rep("CCOD", length(t)*M), rep("FED", length(t)*M), rep("FED (BOOT)", length(t)*M)), 10), 
                 levels=c("FED", "FED (BOOT)", "TGT", "IVW", "POOL", "CCOD"))
Treatment <- rep(c(rep("Treated", n/10), rep("Control", n/10)), 5)
df <- data.frame(Bias=all.bias, Case=Case, Method=Method, time=time, Treatment=Treatment)

colors <- brewer.pal(11, "Paired")[c(1,2,3,7,8,10)]
ggplot(df, aes(x=factor(time), y=Bias, fill=Method)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.6, color="gray50") + scale_fill_manual(values=colors) + 
  geom_hline(yintercept=0, col="hotpink2", lty=2, linewidth=0.6) + 
  facet_grid(Treatment~Case) + labs(x="Time (day)", y="Bias") + theme_bw() -> p.bias

pdf(file="bias_s.pdf", width=12, height=4)
p.bias
dev.off()



load("Res_homo_l.Rdata")
load("Res_diffX_l.Rdata")
load("Res_diffT_l.Rdata")
load("Res_diffC_l.Rdata")
load("Res_diffAll_l.Rdata")

homo <- all.bias.list(result.homo)
diffX <- all.bias.list(result.diffX)
diffT <- all.bias.list(result.diffT)
diffC <- all.bias.list(result.diffC)
diffAll <- all.bias.list(result.diffAll)

all.bias <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.bias)
df <- data.frame(Bias=all.bias, Case=Case, Method=Method, time=time, Treatment=Treatment)

ggplot(df, aes(x=factor(time), y=Bias, fill=Method)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.6, color="gray50") + scale_fill_manual(values=colors) + 
  geom_hline(yintercept=0, col="hotpink2", lty=2, linewidth=0.6) + 
  facet_grid(Treatment~Case) + labs(x="Time (day)", y="Bias") + theme_bw() -> p.bias

pdf(file="bias_l.pdf", width=12, height=4)
p.bias
dev.off()

### main results in the paper
df.main <- df %>% filter(time==90)
ggplot(df.main, aes(x=Case, y=Bias, fill=Method)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.4, color="gray50") + scale_fill_manual(values=colors) + 
  geom_hline(yintercept=0, col="hotpink2", lty=2, linewidth=0.6) + 
  facet_grid(Treatment~.) + labs(x="", y="Bias") + 
  theme_bw() + theme(axis.text.x=element_text(angle=30, hjust=1)) -> p.bias.main

pdf(file="bias_main.pdf", width=7, height=3.2)
p.bias.main
dev.off()



load("Res_homo_l2.Rdata")
load("Res_diffX_l2.Rdata")
load("Res_diffT_l2.Rdata")
load("Res_diffC_l2.Rdata")
load("Res_diffAll_l2.Rdata")

homo <- all.bias.list(result.homo)
diffX <- all.bias.list(result.diffX)
diffT <- all.bias.list(result.diffT)
diffC <- all.bias.list(result.diffC)
diffAll <- all.bias.list(result.diffAll)

all.bias <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.bias)
df <- data.frame(Bias=all.bias, Case=Case, Method=Method, time=time, Treatment=Treatment)

ggplot(df, aes(x=factor(time), y=Bias, fill=Method)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.6, color="gray50") + scale_fill_manual(values=colors) + 
  geom_hline(yintercept=0, col="hotpink2", lty=2, linewidth=0.6) +
  facet_grid(Treatment~Case) + labs(x="Time (day)", y="Bias") + theme_bw() -> p.bias

pdf(file="bias_l2.pdf", width=12, height=4)
p.bias
dev.off()
