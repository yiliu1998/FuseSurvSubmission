CIwidth.results <- function(datlist) {
  M <- 500
  ciw.TGT.1 <- ciw.IVW.1 <- ciw.POOL.1 <- ciw.CCOD.1 <- ciw.FED.1 <- ciw.FED.bt.1 <- c()
  ciw.TGT.0 <- ciw.IVW.0 <- ciw.POOL.0 <- ciw.CCOD.0 <- ciw.FED.0 <- ciw.FED.bt.0 <- c()
  quant <- 1.96
  for(i in 1:length(t)) {
    TGT.1 <- IVW.1 <- POOL.1 <- CCOD.1 <- FED.1 <- FED.bt.1 <- matrix(NA, ncol=2, nrow=M)
    TGT.0 <- IVW.0 <- POOL.0 <- CCOD.0 <- FED.0 <- FED.bt.0 <- matrix(NA, ncol=2, nrow=M)
    for(j in 1:M) {
      TGT.1[j,] <- datlist[[j]]$df.TGT[,c("surv1", "surv1.sd")][i,] %>% as.numeric()
      IVW.1[j,] <- datlist[[j]]$df.IVW[,c("surv1", "surv1.sd")][i,] %>% as.numeric()
      POOL.1[j,] <- datlist[[j]]$df.POOL[,c("surv1", "surv1.sd")][i,] %>% as.numeric()
      CCOD.1[j,] <- datlist[[j]]$df.CCOD[,c("surv1", "surv1.sd")][i,] %>% as.numeric()
      FED.1[j,] <- datlist[[j]]$df.FED[,c("surv1", "surv1.sd")][i,] %>% as.numeric()
      FED.bt.1[j,] <- datlist[[j]]$df.FED.bt[,c("surv1", "surv1.sd")][i,] %>% as.numeric()
      
      TGT.0[j,] <- datlist[[j]]$df.TGT[,c("surv0", "surv0.sd")][i,] %>% as.numeric()
      IVW.0[j,] <- datlist[[j]]$df.IVW[,c("surv0", "surv0.sd")][i,] %>% as.numeric()
      POOL.0[j,] <- datlist[[j]]$df.POOL[,c("surv0", "surv0.sd")][i,] %>% as.numeric()
      CCOD.0[j,] <- datlist[[j]]$df.CCOD[,c("surv0", "surv0.sd")][i,] %>% as.numeric()
      FED.0[j,] <- datlist[[j]]$df.FED[,c("surv0", "surv0.sd")][i,] %>% as.numeric()
      FED.bt.0[j,] <- datlist[[j]]$df.FED.bt[,c("surv0", "surv0.sd")][i,] %>% as.numeric()
    }
    
    ciw.TGT.1[i] <- mean(2*quant*TGT.1[,2], na.rm=T)
    ciw.IVW.1[i] <- mean(2*quant*IVW.1[,2], na.rm=T)
    ciw.POOL.1[i] <- mean(2*quant*POOL.1[,2], na.rm=T)
    ciw.CCOD.1[i] <- mean(2*quant*CCOD.1[,2], na.rm=T)
    ciw.FED.1[i] <- mean(2*quant*FED.1[,2], na.rm=T)
    ciw.FED.bt.1[i] <- mean(2*quant*FED.bt.1[,2], na.rm=T)
  
    ciw.TGT.0[i] <- mean(2*quant*TGT.0[,2], na.rm=T)
    ciw.IVW.0[i] <- mean(2*quant*IVW.0[,2], na.rm=T)
    ciw.POOL.0[i] <- mean(2*quant*POOL.0[,2], na.rm=T)
    ciw.CCOD.0[i] <- mean(2*quant*CCOD.0[,2], na.rm=T)
    ciw.FED.0[i] <- mean(2*quant*FED.0[,2], na.rm=T)
    ciw.FED.bt.0[i] <- mean(2*quant*FED.bt.0[,2], na.rm=T)
    }
  return(list(ciw.TGT.1, ciw.IVW.1, ciw.POOL.1, ciw.CCOD.1, ciw.FED.1, ciw.FED.bt.1,
              ciw.TGT.0, ciw.IVW.0, ciw.POOL.0, ciw.CCOD.0, ciw.FED.0, ciw.FED.bt.0))
}

library(ggplot2)
library(dplyr)
library(reshape2)
load("truth.Rdata")

load("Res_homo_s.Rdata")
load("Res_diffX_s.Rdata")
load("Res_diffT_s.Rdata")
load("Res_diffC_s.Rdata")
load("Res_diffAll_s.Rdata")

homo <- CIwidth.results(result.homo)
diffX <- CIwidth.results(result.diffX)
diffT <- CIwidth.results(result.diffT)
diffC <- CIwidth.results(result.diffC)
diffAll <- CIwidth.results(result.diffAll)

all.ciw <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.ciw)
time <- rep(t, n/length(t))

Case <- factor(c(rep("Homogeneous", n/5), rep("Covariate Shift", n/5), rep("Outcome Shift", n/5), 
                 rep("Censoring Shift", n/5), rep("All Shift", n/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift"))
Method <- factor(rep(c(rep("TGT", length(t)), rep("IVW", length(t)), rep("POOL", length(t)), 
                       rep("CCOD", length(t)), rep("FED", length(t)), rep("FED (BOOT)", length(t))), 10), 
                 levels=c("FED", "FED (BOOT)", "TGT", "IVW", "POOL", "CCOD"))
Treatment <- rep(c(rep("Treated", n/10), rep("Control", n/10)), 5)
df <- data.frame(CIW=all.ciw, Case=Case, Method=Method, time=time, Treatment=Treatment)
df_heatmap <- df %>%
  group_by(Case, Method, Treatment, time) %>% summarize(CIW=median(CIW), .groups='drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=CIW)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(CIW, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white", "slateblue2"),
    name="CI width"
  ) +
  facet_grid(Treatment~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.ciw

pdf(file="ciw_s.pdf", width=10, height=3.25)
p.ciw
dev.off()


load("Res_homo_l.Rdata")
load("Res_diffX_l.Rdata")
load("Res_diffT_l.Rdata")
load("Res_diffC_l.Rdata")
load("Res_diffAll_l.Rdata")

homo <- CIwidth.results(result.homo)
diffX <- CIwidth.results(result.diffX)
diffT <- CIwidth.results(result.diffT)
diffC <- CIwidth.results(result.diffC)
diffAll <- CIwidth.results(result.diffAll)

all.ciw <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
df <- data.frame(CIW=all.ciw, Case=Case, Method=Method, time=time, Treatment=Treatment)
df_heatmap <- df %>%
  group_by(Case, Method, Treatment, time) %>% summarize(CIW=median(CIW), .groups='drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=CIW)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(CIW, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white", "slateblue2"),
    name="CI width"
  ) +
  facet_grid(Treatment~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.ciw

pdf(file="ciw_l.pdf", width=10, height=3.25)
p.ciw
dev.off()

### main results in the paper
df.main <- df %>% filter(time==90)
df_heatmap <- df.main %>%
  group_by(Case, Method, Treatment) %>% summarize(CIW=median(CIW), .groups='drop')

ggplot(df_heatmap, aes(x=factor(Case), y=Method, fill=CIW)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(CIW, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white", "slateblue2"),
    name="CI width"
  ) +
  facet_grid(.~Treatment) + labs(x="", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank(), axis.text.x=element_text(angle=35, hjust=1))  -> p.ciw.main

pdf(file="ciw_main.pdf", width=6, height=2.5)
p.ciw.main
dev.off()


load("Res_homo_l2.Rdata")
load("Res_diffX_l2.Rdata")
load("Res_diffT_l2.Rdata")
load("Res_diffC_l2.Rdata")
load("Res_diffAll_l2.Rdata")

homo <- CIwidth.results(result.homo)
diffX <- CIwidth.results(result.diffX)
diffT <- CIwidth.results(result.diffT)
diffC <- CIwidth.results(result.diffC)
diffAll <- CIwidth.results(result.diffAll)

all.ciw <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
df <- data.frame(CIW=all.ciw, Case=Case, Method=Method, time=time, Treatment=Treatment)
df_heatmap <- df %>%
  group_by(Case, Method, Treatment, time) %>% summarize(CIW=median(CIW), .groups='drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=CIW)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(CIW, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white", "slateblue2"),
    name="CI width"
  ) +
  facet_grid(Treatment~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.ciw

pdf(file="ciw_l2.pdf", width=10, height=3.25)
p.ciw
dev.off()