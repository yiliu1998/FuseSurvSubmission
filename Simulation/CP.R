CP.results <- function(datlist) {
  M <- 500
  cp.TGT.1 <- cp.IVW.1 <- cp.POOL.1 <- cp.CCOD.1 <- cp.FED.1 <- cp.FED.bt.1 <- c()
  cp.TGT.0 <- cp.IVW.0 <- cp.POOL.0 <- cp.CCOD.0 <- cp.FED.0 <- cp.FED.bt.0 <- c()
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
    
    cp.TGT.1[i] <- mean((S1.true[i]<TGT.1[,1] + quant*TGT.1[,2]) & (S1.true[i]>TGT.1[,1] - quant*TGT.1[,2]), na.rm=T)
    cp.IVW.1[i] <- mean((S1.true[i]<IVW.1[,1] + quant*IVW.1[,2]) & (S1.true[i]>IVW.1[,1] - quant*IVW.1[,2]), na.rm=T)
    cp.POOL.1[i] <- mean((S1.true[i]<POOL.1[,1] + quant*POOL.1[,2]) & (S1.true[i]>POOL.1[,1] - quant*POOL.1[,2]), na.rm=T)
    cp.CCOD.1[i] <- mean((S1.true[i]<CCOD.1[,1] + quant*CCOD.1[,2]) & (S1.true[i]>CCOD.1[,1] - quant*CCOD.1[,2]), na.rm=T)
    cp.FED.1[i] <- mean((S1.true[i]<FED.1[,1] + quant*FED.1[,2]) & (S1.true[i]>FED.1[,1] - quant*FED.1[,2]), na.rm=T)
    cp.FED.bt.1[i] <- mean((S1.true[i]<FED.bt.1[,1] + quant*FED.bt.1[,2]) & (S1.true[i]>FED.bt.1[,1] - quant*FED.bt.1[,2]), na.rm=T)
  
    cp.TGT.0[i] <- mean((S0.true[i]<TGT.0[,1] + quant*TGT.0[,2]) & (S0.true[i]>TGT.0[,1] - quant*TGT.0[,2]), na.rm=T)
    cp.IVW.0[i] <- mean((S0.true[i]<IVW.0[,1] + quant*IVW.0[,2]) & (S0.true[i]>IVW.0[,1] - quant*IVW.0[,2]), na.rm=T)
    cp.POOL.0[i] <- mean((S0.true[i]<POOL.0[,1] + quant*POOL.0[,2]) & (S0.true[i]>POOL.0[,1] - quant*POOL.0[,2]), na.rm=T)
    cp.CCOD.0[i] <- mean((S0.true[i]<CCOD.0[,1] + quant*CCOD.0[,2]) & (S0.true[i]>CCOD.0[,1] - quant*CCOD.0[,2]), na.rm=T)
    cp.FED.0[i] <- mean((S0.true[i]<FED.0[,1] + quant*FED.0[,2]) & (S0.true[i]>FED.0[,1] - quant*FED.0[,2]), na.rm=T)
    cp.FED.bt.0[i] <- mean((S0.true[i]<FED.bt.0[,1] + quant*FED.bt.0[,2]) & (S0.true[i]>FED.bt.0[,1] - quant*FED.bt.0[,2]), na.rm=T)
    
    }
  return(list(cp.TGT.1, cp.IVW.1, cp.POOL.1, cp.CCOD.1, cp.FED.1, cp.FED.bt.1,
              cp.TGT.0, cp.IVW.0, cp.POOL.0, cp.CCOD.0, cp.FED.0, cp.FED.bt.0))
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

homo <- CP.results(result.homo)
diffX <- CP.results(result.diffX)
diffT <- CP.results(result.diffT)
diffC <- CP.results(result.diffC)
diffAll <- CP.results(result.diffAll)

all.cp <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.cp)
time <- rep(t, n/length(t))

Case <- factor(c(rep("Homogeneous", n/5), rep("Covariate Shift", n/5), rep("Outcome Shift", n/5), 
                 rep("Censoring Shift", n/5), rep("All Shift", n/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift"))
Method <- factor(rep(c(rep("TGT", length(t)), rep("IVW", length(t)), rep("POOL", length(t)), 
                       rep("CCOD", length(t)), rep("FED", length(t)), rep("FED (BOOT)", length(t))), 10), 
                 levels=c("FED", "FED (BOOT)", "TGT", "IVW", "POOL", "CCOD"))
Treatment <- rep(c(rep("Treated", n/10), rep("Control", n/10)), 5)
df <- data.frame(CP=all.cp, Case=Case, Method=Method, time=time, Treatment=Treatment)

df_heatmap <- df %>%
  group_by(Case, Method, Treatment, time) %>% summarize(CP=mean(CP * 100), .groups='drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=CP)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(CP, 1)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("mediumvioletred", "white","mediumvioletred"), 
    values=scales::rescale(c(0, 87, 95, 99, 100)), # Soft transition around 95%
    name="CP%",
    limits=c(0, 100)
  ) +
  facet_grid(Treatment~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.cp

pdf(file="cp_s.pdf", width=10, height=3.25)
p.cp
dev.off()


load("Res_homo_l.Rdata")
load("Res_diffX_l.Rdata")
load("Res_diffT_l.Rdata")
load("Res_diffC_l.Rdata")
load("Res_diffAll_l.Rdata")

homo <- CP.results(result.homo)
diffX <- CP.results(result.diffX)
diffT <- CP.results(result.diffT)
diffC <- CP.results(result.diffC)
diffAll <- CP.results(result.diffAll)

all.cp <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.cp)
time <- rep(t, n/length(t))
df <- data.frame(CP=all.cp, Case=Case, Method=Method, time=time, Treatment=Treatment)

df_heatmap <- df %>%
  group_by(Case, Method, Treatment, time) %>% summarize(CP=mean(CP * 100), .groups='drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=CP)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(CP, 1)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("mediumvioletred", "white","mediumvioletred"), 
    values=scales::rescale(c(0, 87, 95, 99, 100)), # Soft transition around 95%
    name="CP%",
    limits=c(0, 100)
  ) +
  facet_grid(Treatment~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.cp

pdf(file="cp_l.pdf", width=10, height=3.25)
p.cp
dev.off()

### main results in the paper
df.main <- df %>% filter(time==90)
df_heatmap <- df.main %>%
  group_by(Case, Method, Treatment) %>% summarize(CP=mean(CP * 100), .groups='drop')

ggplot(df_heatmap, aes(x=Case, y=Method, fill=CP)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(CP, 1)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("mediumvioletred", "white","mediumvioletred"),
    values=scales::rescale(c(0, 87, 95, 99, 100)), # Soft transition around 95%
    name="CP%",
    limits=c(0, 100)
  ) +
  facet_grid(.~Treatment) + labs(x="", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank(), axis.text.x=element_text(angle=35, hjust=1)) -> p.cp.main

pdf(file="cp_main.pdf", width=6, height=2.5)
p.cp.main
dev.off()



load("Res_homo_l2.Rdata")
load("Res_diffX_l2.Rdata")
load("Res_diffT_l2.Rdata")
load("Res_diffC_l2.Rdata")
load("Res_diffAll_l2.Rdata")

homo <- CP.results(result.homo)
diffX <- CP.results(result.diffX)
diffT <- CP.results(result.diffT)
diffC <- CP.results(result.diffC)
diffAll <- CP.results(result.diffAll)

all.cp <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.cp)
time <- rep(t, n/length(t))
df <- data.frame(CP=all.cp, Case=Case, Method=Method, time=time, Treatment=Treatment)

df_heatmap <- df %>%
  group_by(Case, Method, Treatment, time) %>% summarize(CP=mean(CP * 100), .groups='drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=CP)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(CP, 1)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("mediumvioletred", "white","mediumvioletred"), 
    values=scales::rescale(c(0, 87, 95, 99, 100)), # Soft transition around 95%
    name="CP%",
    limits=c(0, 100)
  ) +
  facet_grid(Treatment~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.cp

pdf(file="cp_l2.pdf", width=10, height=3.25)
p.cp
dev.off()

