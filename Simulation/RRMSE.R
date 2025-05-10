RRMSE.results <- function(datlist) {
  M <- 500
  TGT.1 <- IVW.1 <- POOL.1 <- CCOD.1 <- FED.1 <- FED.bt.1 <- c()
  TGT.0 <- IVW.0 <- POOL.0 <- CCOD.0 <- FED.0 <- FED.bt.0 <- c()
  for(j in 1:M) {
    TGT.1 <- c(TGT.1, (datlist[[j]]$df.TGT[,"surv1"] - S1.true)^2)
    IVW.1 <- c(IVW.1, (datlist[[j]]$df.IVW[,"surv1"] - S1.true)^2)
    POOL.1 <- c(POOL.1, (datlist[[j]]$df.POOL[,"surv1"] - S1.true)^2)
    CCOD.1 <- c(CCOD.1, (datlist[[j]]$df.CCOD[,"surv1"] - S1.true)^2)
    FED.1 <- c(FED.1, (datlist[[j]]$df.FED[,"surv1"] - S1.true)^2)
    FED.bt.1 <- c(FED.bt.1, (datlist[[j]]$df.FED.bt[,"surv1"] - S1.true)^2)
    
    TGT.0 <- c(TGT.0, (datlist[[j]]$df.TGT[,"surv0"] - S0.true)^2)
    IVW.0 <- c(IVW.0, (datlist[[j]]$df.IVW[,"surv0"] - S0.true)^2)
    POOL.0 <- c(POOL.0, (datlist[[j]]$df.POOL[,"surv0"] - S0.true)^2)
    CCOD.0 <- c(CCOD.0, (datlist[[j]]$df.CCOD[,"surv0"] - S0.true)^2)
    FED.0 <- c(FED.0, (datlist[[j]]$df.FED[,"surv0"] - S0.true)^2)
    FED.bt.0 <- c(FED.bt.0, (datlist[[j]]$df.FED.bt[,"surv0"] - S0.true)^2)
  }
  rrmse.TGT.1 <- TGT.1/TGT.1
  rrmse.IVW.1 <- IVW.1/TGT.1
  rrmse.POOL.1 <- POOL.1/TGT.1
  rrmse.CCOD.1 <- CCOD.1/TGT.1 
  rrmse.FED.1 <- FED.1/TGT.1
  rrmse.FED.bt.1 <- FED.bt.1/TGT.1 
  
  rrmse.TGT.0 <- TGT.0/TGT.0
  rrmse.IVW.0 <- IVW.0/TGT.0
  rrmse.POOL.0 <- POOL.0/TGT.0 
  rrmse.CCOD.0 <- CCOD.0/TGT.0 
  rrmse.FED.0 <- FED.0/TGT.0
  rrmse.FED.bt.0 <- FED.bt.0/TGT.0
  
  return(list(rrmse.TGT.1, rrmse.IVW.1, rrmse.POOL.1, rrmse.CCOD.1, rrmse.FED.1, rrmse.FED.bt.1,
              rrmse.TGT.0, rrmse.IVW.0, rrmse.POOL.0, rrmse.CCOD.0, rrmse.FED.0, rrmse.FED.bt.0))
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

homo <- RRMSE.results(result.homo)
diffX <- RRMSE.results(result.diffX)
diffT <- RRMSE.results(result.diffT)
diffC <- RRMSE.results(result.diffC)
diffAll <- RRMSE.results(result.diffAll)

all.rrmse <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.rrmse)
time <- rep(t, n/length(t))

Case <- factor(c(rep("Homogeneous", n/5), rep("Covariate Shift", n/5), rep("Outcome Shift", n/5), 
                 rep("Censoring Shift", n/5), rep("All Shift", n/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift"))
Method <- factor(rep(c(rep("TGT", length(t)*500), rep("IVW", length(t)*500), rep("POOL", length(t)*500), 
                       rep("CCOD", length(t)*500), rep("FED", length(t)*500), rep("FED (BOOT)", length(t)*500)), 10), 
                 levels=c("FED", "FED (BOOT)", "TGT", "IVW", "POOL", "CCOD"))
Treatment <- rep(c(rep("Treated", n/10), rep("Control", n/10)), 5)

df <- data.frame(RRMSE=all.rrmse, Case=Case, Method=Method, time=time, Treatment=Treatment)
df_heatmap <- df %>%
  group_by(Case, Method, Treatment, time) %>% summarize(RRMSE=median(RRMSE), .groups='drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=RRMSE)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RRMSE, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white","white", "orange"), 
    values=scales::rescale(c(0, 1, 100)),
    name="RRMSE"
  ) +
  facet_grid(Treatment~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.rrmse

pdf(file="rrmse_s.pdf", width=10, height=3.25)
p.rrmse
dev.off()


load("Res_homo_l.Rdata")
load("Res_diffX_l.Rdata")
load("Res_diffT_l.Rdata")
load("Res_diffC_l.Rdata")
load("Res_diffAll_l.Rdata")

homo <- RRMSE.results(result.homo)
diffX <- RRMSE.results(result.diffX)
diffT <- RRMSE.results(result.diffT)
diffC <- RRMSE.results(result.diffC)
diffAll <- RRMSE.results(result.diffAll)

all.rrmse <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
df <- data.frame(RRMSE=all.rrmse, Case=Case, Method=Method, time=time, Treatment=Treatment)
df_heatmap <- df %>%
  group_by(Case, Method, Treatment, time) %>% summarize(RRMSE=median(RRMSE), .groups='drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=RRMSE)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RRMSE, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white","white", "orange"), 
    values=scales::rescale(c(0, 1, 100)),
    name="RRMSE"
  ) +
  facet_grid(Treatment~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.rrmse

pdf(file="rrmse_l.pdf", width=10, height=3.25)
p.rrmse
dev.off()

### main results in the paper
df.main <- df %>% filter(time==90)
df_heatmap <- df.main %>%
  group_by(Case, Method, Treatment) %>% summarize(RRMSE=median(RRMSE), .groups='drop')

ggplot(df_heatmap, aes(x=factor(Case), y=Method, fill=RRMSE)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RRMSE, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white","white", "orange"), 
    values=scales::rescale(c(0, 1, 100)),
    name="RRMSE"
  ) +
  facet_grid(.~Treatment) + labs(x="", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank(), axis.text.x=element_text(angle=35, hjust=1))  -> p.rrmse.main

pdf(file="rrmse_main.pdf", width=6, height=2.5)
p.rrmse.main
dev.off()



load("Res_homo_l2.Rdata")
load("Res_diffX_l2.Rdata")
load("Res_diffT_l2.Rdata")
load("Res_diffC_l2.Rdata")
load("Res_diffAll_l2.Rdata")

homo <- RRMSE.results(result.homo)
diffX <- RRMSE.results(result.diffX)
diffT <- RRMSE.results(result.diffT)
diffC <- RRMSE.results(result.diffC)
diffAll <- RRMSE.results(result.diffAll)

all.rrmse <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
df <- data.frame(RRMSE=all.rrmse, Case=Case, Method=Method, time=time, Treatment=Treatment)
df_heatmap <- df %>%
  group_by(Case, Method, Treatment, time) %>% summarize(RRMSE=median(RRMSE), .groups='drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=RRMSE)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RRMSE, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white","white", "orange"), 
    values=scales::rescale(c(0, 1, 100)),
    name="RRMSE"
  ) +
  facet_grid(Treatment~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.rrmse

pdf(file="rrmse_l2.pdf", width=10, height=3.25)
p.rrmse
dev.off()
