# data pre-processing
dat <- read.csv("amp_survival.csv")
head(dat)
mean(dat$hiv1event)
dat$hiv1event
library(dplyr)
library(ggplot2)
library(gridExtra)
library(SuperLearner)
library(survSuperLearner)
library(CFsurvival)
library(glmnet)
library(caret)
library(tidyr)

Delta <- dat$hiv1event
mean(Delta)
Y <- dat$hiv1survday
mean(Y<601)
sum(Y==601)
mean(Y==601)
A <- as.numeric(dat$rx_pool=="T1+T2")
sum(Y==601 & Delta==1)
X <- data.frame(bweight=dat$bweight, score=dat$standardized_risk_score, age=dat$bbmi)

unique(dat$country)
site <- dat$country
site[site%in%c("Tanzania, Mozambique, Kenya", "Zimbabwe", "Botswana", "Malawi")] <- "African country other than South Africa"
site[site%in%c("Peru", "Brazil")] <- "Brazil or Peru"
site[site%in%c("United States", "Switzerland")] <- "United States or Switzerland"
site <- factor(site, levels=c("South Africa", "African country other than South Africa", "Brazil or Peru", "United States or Switzerland"))
unique(site)

dat.hiv <- data.frame(cbind(A, Y, Delta, X, site))
head(dat.hiv)
range(dat.hiv$Y)

### Summary statistics by region
table(dat.hiv$site)

mean(dat.hiv$age)
sd(dat.hiv$age)

dat.hiv %>% group_by(A) %>% summarize(mean(age), sd(age))

mean(dat.hiv$score)
sd(dat.hiv$score)

mean(dat.hiv$bweight)
sd(dat.hiv$bweight)

sum(dat.hiv$Delta)
mean(dat.hiv$Delta)

dat.sum <- dat.hiv %>% group_by(A) %>% 
  summarize(n=n(),
            age.mean=mean(age), age.sd=sd(age),
            score.mean=mean(score), score.sd=sd(score),
            bweight.mean=mean(bweight), bweight.sd=sd(bweight),
            event.count=sum(Delta), event.rate=mean(Delta))
View(dat.sum)

dat.sum.site <- dat.hiv %>% group_by(A, site) %>% 
  summarize(n=n(),
            age.mean=mean(age), age.sd=sd(age),
            score.mean=mean(score), score.sd=sd(score),
            bweight.mean=mean(bweight), bweight.sd=sd(bweight),
            event.count=sum(Delta), event.rate=mean(Delta))
View(dat.sum.site)

### Sensitivity analysis
names(dat.hiv)
library(survival)
library(EValue)
model <- coxph(Surv(Y, Delta)~A+age+score+bweight, data=dat.hiv)
summary_fit <- summary(model)$conf.int
HR <- as.numeric(summary_fit["A", "exp(coef)"])
lo_CI <- as.numeric(summary_fit["A", "lower .95"])
hi_CI <- as.numeric(summary_fit["A", "upper .95"])
evalues.RR(est = HR, lo = lo_CI, hi = hi_CI, true=1)

### Analysis 
# source("FuseSurv.R")
# result.SA <- FuseSurv(data=dat.hiv, n.folds=2, tgt.name="South Africa", s=1)
# save(file="result_SA.Rdata", result.SA)
# result.OA <- FuseSurv(data=dat.hiv, n.folds=2, tgt.name="African country other than South Africa", s=2444)
# save(file="result_OA.Rdata", result.OA)
# result.BP <- FuseSurv(data=dat.hiv, n.folds=2, tgt.name="Brazil or Peru", s=2444)
# save(file="result_BP.Rdata", result.BP)
# result.US <- FuseSurv(data=dat.hiv, n.folds=2, tgt.name="United States or Switzerland", s=2333)
# save(file="result_US.Rdata", result.US)

### Results plotting and summary
load("result_SA.Rdata")
load("result_OA.Rdata")
load("result_BP.Rdata")
load("result_US.Rdata")

plot_survival_CI <- function(
    df, 
    time_col           ="time",
    surv_treated_col   ="surv1",
    surv_treated_sd_col="surv1.sd",
    surv_control_col   ="surv0",
    surv_control_sd_col="surv0.sd",
    ci_multiplier      =1.96,
    color_treated      ="blue",
    color_control      ="red",
    alpha_fill         =0.2,
    line_size          =0.8,
    xlab               ="Time (day)",
    ylab               ="Survival Probability",
    fig.title          ="CCOD"
){
  ggplot(df, aes_string(x=time_col)) +
    geom_line(aes_string(y=surv_treated_col, color="'Treated'"), size=line_size) +
    geom_ribbon(aes_string(
      ymin=paste0(surv_treated_col, " - ", ci_multiplier, " * ", surv_treated_sd_col),
      ymax=paste0(surv_treated_col, " + ", ci_multiplier, " * ", surv_treated_sd_col),
      fill="'Treated'"), alpha=alpha_fill) +
    geom_line(aes_string(y=surv_control_col, color="'Control'"), size=line_size) +
    geom_ribbon(aes_string(
      ymin=paste0(surv_control_col, " - ", ci_multiplier, " * ", surv_control_sd_col),
      ymax=paste0(surv_control_col, " + ", ci_multiplier, " * ", surv_control_sd_col),
      fill="'Control'"), alpha=alpha_fill) +
    scale_color_manual(name ="Group", values=c("Treated"=color_treated, "Control"=color_control)) +
    scale_fill_manual(name  ="Group", values=c("Treated"=color_treated, "Control"=color_control)) +
    labs(x=xlab, y=ylab, title=fig.title) + ylim(0.92,1) + 
    theme_minimal()
}

plot_fedweights <- function(results.list, site.names=c("SA", "OA", "BP", "US")) {
  df1 <- results.list$weights.bt[,1:4] %>% as.data.frame()
  df0 <- results.list$weights.bt[,5:8] %>% as.data.frame()
  df1$time <- seq_len(nrow(df1)) * 7
  df0$time <- seq_len(nrow(df0)) * 7
  colnames(df1) <- colnames(df0) <- c(site.names, "time")
  
  df_long1 <- df1 %>%
    pivot_longer(cols=-time, names_to="Region", values_to="Weight")
  df_long0 <- df0 %>%
    pivot_longer(cols=-time, names_to="Region", values_to="Weight")
  df_long1$Region <- factor(df_long1$Region, levels=c("SA", "OA", "BP", "US"))
  df_long0$Region <- factor(df_long0$Region, levels=c("SA", "OA", "BP", "US"))
  
  ggplot(df_long1, aes(x=time, y=Weight, color=Region)) +
    geom_line(alpha=0.5) + ylim(0,1) +            
    geom_smooth(se=FALSE, size=1) +
    labs(
      title="Treated group",
      x="Time (day)",
      y="Federated weight"
    ) +
    theme_minimal(base_size=14) +
    theme(legend.position="top") -> p.wt.1   
  
  ggplot(df_long0, aes(x=time, y=Weight, color=Region)) +
    geom_line(alpha=0.5) + ylim(0,1) +                      
    geom_smooth(se=FALSE, size=1) +
    labs(
      title="Control group",
      x="Time (day)",
      y="Federated weight"
    ) +
    theme_minimal(base_size=14) +
    theme(legend.position="top") -> p.wt.0
  
  return(list(p.wt.1=p.wt.1, p.wt.0=p.wt.0))
}

# target: South Africa
plot_survival_CI(df=result.SA$df.TGT, color_treated="springgreen4", color_control="springgreen", fig.title="TGT") -> p.tgt.SA
plot_survival_CI(df=result.SA$df.IVW, color_treated="darkorange3", color_control="orange", fig.title="IVW") -> p.ivw.SA
plot_survival_CI(df=result.SA$df.POOL, color_treated="firebrick4", color_control="tomato1", fig.title="POOL") -> p.pool.SA
plot_survival_CI(df=result.SA$df.CCOD, color_treated="darkslateblue", color_control="darkorchid1", fig.title="CCOD") -> p.ccod.SA
plot_survival_CI(df=result.SA$df.FED, color_treated="cyan4", color_control="cyan", fig.title="FED") -> p.fed.SA
plot_survival_CI(df=result.SA$df.FED.bt, color_treated="blue", color_control="cornflowerblue", fig.title="FED (BOOT)") -> p.fedbt.SA

pdf(file="AMP_SA.pdf", width=13, height=5.5)
grid.arrange(p.fed.SA, p.fedbt.SA, p.ccod.SA, p.tgt.SA, p.ivw.SA, p.pool.SA, ncol=3)
dev.off()

p.wt.SA0 <- plot_fedweights(result.SA)$p.wt.0
p.wt.SA1 <- plot_fedweights(result.SA)$p.wt.1

pdf(file="wts_SA1.pdf", width=5.2, height=3.5)
p.wt.SA1
dev.off()

pdf(file="wts_SA0.pdf", width=5.2, height=3.5)
p.wt.SA0
dev.off()

time.RE <- c(148, 330, 512)
TGT.3time <- result.SA$df.TGT[result.SA$df.TGT$time%in%time.RE,]
FED.3time <- result.SA$df.FED[result.SA$df.FED$time%in%time.RE,]
FEDbt.3time <- result.SA$df.FED.bt[result.SA$df.FED.bt$time%in%time.RE,]
CCOD.3time <- result.SA$df.CCOD[result.SA$df.CCOD$time%in%time.RE,]
IVW.3time <- result.SA$df.IVW[result.SA$df.IVW$time%in%time.RE,]
POOL.3time <- result.SA$df.POOL[result.SA$df.POOL$time%in%time.RE,]

FED.3time$surv1.sd / TGT.3time$surv1.sd
FED.3time$surv0.sd / TGT.3time$surv0.sd

FEDbt.3time$surv1.sd / TGT.3time$surv1.sd
FEDbt.3time$surv0.sd / TGT.3time$surv0.sd

CCOD.3time$surv1.sd / TGT.3time$surv1.sd
CCOD.3time$surv0.sd / TGT.3time$surv0.sd

IVW.3time$surv1.sd / TGT.3time$surv1.sd
IVW.3time$surv0.sd / TGT.3time$surv0.sd

POOL.3time$surv1.sd / TGT.3time$surv1.sd
POOL.3time$surv0.sd / TGT.3time$surv0.sd

# target: Other African countries
plot_survival_CI(df=result.OA$df.TGT, color_treated="springgreen4", color_control="springgreen", fig.title="TGT") -> p.tgt.OA
plot_survival_CI(df=result.OA$df.CCOD, color_treated="darkslateblue", color_control="darkorchid1", fig.title="CCOD") -> p.ccod.OA
plot_survival_CI(df=result.OA$df.FED, color_treated="cyan4", color_control="cyan", fig.title="FED") -> p.fed.OA
plot_survival_CI(df=result.OA$df.FED.bt, color_treated="blue", color_control="cornflowerblue", fig.title="FED (BOOT)") -> p.fedbt.OA
plot_survival_CI(df=result.OA$df.IVW, color_treated="darkorange3", color_control="orange", fig.title="IVW") -> p.ivw.OA
plot_survival_CI(df=result.OA$df.POOL, color_treated="firebrick4", color_control="tomato1", fig.title="POOL") -> p.pool.OA

pdf(file="AMP_OA.pdf", width=12, height=5)
grid.arrange(p.fed.OA, p.fedbt.OA, p.ccod.OA, p.tgt.OA, p.ivw.OA, p.pool.OA, ncol=3)
dev.off()

p.wt.OA0 <- plot_fedweights(result.OA, site.names=c("OA", "SA", "BP", "US"))$p.wt.0
p.wt.OA1 <- plot_fedweights(result.OA, site.names=c("OA", "SA", "BP", "US"))$p.wt.1

pdf(file="wts_OA1.pdf", width=6, height=3.5)
p.wt.OA1
dev.off()

pdf(file="wts_OA0.pdf", width=6, height=3.5)
p.wt.OA0
dev.off()

# target: Brazil or Peru
plot_survival_CI(df=result.BP$df.TGT, color_treated="springgreen4", color_control="springgreen", fig.title="TGT") -> p.tgt.BP
plot_survival_CI(df=result.BP$df.CCOD, color_treated="darkslateblue", color_control="darkorchid1", fig.title="CCOD") -> p.ccod.BP
plot_survival_CI(df=result.BP$df.FED, color_treated="cyan4", color_control="cyan", fig.title="FED") -> p.fed.BP
plot_survival_CI(df=result.BP$df.FED.bt, color_treated="blue", color_control="cornflowerblue", fig.title="FED (BOOT)") -> p.fedbt.BP
plot_survival_CI(df=result.BP$df.IVW, color_treated="darkorange3", color_control="orange", fig.title="IVW") -> p.ivw.BP
plot_survival_CI(df=result.BP$df.POOL, color_treated="firebrick4", color_control="tomato1", fig.title="POOL") -> p.pool.BP

pdf(file="AMP_BP.pdf", width=12, height=5)
grid.arrange(p.fed.BP, p.fedbt.BP, p.ccod.BP, p.tgt.BP, p.ivw.BP, p.pool.BP, ncol=3)
dev.off()

p.wt.BP0 <- plot_fedweights(result.BP, site.names=c("BP", "OA", "SA", "US"))$p.wt.0
p.wt.BP1 <- plot_fedweights(result.BP, site.names=c("BP", "OA", "SA", "US"))$p.wt.1

pdf(file="wts_BP1.pdf", width=6, height=3.5)
p.wt.BP1
dev.off()

pdf(file="wts_BP0.pdf", width=6, height=3.5)
p.wt.BP0
dev.off()

# target: United States or Switzerland
plot_survival_CI(df=result.US$df.TGT, color_treated="springgreen4", color_control="springgreen", fig.title="TGT") -> p.tgt.US
plot_survival_CI(df=result.US$df.CCOD, color_treated="darkslateblue", color_control="darkorchid1", fig.title="CCOD") -> p.ccod.US
plot_survival_CI(df=result.US$df.FED, color_treated="cyan4", color_control="cyan", fig.title="FED") -> p.fed.US
plot_survival_CI(df=result.US$df.FED.bt, color_treated="blue", color_control="cornflowerblue", fig.title="FED (BOOT)") -> p.fedbt.US
plot_survival_CI(df=result.US$df.IVW, color_treated="darkorange3", color_control="orange", fig.title="IVW") -> p.ivw.US
plot_survival_CI(df=result.US$df.POOL, color_treated="firebrick4", color_control="tomato1", fig.title="POOL") -> p.pool.US

pdf(file="AMP_US.pdf", width=12, height=5)
grid.arrange(p.fed.US, p.fedbt.US, p.ccod.US, p.tgt.US, p.ivw.US, p.pool.US, ncol=3)
dev.off()

p.wt.US0 <- plot_fedweights(result.US, site.names=c("US", "OA", "SA", "BP"))$p.wt.0
p.wt.US1 <- plot_fedweights(result.BP, site.names=c("US", "OA", "SA", "BP"))$p.wt.1

pdf(file="wts_US1.pdf", width=6, height=3.5)
p.wt.US1
dev.off()

pdf(file="wts_US0.pdf", width=6, height=3.5)
p.wt.US0
dev.off()

# four site-specific survival curves
pdf(file="AMP_siteSurvs.pdf", width=10, height=5)
grid.arrange(p.tgt.SA, p.tgt.OA, p.tgt.BP, p.tgt.US, ncol=2)
dev.off()


# data=dat.hiv
# covar.name=c("age","score","bweight")
# site.var="site"
# tgt.name="Brazil or Peru"
# trt.name="A"
# time.var="Y"
# event="Delta"
# fit.times=0:601
# eval.times=seq(7, 601, by=7)
# prop.SL.library=c("SL.mean", "SL.glm")
# event.SL.library=c("survSL.km", "survSL.coxph", "survSL.rfsrc")
# cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.rfsrc")
# n.folds=2
# n.fed.boot=200
# s=33333



