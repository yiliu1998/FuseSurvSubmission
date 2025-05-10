weights.results <- function(datlist) {
  M <- 500
  N.site <- ncol(datlist[[1]]$weights)/2
  weights.1 <- weights.0 <- matrix(NA, nrow=length(t), ncol=N.site)
  for(i in 1:length(t)) {
    weights.1.i <- weights.0.i <- matrix(NA, nrow=M, ncol=N.site)
    for(j in 1:M) {
      weights.0.i[j,] <- datlist[[j]]$weights.bt[i,1:N.site]
      weights.1.i[j,] <- datlist[[j]]$weights.bt[i,(N.site+1):(2*N.site)]
    }
    weights.1[i,] <- apply(weights.1.i, 2, mean)
    weights.0[i,] <- apply(weights.0.i, 2, mean)
  }
  return(list(weights.1 = weights.1, weights.0 = weights.0))
}

library(tidyverse)
library(RColorBrewer)
load("Res_homo_s.Rdata")
load("Res_diffX_s.Rdata")
load("Res_diffT_s.Rdata")
load("Res_diffC_s.Rdata")
load("Res_diffAll_s.Rdata")
load("truth.Rdata")

n.time <- length(t)
homo.weights <- weights.results(result.homo)
diffX.weights <- weights.results(result.diffX)
diffT.weights <- weights.results(result.diffT)
diffC.weights <- weights.results(result.diffC)
diffAll.weights <- weights.results(result.diffAll)

weights <- c(homo.weights$weights.1, homo.weights$weights.0,
             diffX.weights$weights.1, diffX.weights$weights.0,
             diffT.weights$weights.1, diffT.weights$weights.0,
             diffC.weights$weights.1, diffC.weights$weights.0,
             diffAll.weights$weights.1, diffAll.weights$weights.0)
n <- length(weights)
time <- rep(t, n/n.time)
Case <- factor(c(rep("Homogeneous", n/5), rep("Covariate Shift", n/5), rep("Outcome Shift", n/5), 
                 rep("Censoring Shift", n/5), rep("All Shift", n/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift"))
Site <- factor(rep(c(rep("target", n.time), rep("source 1", n.time), rep("source 2", n.time), 
                     rep("source 3", n.time), rep("source 4", n.time)), n/(n.time*5)), 
               levels=c("target", "source 1", "source 2", "source 3", "source 4"))
Treatment <- rep(c(rep("Treated", n/10), rep("Control", n/10)), 5)

df.weights <- data.frame(weights=weights, Case=Case, Site=Site, time=time, Treatment=Treatment)
colors <- brewer.pal(11, name="Paired")[c(2,5,7,8,9)]
ggplot(df.weights, aes(x=factor(time), y=weights, color=Site, shape=Site, group=Site)) + 
  geom_point(size=3, alpha=0.7) + geom_line() + geom_hline(yintercept=0.2, col="hotpink2", lty=2, linewidth=0.6) + 
  scale_color_manual(values = colors) + ylim(0,0.7) +   
  facet_grid(Treatment~Case) + labs(x="time", y="Federated weights") + theme_bw() -> p.weights

pdf(file="weights.pdf", width=8, height=3)
p.weights
dev.off()

### weights vs chisq values
weights.vs.chi <- function(datlist, t=c(30,60,90)) {
  M <- 500
  N.site <- N.site <- ncol(datlist[[1]]$weights)/2
  n.time <- length(t)
  weights <- chi <- site <- time <- trt <- c()
  for(j in 1:M) {
    chi <- c(chi, datlist[[j]]$chi)
    weights <- c(weights, datlist[[j]]$weights[,-c(1,N.site+1)])
    site <- c(site, rep(paste("source", sort(rep(1:(N.site-1), n.time))), 2))
    time <- c(time, rep(t, 2*(N.site-1)))
    trt <- c(trt, rep("Control", n.time*(N.site-1)), rep("Treated", n.time*(N.site-1)))
  }
  df <- data.frame(weights=weights, chisq=chi^2, site=site, time=time, trt=trt)
  return(df=df)
}

homo.weights.df <- weights.vs.chi(result.homo)
diffX.weights.df <- weights.vs.chi(result.diffX)
diffT.weights.df <- weights.vs.chi(result.diffT)
diffC.weights.df <- weights.vs.chi(result.diffC)
diffAll.weights.df <- weights.vs.chi(result.diffAll)
n.all <- nrow(homo.weights.df)

all.weights.df <- rbind(homo.weights.df, diffX.weights.df, diffT.weights.df,
                        diffC.weights.df, diffAll.weights.df)
Case <- factor(c(rep("Homogeneous", n.all), rep("Covariate Shift", n.all), rep("Outcome Shift", n.all),
                 rep("Censoring Shift", n.all), rep("All Shift", n.all)),
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift"))
all.weights.df <- all.weights.df %>% mutate(Case=Case)
all.weights.df$time <- as.factor(all.weights.df$time)

colors <- brewer.pal(11, "Paired")[c(2,4,7)]
png("wts_chi.png", width=1600, height=900, res=144)
ggplot(data=all.weights.df, aes(x=chisq, y=weights, color=time)) +
  geom_point(alpha=0.4) + geom_hline(yintercept=0.2, col="hotpink2", lty=2, linewidth=0.6) +
  labs(x=expression(paste(chi^2, " value"))) + scale_color_manual(values=colors) +
  scale_x_continuous(limits=c(0, 0.3)) + scale_y_continuous(limits=c(0,1)) + 
  facet_grid(Case~site) + theme_minimal()
dev.off()

