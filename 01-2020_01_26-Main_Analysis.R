packages <- c("data.table","broom","gridExtra",
              "here","VIM","haven","skimr","beepr",
              "GGally","ranger","tmle")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN') 
  }
}

for (package in packages) {
  library(package, character.only=T)
}

devtools::install_github("hadley/tidyverse")
devtools::install_github("hadley/ggplot2")
devtools::install_github("ehkennedy/npcausal")

for (package in c("tidyverse","ggplot2","npcausal")) {
  library(package, character.only=T)
}

thm <- theme_classic() +
  theme(
    legend.position = "top",
    #legend.title=element_blank(),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

set.seed(123)
#' Load data, create first and last id variables, and explore:
main <- read_csv("data/numomData_small_12.04.2019.csv") %>% 
  mutate(outcome=pree_acog,
         exposure=as.numeric(v_totdens>quantile(v_totdens,.2)),
         age=as.numeric(age),
         age=if_else(is.na(age),mean(age,na.rm=T),age),
         education=as.numeric(education),
         married=married-1,
         prepregbmi=as.numeric(scale(prepregbmi)),
         pct_emptyc=pct_emptyc)

table(main$education)

# main %>%
#   select(prepregbmi,age,dt_kcal,
#          pct_emptyc,f_totdens,v_totdens,
#          d_totdens,p_totdens,g_whldens,
#          g_nwhldens,fa_dens,
#          seaplant_dens,sodium_dens) %>%  ggpairs()
  
outcome <- main$outcome
exposure <- main$exposure
covs <- main %>% select(age,black,education,married,smokerpre,prepregbmi,insurpub,
                        pct_emptyc,d_totdens,p_totdens,f_totdens,sodium_dens,fa_dens,g_nwhldens,g_whldens)
id <- 1:nrow(main)
time <- rep(1,nrow(main))

#' Exposure Model (GLM)
gFit <- glm(exposure~.,data=data.frame(covs),family=binomial("logit"))
glm_dat <- tidy(gFit)[-1,] #%>% select(term,estimate) %>% rename(names=term)
glm_dat

#' Variable Importance for Exposure(Ranger)
# rFit <- ranger(factor(exposure)~.,
#                data=data.frame(covs),
#                probability = T,
#                importance = "impurity")
# 
# plotDat <- data.frame(rFit$variable.importance)
# plotDat$names <- row.names(plotDat)
# names(plotDat)[1] <- "importance"
# row.names(plotDat) <- NULL
# plotDat[order(plotDat$importance),]
# 
# pdf("./figures/variable_importance.pdf",width=5,height=4)
# ggplot(plotDat) + 
#   geom_bar(aes(x=reorder(names,-importance), weight=importance)) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   ylab("Relative Importance") + xlab("Variable")
# dev.off()
# 
# scatter_dat <- left_join(plotDat,glm_dat,by="names")
# ggplot(scatter_dat) + geom_point(aes(x=abs(estimate),y=importance))

ex_sl <- as.numeric(exposure)

SL.fit <- SuperLearner(Y=exposure,X=covs,family = "binomial",
                       SL.library = c("SL.ranger","SL.glm","SL.gam"))
ps <- predict(SL.fit,onlySL=T)$pred

overlap_dat <- data.frame(exposure=exposure,ps=ps)

f1 <- ggplot(overlap_dat) + 
  geom_density(aes(x=ps,
                   group=factor(exposure),
                   fill=factor(exposure)),
               alpha=.75) +
  scale_fill_grey() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Propensity Score") +
  ylab("Estimated Density") +
  theme(text = element_text(size=17.5)) +
  guides(fill=guide_legend(title="Observed Exposure"))

pdf(here("figures","2019_1_26-PS_Overlap_NuMom_Figure1.pdf"),width = 5,height = 5)
f1
dev.off()

ggplot(overlap_dat) + 
  geom_boxplot(aes(x=ps,
                   group=factor(exposure),
                   fill=factor(exposure))) +
  scale_fill_grey()

ate(y=outcome, a=exposure, x=covs, nsplits=2, sl.lib=c("SL.ranger","SL.glm","SL.gam"))
att(y=outcome, a=exposure, x=covs, nsplits=2, sl.lib=c("SL.ranger","SL.glm","SL.gam"))

tml.fit <- tmle(Y=outcome,
                A=exposure,
                W=data.frame(covs),
                family="binomial",
                g.SL.library = c("SL.ranger","SL.glm","SL.gam"),
                Q.SL.library = c("SL.ranger","SL.glm","SL.gam"))
tml.fit

res <- ipsi(y = outcome,
            a = exposure,
            x.trt = covs,
            x.out = covs, 
            nsplits = 5, 
            delta.seq = seq(0.2,5,length.out=15), 
            id=id, 
            time=time)

plotDat <- res$res

plot_obj <- ggplot(plotDat) + scale_x_continuous(trans='log2') +
  geom_ribbon(aes(x=increment,ymin=ci.ll,ymax=ci.ul),
              fill="lightblue",color="lightblue",alpha=.5) + 
  geom_line(aes(x=increment,y=est)) + 
  geom_vline(xintercept=1,linetype="dashed") +
  ylab("Risk of Preeclampsia") +
  xlab("Change in Exposure Propensity (Odds Ratio)")

plot_obj

saveRDS(plot_obj,file="data/sga_g_whldens.rds")
beep(sound=8)
