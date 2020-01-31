packages <- c("data.table","broom","gridExtra",
              "here","VIM","haven","skimr","beepr",
              "GGally","ranger","tmle","xtable")

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
         prepregbmi=as.numeric(prepregbmi),
         pct_emptyc=pct_emptyc)

table(main$education)

# main %>%
#   select(prepregbmi,age,dt_kcal,
#          pct_emptyc,f_totdens,v_totdens,
#          d_totdens,p_totdens,g_whldens,
#          g_nwhldens,fa_dens,
#          seaplant_dens,sodium_dens) %>%  ggpairs()

#' Define the exposure and outcome variables,
#'  create confounder vector
outcome <- main$outcome
exposure <- main$exposure
covs <- main %>% select(age,black,education,married,smokerpre,prepregbmi,insurpub,
                        pct_emptyc,d_totdens,p_totdens,f_totdens,sodium_dens,fa_dens,g_nwhldens,g_whldens)
#' incremental PS estimator requires 
#' an ID and time variable. For time-fixed 
#' exposure, set time to 1 for all 
id <- 1:nrow(main)
time <- rep(1,nrow(main))

#' Begin exploring exposure propensity
#'  via Exposure Model (GLM)
gFit <- glm(exposure~.,data=data.frame(covs),family=binomial("logit"))
glm_dat <- tidy(gFit)[-1,] #%>% select(term,estimate) %>% rename(names=term)
glm_dat

## Descriptives Tables
sum(exposure)
sum(1-exposure)
length(exposure)

t1 <- covs %>% 
  mutate(outcome=outcome,exposure=exposure) %>% 
  group_by(exposure) %>% 
  summarize_all(.,mean)

t2 <- covs %>% 
  mutate(outcome=outcome) %>% 
  summarize_all(.,mean) %>% 
  mutate(exposure=99)

xtable(t(rbind(t1,t2)))

#' Firt stacking algo for obtain PS 
#' for shift and overlap plots
ex_sl <- as.numeric(exposure)
SL.fit <- SuperLearner(Y=exposure,X=covs,family = "binomial",
                       SL.library = c("SL.ranger","SL.glm","SL.gam"))
ps <- predict(SL.fit,onlySL=T)$pred
colnames(ps) <- "pi_hat"

#' how many women woth PS = .5?
ps %>% as_tibble(ps) %>% filter(round(pi_hat,2)==.50) %>% count()

#' restrict to only unique PS values and add
#'  column for delta shift parameter 
ps <- as_tibble(round(ps,2)) %>% 
  distinct() 

ps <- cbind(rep(seq(0.2,5,length.out=15),each=nrow(ps)),
            unlist(t(rep(ps,15))))
ps <- data.frame(ps)
names(ps) <- c("delta","pi_hat")

#' construct shifted PS values
ps <- as_tibble(ps) %>% 
  mutate(pi_delta=(delta*pi_hat)/(delta*pi_hat + 1 - pi_hat))

#' generate shift plot 1, which tracks how each unique PS value
#' changes as a function of delta
shift_plot1 <- ggplot() + 
  geom_line(data=ps,aes(x=delta,y=pi_delta,group=pi_hat),
            color="darkgray",size=.2) +
  geom_vline(xintercept=1,color="black",size=.5) +
  scale_x_continuous(expand=c(0,0),trans="log2") +
  scale_y_continuous(expand=c(0,0)) +
  theme(text = element_text(size=17.5)) +
  ylab(expression(paste("Shifted PS, ", pi[delta]))) +
  xlab(expression(paste("Odds Ratio Change, ", delta)))

#' output plot to folder
pdf(here("figures","2019_1_26-Shift_Plot1-NuMom.pdf"),width = 5,height = 5)
shift_plot1
dev.off()

#' generate shiftplot 2, which shows how entire
#'  PS distribution changes under the extreme delta 
#'  values
ps <- predict(SL.fit,onlySL=T)$pred
colnames(ps) <- "pi_hat"
ps <- cbind(rep(c(0.2,5),each=nrow(ps)),
            rep(ps,2))
ps <- data.frame(ps)
names(ps) <- c("delta","pi_hat")
ps <- as_tibble(ps) %>% 
  mutate(pi_delta=(delta*pi_hat)/(delta*pi_hat + 1 - pi_hat))

#' create dataset with shifted PS values under
#' delta = 0.2 and 5
p1 <- ps %>%  
  select(delta,pi_delta)
#' add on observed PS values
p2 <- ps %>% 
  filter(delta == .2) %>% 
  select(delta,pi_hat) %>%
  mutate(delta=1) %>% 
  rename(pi_delta=pi_hat)

p <- rbind(p1,p2)
shift_plot2 <- ggplot(p) + 
  geom_density(aes(pi_delta,
                   group=delta,
                   fill=factor(delta)),
                  alpha=.5,bw=.0125) +
  scale_fill_grey() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(text = element_text(size=17.5)) +
  guides(fill=guide_legend(title=paste(expression("Shift OR")))) +
  ylab("Shifted PS Density") +
  xlab(expression(paste("Shifted PS, ",pi[delta])))

#' output plot to folder
pdf(here("figures","2019_1_26-Shift_Plot2-NuMom.pdf"),width = 5,height = 5)
shift_plot2
dev.off()

#' create combined two panel plot
pdf(here("figures","2019_1_26-Shift_Plot_combined-NuMom.pdf"),
    width = 10,
    height = 5)
  grid.arrange(shift_plot1,shift_plot2,ncol=2)
dev.off()

#' create PS overlap plot comparing exposed and unexposed
ps <- predict(SL.fit,onlySL=T)$pred
overlap_dat <- data.frame(exposure=exposure,ps=ps)

#' how low does PS get and how many ppl are there?
summary(overlap_dat$ps)
overlap_dat %>% filter(ps<.05)

#' PS overlap plot
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

#' ouput plot to folder
pdf(here("figures","2019_1_26-PS_Overlap-NuMom.pdf"),width = 5,height = 5)
f1
dev.off()


#' estimate ATE using AIPW
ate(y=outcome, a=exposure, x=covs, nsplits=2, sl.lib=c("SL.ranger","SL.glm","SL.gam"))

#' estimate ATE using TMLE
tml.fit <- tmle(Y=outcome,
                A=exposure,
                W=data.frame(covs),
                family="binomial",
                g.SL.library = c("SL.ranger","SL.glm","SL.gam"),
                Q.SL.library = c("SL.ranger","SL.glm","SL.gam"))
tml.fit

#' estimate IPS using ipsi function
res <- ipsi(y = outcome,
            a = exposure,
            x.trt = covs,
            x.out = covs, 
            nsplits = 5, 
            delta.seq = seq(0.2,5,length.out=15), 
            id=id, 
            time=time,
            sl.lib=c("SL.ranger","SL.glm","SL.gam"))

#' generate risk differences for Table 3

#' generate figure depicting preeclampsia risk over range of delta
plotDat <- res$res
plot_obj <- ggplot(plotDat) + 
  geom_ribbon(aes(x=increment,ymin=ci.ll,ymax=ci.ul),
              fill="lightgrey",color="lightgrey",alpha=.5) + 
  geom_line(aes(x=increment,y=est)) + 
  geom_vline(xintercept=1,linetype="dashed") +
  scale_x_continuous(expand = c(0,0),trans='log2',limits = c(.2,5)) +
  scale_y_continuous(expand = c(0,0),limits=c(.05,.1)) +
  theme(text = element_text(size=17.5)) +
  ylab("Risk of Preeclampsia") +
  xlab(expression(paste("Odds Ratio Change, ", delta)))

#' output figure to file
pdf(here("figures","2019_1_26-IPS_Estimate-NuMom.pdf"),width = 5,height = 5)
plot_obj
dev.off()



