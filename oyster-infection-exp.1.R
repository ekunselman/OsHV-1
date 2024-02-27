# Oyster infection experiment 1

# kaplan meier curves --------------------

# This code will plot survival data and test survival difference between two or more groups
# data sheet contains timepoints at which individual was "censored" or checked for mortality
# and then a column that specifies whether the individual was dead or alive
# my data contains 0 for alive and 1 for death
# the survivor function asses the probability that an individual survives from the time of origin (infection) until some future specific time

# model temperatures (all tanks)
# model antibiotics versus non (by temperature)
# model replicates by condition and temp to see how closely they follow

#transition to this tutorial http://www.sthda.com/english/wiki/survival-analysis-basics

# load data
library(readr)
oshv1_km<- read_csv("/Users/emilykunselman/Documents/3 - OsHV-1 Spat SDB Project/Data/Infection experiment #1/OSHV-1 infection shared data sheets - kaplan meier mortality.csv")

# install packages! 

# library packages
library(survminer)
library(survival)

# the function survfit computes the KM survivial estimate
# time_days is time at which indiviudal was censored
# death is 0 or 1 - was individual alive or not
# temperature is my grouping of interest
km.model.temp<-survfit(Surv(time_days, death) ~ temperature, data = oshv1_km)
print(km.model.temp)
# n is number of occurrences (we miscounted one 18 degree tank at the beginning)
# events is number of deaths recorded for that group
# median is median survival(time at which 50% have died) for that group and then you have confidence intervals around median
summary(km.model.temp)
# summarizes for each timepoint of censoring: number of individuals at risk for each temp, number of deaths, total survival, and CIs

# plot

ggsurvplot(km.model.temp,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           surv.median.line = "hv",
           ggtheme = theme_dark(base_size = 15),
           palette = "YlOrRd",
           font.x = "black",
           font.tickslab = "black")
#colorblind friendly
#other option below

cbPalette = c("#CC79A7", "#E69F00", "#E69F00", "#009E73")

# note median black dotted line shows median for each temperature and the time at whcih it was achieved

# Temp x Time
km.model.temp.2<-survfit(Surv(tempXtime, death) ~ temperature, data = oshv1_km)


ggsurvplot(km.model.temp.2,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           surv.median.line = "hv",
           ggtheme = theme_dark(base_size = 15),
           palette = "YlOrRd")
# test pairwise differences in temperature

# Log Rank test comparing survival
surv_diff <- survdiff(Surv(time_days, death) ~ temperature, data = oshv1_km)
surv_diff

# pairwise test (24 and 21 not significantly different)
# temperature
p_surv_diff<-pairwise_survdiff(Surv(time_days, death) ~ temperature, data = oshv1_km)
p_surv_diff

#temperature x time
p_surv_diff<-pairwise_survdiff(Surv(tempXtime, death) ~ temperature, data = oshv1_km)
p_surv_diff

# antibiotics
# compare by temperature

# filter to 

# complex curves: Temperature and antibiotics

km.model.temp.AB <- survfit( Surv(time_days, death) ~ temperature + antibiotics,
                 data = oshv1_km )
# plot with fecat by temperature

ggsurv <- ggsurvplot(km.model.temp.AB, conf.int = TRUE,
                     ggtheme = theme_dark(base_size = 15), palette = "YlOrRd")
ggsurv$plot +theme_dark(base_size = 15) + 
  theme (legend.position = "right")+
  facet_grid(temperature ~ .)

# is antibiotic exposure impacting death?
temp18<-subset(oshv1_km, temperature == "18")
p_surv_diff<-pairwise_survdiff(Surv(time_days, death) ~ antibiotics, data = temp18)
p_surv_diff
 # no

temp21<-subset(oshv1_km, temperature == "21")
p_surv_diff<-pairwise_survdiff(Surv(time_days, death) ~ antibiotics, data = temp21)
p_surv_diff
  # no

temp24<-subset(oshv1_km, temperature == "24")
p_surv_diff<-pairwise_survdiff(Surv(time_days, death) ~ antibiotics, data = temp24)
p_surv_diff
  # no

# add in replicate variable
km.model.temp.AB.rep <- survfit( Surv(time_days, death) ~ temperature + antibiotics + replicate,
                             data = oshv1_km )

ggsurv <- ggsurvplot(km.model.temp.AB.rep, conf.int = TRUE,
                     ggtheme = theme_bw())
ggsurv$plot +theme_bw() + 
  theme (legend.position = "right")+
  facet_grid(temperature ~ antibiotics)

# Daysi Senior Thesis Data ---------------------
# 18C antibiotic vs non
# load data
library(readr)
daysi_km<- read_csv("/Users/emilykunselman/Documents/3 - OsHV-1 Spat SDB Project/Data/Daysi_ABX_18_mortality.csv")

# library packages
library(survminer)
library(survival)

# the function survfit computes the KM survivial estimate
# time_days is time at which indiviudal was censored
# death is 0 or 1 - was individual alive or not
# temperature is my grouping of interest
km.model.ABX<-survfit(Surv(Days_since_exposure, death) ~ antibiotics, data = daysi_km)
print(km.model.ABX)
# n is number of occurrences (we miscounted one 18 degree tank at the beginning)
# events is number of deaths recorded for that group
# median is median survival(time at which 50% have died) for that group and then you have confidence intervals around median
summary(km.model.ABX)
# summarizes for each timepoint of censoring: number of individuals at risk for each temp, number of deaths, total survival, and CIs

# plot

ggsurvplot(km.model.ABX,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme = theme_dark(base_size = 15),
           palette = c("#FED976", "#FEB24C"))
library(RColorBrewer)
brewer.pal(n = 8, name = "YlOrRd")
"#FFFFCC" "#FFEDA0" "#FED976" "#FEB24C" "#FD8D3C" "#FC4E2A" "#E31A1C" "#B10026"

# qPCR OsHV-1 Load figures -----------------

setwd("~/Documents/3 - OsHV-1 Spat SDB Project/Data/Infection experiment #1")
library(readr)
oys_qPCR<- read_csv("oyster_only_metadata_qpcr.csv")

# filter to OsHV-1 only
oys_qPCR_oshv1<-subset(oys_qPCR, condition=='oshv1')
# plot dead vs alive
library(ggplot2)
library(RColorBrewer)
ggplot(oys_qPCR_oshv1, aes(x=dead, y=log_oshv1_copies)) +
  geom_violin()+
  ggtitle("OsHV-1 load in dead oysters vs remaining live oysters")

# are copy numbers different in live verus dead?
kruskal.test(log_oshv1_copies ~ dead,
             data = oys_qPCR_oshv1
)
  #chi-squared = 104.18, df = 1, p-value < 2.2e-16

wilcox.test(log_oshv1_copies ~ dead,
            data = oys_qPCR_oshv1)

# filter to dead oysters in oshv1 only
oys_qPCR_oshv1_dead<-subset(oys_qPCR, dead=='dead')
oys_qPCR_oshv1_dead$temperature <- as.factor(oys_qPCR_oshv1_dead$temperature)
# plot temperature comparison
p1<-ggplot(oys_qPCR_oshv1_dead, aes(x=temperature, y=log_oshv1_copies, fill = temperature)) +
  geom_violin()+
  geom_jitter(width = 0.25)+
  ggtitle("OsHV-1 load, dead oysters, by temperature")+
  ylab("log(OsHV-1 copies/g)")+
  ylim(-1,8)+
  theme_dark(base_size = 15)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#FFCC00", "Orange", "Red"))
p3<-ggplot(oys_qPCR_oshv1_dead, aes(x=temperature, y=log_oshv1_copies)) +
  geom_boxplot()+
  ggtitle("OsHV-1 load, dead oysters, by temperature")+
  ylim(-1,8)

# filter to live oysters in oshv1 only
oys_qPCR_oshv1_alive<-subset(oys_qPCR, dead=='alive')
oys_qPCR_oshv1_alive$temperature <- as.factor(oys_qPCR_oshv1_alive$temperature)
# plot temperature comparison
p2<-ggplot(oys_qPCR_oshv1_alive, aes(x=temperature, y=log_oshv1_copies, fill = temperature)) +
  geom_violin()+
  geom_jitter(width = 0.25)+
  ggtitle("OsHV-1 load, live oysters, by temperature")+
  ylim(-1,8)+
  theme_dark(base_size = 15)+
  theme(axis.title.y = element_blank(), legend.position = "none")+
  scale_fill_brewer(palette = "YlOrRd")
p4<-ggplot(oys_qPCR_oshv1_alive, aes(x=temperature, y=log_oshv1_copies)) +
  geom_boxplot()+
  ggtitle("OsHV-1 load, live oysters, by temperature")+
  ylim(-1,8)+
  theme(axis.title.y = element_blank())

library(gridExtra)
grid.arrange(p1,p2, ncol = 2)
grid.arrange(p3,p4, ncol = 2)

# get legend from p2
ggplot(oys_qPCR_oshv1_alive, aes(x=temperature, y=log_oshv1_copies)) +
  geom_violin(aes(color = temperature))+
  ggtitle("OsHV-1 load, live oysters, by temperature")+
  ylim(-1,8)+
  theme(axis.title.y = element_blank())+
  scale_color_manual(values = c("Green", "Purple", "Orange", "Red"))

# Run anovas to determine significant differences in pathogen load between temps

# are data normally distributed?
hist(oys_qPCR_oshv1_dead$log_oshv1_copies) # no
library(car)
qqPlot(oys_qPCR_oshv1_dead$log_oshv1_copies) # a few outliers
hist(oys_qPCR_oshv1_alive$log_oshv1_copies) # no
qqPlot(oys_qPCR_oshv1_alive$log_oshv1_copies) # many outliers

# need to run Kruskal Wallis rather than ANOVA

kruskal.test(log_oshv1_copies ~ temperature,
             data = oys_qPCR_oshv1_dead
)
  # chi-squared = 3.5675, df = 2, p-value = 0.168

kruskal.test(log_oshv1_copies ~ temperature,
             data = oys_qPCR_oshv1_alive
)

  # chi-squared = 1.3759, df = 3, p-value = 0.7112


### RESULT: no difference in oshv-1 load by temperature, in either alive or dead oysters

# clearly though, viral loads are vary high in dead oysters, but lower and variables in live oysters
