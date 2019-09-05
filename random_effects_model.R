library(tidyverse)
library(lmerTest)
library(lsmeans)
library(multcomp)


setwd("/Volumes/GoogleDrive/My Drive/Currie_lab/Garden_Bacteria/Plant Defense Compound/GC_MS/headspace_sampling")

## load data ###

apinene <- read.csv("apinene_colony.csv")
linalool <- read.csv("linalool_by_colony.csv")

apinene$value <- apinene$value / 1000000
linalool$value <- linalool$value / 1000000

##random effects with colony as random effect##

ap_col <- lmer(value~time + treatment + (1|colony),data=apinene)
summary(ap_col)

lina <- lmer(value~time + treatment + (1|colony), data=linalool)
summary(lina)


##emmeans

ap.lm <-lm(value~time*treatment, data=apinene)
lina.lm <- lm(value~time*treatment, data=linalool)

apemm <- emmeans(ap.lm, "treatment", by="time")
apemm2 <-emmeans(ap.lm, "time", by="treatment")

summary(apemm)
pairs(apemm)
pairs(apemm2)

linaemm<- emmeans(lina.lm, "treatment", by="time")
linaemm2<- emmeans(lina.lm, "time", by="treatment")

summary(linaemm)
pairs(linaemm)
pairs(linaemm2)

summary(lina.lm)
anova(lina.lm)

## checking assumptions ##
par(mfrow=c(2,2))
plot(ap.lm)

par(mfrow=c(2,2))
plot(lina.lm)



