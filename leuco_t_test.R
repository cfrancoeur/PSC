library(tidyverse)
library(lmerTest)
library(lsmeans)
library(multcomp)


setwd("/Volumes/GoogleDrive/My Drive/Currie_lab/Garden_Bacteria/Plant Defense Compound/GC_MS/headspace_sampling/leuco_june_2019")

## load data ###

leuco_apinene <- read.csv("ap_leuco.csv")
leuco_linalool <- read.csv("linalool_leuco.csv")

t.test(ap ~ sample, data=leuco_apinene)
t.test(lina ~ sample, data=leuco_linalool)
