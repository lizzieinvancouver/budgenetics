### Started 13 March 2016 ###
### By Lizzie ### 

## Happy pi day! ##
## Quick look at distance genotypes and budburst data ## 

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("~/Documents/git/projects/treegarden/genetics/analyses")

## libraries
library(ape)
library(caper)


alldater <- read.csv("output/indforGBS.csv", header=TRUE)
dater <- read.csv("output/budsummary.csv", header=TRUE) # see phentissue.R for where this file is created
mytree <- read.nexus("input/tree_intra.tre")

treat20 <- subset(dater, main.effect=="20" | main.effect=="15")
treat20.sm <- subset(treat20, select=c("jolyid", "mean", "sp", "main.effect"))
treat20.sm$main.effect <- as.numeric(treat20.sm$main.effect)

# wait, check for what's missing
treat20[which(!treat20$jolyid %in% mytree$tip.label),]
length(unique(treat20$jolyid))
mytree$tip.label[which(!mytree$tip.label %in% treat20$jolyid)]

# and check what's missing when you consider all treatments
mytree$tip.label[which(!mytree$tip.label %in% alldater$jolyid)]

compdat <- comparative.data(mytree, treat20.sm, jolyid) ## hmm, can't run with multiple values per species ....

# randomize data within species
