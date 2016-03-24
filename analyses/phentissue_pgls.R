### Started 13 March 2016 ###
### By Lizzie ### 

## Happy pi day! ##
## Quick look at distance genotypes and budburst data ## 

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("~/Documents/git/projects/treegarden/genetics/analyses")

if(length(grep("danflynn", getwd()))>0){ setwd("~/Documents/git/budgenetics/analyses") }

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

# same, just checking
( tissue.no.exp <- data.frame(id = mytree$tip.label[is.na(match(mytree$tip.label, alldater$jolyid))]) )

notes <- c(
  "Tissue individual, but not enough material, no terminal bud. Replaced with another one.",  # vacmyr02_HF, from common garden individuals google doc
  "Was not on common garden indiv sheet; note from HF_TISSUE: 'On path past nature trail, between 2 swamps, also collected fruit'",
  "Notes from Winter sampling 15: 'Tim and Harold searched for an hour, could not locate.'",
  "Mostly dead. only 1 semi-alive twig. Didn't clip. 12' tall, sick.",
  "Could not find in Winter 15, replaced the a new one",
  "Could not find in Winter 15, replaced the a new one",
  "Could not find in Winter 15, replaced the a new one",
  "Could not find in Winter 15, replaced the a new one",
  "Could not find in Winter 15, replaced the a new one")

tissue.no.exp <- data.frame(tissue.no.exp, notes)

# data.frame(dx$id[grep("POPGRA02_HF", dx$id)]) # present in dx dataframe, now fixed manually. 

compdat <- comparative.data(mytree, treat20.sm, jolyid) ## hmm, can't run with multiple values per species ....

# randomize data within species
