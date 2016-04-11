### Started 13 March 2016 ###
### By Lizzie ### 

## Happy pi day! ##
## Quick look at distance genotypes and budburst data ##

## Last updated on 26 March 2016: Got the PGLS to run but probably need to think of better analysis ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("~/Documents/git/projects/treegarden/genetics/analyses")

# if(length(grep("danflynn", getwd()))>0){ setwd("~/Documents/git/budgenetics/analyses") }

## source
source("duplicate.tips.R")

## libraries
library(caper)
library(nlme)
library(lme4)
library(ape) # for varcomp for nlme
# library(HLMdiag) # for varcomp for lme4, which is basically already in lme4 output

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

mytree$tip.label[which(!mytree$tip.label %in% alldater$ind)]

# same, just checking
(tissue.no.exp <- data.frame(id = mytree$tip.label[is.na(match(mytree$tip.label, alldater$ind))]) )

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
write.csv(tissue.no.exp, "output/tissue.no.expnotes.csv", row.names=FALSE)


# subset the data 
daterch0 <- subset(alldater, chill=="chill0")

##
## fit simple mixed-effects models (ignoring the phylogeny)
# with and without individual
##

# first I tried nlme package but with this we cannot have site/ind and sp crossed ...
# so I just did models with site OR ind (and both with species)
lme.mod.sp.ind <- lme(lday~warm*photo, random = list(ind=~1, sp=~1),
    data=daterch0, na.action=na.exclude)
lme.mod.sp.site <- lme(lday~warm*photo, random = list(site=~1, sp=~1),
    data=daterch0, na.action=na.exclude)

# interestingly the model comparisons show the model with site better
# this is probably because of DF penalties
anova(lme.mod.sp.ind, lme.mod.sp.site)

# when you look at the varcomp though ...
# the model with ind takes half the variance away from sp and gives it to ind
varcomp(lme.mod.sp.site, scale=TRUE)
varcomp(lme.mod.sp.ind, scale=TRUE) 

# let's try lme4 to get the exact model I think we have
mod.site.sp <- lmer(lday~warm*photo + (1|sp) + (1|site), data=daterch0,
    na.action=na.exclude)
mod.sp.ind <- lmer(lday~warm*photo + (1|sp) + (1|ind), data=daterch0,
    na.action=na.exclude)
mod.site.sp.ind <- lmer(lday~warm*photo + (1|sp) + (1|site/ind), data=daterch0,
    na.action=na.exclude)

summary(mod.site.sp) # site explains nada!
# so the below models are identical (wow)
summary(mod.sp.ind)
summary(mod.site.sp.ind)

varcomp.mer(mod.site.sp.ind)

##
## back to the phylogeny, get down the data needed and duplicate tips
##
daterch0$indX <- paste(daterch0$ind, rep(1:4, 102), sep="")
daterch0sm <- subset(daterch0, select=c("site", "sp", "ind", "treatcode", "warm",
    "photo", "lday", "indX"))

mytree.dup <- duplicate.tips.jd(mytree, 4, sep.char="")

# now try to run PGLS
compdat <- comparative.data(mytree.dup, daterch0sm, indX) 
pgls(lday~warm*photo, lambda="ML", data=compdat)

# why won't it run? Zero branch lengths....
# cheap trick: make all branch lengths 1
whee <- mytree.dup
whee$edge.length <- rep(1, length(whee$edge.length))
compdat.whee <-  comparative.data(whee, daterch0sm, indX) 
mod1 <- pgls(lday~warm*photo, lambda="ML", data=compdat.whee)

# better idea ...
# resolve tree and set min branch length to depth of tree/1000
whee2 <- multi2di(mytree.dup)
whee2$edge.length[whee2$edge.length==0] <- max(cophenetic(whee2))/2000
compdat2 <- comparative.data(whee2, daterch0sm, indX)
mod2 <- pgls(lday~warm*photo, lambda="ML", data=compdat2)
summary(mod2)

# randomize data within species
randat <- c()
sphere <- unique(daterch0sm$sp)
for (species in c(1:length(sphere))) {
    subby <- subset(daterch0sm, sp==sphere[species])
    randat <- c(randat, sample(subby$indX, length(subby$indX)))
}

daterch0sm$randind <- randat

# arghh! why on earth is caper dropping more rows (and 2 ind)
# in compat2 vs. compdat3 ??
compdat2 <- comparative.data(whee2, daterch0sm, indX) 
compdat3 <- comparative.data(whee2, daterch0sm, randind)

compdat2$data$ind[which(!unique(!compdat3$data$ind)
    %in% unique(compdat2$data$ind))]
length(compdat3$dropped$unmatched.rows)
length(compdat2$dropped$unmatched.rows)

# not useful yet but here's the randomized tips model
# JD pointed out since the treatments produce big effects
# we're forcing a ton of evolution at the tips ....
mod2rand <- pgls(lday~warm*photo, lambda="ML", data=compdat3)
summary(mod2rand)
