### Started 2 February 2016 ###
### By Lizzie ### 

## Quick look at responses of species we collected tissue for genotyping from ## 

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## libraries
library(plyr)
library(dplyr)

setwd("~/Documents/git/projects/treegarden/budexperiments/analyses")

if(length(grep("danflynn", getwd()))>0){ setwd("~/Documents/git/buds/analyses") }


# get latest data
print(toload <- sort(dir("./input")[grep("Budburst Data", dir('./input'))], T)[1])

load(file.path("input", toload))

source("source/simpleplot.R")

# Analysis of where the leafout cuttings were
lx <- dx[dx$nl == 1,]

summary(lx)

##
setwd("~/Documents/git/projects/treegarden/bud-genetics/analyses")

if(length(grep("danflynn", getwd()))>0){ setwd("~/Documents/git/budgenetics/analyses") }

## get the tissue individuals
hftissue <- read.csv("input/HF_TISSUE_data - Tissue individuals, HF.csv", header=TRUE)
shtissue <- read.csv("input/St.Hip_TISSUE_data - Tissue individuals.csv", header=TRUE)

hftissvector <- paste(hftissue$Individual, "HF", sep="_")
shtissvector <- paste(shtissue$Individual, "SH", sep="_")

# Manually adding in POPGRA02_HF, not in HF_TISSUE_data sheet for some reason


tissues <- c(hftissvector, shtissvector, "POPGRA02_HF")

budtiss <- dx[which(dx$ind %in% tissues),]

howmanychill <- aggregate(budtiss["bday"], budtiss[c("chill","sp" )], FUN=length) # ACEPEN, FAGGRA, POPGRA, QUERUB, VIBLAN also have full chilling (so 5 species do and 5 species don't .... sigh)

# do a little formatting to get names to match Simon's naming
# my unname and lookup did not working, so doing it the hard way...
budtiss$jolyname <- NA
budtiss$jolyname[budtiss$sp=="SPIALB"] <- "Spiraea_alba"
budtiss$jolyname[budtiss$sp=="LONCAN"] <- "Lonicera_canadensis"
budtiss$jolyname[budtiss$sp=="PRUPEN"] <- "Prunus_pensylvanica"
budtiss$jolyname[budtiss$sp=="FAGGRA"] <- "Fagus_grandifolia"
budtiss$jolyname[budtiss$sp=="POPGRA"] <- "Populus_grandidentata"
budtiss$jolyname[budtiss$sp=="ACEPEN"] <- "Acer_pensylvanicum"
budtiss$jolyname[budtiss$sp=="QUERUB"] <- "Quercus_rubra"
budtiss$jolyname[budtiss$sp=="VIBLAN"] <- "Viburnum_lantanoides"
budtiss$jolyname[budtiss$sp=="VACMYR"] <- "Vaccinium_myrtilloides"
budtiss$jolyname[budtiss$sp=="ALNINC"] <- "Alnus_incana"
# budtiss$jolyname <- unname(lookupnames[budtiss$sp])
budtiss$newsite <- NA
budtiss$newsite[budtiss$site=="HF"] <- "MA"
budtiss$newsite[budtiss$site=="SH"] <- "QC"
budtiss$idnumber <- substr(budtiss$ind, 7, 8)

budtiss$jolyid <- paste(budtiss$jolyname, budtiss$newsite, budtiss$idnumber, sep="_")
changelabelsall <- subset(budtiss, select=c("ind", "sp", "jolyid")) # we'll use this later
changelabels <- changelabelsall[!duplicated(changelabelsall),]

foo <- lm(bday~chill+warm+photo, data=subset(budtiss, ind=="ACEPEN01_HF"))

## summarize the budburst data
budsummary.chill <-
      ddply(budtiss, c("ind", "chill"), summarise,
      mean = mean(bday, na.rm=TRUE), sd = sd(bday, na.rm=TRUE),
      sem = sd(bday, na.rm=TRUE)/sqrt(length(bday)))

budsummary.warm <-
      ddply(budtiss, c("ind", "warm"), summarise,
      mean = mean(bday, na.rm=TRUE), sd = sd(bday, na.rm=TRUE),
      sem = sd(bday, na.rm=TRUE)/sqrt(length(bday)))

budsummary.photo <-
      ddply(budtiss, c("ind", "photo"), summarise,
      mean = mean(bday, na.rm=TRUE), sd = sd(bday, na.rm=TRUE),
      sem = sd(bday, na.rm=TRUE)/sqrt(length(bday)))

names(budsummary.chill)[names(budsummary.chill)=="chill"] <- "main.effect"
names(budsummary.warm)[names(budsummary.warm)=="warm"] <- "main.effect"
names(budsummary.photo)[names(budsummary.photo)=="photo"] <- "main.effect"

budsummary1 <- rbind(budsummary.chill, budsummary.warm, budsummary.photo)
budsummary <- merge(budsummary1, changelabels, by="ind", all.x=TRUE, all.y=FALSE)

## summarize site effects
whichspp <- unique(budtiss$sp)
siteffects <- c()
for (i in c(1:length(whichspp))){
    sp.model <- lm(bday~site*photo*warm, data=subset(budtiss, sp==whichspp[i]))
    siteffects[i] <- summary(sp.model)$coef[2,1] # site effect
}
# SPIALB has strong site effect
siteffects.out <- data.frame(spp=unique(budtiss$sp), siteeffect=siteffects)

write.csv(budtiss, "output/indforGBS.csv", row.names=FALSE)
write.csv(budsummary, "output/budsummary.csv", row.names=FALSE)
write.csv(siteffects.out, "output/siteffects.csv", row.names=FALSE)



