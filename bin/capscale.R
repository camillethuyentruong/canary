# install_github("vqv/ggbiplot")
# install.packages("raster")
# update.packages("ggplot2")

library(vegan)
library(mvabund)
library(VennDiagram)
library(venneuler)
library(indicspecies)
library(ggplot2)
library(raster)
library(grid)
library(stats)
library(devtools)
library(reshape2)
library(car)
library(MASS)
library(PerformanceAnalytics)
library(ggbiplot)
library(devtools)
library(iNEXT)
library(RColorBrewer)
library(gridExtra)
library(lme4)
library(lmerTest)
library(lsmeans)
library(GAD)

display.brewer.all(type="seq")
display.brewer.all(type="div")
display.brewer.all(type="qual")


setwd("/Users/Camille/Documents/Luis_Quijada/R")


# IMPORT DATA
species_table<-read.csv("/Users/Camille/Documents/Luis_Quijada/R/species_table.csv", row.names=1, check.names=FALSE)
head(species_table) # check
row.names(species_table)
colnames(species_table)

meta_table<-read.csv("/Users/Camille/Documents/Luis_Quijada/R/meta_table.csv",row.names=1,check.names=FALSE)
head(meta_table) # check
row.names(meta_table)
colnames(meta_table)

# Check to ensure that the samples in meta_table are in the same order as in species_table
meta_table<-meta_table[rownames(species_table),]

# Transform numerical values
num <- c("BranchLenght","BranchPerimeter","DecayRange","PercBarkAttached","Elevation","CanopyHeight","PercCanopyCover","PerLeafLitter")
num_table <- subset(meta_table, select = num)
num_table
scale_table=scale(num_table, center=TRUE, scale = TRUE)
scale_table

cat <- c("HostPlant","VegetationSubtype","VegetationType","LocName","Exposure","CollectionMonth","CollectionYear")
cat_table <- subset(meta_table, select = cat)
cat_table

##########################

# DISTANCE-BASED REDUDANCY ANALYSIS

# 2 tables: species matrix (raup) and envir variables (data)

raup<-raupcrick(species_table, null="r1", nsimul=999, chase=FALSE)

data=data.frame(meta_table, row.names=rownames(meta_table))
data

raup.min<-capscale(raup ~ 1 , data)


# Model fitting with all factors
# PercBarkAttached removed because of missing data
raup.all<-capscale(raup ~ BranchLenght+BranchPerimeter+DecayRange+Elevation+CanopyHeight+PercCanopyCover+PerLeafLitter+HostPlant+VegetationSubtype+VegetationType+LocName+Exposure+CollectionMonth+CollectionYear, data)
RsquareAdj(raup.all)$adj.r.squared # 0.3393267

# Test the global model (including all explanatory variables from which you plan to do selection) to see whether it is significant
# if it is not, the forward selection should not be done.
plot(raup.all)
RsquareAdj (raup.all)$adj.r.squared
anova (raup.all)

# check for colinearity of soil/elevation variables (so we can use them independently)
# values over 10 indicate redundant constraints, otherwise can be used independently
vif.cca (raup.all) 



# selection based on adjR2 and P value (Pin=0.05)
raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction = "both", Pin = 0.05, R2scope = TRUE, permutations = how(nperm = 999), 
                    trace = TRUE, R2permutations = 999)

RsquareAdj(raup.db)$adj.r.squared

# Assumptions in linear models:
stressplot(raup.db)
qqnorm(residuals(raup.db))
summary(raup.db)
plot(raup.db)
anova (raup.db)
anova(raup.db, by = "margin")
anova(raup.db, by = "terms")

# Host plant, F=6.6287***
# Loc Name, F=2.2876***
# CollectionMonth F=1.4436***
