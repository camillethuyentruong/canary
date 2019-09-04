# MODELIZATION

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
library(dplyr)

display.brewer.all(type="seq")
display.brewer.all(type="div")
display.brewer.all(type="qual")


setwd("/Users/Camille/Documents/Luis_Quijada/Canary/Git")


# IMPORT DATA
species_table<-read.csv("data/species_table.csv", row.names=1, check.names=FALSE)
head(species_table) # check
row.names(species_table)
colnames(species_table)

meta_table<-read.csv("data/meta_table.csv",row.names=1,check.names=FALSE)
head(meta_table) # check
row.names(meta_table)
colnames(meta_table)

# Check to ensure that the samples in meta_table are in the same order as in species_table
meta_table<-meta_table[rownames(species_table),]


# TRANSFORM DATA
num <- c("BranchLenght","BranchPerimeter","DecayRange","PercBarkAttached","Elevation","CanopyHeight","PercCanopyCover","PerLeafLitter")
num_table <- subset(meta_table, select = num)
num_table
scale_table=scale(num_table, center=TRUE, scale = TRUE)
scale_table


cat <- c("HostPlant","VegetationSubtype","VegetationType","LocName","Exposure","CollectionMonth","CollectionYear")
cat_table <- subset(meta_table, select = cat)
cat_table


species_table_merge <- species_table %>% 
  group_by(meta_table$HostPlant, meta_table$LocName) %>% 
  summarise_all(funs(sum))
species_table_merge
species_table_merge <- species_table_merge[ -c(1:2) ]
species_table_merge

meta_table_merge <- meta_table %>% 
  group_by(HostPlant, LocName, VegetationType, VegetationSubtype, Exposure) %>% 
  summarise_if(is.numeric, mean, keepFACTORS = TRUE)
meta_table_merge


# MAKE MATRICES
# 2 tables: species matrix (raup) and envir variables (data)

raup<-raupcrick(species_table, null="r1", nsimul=999, chase=FALSE)

data=data.frame(meta_table, row.names=rownames(meta_table))
data


binary_table=decostand(species_table_merge,method="pa") 
binary_table
raup_merge<-raupcrick(binary_table, null="r1", nsimul=999, chase=FALSE)

data_merge=data.frame(meta_table_merge, row.names=rownames(meta_table_merge))
data_merge



#####################

# VARIATION PARTITIONING
varp <- varpart(raup, ~ data$HostPlant, ~data$VegetationSubtype, ~data$CollectionMonth)
plot(varp, Xnames=c("HostPlant", "VegSubtype", "CollMonth"))



# DISTANCE-BASED REDUDANCY ANALYSIS

raup.min<-capscale(raup ~ 1 , data)

# Model fitting with all factors
# PercBarkAttached removed because of missing data
raup.all<-capscale(raup ~ BranchLenght+BranchPerimeter+DecayRange+Elevation+CanopyHeight+PercCanopyCover+PerLeafLitter+HostPlant+VegetationSubtype+VegetationType+LocName+Exposure+CollectionMonth+CollectionYear, data)
RsquareAdj(raup.all)$adj.r.squared # 0.3393267

# selection based on adjR2 and P value (Pin=0.05)
raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction = "both", Pin = 0.05, R2scope = TRUE, permutations = how(nperm = 999), 
                    trace = TRUE, R2permutations = 999)

# Host plant, F=6.6287***
# Loc Name, F=2.2876***
# CollectionMonth F=1.4436***
# DecayRange F=3.0478***
# CollectionYear F=1.6228**
# BranchPerimeter F=1.717*



# Removed LocName, VegetationType, CollYear
raup.all<-capscale(raup ~ HostPlant+VegetationSubtype+Exposure+CollectionMonth+BranchLenght+BranchPerimeter+DecayRange+Elevation+CanopyHeight+PercCanopyCover+PerLeafLitter, data)
RsquareAdj(raup.all)$adj.r.squared # 0.3113949

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction = "both", Pin = 0.05, R2scope = TRUE, permutations = how(nperm = 999), 
                    trace = TRUE, R2permutations = 999)

# HostPlant F=6.6279***
# VegetationSubtype F=2.8397***
# CollectionMonth F=1.7901***
# DecayRange F=3.2493***
# PercCanopyCover F=2.9168***
# Elevation F=2.3192**
# Exposure  F=2.7415***
# PerLeafLitter F=1.809*
# BranchLenght F=1.7505*
# CanopyHeight F=1.6488*


# Only envir variables
raup.db<-capscale(raup ~ Elevation+HostPlant+VegetationSubtype+Exposure+CollectionMonth,data)
RsquareAdj(raup.all)$adj.r.squared # 0.3055888

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction = "both", Pin = 0.05, R2scope = TRUE, permutations = how(nperm = 999), 
                        + trace = TRUE, R2permutations = 999)

# HostPlant F=6.6289***
# VegetationSubtype F=2.8399***
# CollectionMonth F=1.7902***
# Exposure F=2.7284***
# Elevation F=3.01***


#####################

# MERGED DATA

# VARIATION PARTITIONING
varp <- varpart(raup_merge, ~ data_merge$HostPlant, ~data_merge$VegetationSubtype, ~data_merge$Exposure) 
plot(varp, Xnames=c("HostPlant", "VegSubtype", "Exposure"))



# DISTANCE-BASED REDUDANCY ANALYSIS

raup.min<-capscale(raup_merge ~ 1 , data_merge)

# Model fitting with all factors
# PercBarkAttached removed because of missing data
# No more CollMonth and CollYear
raup.all<-capscale(raup_merge ~ BranchLenght+BranchPerimeter+DecayRange+Elevation+CanopyHeight+PercCanopyCover+PerLeafLitter+HostPlant+VegetationSubtype+VegetationType+LocName+Exposure, data_merge)
RsquareAdj(raup.all)$adj.r.squared # 0.6308713

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction = "both", Pin = 0.05, R2scope = TRUE, permutations = how(nperm = 999), 
                    trace = TRUE, R2permutations = 999)

# Host plant, F=5.7918***


# Partialling variation from Site
raup.all<-capscale(raup_merge ~ BranchLenght+BranchPerimeter+DecayRange+Elevation+CanopyHeight+PercCanopyCover+PerLeafLitter+HostPlant+VegetationSubtype+VegetationType+Exposure+Condition(LocName), data_merge)
RsquareAdj(raup.all)$adj.r.squared # 0.3254982

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction = "both", Pin = 0.05, R2scope = TRUE, permutations = how(nperm = 999), 
                    trace = TRUE, R2permutations = 999)

# none

# Removed LocName, VegetationType
raup.all<-capscale(raup_merge ~ HostPlant+VegetationSubtype+Exposure+BranchLenght+BranchPerimeter+DecayRange+Elevation+CanopyHeight+PercCanopyCover+PerLeafLitter, data_merge)
RsquareAdj(raup.all)$adj.r.squared # 0.5796248

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction = "both", Pin = 0.05, R2scope = TRUE, permutations = how(nperm = 999), 
                    trace = TRUE, R2permutations = 999)

# HostPlant F=5.7918***
# VegetationSubtype F=2.3375***
# Exposure F=4.0472***
# Elevation F=2.7589**
# PerLeafLitter F=2.6816*
# CanopyHeight F=2.5644*


# Only envir variables
raup.all<-capscale(raup_merge ~ Elevation+HostPlant+VegetationSubtype+Exposure,data_merge)
RsquareAdj(raup.all)$adj.r.squared # 0.5685693

raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction = "both", Pin = 0.05, R2scope = TRUE, permutations = how(nperm = 999), 
                    trace = TRUE, R2permutations = 999)

# HostPlant F=5.7918***
# VegetationSubtype F=2.3375***
# Exposure F=4.0472**



# Assumptions in linear models:
stressplot(raup.db)
qqnorm(residuals(raup.db))
summary(raup.db)
plot(raup.db)
anova (raup.db)
anova(raup.db, by = "margin")
anova(raup.db, by = "terms")

# check for colinearity of soil/elevation variables (so we can use them independently)
# values over 10 indicate redundant constraints, otherwise can be used independently
vif.cca (raup.all) 



#####################

# PLOTS

# Dataframe of CAP scores
CAP=data.frame(scores(raup.db,display=c("sites")), HostPlant=data_merge$HostPlant, VegSubtype=data_merge$VegetationSubtype, Exposure=data_merge$Exposure)
CAP


p <- ggplot(CAP)

p + geom_point(mapping = aes(x = CAP1, y = CAP2, colour=VegSubtype, shape = Exposure), size=3) +
  coord_fixed() +
  theme_bw()

  labs(x = "CAP1 (71.78% of variance explained)") +
  labs(y = "CAP1 (12.25% of variance explained)") +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12)) +
  theme(legend.title=element_text(size=12) , legend.text=element_text(size=8)) +
  theme(legend.key = element_blank(),  #removes the box around each legend item
        legend.position = "right", #legend at the bottom
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.box.just = "centre") +
  theme(panel.grid.major = element_blank(), # for paper
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
