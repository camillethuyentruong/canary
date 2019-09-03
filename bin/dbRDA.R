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

RsquareAdj(raup.db)$adj.r.squared # 0.3392523

# Host plant, F=6.6287***
# Loc Name, F=2.2876***
# CollectionMonth F=1.4436***
# DecayRange F=3.0478***
# CollectionYear F=1.6228**
# BranchPerimeter F=1.717*

# Assumptions in linear models:
stressplot(raup.db)
qqnorm(residuals(raup.db))
summary(raup.db)
plot(raup.db)
anova (raup.db)
anova(raup.db, by = "margin")
anova(raup.db, by = "terms")





# Model fitting with best variables
# Controlling variation explained by site
# Veg subtype nested within veg type
raup.all<-capscale(raup ~ Elevation+CanopyHeight+PercCanopyCover+PerLeafLitter+HostPlant+VegetationType/VegetationSubtype+Exposure+CollectionMonth+CollectionYear+Condition(LocName),data)
RsquareAdj(raup.all)$adj.r.squared # 0.05829149

# selection based on adjR2 and P value (Pin=0.05)
raup.db<-ordiR2step(raup.min, scope = formula (raup.all), direction = "both", Pin = 0.05, R2scope = TRUE, permutations = how(nperm = 999), 
                    trace = TRUE, R2permutations = 999)

RsquareAdj(raup.db)$adj.r.squared


# model fitting
raup.all<-capscale(raup ~ SoilMoist+pH+C+N+Nav+P+Elevation+Exposure+Condition(Transect), data)
RsquareAdj (raup.all)$adj.r.squared # 0.6040367, good!!!

raup.all2<-capscale(raup ~ (SoilMoist+C+N+Nav+P+pH)*Elevation+Exposure+Condition(Transect), data)
RsquareAdj (raup.all2)$adj.r.squared # 0.6823664

raup.all3<-capscale(raup ~ (SoilMoist+C+N+Nav+P+pH)*Elevation*Exposure+Condition(Transect), data)
RsquareAdj (raup.all3)$adj.r.squared # 0.723

anova(raup.all,raup.all2,raup.all3) # better with interactions of elevation and exposure

plot(raup.all3)
qqnorm(resid(raup.all3))
anova(raup.all3)
anova(raup.all3, by = "margin")

raup.db<-ordiR2step(raup.min, scope = formula (raup.all),direction =  "both", Pin = 0.05, R2scope = TRUE, permutations = how(nperm = 999), 
                    trace = TRUE, R2permutations = 999)



#####################

# VARIATION PARTITIONING
varp <- varpart(raup, ~ data$pH, ~data$Site) 
plot(varp)


#####################

# PLOTS

# Dataframe of CAP scores
CAP=data.frame(scores(raup.db,display=c("sites")),Elevation=data$Elevation,Transect=data$Transect, pH=meta_table$pH)
CAP

# Generate Elipse (for elevation)
CAPscores=data.frame(scores(raup.all,display=c("sites")))
CAPscores

ord<-ordiellipse(CAPscores, as.factor(data$Elevation) ,display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# Generate ellipse points
df_ell <- data.frame()
for(g in levels(data$Elevation)){
  if(g!="" && (g %in% names(ord))){
    
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(data[data$Elevation==g,],
                                                     veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                  ,Elevation=g))
  }
}

head(df_ell)


# Create dataframe for vectors (pH)
vector=scores(raup.all,display=c("bp"))

rownames(vector)[2:3]->remove # select only pH
vector=vector[!rownames(vector) %in% remove, ]
vector=data.frame(vector)
vector

# Plot
plot(raup.all2)

p <- ggplot(CAP)

p + geom_point(mapping = aes(x = CAP1, y = CAP2, colour=Elevation, shape = Transect), size=3) +
  geom_segment(aes(x = 0, y = 0, xend = CAP1, yend = CAP2, colour = "segment"), data = vector,
               arrow = arrow(length = unit(0.25, "cm")), colour = "black", size = 1) +
  coord_fixed() +
  theme_bw()

p + geom_point(mapping = aes(x = CAP1, y = CAP2, colour=pH, shape = Transect), size=3) +
  geom_path(data=df_ell, aes(x=CAP1, y=CAP2,fill=Elevation), size=0.75, linetype=2,color="darkgrey") +
  coord_fixed() +
  theme_bw() +
  scale_colour_gradient(low = "gold", high = "dark blue", space = "Lab", guide = "colourbar") +
  labs(colour = "pH") #another way to set the labels, in this case, for the colour legend

p + geom_point(mapping = aes(x = CAP1, y = CAP2, colour=pH, shape = Elevation, size=Transect)) +
  coord_fixed() +
  theme_bw() +
  scale_colour_gradient(low = "gold", high = "dark blue", space = "Lab", guide = "colourbar") +
  labs(colour = "pH") + #another way to set the labels, in this case, for the colour legend
  
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










##################################
# EFFECT OF ELEVATION ON Fungal communities
# permanova

d = vegdist(species_table_merge_ecm, "raup")

elevation_Adonis = adonis(d ~ X1, grouping_info_merge) 
elevation_Adonis

exposure_Adonis = adonis(d ~ X2, grouping_info_merge)  
exposure_Adonis

transect_Adonis = adonis(d ~ X3, grouping_info_merge)  
transect_Adonis


# differences within/among sites
row.names(species_table)
mid <- species_table[1:48,]
treeline <- species_table[49:96,]
low <- species_table[97:156,]

meta_mid <- meta_table[1:48,]
meta_treeline <- meta_table[49:96,]
meta_low <- meta_table[97:156,]

site_anosim = anosim(mid, meta_mid$SiteName,permutations = 999, distance = "raup")
site_anosim

# Significatif: dissimilarities among sites greater than within sites!!





# CORRELATION AMONG ENV/ENZ FACTORS

# Most highly correlated factors
mosthighlycorrelated <- function(mydataframe,numtoreport)
{
  # find the correlations
  cormatrix <- cor(mydataframe)
  # set the correlations on the diagonal or lower triangle to zero,
  # so they will not be reported as the highest ones:
  diag(cormatrix) <- 0
  cormatrix[lower.tri(cormatrix)] <- 0
  # flatten the matrix into a dataframe for easy sorting
  fm <- as.data.frame(as.table(cormatrix))
  # assign human-friendly names
  names(fm) <- c("First.Variable", "Second.Variable","Correlation")
  # sort and print the top n correlations
  head(fm[order(abs(fm$Correlation),decreasing=T),],n=numtoreport)
}

mosthighlycorrelated(scale(log(env_all_table_merge), center = TRUE, scale = TRUE), 20)



# Pairwise comparisons
pairs(data=scale_all_merge,~Nav+NAG+LEU, 
      main="Simple Scatterplot Matrix")


# FOR PAPER
# SCATTERPLOT MATRIX

colnames(select_table_merge) # change names of variables
colnames(select_table_merge) <- c("NAG", "GLUCU", "BGLU", "PHOS", "CEL", "XIL", "AGLU", "LEU","elevation") # change names of taxonomic table
colnames(select_table_merge)

c=scale(log(env_table_merge), center = TRUE, scale = TRUE)
scatterplotMatrix(c) # groups = (grouping_info[,1])

cor.prob <- function (X, dfr = nrow(X) - 2) {  # Adding corrrelation coefficients and p-values
  R <- cor(X, use="pairwise.complete.obs",method = c("spearman"))
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr/(1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R[row(R) == col(R)] <- NA
  R}

flattenSquareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(i = rownames(m)[row(m)[ut]],
             j = rownames(m)[col(m)[ut]],
             cor=t(m)[ut],
             p=m[ut])}

cor(c) # correlation matrix
cor.prob(c) # correlation matrix with p-values
flattenSquareMatrix(cor.prob(c)) # "flatten" that table
chart.Correlation(c) # plot the data



#PLOT CORRELATIONS ENV/ENZ FACTORS WITH ELEVATION
# DBH
p = ggplot(meta_table) + geom_boxplot(mapping = aes(x = (grouping_info[,1]), y = DBH, fill=(grouping_info[,1])))

p + labs(y ="DBH [cm]", x = "", size=20) +
  theme_bw() +
  scale_fill_discrete(breaks=c("L","M","T"), labels=c("Lowland", "Mid-elevation", "Treeline")) +
  guides(fill=guide_legend(title=NULL)) + # title=NULL to remove legend title
  theme(legend.key = element_blank(),  #removes the box around each legend item
        legend.position = c(0.5,0.99), 
        legend.direction = "horizontal",
        text = element_text(size = 25),
        legend.text = element_text(size = 20, face="bold"))


par(bg=NA)
p=plot(pH ~ Elevation, meta_table, ylab="pH", xlab="")
p+title("")


# Significant variables (from scatterplot)
aov <- aov(pH ~ Elevation, meta_table) # (one-way anova)
summary(aov)
TukeyHSD(aov) #pairwise

a<-lmer(SoilMoist ~ Elevation + (1|TransectName), data=meta_table) # effect of elevation after after controlling for the specific variation of each transect
lmerTest::summary(a)
lmerTest::anova(a) # type III anova
lmerTest::rand(a) # analysis of the random part of a mixed model, LRT (likelihood ratio test)
st = lmerTest::step(a) # pairewise comparison
plot(st)
st

par (mfrow = c(1,3)) # this function divides plotting window into columns

p=plot(SoilMoist ~ Elevation, meta_table, ylab="% Soil Moisture", xlab="")
p+title("")

p=plot(pH ~ Elevation, meta_table, ylab="pH", xlab="")
p+title("")

p=plot(LEU ~ Elevation, meta_table, ylab="Leucin-aminopeptidase", xlab="")
p+title("")


# Regressions (true elevation values)
g1=ggplot(meta_table_merge, aes(x=ElevNum, y=SoilMoist)) +
  geom_point(size=3) +
  theme_classic() +
  geom_smooth(method=lm, col="red", se=FALSE) + # Add linear regression line  (by default includes 95% confidence region, use se to remove)
  labs(y ="Soil Moisture [%]", x = "Elevation [m]") +
  ggtitle("") +
  theme(plot.title = element_text(size=25, face="bold", margin = margin(0, 0, 0, 0)))

g2=ggplot(meta_table_merge, aes(x=ElevNum, y=pH)) +
  geom_point(size=3) +
  theme_classic() +
  geom_smooth(method=lm, col="red", se=FALSE) + # Add linear regression line  (by default includes 95% confidence region, use se to remove)
  labs(y ="pH", x = "Elevation [m]") +
  ggtitle("") +
  theme(plot.title = element_text(size=25, face="bold", margin = margin(0, 0, 0, 0)))

g3=ggplot(meta_table_merge, aes(x=ElevNum, y=LEU)) +
  geom_point(size=3) +
  theme_classic() +
  geom_smooth(method=lm, col="red", se=FALSE) + # Add linear regression line  (by default includes 95% confidence region, use se to remove)
  labs(y ="Leucine-aminopeptidase", x = "Elevation [m]") +
  ggtitle("") +
  theme(plot.title = element_text(size=25, face="bold", margin = margin(0, 0, 0, 0)))


grid.arrange(g1,g2,g3, nrow=1, ncol=3)


# PLOT ENV/ENZ FACTORS WITH ELEVATION AND EXPOSITION
p = ggplot(meta_table) + geom_boxplot(mapping = aes(x = Elevation, y = pH, fill=Exposure))
p + theme_bw() + labs(y ="pH", x = "", size=20) +
  guides(fill=guide_legend(title="Exposure")) + # title=NULL to remove legend title
  scale_fill_brewer(palette="Greys") +
  theme(legend.key = element_blank(),  #removes the box around each legend item
        legend.position = c(0.5,0.99), 
        legend.direction = "horizontal",
        text = element_text(size = 20),
        legend.text = element_text(size = 20, face="bold"))






# For PAPER 
# PCA OF ENV/ENZ FACTORS
pca = prcomp(scale_enz_merge)
print(pca)
summary(pca)
str(pca)
print(pca$x)

# write.csv(pca$x, "env.pca.merge.csv") # export to Excel

plot(pca, type = "l") # plot of the variances (y) associated with the PCs (x)

ggbiplot(pca, obs.scale = 1, var.scale = 1, # plot the PCA
         groups = (grouping_info_merge[,1]), ellipse = TRUE, 
         circle = FALSE,var.axes=FALSE)

# Plot in ggplot
enz=data.frame(scores(pca),Elevation=as.factor(grouping_info_merge[,1]),Exposure=as.factor(grouping_info_merge[,2]),Transect=as.factor(grouping_info_merge[,3]), Site=as.factor(meta_table_merge[,21]),Ctot=meta_table_merge$C)
enz

p <- ggplot(enz) +
  geom_point(mapping = aes(x = PC1, y = PC2, colour = Elevation), size=3) +
  coord_fixed()
p + theme_bw()

# Add Elipse
ord<-ordiellipse(pca, as.factor(grouping_info_merge[,1]) ,display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# Generate ellipse points
df_ell <- data.frame()
for(g in levels(enz$Elevation)){
  if(g!="" && (g %in% names(ord))){
    
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(enz[enz$Elevation==g,],
                                                     veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                  ,Elevation=g))
  }
}

head(df_ell)


p <- ggplot(enz)
p1=p + geom_point(mapping = aes(x = PC1, y = PC2, colour=Ctot, shape = Elevation), size=7) +
  geom_path(data=df_ell, aes(x=PC1, y=PC2,fill=Elevation), size=0.75, linetype=2,color="darkgrey") +
  coord_fixed() +
  theme_bw() +
  scale_colour_gradient(low = "orange", high = "dark blue", space = "Lab", guide = "colourbar") +
  labs(x="PC1 (59.2%)", y="PC2 (17.0%)", colour = "total C") + #another way to set the labels, in this case, for the colour legend
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12)) +
  theme(legend.title=element_text(size=15) , legend.text=element_text(size=10)) +
  theme(legend.key = element_blank(),  #removes the box around each legend item
        legend.position = "bottom", #legend at the bottom
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.box.just = "centre") +
  theme(panel.grid.major = element_blank(), # for paper
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
p1 


# Plot env as factors in PCA
set.seed(123)
envf <- envfit(pca, scale_env_merge, na.rm = TRUE, perm = 999)
envf

# Create dataframe with only vectors with significant pvalues
A <- as.list(envf$vectors)
head(A)
pvals<-as.data.frame(A$pvals)
head(pvals)
arrows<-as.data.frame(A$arrows*sqrt(A$r))
head(arrows)
C<-cbind(arrows, pvals)
head(C)
Cred<-subset(C,pvals<0.051)
head(Cred)

# Link vector coordinates to their names
env.scrs <- cbind(Cred, Envir = rownames(Cred))
head (env.scrs)

# plot with Ctotal as factor # FOR PAPER
p1 + geom_segment(data = env.scrs,
                  aes(x = 0, xend = PC1, y = 0, yend =PC2),
                  arrow = arrow(length = unit(0.25, "cm")), colour = "black", size = 1.5) +
  geom_text(data = env.scrs, aes(x = PC1, y = PC2, label = Envir), size = 5, ,hjust=-0.2, vjust=0)


svg(filename="Figure_4.svg", 
    width=5, 
    height=6, 
    pointsize=10)
p1 + geom_segment(data = env.scrs,
                  aes(x = 0, xend = PC1, y = 0, yend =PC2),
                  arrow = arrow(length = unit(0.25, "cm")), colour = "black", size = 1.5) +
  geom_text(data = env.scrs, aes(x = PC1, y = PC2, label = Envir), size = 5, ,hjust=-0.2, vjust=0)
dev.off()





# generate species vectors on enzyme plot
vf <- envfit(pca, species_table_merge_top50_ecm, perm = 999) # 50 most common OTU
vf

# with best R values
A <- as.list(vf$vectors)
head(A)
r<-as.data.frame(A$r)
head(r)
arrows<-as.data.frame(A$arrows*sqrt(A$r))
head(arrows)
C<-cbind(arrows, r)
head(C)
Cred<-subset(C,A$pvals<0.051)
head(Cred)
print(Cred)

# Link vector coordinates to OTU names
spp.scrs <- cbind(Cred, Species = rownames(Cred))
head (spp.scrs)

# Plotting
p1 + geom_segment(data = spp.scrs,
                  aes(x = 0, xend = PC1, y = 0, yend =PC2),
                  arrow = arrow(length = unit(0.25, "cm")), colour = "black", size = 1.5) +
  geom_text(data = spp.scrs, aes(x = PC1, y = PC2,label = Species), size = 5, ,hjust=-0.2, vjust=0)






# FOR PAPER
# TEST EFFECT OF ELEVATION/TRANSECT ON PCA OF ENZYMES
pca = prcomp(scale_enz_merge)
summary(pca)

enz=data.frame(scores(pca),Elevation=as.factor(grouping_info_merge[,1]),Exposure=as.factor(grouping_info_merge[,2]),Transect=as.factor(grouping_info_merge[,3]), Site=as.factor(meta_table_merge[,21]))
enz

# (one-way anova)
transect.anova <- aov(PC3 ~ Transect, enz)
summary(transect.anova)
TukeyHSD(transect.anova) #pairwise

# Test for homogeneity of variances/normality
resid <- residuals(transect.anova.anova)
shapiro.test(resid) # normality test
leveneTest(PC3 ~ Transect, enz) # heteroscedasticity, homogeneity of variances

model <- lm(PC1 ~ Elevation, enz)
C.test(model) # Cochran test of homogeneity of variances


a<-lmer(PC1 ~ Elevation + (1|Transect), data=enz) # effect of elevation after after controlling for the specific variation of each transect
lmerTest::anova(a) # type III anova


# Plot PC axis per elevation/transect
svg(filename="Figure_S7.svg", 
    width=7, 
    height=6, 
    pointsize=10)

par (mfrow = c(2,3))

p1=plot(PC1 ~ Transect, enz, ylab="PC1", xlab="",cex.axis=1.25,cex.lab=1.5)
p2=plot(PC2 ~ Transect, enz, ylab="PC2", xlab="",cex.axis=1.25,cex.lab=1.5)
p3=plot(PC3 ~ Transect, enz, ylab="PC3", xlab="",cex.axis=1.25,cex.lab=1.5)

p1=plot(PC1 ~ Elevation, env, ylab="PC1", xlab="",cex.axis=1.25,cex.lab=1.5)
p2=plot(PC2 ~ Elevation, env, ylab="PC2", xlab="",cex.axis=1.25,cex.lab=1.5)
p3=plot(PC3 ~ Elevation, env, ylab="PC3", xlab="",cex.axis=1.25,cex.lab=1.5)

dev.off()




# FOR PAPER
# TEST EFFECT OF ELEVATION/TRANSECT ON PCA OF SOIL CHEMSTRY
pca = prcomp(scale_env_merge)
summary(pca)

env=data.frame(scores(pca),Elevation=as.factor(grouping_info_merge[,1]),Exposure=as.factor(grouping_info_merge[,2]),Transect=as.factor(grouping_info_merge[,3]), Site=as.factor(meta_table_merge[,21]))
env

a<-lmer(PC2 ~ Elevation + (1|Transect), data=env) # effect of elevation after controlling for the specific variation of each transect
lmerTest::summary(a)
lmerTest::anova(a) # type III anova

st = lmerTest::step(a) 
plot(st)
st


# (one-way anova)
transect.anova <- aov(PC3 ~ Transect, env)
summary(transect.anova)
TukeyHSD(transect.anova) #pairwise

resid <- residuals(transect.anova)
shapiro.test(resid) # normality test
leveneTest(PC3 ~ Transect, env) # heteroscedasticity, homogeneity of variances

model <- lm(PC3 ~ Transect, env)
C.test(model) # Cochran test of homogeneity of variances

# Plot PC axis per elevation/transect
svg(filename="Figure_S6.svg", 
    width=7, 
    height=7, 
    pointsize=10)

par (mfrow = c(2,3))

p1=plot(PC1 ~ Transect, env, ylab="PC1", xlab="",cex.axis=1.25,cex.lab=1.5)
p2=plot(PC2 ~ Transect, env, ylab="PC2", xlab="",cex.axis=1.25,cex.lab=1.5)
p3=plot(PC3 ~ Transect, env, ylab="PC3", xlab="",cex.axis=1.25,cex.lab=1.5)

p1=plot(PC1 ~ Elevation, env, ylab="PC1", xlab="",cex.axis=1.25,cex.lab=1.5)
p2=plot(PC2 ~ Elevation, env, ylab="PC2", xlab="",cex.axis=1.25,cex.lab=1.5)
p3=plot(PC3 ~ Elevation, env, ylab="PC3", xlab="",cex.axis=1.25,cex.lab=1.5)

dev.off()


























# EFFECT OF ELEVATION ON NMDS
elevation_anosim = anosim(species_table_merge_ecm, grouping_info_merge$X1,permutations = 999, distance = "raup")
elevation_anosim$signif
elevation_anosim$statistic

exposure_anosim = anosim(species_table_merge_ecm, grouping_info_merge$X2,permutations = 999, distance = "raup")
exposure_anosim$signif
exposure_anosim$statistic

site_anosim = anosim(species_table_merge_ecm, grouping_info_merge$X3,permutations = 999, distance = "raup")
site_anosim$signif
site_anosim$statistic



# INTERACTION ELEVATION WITH pH
d = vegdist(species_table_merge_ecm, "raup")
d

# permanova: depends on order of variables in formula! Not good
elevation_Adonis = adonis(d ~ ElevNum + pH, meta_table_merge) 
elevation_Adonis

elevation_Adonis = adonis(d ~ pH + ElevNum, meta_table_merge)
elevation_Adonis

Adonis.data=data.frame(elevation_Adonis$aov.tab)
Adonis.data
Adonis.sub = subset(Adonis.data, Pr..F. < 0.051)
Adonis.sub
Adonis.sort = Adonis.sub[order(-Adonis.sub$R2),] 
Adonis.sort
barplot(Adonis.sort$R2,names.arg=rownames(Adonis.sort), ylab="R2", col=terrain.colors(5),cex.axis=1.5, cex.names=1.5, cex.lab=1.5)


# GLM simple
sol=metaMDS(species_table_merge_ecm, distance = "raup", k = 2, trymax = 1000)
NMDS=data.frame(scores(sol),Elevation=as.factor(grouping_info_merge[,1]),Exposure=as.factor(grouping_info_merge[,2]),Transect=as.factor(grouping_info_merge[,3]),pH=as.factor(meta_table_merge$pH))
NMDS
colnames(NMDS)

glm<-glm(NMDS1 ~ as.numeric(pH) * Elevation, data=NMDS)
summary(glm)

a<-lm(NMDS1 ~ as.numeric(pH) * Elevation, data=NMDS) # order matter
anova(a)

# GLM with random effect
a<-lmer(NMDS1 ~ Elevation * as.numeric(pH) + (1|Transect), data=NMDS) # effect of elevation after controlling for the specific variation of each transect
lmerTest::summary(a)
lmerTest::anova(a) # type III anova

a<-lmer(NMDS1 ~ as.numeric(pH) + (1|Transect), data=NMDS)
lmerTest::summary(a)
lmerTest::anova(a) # type III anova







# MULTIVARIATE + SPECIES INTERRACTION
# greater power to detect patterns when analysing all species simultaneously than when looking for a pattern separately in each species

### Global
global=data.frame(species_table_merge_ecm, meta_table_merge, Elevation=as.factor(grouping_info_merge[,1]), Exposure=as.factor(grouping_info_merge[,2]), Transect=as.factor(grouping_info_merge[,3]))
colnames(global)
attach(global)

glm_global <- mvabund(global[,1:487])
glm_global

boxplot(global[,1:487],horizontal = TRUE,las=2, main="Abundance")
plot(glm_global~Elevation, cex.axis=0.8, cex=0.8)

manyglm_global <- manyglm(glm_global ~ pH * factor(Elevation) + factor(Transect), family="negative.binomial")
plot(manyglm_global) # little pattern suggests the assumption is plausible.
meanvar.plot(glm_global ~ Elevation)

aov.global=anova.manyglm(manyglm_global, show.time="all", nBoot=199) # no univariate here
aov.global
str(aov.global)




#### FOR PAPER
# Top50 OTU (all and ecm)

top50=data.frame(species_table_merge_top50, meta_table_merge, Elevation=as.factor(grouping_info_merge[,1]), Exposure=as.factor(grouping_info_merge[,2]), Transect=as.factor(grouping_info_merge[,3]))
colnames(top50)
attach(top50)

glm_top50 <- mvabund(top50[,1:50])
glm_top50

boxplot(top50[,1:50],horizontal = TRUE,las=2, main="Abundance")
plot(glm_top50~Elevation, cex.axis=0.8, cex=0.8)

manyglm_top50 <- manyglm(glm_top50 ~ pH * Elevation, family="negative.binomial")
plot(manyglm_top50) # little pattern suggests the assumption is plausible.

aov.top50=anova.manyglm(manyglm_top50, p.uni="adjusted",show.time="all") # univariate (per species) and multivariate (all species) test
aov.top50

aov.top50$table
write.csv(aov.top50$table, "top50_multivar3.csv")
manyglm_top50 = data.frame(t(aov.top50$uni.test), t(aov.top50$uni.p))
write.csv(manyglm_top50, "top50_univar3.csv")


# Top 50 ecm
top50=data.frame(species_table_merge_top50_ecm, meta_table_merge, Elevation=as.factor(grouping_info_merge[,1]), Exposure=as.factor(grouping_info_merge[,2]), Transect=as.factor(grouping_info_merge[,3]))
colnames(top50)
attach(top50)

glm_top50 <- mvabund(top50[,1:50])
glm_top50

boxplot(top50[,1:50],horizontal = TRUE,las=2, main="Abundance")
plot(glm_top50~Elevation, cex.axis=0.8, cex=0.8)

manyglm_top50 <- manyglm(glm_top50 ~ pH * Elevation, family="negative.binomial")
plot(manyglm_top50) # little pattern suggests the assumption is plausible.

aov.top50=anova.manyglm(manyglm_top50, p.uni="adjusted",show.time="all") # univariate (per species) and multivariate (all species) test
aov.top50
str(aov.top50)

aov.top50$table
write.csv(aov.top50$table, "top50_ecm_multivar3.csv")
manyglm_top50 = data.frame(t(aov.top50$uni.test), t(aov.top50$uni.p))
write.csv(manyglm_top50, "top50_ecm_univar3.csv")


# otu merged per genus
genus=data.frame(species_table_merge_genus, meta_table_merge, Elevation=as.factor(grouping_info_merge[,1]), Exposure=as.factor(grouping_info_merge[,2]), Transect=as.factor(grouping_info_merge[,3]))
colnames(genus)
attach(genus)

glm_genus <- mvabund(genus[,1:129]) # eliminate OTU in only one site or sample = the same!

boxplot(genus[,1:129],horizontal = TRUE,las=2, main="Abundance")
plot(glm_genus~Elevation, cex.axis=0.8, cex=0.8)

manyglm_genus <- manyglm(glm_genus ~ pH * Elevation, family="negative.binomial")
plot(manyglm_genus) # little pattern suggests the assumption is plausible.

aov.genus=anova.manyglm(manyglm_genus, p.uni="adjusted",show.time="all") # univariate (per species) and multivariate (all species) test
aov.genus$table

write.csv(aov.genus$table, "genus_multivar2.csv")
manyglm_genus = data.frame(t(aov.genus$uni.test), t(aov.genus$uni.p))
write.csv(manyglm_genus, "genus_univar2.csv")


# ecm genus
genus=data.frame(species_table_merge_genus_ecm, meta_table_merge, Elevation=as.factor(grouping_info_merge[,1]), Exposure=as.factor(grouping_info_merge[,2]), Transect=as.factor(grouping_info_merge[,3]))
colnames(genus)
attach(genus)

glm_genus <- mvabund(genus[,1:34]) # eliminate OTU in only one site or sample = the same!

boxplot(genus[,1:34],horizontal = TRUE,las=2, main="Abundance")
plot(glm_genus~Elevation, cex.axis=0.8, cex=0.8)

manyglm_genus <- manyglm(glm_genus ~ pH * Elevation, family="negative.binomial") # factor(Elevation)=more or less the same!
plot(manyglm_genus) # little pattern suggests the assumption is plausible.

aov.genus=anova.manyglm(manyglm_genus, p.uni="adjusted",show.time="all") # univariate (per species) and multivariate (all species) test
aov.genus$table

write.csv(aov.genus$table, "genus_ecm_multivar2.csv")
manyglm_genus = data.frame(t(aov.genus$uni.test), t(aov.genus$uni.p))
write.csv(manyglm_genus, "genus_ecm_univar2.csv")



# Each family separatedely: Mortierellaceae file does not work anymore, file corrupted???
Mortierellaceae.mergesite<-read.csv("/Users/Camille/Documents/NOTHOFAGUS/R/elevation/merge.site.Mortierellaceae(mc2).csv",row.names=1,check.names=FALSE) # re-import table
Mortierellaceae.mergesite
Mortierellaceae.mergesite<-Mortierellaceae.mergesite[rownames(meta_table_merge),]
Mortierellaceae.mergesite

Mortierellaceae=data.frame(Mortierellaceae.mergesite, meta_table_merge, Elevation=as.factor(grouping_info_merge[,1]), Exposure=as.factor(grouping_info_merge[,2]), Transect=as.factor(grouping_info_merge[,3]))
colnames(Mortierellaceae)
attach(Mortierellaceae)

glm_Mortierellaceae <- mvabund(Mortierellaceae[,1:50]) # eliminate OTU prsent only in one site

boxplot(Mortierellaceae[,1:50],horizontal = TRUE,las=2, main="Abundance")
plot(glm_Mortierellaceae~Elevation, cex.axis=0.8, cex=0.8)

manyglm_Mortierellaceae <- manyglm(glm_Mortierellaceae ~ pH * Elevation, family="negative.binomial")
plot(manyglm_Mortierellaceae) # little pattern suggests the assumption is plausible.

aov.Mortierellaceae=anova.manyglm(manyglm_Mortierellaceae, p.uni="adjusted",show.time="all") # univariate (per species) and multivariate (all species) test
aov.Mortierellaceae$table

write.csv(aov.Mortierellaceae$table, "Mortierellaceae_multivar2.csv")
manyglm_Mortierellaceae = data.frame(t(aov.Mortierellaceae$uni.test), t(aov.Mortierellaceae$uni.p))
write.csv(manyglm_Mortierellaceae, "Mortierellaceae_univar2.csv")



Cortinariaceae.mergesite<-read.csv("/Users/Camille/Documents/NOTHOFAGUS/R/elevation/merge.site.Cortinariaceae(mc2).csv",row.names=1,check.names=FALSE) # re-import table
Cortinariaceae.mergesite
Cortinariaceae.mergesite<-Cortinariaceae.mergesite[rownames(meta_table_merge),]
Cortinariaceae.mergesite

Cortinariaceae=data.frame(Cortinariaceae.mergesite, meta_table_merge, Elevation=as.factor(grouping_info_merge[,1]), Exposure=as.factor(grouping_info_merge[,2]), Transect=as.factor(grouping_info_merge[,3]))
colnames(Cortinariaceae)
attach(Cortinariaceae)

glm_Cortinariaceae <- mvabund(Cortinariaceae[,1:241]) # eliminate OTU prsent only in one sample

boxplot(Cortinariaceae[,1:241],horizontal = TRUE,las=2, main="Abundance")
plot(glm_Cortinariaceae~Elevation, cex.axis=0.8, cex=0.8)

manyglm_Cortinariaceae <- manyglm(glm_Cortinariaceae ~ pH + Elevation + Transect, family="negative.binomial")
plot(manyglm_Cortinariaceae) # little pattern suggests the assumption is plausible.

aov.Cortinariaceae=anova.manyglm(manyglm_Cortinariaceae, p.uni="adjusted",show.time="all") # univariate (per species) and multivariate (all species) test
aov.Cortinariaceae$table

write.csv(aov.Cortinariaceae$table, "Cortinariaceae_multivar3.csv")
manyglm_Cortinariaceae = data.frame(t(aov.Cortinariaceae$uni.test), t(aov.Cortinariaceae$uni.p))
write.csv(manyglm_Cortinariaceae, "Cortinariaceae_univar3.csv")







# species indicator (ELEVATION and pH)
indval = multipatt(species_table_merge_family, grouping_info_merge$X1, control = how(nperm=999))
summary(indval)
summary(indval, indvalcomp=TRUE)
str(indval)
write.csv(data.frame(indval$str,indval$sign), "multipatt_family_elevation.csv")

prefsign = signassoc(species_table_merge_genus_ecm, cluster=meta_table_merge$pH.cat,  alternative = "two.sided", control = how(nperm=199))
head(prefsign)
print(prefsign)
write.csv(prefsign, "signassoc_genus_ecm_pH.csv")






### Species accumulation curves

otu_table(elevation.mc1..rel)
data=t(otu_table(elevation.mc1.rel))
data
accum=specaccum(data, method = "random", permutations=1000) # default method = "exact"
plot(accum)
plot(accum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue") #males a prettier plot

# per elevation belt
# all fungi
belt=subset_samples(elevation.mc1.rel, sample_data(elevation.mc1.rel)$Elevation == "L")
belt
data=t(otu_table(belt))
accum=specaccum(data, method = "random", permutations=1000) # default method = "exact"
plot(accum)

belt=subset_samples(elevation.mc1.rel, sample_data(elevation.mc1.rel)$Elevation == "M")
belt
data=t(otu_table(belt))
accumM=specaccum(data, method = "random", permutations=1000) # default method = "exact"
plot(accumM, add=TRUE)

belt=subset_samples(elevation.mc1, sample_data(elevation.mc1)$Elevation == "T")
belt
data=t(otu_table(belt))
accumT=specaccum(data, method = "random", permutations=1000) # default method = "exact"
plot(accumT, add=TRUE)

# ecm fungi
belt=subset_samples(elevation.mc1.ecm.rel, sample_data(elevation.mc1.ecm.rel)$Elevation == "L")
belt
data=t(otu_table(belt))
accumL=specaccum(data, method = "random", permutations=1000) # default method = "exact"
plot(accumL)

belt=subset_samples(elevation.mc1.ecm.rel, sample_data(elevation.mc1.ecm.rel)$Elevation == "M")
belt
data=t(otu_table(belt))
accumM=specaccum(data, method = "random", permutations=1000) # default method = "exact"
plot(accumM, add=TRUE)

belt=subset_samples(elevation.mc1.ecm.rel, sample_data(elevation.mc1.ecm.rel)$Elevation == "T")
belt
data=t(otu_table(belt))
accumT=specaccum(data, method = "random", permutations=1000) # default method = "exact"
plot(accumT, add=TRUE)




# idem per site
site=subset_samples(elevation.mc1, sample_data(elevation.mc1)$SiteName == "LV1")
site
data=t(otu_table(site))
accum=specaccum(data, method = "random", permutations=1000) # default method = "exact"
plot(accum)

site=subset_samples(elevation.mc1, sample_data(elevation.mc1)$SiteName == "LM1")
site
data=t(otu_table(site))
accum=specaccum(data, method = "random", permutations=1000) # default method = "exact"
plot(accum, add=TRUE)



sample_data(elevation.mc1)

plot(specaccum(data,method="random"))

plot(random)
plot(random, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", 
     main = "Default: Prettier CI")

raref<-specaccum(data, method="rarefaction", permutations=1000,xvar=180)
p=plot(raref$sites, raref$richness,
       xlab="Number of Sites",
       ylab="Number of OTUs",
       main="")


# observed abundances: rarefy!
otu_table(elevation.mc1)
s=sample_sums(otu_table(elevation.mc1))
s
raremax <- min(rowSums(s)) # Number of sequences per samles ; min=100001
raremax

data=t(otu_table(elevation.mc1.ecm))
raremax <- min(rowSums(data)) # Number of sequences per samles ; min=100001
raremax

data

S <- specnumber(data) #total number of species in each sample (row of data)
S
raremax <- min(rowSums(data)) # Number of sequences per samles ; min=100001
raremax
Srare <- rarefy(data, raremax) # rarefy, with raremax as input
Srare

#Plot rarefaction results
par(mfrow = c(1,2))
plot(S, Srare, xlab = "Observed No. of Species", 
     ylab = "Rarefied No. of Species",
     main = " plot(rarefy(BCI, raremax))")
abline(0, 1)
rarecurve(data, step = 20, 
          sample = raremax, 
          col = "blue", 
          cex = 0.6,
          main = "rarecurve()")

rarecurve(data, step = 20, 
          sample = raremax, 
          col = "blue", 
          cex = 0.6,
          main = "rarecurve()")





sp1 <- specaccum(otu_table(elevation.mc1.rel))
sp2 <- specaccum(otu_table(elevation.mc1.rel), "random")
summary(sp2)
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sp2, col="yellow", add=TRUE, pch="+")

## Fit Lomolino model to the exact accumulation
mod1 <- fitspecaccum(sp1, "lomolino")
coef(mod1)
fitted(mod1)
plot(sp1)

## Add Lomolino model using argument 'add'
plot(mod1, add = TRUE, col=2, lwd=2)

## Fit Arrhenius models to all random accumulations
mods <- fitspecaccum(sp2, "arrh")
plot(mods, col="hotpink")
boxplot(sp2, col = "yellow", border = "blue", lty=1, cex=0.3, add= TRUE)

## Use nls() methods to the list of models
sapply(mods$models, AIC)






# Rarefaction curves were generated using the rarefy function in vegan by repeatedly sub-sampling the OTU table at steps of 1000
# A mean number of OTUs and the standard deviationfor each sample per step was calculated using the R package plyr and plotted using the R package ggplot2

S <- specnumber(t(otu_table(elevation.mc1))) # observed number of species
S
sp.abund <- rowSums(t(otu_table(elevation.mc1)))  #gives the number of sequences found in each samples
sp.abund
raremax <- min(rowSums(t(otu_table(elevation.mc1))))  # smallest number of observations (reads) per sample to extrapolate the expected number if all other samples only had that number of observations
raremax

Srare <- rarefy(t(otu_table(elevation.mc1)), raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(t(otu_table(elevation.mc1)), step = 20, sample = raremax, col = "blue", cex = 0.6,xlab = "Sample Size", ylab = "Species",
          label = TRUE)






### ordir2step: TO DO

# raup crick distance matrix: for binary data
binary_table_ecm=decostand(species_table_merge_ecm,method="pa") # use all OTU (mc1, with OTU present in one sample)
# do not merge per site

raup.ecm<-raupcrick(binary_table_ecm, null="r1", nsimul=999, chase=FALSE)

# capscale for binary data (pcoa)
# dbrda for abundance data
ecm.full<-capscale(raup.ecm~soil.PC1+soil.PC2+soil.PC3,soil.PCA_table_merge) # add all PC axis
ecm.min<-capscale(raup.ecm~1,soil.PCA_table_merge)

ecm.raup<-ordiR2step(ecm.min, scope = formula (ecm.full))

plot(ecm.raup)
summary(ecm.raup)
anova(ecm.raup, by="term")
anova(ecm.raup, by="axis")


# soil variable independently:
# check for colinearity of soil variables
# if they are not you can put them in the model as they are

ecm.soil<-capscale(raup.ecm~SoilMoist+pH+C+N+NO3+NH4+Nav+P,env_all_table_merge)
vif.cca (ecm.soil) # values over 10 indicate redundant constraints, otherwise can be used independently

ecm.soil<-capscale(raup.ecm~SoilMoist+pH+C+N+NO3+NH4+Nav+P,env_all_table_merge) # instead of ecm.full
ecm.min<-capscale(raup.ecm~1,env_all_table_merge)

ecm.raup<-ordiR2step(ecm.min, scope = formula (ecm.soil))

plot(ecm.raup)
summary(ecm.raup)
anova(ecm.raup, by="term")
anova(ecm.raup, by="axis")


