rm(list=ls())

setwd("~/Library/CloudStorage/OneDrive-ImperialCollegeLondon/MSc EEC/Mini Project")

##Load packages

library(tidyverse)
library(ape)
library(caper)
require(MASS)

#Read in coursework_trait_data for EcoEDGE calculations

coursework_trait_data <- read.csv("coursework_trait_data.csv")

#Remove extinct species

trait_data <- coursework_trait_data%>% filter(!(Redlist_cat == "EX" | Redlist_cat == "EW"))

##Calculate EDGE scores##

##Calculate ED Scores

#Load in and plot the phylogenetic tree

bird_tree <- read.tree("all_birds.tre")
plot(bird_tree)

#Transform our tree into a matrix of distances from each tip to tip. (Takes a while to run)
bird_matrix <- clade.matrix(bird_tree)

#Calculate ED scores for each species. 
ED <- ed.calc(bird_matrix)$spp
head(ED)

#Log transform and normalise scores.

ED$EDlog <- log(1+ED$ED)

#normalise scores so they're scaled between 0 and 1
ED$EDn <- (ED$EDlog - min(ED$EDlog)) / (max(ED$EDlog) - min(ED$EDlog))
head(ED)

# Find the highest ED score
ED[ED$EDn == max(ED$EDn),] #Steatornis_caripensis - oilbird

##Calculate GE Scores

# Create an empty column to store GE scores.
trait_data$GE <- NA

# Create a vector to increase with each new ranking, starting at 0 for least concern.
i <- 0

# Create a list to loop through in the order of GE scores.
redlist_cats <- c("LC", "NT", "DD", "VU", "EN", "CR")

# Loop through each different category in the redlist categories.
for (category in redlist_cats){
  
  # Add the GE score for that category.
  trait_data[trait_data$Redlist_cat == category, "GE"] <- i
  
  # Because DD comes after NT, and both are scored as 1, don't want to change i if the category is NT.
  if (category != "NT"){
    i = i + 1
  }
}

unique(trait_data$Redlist_cat)
unique(trait_data$GE)

#Merge GE scores with ED scores in one dataframe.

# Join the last two columns of UK_Jetz to ED scores. 
Bird_EDGE <- left_join(trait_data, ED,  by = c("Jetz_Name" = "species"))
head(Bird_EDGE)

#Calculate EDGE scores - 𝐸𝐷𝐺𝐸=𝑙𝑛(1 +𝐸𝐷)+𝐺𝐸×𝑙𝑛(2)

Bird_EDGE$EDGE <- Bird_EDGE$EDlog + Bird_EDGE$GE * log(2)
head(Bird_EDGE)

#Find the highest EDGE score.
Bird_EDGE[Bird_EDGE$EDGE == max(Bird_EDGE$EDGE),] #Thaumatibis gigantea - Giant Ibis

##Calculate FD Scores##

#change row names to species names and remove all the columns except traits. Then normalise trait data so that body_mass and beak have the same scale (the same variance).

#Make a copy of bird data
Bird_traits <- Bird_EDGE

#Change row names and keep just trait data.
rownames(Bird_traits) <- Bird_traits$Jetz_Name
Bird_traits <- Bird_traits[,6:7]

#Make each column have the same scale.
Bird_traits <- scale(Bird_traits, scale=T)
head(Bird_traits)

summary(Bird_traits)

#Remove any NAs
Bird_traits <- na.omit(Bird_traits)

#Create a matrix
traits_matrix <- as.matrix(Bird_traits)

#Convert traits into 'distance' in trait space. Species with similar traits will have smaller ‘distances’.
distance_matrix <- dist(traits_matrix)

#Create a new tree using the neighbour-joining method

#Create the tree (takes a long time to run ~ 30 mins)
trait_tree <- nj(distance_matrix)

#Create a matrix of distance from tip to tip (takes a very long time to run ~ 30 mins)
tree_matrix <- clade.matrix(trait_tree)

##Calculate FD scores
FD <- ed.calc(tree_matrix)$spp

#Change the name to FD
colnames(FD)[2] <- "FD"
head(FD)

#Log and normalise the data
FD$FDlog <- log(1+FD$FD)
FD$FDn <- (FD$FDlog - min(FD$FDlog)) / (max(FD$FDlog) - min(FD$FDlog))

# Find the highest FD score
FD[FD$FDn == max(FD$FDn),] #Struthio_camelus

#Combine GE scores. Use the same formula as before:𝐹𝑈𝐷𝐺𝐸=𝑙𝑛(1+𝐹𝐷)+𝐺𝐸×𝑙𝑛(2)

# Join FD and GE scores
Bird_FUDGE <- left_join(trait_data, FD, by = c("Jetz_Name" = "species"))

# Calculate FUDGE scores
Bird_FUDGE$FUDGE <- Bird_FUDGE$FDlog + Bird_FUDGE$GE * log(2)
head(Bird_FUDGE)

#Find the highest FUDGE score
Bird_FUDGE[Bird_FUDGE$FUDGE == max(Bird_FUDGE$FUDGE),] #Struthio camelus - Common Ostrich 

#Omit NAs in both datasets
Bird_FUDGE <- na.omit(Bird_FUDGE)
Bird_EDGE <- na.omit(Bird_EDGE)

##Calculating EcoEDGE Scores##

#We give ED and FD scores equal weighting: 𝐸𝑐𝑜𝐸𝐷𝐺𝐸=(0.5×𝐸𝐷𝑛+0.5×𝐹𝐷𝑛)+𝐺𝐸×𝑙𝑛(2)

#Merge FD and ED scores.
Bird_EcoEDGE <- left_join(Bird_EDGE, Bird_FUDGE)

#Calculate EcoEDGE scores
Bird_EcoEDGE$EcoEDGE <- (0.5*Bird_EcoEDGE$EDn + 0.5*Bird_EcoEDGE$FDn) + Bird_EcoEDGE$GE*log(2)
head(Bird_EcoEDGE)

#Get the highest scoring species
Bird_EcoEDGE[Bird_EcoEDGE$EcoEDGE == max(Bird_EcoEDGE$EcoEDGE),] #Thaumatibis gigantea - Giant Ibis

#See the spread
hist(Bird_EcoEDGE$EcoEDGE, breaks = 20)


###Data Analysis###


#Read in HWI dataset

Dataset_HWI <- read.csv("Dataset_HWI.csv")

#Join the EcoEDGE scores to the new dataset
All_Data <- left_join(Bird_EcoEDGE, Dataset_HWI, by = c("Jetz_Name" = "Tree.name"))

summary(All_Data)
str(All_Data)

#Drop unwanted columns

colnames(All_Data)

All_Data <- All_Data%>%dplyr::select(-c(19,20,21,22,24,25,27,29,30,31,33,39,40))

summary(All_Data)

All_Data <- All_Data%>%dplyr::select(-c(20,21,24,25,26,27))

#Remove NAs

All_Data <- na.omit(All_Data)

summary(All_Data)

write.csv(All_Data, "All_Data.csv")

#Preliminary analysis - trying to see what predicts EcoEDGE/other conservation scores

hist(All_Data$HWI)
hist(log(All_Data$HWI))
hist(scale(All_Data$HWI))
hist(All_Data$EcoEDGE)
hist(log(All_Data$EcoEDGE))

plot(EcoEDGE ~ HWI, All_Data)

boxplot(EcoEDGE ~ Diet, All_Data)

dev.off()

###Model selection###

require(usdm)
require(psych)
library(lmtest)
require(lme4)
require(sjPlot)
require(lmerTest) #This gives p values for LMMs

#Test backwards from the maximum optimal model

model1 <- lmer(log(EcoEDGE) ~ scale(HWI)*Diet + (1|Jetz_order/Jetz_family), data=All_Data)
summary(model1)

plot(model1)
anova(model1)

#Information Criteria (AIC)

Test_model <- step(model1, direction = "backward", scope = list(lower=~1, upper=~scale(HWI)*Diet + (1|Jetz_order/Jetz_family)), data=All_Data) # scope here is indicating the lower (null model) and upper (maximal model)
Test_model #Take out the interaction

model2 <- lmer(log(EcoEDGE) ~ scale(HWI) + Diet + (1 | Jetz_order/Jetz_family), data = All_Data)
summary(model2)

plot(model2)
hist(resid(model2))

plot_model(model2, show.values = TRUE, show.intercept = TRUE) #Scavengers have a mean EcoEDGE of 0.55 higher than the intercept

#Pairwise ad hoc test to look at differences in diets
library(emmeans)

emmeans(model2, pairwise ~ Diet)

#Calculate variance explained by random effects
0.07+0.1+0.76
0.07/0.93
0.1/0.93
0.17/0.93


##Barplot of mean EcoEDGE scores of diets 

#Standard deviation of the mean as error bar

meanEcoEDGE <- tapply(All_Data$EcoEDGE, list(All_Data$Diet), mean)
se <- function(x) sd(x)/sqrt(length(x))
seEcoEDGE <- tapply(All_Data$EcoEDGE, list(All_Data$Diet), se)

Diet <- c("fruit", "invert", "nectar", "omni", "plants", "scav",  "seeds", "vertebrates")

Diet_plot <- data_frame(Diet, meanEcoEDGE, seEcoEDGE)

p <- ggplot(Diet_plot, aes(x=Diet, y=meanEcoEDGE)) +  
  geom_bar(stat="identity", fill="#70a858") + 
  theme_minimal() +
  xlab("Diet") + ylab("EcoEDGE Score") +
  theme(axis.text=element_text(size=13, color="black")) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.title.x = element_text(size = 15)) +
  geom_errorbar(aes(ymin=meanEcoEDGE-seEcoEDGE, ymax=meanEcoEDGE+seEcoEDGE), width=0.2, size=0.3, color="black")

p







