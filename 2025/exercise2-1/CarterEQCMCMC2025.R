
# This R script was written by Pat Carter in June 2025 for the Evolutionary Quantitative Genetics Workshop held in June 2025.
# It is a modification of the R script written by Pat Carter in April 2024 for the WSU Biol 521 Quantitative Genetics class.
# It is a modification of the R script written by Pat Carter n July 2022 for the Evolutionary Quantitative Genetics Workshop.
# It is a modification of the R script written by Pat Carter in April 2022 for the WSU Biol 521 Quantitative Genetics class.
# This in turn is a modification of the R script written by Pat Carter in July 2021 for the Friday Harbor Evolutionary Quantitative Genetics class, 
# and is a modification of R Scripts written in June 2019 and June 2018
# for the Friday Harbor Evolutionary Quantitative Genetics class,  April 2017 for WSU Biology 521, and in August 2016 for the 
# NIMBioS Evolutionary Quantitative Genetics course.

# The purpose of this script is to highlight code needed to run quantitative genetic analyses in R using MCMCglmm
# The data to be analyzed are 3 phenotypic traits from a toy data set: Ptype1, Ptype2 and Ptype3 generated in Biol521-2018-1.R.
# Data were generated for only one generation but pedigree information is known for parents and grandparents
# Individuals were measured in one of two batches

#IMPORTANT NOTE TO STUDENTS: You will need to change the path for opening and saving data files to whatever is appropriate to your computer


#######################################################
#Visualize and examine the data 

#Clear memory if needed
rm(list=ls()) 

#load graphics library ggplot2
library(ggplot2)

#open the data file and read it in as an R data file
Toy4 <- read.table ("C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\Toy4SimZ.dat",header = T)



#look at the header of the data file
head(Toy4)
#summarize the data file
summary (Toy4)

# look at distributions of the 3 Ptypes
Toy41.hist = hist(Toy4$Ptype1, breaks = 20)
Toy42.hist = hist(Toy4$Ptype2, breaks = 20)
Toy43.hist = hist(Toy4$Ptype3, breaks = 20)

#plot the phenotypes with each other
qplot(x = Ptype1, y = Ptype2, data = Toy4, geom = "point", color = factor(Batch) )
qplot(x = Ptype1, y = Ptype3, data = Toy4, geom = "point", color = factor(Batch) )
qplot(x = Ptype2, y = Ptype3, data = Toy4, geom = "point", color = factor(Batch) )

plot(density(x=Toy4$Ptype1),main="Density Estimate")

#Look at covariance structure and correlations
#Choose variables for which we want covariances
TempVar <- subset(Toy4, select = c(Ptype1,Ptype2,Ptype3))
summary(TempVar)
#Get covariance matrix  
CV = cov(TempVar)
#Get correlation matrix
CR = cor(TempVar)
#Look at eigen structure just for fun
EigCV <- eigen(CV)
EigCV


#######################################################
# create pedigree file 

# Important information about the structure and construction of pedigree files to be used by MCMCglmm
# pedigree files contain identification information for measured individuals and their mothers (Dam) and fathers (Sire) and other ancestors
# this file may contain many generations of individuals for whom phenotypes were never measured
# note that offspring always must have larger id numbers than their parents, parents must have larger id numbers than grandparents, etc

# make sure missing values for Sire and Dam identification numbers are coded as NA
# make sure that dam and sire id numbers are coded as factors
# the identification variable for the individuals with phenotypic information in the data set MUST be called animal and it should be a factor variable
# the pedigree file should be sorted by the variable animal

# individuals may appear in the pedigree as animal and they may also appear as a Sire or Dam

# frequently you must create the pedigree file from the data file, as we have to do here
# the data file contains one record for each measured individual
# each record contains data on all 3 phenotypes as well as identification information on parents and grandparents

rm(list=ls())

#the nadiv and plyr libraries are needed to make the pedigree

library(plyr)
library(nadiv)

#open the data file
Toy4 <- read.table ("C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\Toy4SimZ.dat",header = T)

#Create Offspring pedigree
Toy4Pedid<-subset(Toy4,select=c(Id,Sireid,Damid))
Toy4Pedid<-rename(Toy4Pedid,c(Id="animal"))

#Create Pedigree of dams
#save one record per dam
Toy4PedDam<- aggregate(Toy4, list(Toy4$Damid), FUN=head, 1)
#choose only identification information for pedigree
Toy4PedDam<-subset(Toy4PedDam,select=c(Damid,MGsireid,MGdamid))
#rename variables
Toy4PedDam<-rename(Toy4PedDam,c(Damid="animal",MGsireid="Sireid",MGdamid="Damid"))

#Create Pedigree of sires
#save one record per sire
Toy4PedSire<- aggregate(Toy4, list(Toy4$Sireid), FUN=head, 1)
#choose only identification information for pedigree
Toy4PedSire<-subset(Toy4PedSire,select=c(Sireid,PGsireid,PGdamid))
#rename variables
Toy4PedSire<-rename(Toy4PedSire,c(Sireid="animal",PGsireid="Sireid",PGdamid="Damid"))

#Combine all 3 files using rbind to make one full pedigree file
Toy4ped<-rbind(Toy4Pedid,Toy4PedDam,Toy4PedSire)
#sort by animal number
Toy4ped<-Toy4ped[order(Toy4ped$animal),]


#replace "." with NA for missing values for Sire and Dam id numbers
Toy4ped$Sireid[Toy4ped$Sireid=="."]<-NA
Toy4ped$Damid[Toy4ped$Damid=="."]<-NA

#identify animal, Dam and Sire as factors
Toy4ped$Damid <- as.factor(Toy4ped$Damid)
Toy4ped$Sireid <- as.factor(Toy4ped$Sireid)
Toy4ped$animal <- as.factor(Toy4ped$animal)

# this command from the nadiv library completes a pedigree with missing information for some sires and dams 
# by adding the generation in which all Dams and Sires were unknown; you will need to do this when making your pedigree file:
# Toy4ped<-insertPed(Toy4ped, founders=NULL) try pedtools?

Toy4pedT<- prepPed(Toy4ped)

#save pedigree file as an R data file
save(Toy4ped,file="C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\Toy4ped.rat")



###############################################
# modify data file

# make sure missing values for Sire and Dam id numbers are coded as NA
# make sure that animal, dam and sire id numbers are coded as factors
# the identification variable for individuals with phenotypic data MUST called animal 
# drop cases of individuals without phenotypic data 
# file should be sorted by the variable animal

Toy4 <- read.table ("C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\Toy4SimZ.dat",header = T)  #dont forget to change the path#

#replace "." with NA for missing values for Sire and Dam id numbers
Toy4$Sireid[Toy4$Sireid=="."]<-NA
Toy4$Damid[Toy4$Damid=="."]<-NA

#identify animal, Dam and Sire as factors
Toy4$Damid <- as.factor(Toy4$Damid)
Toy4$Sireid <- as.factor(Toy4$Sireid)


#rename the variable id as animal
Toy4<-rename(Toy4,c(Id="animal"))

Toy4$animal <- as.factor(Toy4$animal)

#extract only the variables we need
Toy4<-subset(Toy4, select = c(animal,Sireid,Damid,Batch, Ptype1, Ptype2, Ptype3))

#summarize data file (Note: before this point you should already have thoroughly graphed and examined your data)
summary(Toy4)
#save data file as an R data file
save(Toy4,file="C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\Toy4.rat")

######################################################################################################################################################
######################################################################################################################################################

# Model 101: Trait = Ptype1, Batch fixed, Additive effects only random effect

#Clear memory if needed
rm(list=ls()) 

#needed library for running generalized linear mixed models
library(MCMCglmm)
#Using remotes install Matthew Wolak's package wolakR from github which provides better graphics for posterior distributions
remotes::install_github("matthewwolak/wolakR")
#Load wolakR library
library(wolakR)

#load the data file and the pedigree file
load("C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\Toy4.rat")
load("C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\Toy4ped.rat")

# set up parameters for the priors
# Fixed effects are automatically set by MCMCglmm to follow a normal distribution and do not need to be specified
# G is for specified random effects (additive, maternal, etc). 
# Here we set for weak priors that will be used in an inverse gamma distribution automatically set by MCMCglmm
# R is residual effects for each specified random effect and follows same rules as G
prior101 <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=.002)))

# prunePed ensures that ancestors of focal animals are retained in the pedigree
Toy4PrunePed<-prunePed(Toy4ped,keep=1001:4500, make.base=TRUE)

# model statement
model101 <- MCMCglmm(Ptype1 ~ 1 + Batch,                #intercept (the 1) and Batch are the fixed effect
                     random = ~animal,                  #additive effects (animal) the only random effect
                     family = "gaussian",               #phenotype has gaussian distribution
                     prior = prior101,                  #call the priors parameters defined above
                     data = Toy4,                       #call the data file
                     nitt = 10000,                      #number of MCMC iterations 
                     burnin = 2000,                     #number of iterations for burnin
                     thin = 500,                        #sampling interval
                     pedigree = Toy4PrunePed)           #call the pedigree to get the inverse of the NRM

#save model output as an R object so we can access it later; this is very important when the nitt is very high
save(model101, file = "C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\model101.obj") #dont forget to change the path#

#load model output file
load(file = "C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\model101.obj") #dont forget to change the path#

############################################################################Examine results of model101

#########################Diagnostics

#plot trace and density of fixed effects; should be no trend in trace
plot(model101$Sol)
#plot trace and density of random (additive and residual (=environmental)) variances; should be no trend in trace
plot(model101$VCV)

#examine autocorrelation of fixed effects
autocorr.diag(model101$Sol)
#examine autocorrelation of random (additive and residual) variances
autocorr.diag(model101$VCV)

#check effective population size for fixed effects; should be gt 1000
effectiveSize(model101$Sol)
#check effective population size for random effects (additve and residual variances); should be gt 1000
effectiveSize(model101$VCV)

#test of convergence, p should be greater than 0.05 for good convergence
heidel.diag(model101$VCV)

#estimates of additive and residual variances
posterior.mode(model101$VCV)

#summary of model; make sure to check DIC score (smaller is better)
summary(model101)

#########################Additive Genetic effects and Heritability (i.e., variance ratio of additive genetic effects)

#estimate posterior distribution of additive genetic effects and get the 95% CI
posterior.mode(model101$VCV[,"animal"])
HPDinterval(model101$VCV[,"animal"])

#plot posterior distribution of additive genetic effects
plot(model101$VCV[,"animal"])
#plot posterior distribution of heritability with histogram, 95% CI, and mean
postPlot(model101$VCV[,"animal"])
#estimate 5, 20, and 90% quantiles of posterior distribution
quantile(model101$VCV[,"animal"],na.rm = T,probs = c(.05, .20,.90))

#estimate posterior distribution of the heritability (animal variance divided by animal + residual variances)
herit <- model101$VCV[, "animal"]/(model101$VCV[, "animal"] + model101$VCV[, "units"])

#effective sample size for heritability should be gt 1000
effectiveSize(herit)

# get the mean and mode from the posterior distribution of heritability
mean(herit)
posterior.mode(herit)

# get confidence interval for heritability
HPDinterval(herit)

#plot the trace of heritability, should not be any pattern
plot(herit)
#plot posterior distribution of heritability with histogram, 95% CI, and mean
postPlot(herit)
#estimate 5, 20, and 90% quantiles of posterior distribution
quantile(herit,na.rm = T,probs = c(.05, .20,.90))



######################################################################################################################################################
######################################################################################################################################################

# Model 1101: Trait = Ptype1, Batch fixed, Additive effects and maternal environmental effects are random

#needed library for running genearlized linear mixed models
library(MCMCglmm)
#Using remotes install Matthew Wolak's package wolakR from github which provides better graphics for posterior distributions
remotes::install_github("matthewwolak/wolakR")
#Load wolakR library
library(wolakR)

#Clear memory if needed
rm(list=ls()) 

load("C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\Toy4.rat")
load("C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\Toy4ped.rat")

# set up parameters for the priors
# Fixed effects are automatically set by MCMCglmm to follow a normal distribution and do not need to be specified
# G is for specified random effects (additive, maternal, etc). 
# Here we set for weak priors that will be used in an inverse gamma disribution automatically set by MCMCglmm
# R is residual effects for each specified random effect and follows same rules as G
prior1101 <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002), G2 = list(V = 1, nu = 0.002)))

# prunePed ensures that ancestors of focal animals are retained in the pedigree
Toy4PrunePed<-prunePed(Toy4ped,keep=1001:4500, make.base=TRUE)

# model statement
model1101 <- MCMCglmm(Ptype1 ~ 1 + Batch,               #intercept (the 1) and Batch are the fixed effect
                     random = ~animal + Damid,          #additive (animal) and maternal (Damid) are random effects
                     family = "gaussian",               #phenotype has gaussian distribution
                     prior = prior1101,                 #call the priors parameters defined above
                     data = Toy4,                       #call the data file
                     nitt = 1000000,                    #number of MCMC iterations 
                     burnin = 10000,                    #number of iterations for burnin
                     thin = 1000,                       #sampling interval
                     pedigree = Toy4PrunePed)                #call the pedigree to get the inverse of the NRM
                     

#save model output as an R object so we can access it later
save(model1101, file = "C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\model1101.obj") #dont forget to change the path#

#load model output file
load(file = "C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\model1101.obj") #dont forget to change the path#


############################################################################Examine results of model1101

#########################Diagnostics

#plot trace and density of fixed effects; should be no trend in trace
plot(model1101$Sol)
#plot trace and density of random (additive, maternal and residual (=environmental)) variances; should be no trend in trace
plot(model1101$VCV)

#examine autocorrelation of fixed effects
autocorr.diag(model1101$Sol)
#examine autocorrelation of random (additive, maternal and residual) variances
autocorr.diag(model1101$VCV)

#check effective population size for fixed effects; should be gt 1000
effectiveSize(model1101$Sol)
#check effective population size for random effects (additve, maternal and residual variances); should be gt 1000
effectiveSize(model1101$VCV)

#test of convergence, p should be greater than 0.05 for good convergence
heidel.diag(model1101$VCV)

#estimates of additive and residual variances
posterior.mode(model1101$VCV)

#summary of model; make sure to check DIC score (smaller is better)
summary(model1101)

#########################Additive Genetic Effects and Heritability (i.e., variance ratio of additive genetic effects)

#estimate posterior distribution of additive genetic effects and get the 95% CI
posterior.mode(model1101$VCV[,"animal"])
HPDinterval(model1101$VCV[,"animal"])

#plot posterior distribution of additive genetic effects
plot(model1101$VCV[,"animal"])
#plot posterior distribution of heritability with histogram, 95% CI, and mean
postPlot(model1101$VCV[,"animal"])
#estimate 5, 20, and 90% quantiles of posterior distribution
quantile(model1101$VCV[,"animal"],na.rm = T,probs = c(.05, .20,.90))

#estimate posterior distribution of the heritability (animal vairance divided by animal + maternal +  residual variances)
herit <- model1101$VCV[, "animal"]/(model1101$VCV[, "animal"] + model1101$VCV[, "Damid"] + model1101$VCV[, "units"])

#effective sample size for heritability should be gt 1000
effectiveSize(herit)

# get the mean from the posterior distribution of heritability
mean(herit)

# get confidence interval for heritability
HPDinterval(herit)

#plot the trace of heritability, should not be any pattern
plot(herit)
#plot posterior distribution of heritability with histogram, 95% CI, and mean
postPlot(herit)
#estimate 5, 20, and 90% quantiles of posterior distribution
quantile(herit,na.rm = T,probs = c(.05, .20,.90))


#########################Maternal environmental effects and variance ratio of maternal environmental effects

#estimate posterior distribution of additive genetic effects and get the 95% CI
posterior.mode(model1101$VCV[,"Damid"])
HPDinterval(model1101$VCV[,"Damid"])

#plot posterior distribution of additive genetic effects
plot(model1101$VCV[,"Damid"])
#plot posterior distribution of heritability with histogram, 95% CI, and mean
postPlot(model1101$VCV[,"Damid"])
#estimate 5, 20, and 90% quantiles of posterior distribution
quantile(model1101$VCV[,"Damid"],na.rm = T,probs = c(.05, .20,.90))


#get variance component ratio for maternal environmental effects
MatEnvEff1101 <- model1101$VCV[, "Damid"]/(model1101$VCV[, "animal"] + model1101$VCV[, "Damid"] + model1101$VCV[, "units"])

#effective size
effectiveSize(MatEnvEff1101)

#get the postrior mode
posterior.mode(MatEnvEff1101)

#get the posterio mean
mean(MatEnvEff1101)

#get the confidence interval
HPDinterval(MatEnvEff1101)

#plot posterior distribution 
plot(MatEnvEff1101)
#plot posterior distribution with histogram, 95% CI, and mean
postPlot(MatEnvEff1101)
#estimate 5, 20, and 90% quantiles of posterior distribution
quantile(MatEnvEff1101,na.rm = T,probs = c(.05, .20,.90))




######################################################################################################################################################
######################################################################################################################################################

# Model 1102: Trait = Ptype1, Batch fixed, Additive effects and maternal genetic effects are random

#needed library for running generalized linear mixed models
library(MCMCglmm)
#Using remotes install Matthew Wolak's package wolakR from github which provides better graphics for posterior distributions
remotes::install_github("matthewwolak/wolakR")
#Load wolakR library
library(wolakR)

#Clear memory if needed
rm(list=ls()) 

load("C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\Toy4.rat")
load("C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\Toy4ped.rat")

# set up parameters for the priors
# Fixed effects are automatically set by MCMCglmm to follow a normal distribution and do not need to be specified
# G is for specified random effects (additive, maternal, etc). 
# Here we set for weak priors that will be used in an inverse gamma disribution automatically set by MCMCglmm
# R is residual effects for each specified random effect and follows same rules as G
prior1102 <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002), G2 = list(V = 1, nu = 0.002)))

#Set up inverse of Numerator Relationship Matrix using MCMCglmm inverseA command
invA <- inverseA(Toy4ped)$Ainv 

# model statement
model1102 <- MCMCglmm(Ptype1 ~ 1 + Batch,                           #intercept (the 1) and Batch are the fixed effect
                      random = ~animal + Damid,                     #additive (animal) and maternal (Damid) are random effects
                      family = "gaussian",                          #phenotype has gaussian distribution
                      prior = prior1102,                            #call the priors parameters defined above
                      data = Toy4,                                  #call the data file
                      nitt = 1000000,                               #number of MCMC iterations 
                      burnin = 10000,                               #number of iterations for burnin
                      thin = 1000,                                  #sampling interval
                      ginverse=list(animal=invA,Damid=invA))        #call the inverse of the NRM for both additive and mat gen effects


#save model as an R object so we can access it later
save(model1102, file = "C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\model1102.obj") #dont forget to change the path#
#load model file
load(file = "C:\\Users\\pacarter\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Teaching\\Biol 521\\model1102.obj") #dont forget to change the path#

#########################Diagnostics

#plot trace and density of fixed effects; should be no trend in trace
plot(model1102$Sol)
#plot trace and density of random (additive, maternal and residual (=environmental)) variances; should be no trend in trace
plot(model1102$VCV)

#examine autocorrelation of fixed effects
autocorr.diag(model1102$Sol)
#examine autocorrelation of random (additive, maternal and residual) variances
autocorr.diag(model1102$VCV)

#check effective population size for fixed effects; should be gt 1000
effectiveSize(model1102$Sol)
#check effective population size for random effects (additve, maternal and residual variances); should be gt 1000
effectiveSize(model1102$VCV)

#test of convergence, p should be greater than 0.05 for good convergence
heidel.diag(model1102$VCV)

#estimates of additive and residual variances
posterior.mode(model1102$VCV)

#summary of model; make sure to check DIC score (smaller is better)
summary(model1102)

#########################Additive Genetic Effects and Heritability (i.e., variance ratio of additive genetic effects)

#estimate posterior distribution of additive genetic effects and get the 95% CI
posterior.mode(model1102$VCV[,"animal"])
HPDinterval(model1102$VCV[,"animal"])

#plot posterior distribution of additive genetic effects
plot(model1102$VCV[,"animal"])
#plot posterior distribution of heritability with histogram, 95% CI, and mean
postPlot(model1102$VCV[,"animal"])
#estimate 5, 20, and 90% quantiles of posterior distribution
quantile(model1102$VCV[,"animal"],na.rm = T,probs = c(.05, .20,.90))

#estimate posterior distribution of the heritability (animal vairance divided by animal + maternal +  residual variances)
herit <- model1102$VCV[, "animal"]/(model1102$VCV[, "animal"] + model1102$VCV[, "Damid"] + model1102$VCV[, "units"])

#effective sample size for heritability should be gt 1000
effectiveSize(herit)

# get the mean from the posterior disribution of heritability
mean(herit)

# get confidence interval for heritability
HPDinterval(herit)

#plot the trace of heritability, should not be any pattern
plot(herit)
#plot posterior distribution of heritability with histogram, 95% CI, and mean
postPlot(herit)
#estimate 5, 20, and 90% quantiles of posterior distribution
quantile(herit,na.rm = T,probs = c(.05, .20,.90))


#########################Maternal genetic effects and variance ratio of maternal environmental effects

#estimate posterior distribution of additive genetic effects and get the 95% CI
posterior.mode(model1102$VCV[,"Damid"])
HPDinterval(model1102$VCV[,"Damid"])

#plot posterior distribution of additive genetic effects
plot(model1102$VCV[,"Damid"])
#plot posterior distribution of heritability with histogram, 95% CI, and mean
postPlot(model1102$VCV[,"Damid"])
#estimate 5, 20, and 90% quantiles of posterior distribution
quantile(model1102$VCV[,"Damid"],na.rm = T,probs = c(.05, .20,.90))


#get variance component ratio for maternal environmental effects
MatEnvEff1102 <- model1102$VCV[, "Damid"]/(model1102$VCV[, "animal"] + model1102$VCV[, "Damid"] + model1102$VCV[, "units"])

#effective size
effectiveSize(MatEnvEff1102)

#get the posterior mode
posterior.mode(MatEnvEff1102)

#get the posterior mean
mean(MatEnvEff1102)

#get the confidence interval
HPDinterval(MatEnvEff1102)

#plot posterior distribution 
plot(MatEnvEff1102)
#plot posterior distribution with histogram, 95% CI, and mean
postPlot(MatEnvEff1102)
#estimate 5, 20, and 90% quantiles of posterior distribution
quantile(MatEnvEff1102,na.rm = T,probs = c(.05, .20,.90))











