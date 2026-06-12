# Set working directory
setwd(".")

## Required resources for the analysis: 
# (i) Properties per peptide of the aaindex database from "seqinr" package. As stated by Fusaro et al. (2009): "For each peptide and a given property, the constituent amino acid numerical values were averaged to produce a single value. Missing values were ignored. The average (rather than median or sum) was chosen because it is sensitive to outliers and normalizes for peptide length. It was assumed that the average physicochemical property across each peptide was sufficient to capture relevant information about peptide response."  
# (ii) Random Forest classifier of peptide detectability: "rfFit" object developed with "caret" R package.


library(Biostrings)
library(seqinr)
data(aaindex)

load("rfFit.rda")
library(caret)




#################################################
##       Properties per peptide                ##  
#################################################

pep_example <- c("AMGIMNSFVNDIFER","LLYAIEETEGFGQE")

stats <- sapply(pep_example, FUN=function(x) AAstat(s2c(x), plot=FALSE))
prop_mean <- sapply(pep_example, FUN=function(x) c(lapply(aaindex, FUN=function(y) mean(y$I[aaa(s2c(x))], na.rm=TRUE)), Length=getLength(x), PMW=pmw(s2c(x)), PI=computePI(s2c(x)), unlist(stats["Prop",x])))


#################################################
##   Classification of peptide sequences       ##
#################################################

# 106 peptide properties are considered by the classifier of the more than 500 peptide properties available 
prop_mean_rf <- prop_mean[predictors(rfFit),]

# Class definitions: "MObs" = detectable peptides, "LObs" = non-detectable peptides 
rfClasses <- apply(prop_mean_rf, 2, FUN=function(x) predict(rfFit, newdata = x))

# Class probabilities: "MObs" = probability of being a detectable peptide, "LObs" = probability of being a non-detectable peptide 
rfProbs <- apply(prop_mean_rf, 2, FUN=function(x) predict(rfFit, newdata = x, type="prob"))


#-----------------------------------------------------------------------------------------------------------------------------

## We have analyzed millions of peptide sequences using parallel computing in a cluster based on this script.

