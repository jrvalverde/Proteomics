# Set working directory
setwd(".")

## Required resources for the analysis: 
# (i) Properties per peptide of the aaindex database from "seqinr" package. As stated by Fusaro et al. (2009): "For each peptide and a given property, the constituent amino acid numerical values were averaged to produce a single value. Missing values were ignored. The average (rather than median or sum) was chosen because it is sensitive to outliers and normalizes for peptide length. It was assumed that the average physicochemical property across each peptide was sufficient to capture relevant information about peptide response."  
# (ii) Random Forest classifier of peptide detectability: "rfFit" object developed with "caret" R package.


library(Biostrings)
library(seqinr)
data(aaindex)
library(fs)
library(data.table)

load("bin/detectability/rfFit.rda")
library(caret)

#################################################
##      Get input peptides from Fasta          ##
#################################################

indb <- "in.fasta"
oudb <- "ou.fasta"

cat("Processing", indb, '\n')
cat("    Reading", indb, '\n')

# read in the fasta file as a list of sequences with attributes
seqs <- read.fasta(indb, seqtype="AA", as.string=TRUE, set.attributes=TRUE)

cat("    Converting to data.frame\n")
# seqs is a list, convert to data.frame
seq <- c()
name <- c()
annot <- c()
i <- 0
for (s in seqs) {
    seq <- c(seq, s[1])
    name <- c(name, attr(s, 'name'))
    annot <- c(annot, attr(s, 'Annot'))
    #i=i+1
    #if (i % 1000) cat(i, ' ')
}
fa.df <- data.frame(seq=seq, name=name, annot=annot, class=NA, plobs=0, pmobs=0)
    

#################################################
##       Properties per peptide                ##  
#################################################

# peptides
pep_list <- fa.df$seq

# calculate statistics
cat("    Calculating AA statistics (slow)\n")
stats <- sapply(pep_list, FUN=function(x) AAstat(s2c(x), plot=FALSE))
# produces a complex data structure with aa statistics for each peptide
# it has 3 rows and as many columns as peptides, rows are
#   composition
#   properties (a list with physico-chemical properties)
#   pi (theoretical isoelectric point)

# calculate mean properties (
cat("    Calculating mean properties (very slow)\n")
# gets 554 properties for each amino acid (aaindex), calculates 
# peptide molecular weight (pmw) and theoretical isoelectric point
# (computePI)
prop_mean <- sapply(pep_list, FUN=function(x) c(lapply(aaindex, FUN=function(y) mean(y$I[aaa(s2c(x))], na.rm=TRUE)), Length=getLength(x), PMW=pmw(s2c(x)), PI=computePI(s2c(x)), unlist(stats["Prop",x])))
# this generates a HUGELY HUMONGOUS list

#################################################
##   Classification of peptide sequences       ##
#################################################

cat("    Classifying with RF model\n")
# 106 peptide properties are considered by the classifier of the more 
# than 500 peptide properties available in prop_mean, so we will select
# only those properties for the prediction
cat("        selecting actual RF predictors from prop_mean\n")
prop_mean_rf <- prop_mean[predictors(rfFit),]


# Make the predictions
# Class definitions: 
#   "MObs" = detectable peptides, 
#   "LObs" = non-detectable peptides 
cat("        predicting MObs vs LObs (slow)\n")
rfClasses <- apply(prop_mean_rf, 2, FUN=function(x) predict(rfFit, newdata = x))

# at this point we would like to cbind the MObs/LObs status to the
# fa.df dataframe
fa.df$class <- as.character(rfClasses)

# Class probabilities: 
#   "MObs" = probability of being a detectable peptide, 
#   "LObs" = probability of being a non-detectable peptide 
cat("        predicting probabilities\n")
rfProbs <- apply(prop_mean_rf, 2, FUN=function(x) predict(rfFit, newdata = x, type="prob"))

# at this point we would like to cbind the probabilities to the dataframe
# rfProbs is a list with as many components as peptides, each named after
# the corresponding peptide, and each component being in turn a list of
# 2 variables, MObs and Lobs, being their probabilities.
for ( i in 1:dim(fa.df)[1] ) {
    fa.df$plobs[i] <- rfProbs[[ i ]]$LObs
    fa.df$pmobs[i] <- rfProbs[[ i ]]$MObs
}

cat("    Saving to", oudb, "\n")
# once we have the complete data frame, we need to generate the
# output file:
#   for each peptide, output "Annot", type, pMObs, pLObs, and sequence
#
fa.list <- list()
for ( i in 1:dim(fa.df)[1] ) {
    fa.entry <- paste(
        ">", fa.df$name[i],
        " ", gsub('>', ' : ', fa.df$annot[i]), 
        " class=", fa.df$class[i], 
        " pLObs=", fa.df$plobs[i], 
        " pMobs=", fa.df$pmobs[i], 
        "\n", fa.df$seq,
        sep='')
    fa.list[[1]] <- fa.entry
}

file_create(oudb)
fwrite(fa.list, file=oudb, sep='', quote=FALSE)
