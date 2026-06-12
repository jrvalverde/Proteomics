library.cran <- function(pkg) {
    if (!library(pkg, quietly=T, character.only=T, logical.return=T)) {
	install.packages(pkg, dependencies=T)
	if (!library(pkg, logical.return=T)) {
	  cat('\n*** *** ***', pkg, 'failed *** *** ***\n\n')
	  stop()
	}
    }
}

require.cran <- function(pkg) {
    if (!require(pkg, quietly=T, character.only=T)) {
	install.packages(pkg, dependencies=T)
	if (!require(pkg)) {
	  cat('\n*** *** ***', pkg, 'failed *** *** ***\n\n')
	  stop()
	}
    }
}

library.bioc <- function(pkg, ...) {
    # force=T may be used to force re-installation EVERY TIME
    # This may be helpful to force an update of a pkg installed
    # in a read-only system directory into a writable user-directory.
    # Normally, it should be left as FALSE

    if ( ! library(pkg, quietly=T, character.only=T, logical.return=T) )  {
        BiocManager::install(pkg, ...)
	if ( ! library(pkg, quietly=T, character.only=T, logical.return=T)) {
            cat('\n*** *** ***', pkg_i, 'failed *** *** ***\n\n')
            stop()
	}
    }
}

require.bioc <- function(pkg, ...) {
    # force=T may be used to force re-installation EVERY TIME
    # This may be helpful to force an update of a pkg installed
    # in a read-only system directory into a writable user-directory.
    # Normally, it should be left as FALSE
    
    if ( ! require(pkg, quietly=T, character.only=T)) {
        BiocManager::install(pkg, ...)
	if ( ! require(pkg, quietly=T, character.only=T)) {
            cat('\n*** *** ***', pkg_i, 'failed *** *** ***\n\n')
            stop()
	}
    }
}

# note: in the following two, we should check availability of 'remotes'
# or we should have loaded it already (which we now do)

library.github <- function(pkg.src, ...) {
    # match start-of-string (^), anything (.*) and, lazily (?), a '/'
    # substitute by nothing ('') in the package github source
    pkg <- gsub('^.*?/', '', pkg.src)
    if (!library(pkg, quietly=T, character.only=T, logical.return=T)) {
        remotes::install_github(pkg.src, ...)
	if (!library(pkg, quietly=T, character.only=T, logical.return=T)) {
            cat('\n*** *** ***', pkg_i, 'failed *** *** ***\n\n')
            stop()
	}
    }
}

require.github <- function(pkg, ...) {
    # match start-of-string (^), anything (.*) and, lazily (?), a '/'
    # substitute by nothing ('') in the package github source
    pkg <- gsub('^.*?/', '', pkg.src)
    if (!require(pkg, quietly=T, character.only=T)) {
        remotes::install_github(pkg, ...)
	if (!require(pkg, quietly=T, character.only=T)) {
            cat('\n*** *** ***', pkg_i, 'failed *** *** ***\n\n')
            stop()
	}
    }
}

## Install missing packages
packages <- c("remotes", "dplyr", "ggplot2", "remotes", "data.table",
              "reshape2", "pheatmap", "ggpubr", "patchwork", "httr",
	      "wrProteo", "wrMisc", "wrGraph", "limma", "MSnbase",
	      "SummarizedExperiment", "vsn", "MaldiQuant", "MSnSet.utils",
	      "MaldiQuantForeign", "BRAIN", "Rdisop", "OrgMassSpecR",
	      "synapter", "qvalue", "isobar", "DEP", "msmsTests",
	      "clusterProfiler", "pathview", "enrichplot", "DOSE",
	      "UniProt.ws", "KEGGRES", "org.At.eg.db", "knitr", 
	      "wrMisc", "wrProteo", "wrGraph")
	      
for (pkg in packages) {
  if (!require(pkg, quietly = T, character.only = T))
    BiocManager::install(pkg, update=F)	# upd=F to avoid unexpected interactive questions
}

library.github("PNNL-Comp-Mass-Spec/MSnSet.utils")
library.github('steven-shuken/permFDP')

## ------------------------

# load packages
library(dplyr)			## for convenience
library(data.table)		## for renaming duplicates
library("RColorBrewer") 	## Color palettes
library("ggplot2")  		## Convenient and nice plotting
library("reshape2") 		## Flexibly reshape data
library("pheatmap")
library("ggpubr")
library("patchwork")
library("httr")
library('remotes')
library("qvalue")		## for q.value calculation
library("knitr")		## 

library(wrProteo)
library("MSnbase")		## for MSnSet manipulation
library(MSnSet.utils)
library("SummarizedExperiment")
library("vsn")			## for VSN normalization
library("MALDIquant")
library("MALDIquantForeign")
library("BRAIN")
library("Rdisop")
library("OrgMassSpecR")

library("synapter")		## for synapter processing

#library("isobar")
library("DEP")			## for DEP processing

library("msmsTests")		## for edgeR processing

library(clusterProfiler)	## for annotation
library(pathview)
library(enrichplot)
library(DOSE)
library(UniProt.ws)
library(org.At.eg.db)
library(KEGGREST)

library("wrMisc")
library("wrProteo")
library("wrGraph")

library("limma")


# for multifile
#abundanceFile <- "Abundances_Scaled.tab"
#abundanceFile <- "Abundances_Normalized.tab"

# for single file
#exprColsPrefix="SclAbun"
exprColsPrefix="NormAbun."
#exprColsPrefix="No.PSMs"

# number of replicas per sample	(used to calculate offsets)
n.replicas <- 4

log2cutoff <- 1
p.cutoff <- 0.05

by.row <- 1
by.col <- 2


pos.left <- 1
pos.above <- 2
pos.below <- 3
pos.right <- 4
out.text <- function( labels, ... ) {
    n.lines <- length(labels)
    plot.new()
    text(x=0, y=seq(1, 1-(0.01*(n.lines-1)), -0.01), 
         labels=labels, 
         pos=pos.right, 
         ...)
}

out <- paste('Analysis', exprColsPrefix, 'pdf', sep='.')
pdf(out, paper="a4")


#
#
# https://www.bioconductor.org/packages/release/data/experiment/vignettes/RforProteomics/inst/doc/RforProteomics.html#4_Quantitative_proteomics
#
#

# read in the data
# using separate files
# (left here for reference just in case we need it in the future)
#abun <- readMSnSet(exprsFile=abundanceFile, 
#		  featureDataFile="Features.tab", 
#                  phenoDataFile="Phenotypes.tab", 
#                  sep='\t' )
                  
# using a single file (plus a metadata file)
# data.tab should have been saved as a TAB file with " string delimiter
#
# 1. get the metadata
phenoData <- read.table('metadata.tab', header=T)
# set column names to use for counts
exprCols <- sub("Line", exprColsPrefix, phenoData$label)
phenoData$label <- sub("Line", exprColsPrefix, exprCols)
rownames(phenoData) <- phenoData$label		#rownames(phenoData) <- rownames(pData(abun))

# 2. get all the data
#    data table will be first read (and col.names converted to R) prior to
#    conversion to MSnSet, so we need to use "R-converted" column names.
fnames <- "pep.mod.id"

abun <- readMSnSet2("data.tab", ecol=exprCols, fnames=fnames, sep='\t')
head(abun)
head(exprs(abun))
head(fData(abun))
head(pData(abun))	# not present in the data file, we need to fill it

# 3. populate pData
#pData(abun) <- cbind(pData(abun), read.table('metadata.tab', header=T)
#	or the equivalent, safer and more general
pData(abun) <- phenoData
pData(abun)


# we can change sample names to something more useful with
# sampleNames(qnt) <- c(sub(exprColsPrefix, '', colnames(s)))

# 4. set row names
# for convenience, we want to have more useful rownames(abun)
#    we cannot use accession as it repeats, so we need a column
#    with "accession.number"
#	This will identify each peptide by the protein it belongs to,
#	and a sequential number (separated by a dot ('.')
fData(abun)$acc.no <- paste(fData(abun)$accession, 
			    rowid(fData(abun)$accession), 
			    sep='.')
rownames(abun) <- fData(abun)$acc.no
rownames(exprs(abun)) <- fData(abun)$acc.no

# now we have the data in suitable a MSnSet structure


# process
# get filtered quantitation data (remove rows with NA (undetected peptides)
qnt <- filterNA(abun)
# process it
processingData(qnt)		# show all data processing history
# get protein quantitation data
#	we will group by protein accession code
protqnt <- combineFeatures(qnt,
                           groupBy = fData(qnt)$accession,
                           method = sum)
#
# When applied to peptide normalized abundances, the last step gives
# us the same abundances as reported in the proteins.data.tab file.
# Thus, further processing should give similar results, although this
# time, shall we use peptide data?
#


# plot N random entries 
n.prots <- 8          

# allocate 1/4 of the page for text
#layout(matrix(c(1, 1, 2, 2, 2, 2, 2, 2), 4, 2, byrow = TRUE))
#layout.show(2)
#out.text(labels=paste("Protein intensity for", n.prots, "last proteins", sep=' '))

set.seed(16160423)		# Shakespeare, Cervantes, Inca Garcilaso die
n.random.prots <- sample( 1:dim(protqnt)[1], n.prots, replace=F )
exprs.to.plot <- exprs(protqnt)[n.random.prots, ]
acc.to.plot <- rownames(protqnt)[n.random.prots]

cls <- brewer.pal(n.prots, "Set1")
matplot(t(exprs.to.plot), type = "b",
        lty = 1, col = cls,
        ylab = "Protein intensity (summed peptides)",
        xlab = "Measure group")
legend("topright", tail(featureNames(protqnt), n=n.prots),
       lty = 1, bty = "n", cex = .8, col = cls)
 
 
# normalize (starting from qnt, the total quantitation data, 
# not from protqnt, the summarized protein data)
qntS <- normalise(qnt, "sum")
qntV <- normalise(qntS, "vsn")
qntV2 <- normalise(qnt, "vsn")

protqntS <- normalise(protqnt, "sum")
protqntV <- normalise(protqntS, "vsn")	# worrysome
protqntV2 <- normalise(protqnt, "vsn")	# worrysome
# plot SD vs means (by row)
#meanSdPlot(qntV)
#meanSdPlot(qntV2)

# Pick a few datums to display


# fData(qnt)$accession gets the accession column from the data in qnt
#	select all peptides for each of the n.prot pre-selected proteins
idx <- sapply(acc.to.plot, grep, fData(qnt)$accession)
#idx <- sapply(acc.to.plot, grep, fData(qnt)$acc.no)

#	select up to the first three peptides for each of the n.prot 
#	pre-selected proteins
idxs <- sapply(idx, head, 3)
# some peptides are on more than one protein because some protein entries
# seem to be duplicated as alone and in a list
small <- qntS[ unique(sort(unlist(idxs))), ]	# unlist idxs and extract them from qntS

#	select up to the first ten peptides for each of the n.prot 
#	pre-selected proteins
idxm <- sapply(idx, head, 10)
medium <- qntV[unique(sort(unlist(idxm))), ]

s <- exprs(small)	# from qntS and head(3)
m <- exprs(medium)	# from qntV and head(10)
head(s)
head(m)
# use more convenient identifiers
colnames(s) <- sub(exprColsPrefix, '', colnames(s))
colnames(m) <- sub(exprColsPrefix, '', colnames(m))

# we must use rownames to identify unique proteins, so we use the accession
rownames(s) <- fData(small)$accession
rownames(m) <- fData(medium)$accession

# change chosen gene names to something more readable
#rownames(m)[grep("CYC", rownames(m))] <- "CYT"
#rownames(m)[grep("ENO", rownames(m))] <- "ENO"
#rownames(m)[grep("ALB", rownames(m))] <- "BSA"
#rownames(m)[grep("PYGM", rownames(m))] <- "PHO"
#rownames(m)[grep("ECA", rownames(m))] <- "Background"


# draw the heatmaps
cls <- c(brewer.pal(length(unique(rownames(m)))-1, "Set1"),
         "grey")
names(cls) <- unique(rownames(m))
wbcol <- colorRampPalette(c("white", "darkblue"))(256)

msg=paste('Heatmap of Sum-normalized', n.prots, 'random proteins')
heatmap(s, col = wbcol, RowSideColors=cls[rownames(s)], main=msg)

msg=paste('Heatmap of Sum+VSN-normalized', n.prots, 'random proteins')
heatmap(m, col = wbcol, RowSideColors=cls[rownames(m)], main=msg)


# draw spikes plots
dfr.s <- data.frame(exprs(small),
                  Protein = as.character(fData(small)$accession),
                  Feature = featureNames(small),
                  stringsAsFactors = FALSE)

colnames(dfr.s) <- c(sub(exprColsPrefix, '', colnames(s)),
                   "Protein", "Feature")
# set easier names
#dfr$Protein[dfr$Protein == "sp|P00924|ENO1_YEAST"] <- "ENO"
#dfr$Protein[dfr$Protein == "sp|P62894|CYC_BOVIN"]  <- "CYT"
#dfr$Protein[dfr$Protein == "sp|P02769|ALBU_BOVIN"] <- "BSA"
#dfr$Protein[dfr$Protein == "sp|P00489|PYGM_RABIT"] <- "PHO"
#dfr$Protein[grep("ECA", dfr$Protein)] <- "Background"
dfr.s.2 <- melt(dfr.s)
## Using Protein, Feature as id variables

ggplot(aes(x = variable, y = value, colour = Protein),
       data = dfr.s.2) +
  geom_point() +
  geom_line(aes(group=as.factor(Feature)), alpha = 0.5) +
  facet_grid(. ~ Protein) + theme(legend.position="none") +
  labs(x = "Reporters", y = "Normalised intensity") +
  ggtitle("Sum-normalized plots")

#ggsave(plot = myplot, filename = "myplot.pdf", device = "pdf")

dfr.m <- data.frame(exprs(medium),
                  Protein = as.character(fData(medium)$accession),
                  Feature = featureNames(medium),
                  stringsAsFactors = FALSE)

colnames(dfr.m) <- c(sub(exprColsPrefix, '', colnames(s)),
                   "Protein", "Feature")
# set easier names
#dfr$Protein[dfr$Protein == "sp|P00924|ENO1_YEAST"] <- "ENO"
#dfr$Protein[dfr$Protein == "sp|P62894|CYC_BOVIN"]  <- "CYT"
#dfr$Protein[dfr$Protein == "sp|P02769|ALBU_BOVIN"] <- "BSA"
#dfr$Protein[dfr$Protein == "sp|P00489|PYGM_RABIT"] <- "PHO"
#dfr$Protein[grep("ECA", dfr$Protein)] <- "Background"
dfr.m.2 <- melt(dfr.m)
## Using Protein, Feature as id variables

ggplot(aes(x = variable, y = value, colour = Protein),
       data = dfr.m.2) +
  geom_point() +
  geom_line(aes(group=as.factor(Feature)), alpha = 0.5) +
  facet_grid(. ~ Protein) + theme(legend.position="none") +
  labs(x = "Reporters", y = "Normalised intensity") +
  ggtitle("Sum+VSN-normalized plots")



#############################################################################
#			ANALYZE DIFFERENTIAL EXPRESSION
#############################################################################

#
#
# DEA
#

# https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/DEA.html

#m <- abun

#m <- filterNA(m)

# ...


#
#
# synapter
#
#

# library(synapter)
# synapterGuide()
#
# Trivial t-test based analysis
#
# we start from an MSnSet like protqnt, qnt, qntS, qntV or qntV2
# ...
#
# the TOP3 approach uses the 3 most intense peptides to compute protein
# intensity.
# As we have protein data already, we can skip it and initial filtering
# and normalisation.

# start from data at the peptide level
if (exprColsPrefix == "No.PSMs.") {
    syn_data.pep <- qntV
} else if (exprColsPrefix == "NormAbun.") {
    syn_data.pep <- qntV2		# data is already sum-normalized
}

# convert expressions to log2
exprs(syn_data.pep) <- log2(exprs(syn_data.pep))



# Return the p.value of a permutation comparison test
# this function is designed to be used with apply, and so,
# it can only take one argument: the two sets must have been
# combined in a single vector cbind(set1, set2)
#
p.perm.test <- function(row) {
    # since we can only take one value in apply, we must FORCE some
    # values and assumptions:
    n.permutations <- 1000	### XXX JR XXX ### HARD CODED!!!
    n.cols <- length(row)
    if (n.cols %% 2 != 0) {	### XXX JR XXX ### ONLY TWO SAME-SIZE GROUPS
       cat('p.perm.test: both subsets must have the same number of columns!\n')
       cat('p.perm.test: returnin p=1 to signal invalid comparison\n')
       return (1)
    }
    
    sss <- n.cols / 2	# subset size
    data <- data.frame(x=row, grp=factor(rep(c("A", "B"), c(sss, sss))))
    
    # Note that oneway_test can check two or more groups, but we
    # are using it for only two groups here, using the Fisher-Pitman
    # permutation test, and that blocking needs not be uniform, but
    # we force both groups to have the same size.
    # Possibilities in coin are
    #    oneway_test	(Fisher-Pitman permutation test)
    #    wilcox_test	(Wilcoxon-Mann-Whitney U test)
    #    kruskal_test	(Kruskal-Wallis test)
    #    normal_test	(van der Waarden test)
    #    median_test    (Brown-Mood median test)
    #    savage_Test	(Savage test)
    pvalue(
	coin::oneway_test(
	  x ~ grp, 
	  data=data, 
	  conf.level=0.95,
	  ties.method='average-scores',  # (def:average_scores) average scores of randomly broken ties
	  # ... independence_test extra options
	  distribution=approximate(nresample=n.permutations),	# use MonteCarlo resampling with n.perms
	  alternative='two.sided'
	)
    )[1]
}

# compute p for all rows using coin::oneway_test() permutation test
compare.perm.test <- function(X, Y) {
    ntests <- nrow(X)
    rslts <- as.data.frame(matrix(NA, nrow<-ntests, ncol<-2))
    names(rslts) <- c("ID", "pvalue")
        test.data <- cbind(X, Y)
    by.row <- 1
    
    rslts$ID <- rownames(X)
    rslts$pvalue <- apply(test.data, by.row, p.perm.test))
    return(rslts)
}

# compute p for all rows using RVAideMemoire::permutation.t.test
compare.perm.t.test <- function(X, Y) {
    # VERY, VERY SLOW
    ntests <- nrow(X)
    rslts <- as.data.frame(matrix(NA, nrow<-ntests, ncol<-2))
    names(rslts) <- c("ID", "pvalue")
        test.data <- cbind(X, Y)
    by.row <- 1
    set1.cols <- 1:ncol(X)
    set2.cols <- set1.cols + ncol(X)	# shift ncol(X)
    test.data <- cbind(X, Y)
    
    rslts$ID <- rownames(X)
    rslts$pvalue <- apply(test.data, by.row, 
      function(x) perm.t.test(x[set1.cols], 
                              x[set2.cols], 
			      nperm=n.permutations, 
			      progress=F
			      )$p.value
	         )
	    )
    return(rslts)
}

## compare all rows using Pearson correlation
compare.pearson <- function(X,Y){
    ntests <- nrow(X)
    rslts <- as.data.frame(matrix(NA, nrow<-ntests, ncol<-2))
    names(rslts) <- c("ID", "pvalue")
    for(i in 1:ntests){
       fit <- cor.test(X[,i],Y[,i],
		      na.action<-"na.exclude",
		      alternative<-"two.sided",        
		      method<-"pearson")       # pearson, kendall, spearman
       rslts[i, "ID"] <- rownames(X)[i]
       rslts[i, "pvalue"] <- fit$p.value
    }
    return(rslts)
} 

# compare all rows using t.test:
compare.t.test <- function(X, Y) {
    ntests <- nrow(X)
    rslts <- as.data.frame(matrix(NA, nrow<-ntests, ncol<-2))
    names(rslts) <- c('ID', 'pvalue')
    for (i in 1:ntests) {
	rslts[i, 'ID'] <- rownames(X)[i]
	rslts[i, 'pvalue'] <- t.test(X[i, ], Y[i, ])$p.value
    }
    return(rslts)
}



is.significant.permFDP <- function(x, y=NULL, 
   design=NULL,  # only used if y == NULL
   threshold=0.05, 
   n.permutations=1000)
{
    #threshold = 0.05
    #n.permutations = 1000
    if (! is.null(y)) {
        test.data <- data.frame(cbind(x, y))
        design = rep(c(1, 2), c(nrow(x), nrow(y)))
    } else {
        test.data <- x
	if ( is.null(design) ) {
	    # in this case we will assume two equal-size groups
	    design <- rep(c(1, 2), rep(ncol(test.data) / 2, 2))
	}
    }

    # we'll do all the comparisons using a MonteCarlo coin::identity_test 
    pt.pv <- apply(test.data, 1, p.perm.test)
    corrThreshold = permFDP::permFDP.adjust.threshold(
	pt.pv, 
	threshold, 
	design, 
	test.data, 
	n.permutations
    )
    return ( p.pv < corrThreshold )
}



# calculate comparison significance using permutations and the 'permFDP' 
# package
#	NOTE: we assume the same threshold is used for p-value
#	significance as for FDR
#	NOTE: permFDP uses internally t-tests for permutation comparisons
#	hence, it makes little sense to use this for something other than
#	straight t-test comparisons
#
# The function will return the corrected rejection threshold to control the FDR
# according to the rejection threshold you supplied. E.g., if you supply an
# uncorrected rejection threshold of 0.05, the estimated FDR will now be 5%
# using the new corrected rejection threshold.
# 
significance.permFDP <- function(x, y=NULL, 
			design=NULL, 		# only used if y == NULL
			threshold=0.05, 	# FP tolerance
			n.permutations=1000,	# seems sensible for now
			init.p.values='t.test', # permFDP uses internally t.test
			verbose=F)
{
    #threshold = 0.05
    #n.permutations = 1000
    
    if (! is.data.frame(x) && ! is.matrix(x)) {
        cat('x must be data frame or a matrix!\n')
	return(NULL)
    }
    # we have raw data and need to make the comparisons
    # do we have a full table or two separate tables?
    if (! is.null(y)) {
        if (! is.data.frame(y) && ! is.matrix(y)) {
            cat('y must be a data frame or a matrix!\n')
	    return(MULL)
	}
        if (verbose) cat('USING x and y\n')
        data <- data.frame(cbind(x, y))
        design = rep(c(1, 2), c(ncol(x), ncol(y)))
    } else {
        if (verbose) cat('USING x only\n')
        data <- x
	if ( ! is.null(design) ) {
	    design <- design
	} else {
	    if (verbose) cat('SPLITTING x evenly\n')
	    # in this case we will assume two equal-size groups
	    design <- rep(c(1, 2), c(ncol(data) / 2, ncol(data) / 2))
	}
    }
    # fix design: there should be only two levels, named 1 and 2
    if (length(levels(as.factor(design))) < 2 ) return(NULL)
    # if the number of levels is > 2, maybe we should consider
    # doing something else, like all pairwise comparisons...
    g1 <- levels(as.factor(design))[1]
    g2 <- levels(as.factor(design))[2]
    # if there are already groups labelled 1, and/or 2 they will be used 
    # because of the way levels() reports levels (numbers first, then
    # chars in alphabetical order, small-cap first for each alpha level
    design[ design == g1 ] <- 1
    design[ design == g2 ] <- 2
    if (verbose) cat('USING design =', design, '\n')
    
    if (verbose) cat('COMPUTING initial', init.p.values, 'p-values\n')
    if (init.p.values == 't.test') {
        # NOTE that design should only have values 1 and 2
	# as required by permFDR
	#
	p.values <- apply(data, 1 ,
            function(r) t.test(r[design == 1], r[design == 2])$p.value)
    }
    # the next shouldn't be used because the permutation function
    # (in C++) uses internally t-test, so the threshold may be not 
    # applicable to other tests
    else if (init.p.values == 'permutation') {
	# we'll do all the comparisons using a MonteCarlo coin::identity_test 
	if (verbose) cat('COMPUTING initial permutation p-values\n')
	p.values <- apply(data, 1, p.perm.test)
    }
    else if (init.p.values == 'permutation.t.test') {
        # use RVAideMemoire::perm.t.test
	if (verbose) cat('COMPUTING initial SLOW permutation t-test p-values\n')
	p.values <- apply(data, 1,
	      function(r) perm.t.test(r[design == "1"], 
                        	      r[design == "2"], 
				      nperm=n.permutations, 
				      progress=F
				      )$p.value
		    )

    }
    
    if (verbose) cat('COMPUTING permutation FDP threshold\n')
    corr.threshold = permFDP::permFDP.adjust.threshold(
    	p.values, 
	threshold, 
	design, 
	data, 
	n.permutations
	)
    p.is.significant <- p.values < corr.threshold

    # compute local-fdr q-values (default is lfdr.out=T)
    if (verbose) cat('COMPUTING q-values\n')
    q.values <- qvalue::qvalue(p.values, fdr.level=threshold, lfdr.out=T)

    pfdp <- list()	# prepare results
    pfdp[['data']]             <- data
    pfdp[['design']]           <- design
    pfdp[['calc.method']] <- paste(init.p.values, 'x t.test')
    pfdp[['threshold']]        <- threshold
    pfdp[['n.permutations']]   <- n.permutations
    pfdp[['p.values']]          <- p.values
    pfdp[['p.threshold']]      <- corr.threshold
    pfdp[['p.is.significant']] <- p.is.significant
    pfdp[['q.values']]          <- q.values$qvalues
    pfdp[['q.threshold']]      <- threshold
    pfdp[['q.is.significant']] <- q.values$significant
    return ( pfdp )
}


# calculate comparison significance using permutations and the 'cit' 
# package
#	NOTE: we assume the same threshold is used for p-value
#	significance as for FDR
# This should be more general as it allows using our own method for
# p-value calculation both initially and in the permutations
#
significance.cit <- function(X, Y, 
			     threshold=0.5, 
			     calc.p.by.row=compare.t.test, 
			     n.permutations=1000,
			     verbose=F) {

    ## calculate permutation-based FDR using 'cit' (see cit::fdr.cit() or fdr.od())
    if (nrow(X) == nrow(Y)) {    # should be the same as for Y
        n.row <- nrow(X)
    } else {
        return(NULL)
    }
    if (ncol(X) == ncol(Y)) {    # should be the same as for Y
        n.col <- ncol(X)
    } else {
        return(NULL)
    }
    
    # we will use a convenience function to calculate p values
    # that we can use to get unpermutated and permutated p-values

    if (verbose) cat('COMPUTING initial p-values\n')
    # Obtain observed results
    obs = calc.p.by.row(X,Y)

    pvalues <- obs$pvalue		# all p-values
    names(pvalues) <- obs$ID
    sigpval <- pvalues < threshold	# significant p-values
    nspv <- sum(sigpval)		# number of significant p-values
    if (nspv == 0) {	# nothing is significant
        cat("No significant p-values (pv <", threshold, ")\n")
	r <- list()
        r$X <- X
	r$Y <- Y
        r$calc.method <- deparse(substitute(calc.p.by.row))
        r$n.permutations <- n.permutations
	r$p.values <- pvalues
	r$p.threshold <- threshold
	r$p.is.significant <- sigpval
	r$fdr <- 1
	r$fdr.lower.conf.limit <- NULL
	r$fdr.upper.conf.limit <- NULL
	r$pi0 <- NULL
	r$overdispersion <- NULL
	r$fdr.observed.positive.tests <- 0
	r$fdr.total.positive.tests.over.all.permutations <- 0
	r$q.values <- NULL		# not computed (we use permutations
	r$q.threshold <- threshold	# to compute q-values, but it doesn't
	r$q.is.significant <- NULL	# make sense if no p-value < 0.05)
	return(r)
    }
    
    
    if (verbose) cat("COMPUTING", n.permutations, "permutations\n")
    ## Generate permuted results
    perml = vector('list',n.permutations)
    for (perm in 1:n.permutations) {
	# permutate rows in X
	Xperm = X[order(runif(n.row)),]
	# and compare to Y
	perml[[perm]] = calc.p.by.row(Xperm,Y)
	#cat('.')
    }
    
    if (verbose)  cat("ESTIMATING threshold for all significant p-values\n")
    
    # check each of the significant p-values one by one
    # start with an empty table
    ftbl <- data.frame(matrix(ncol=8, nrow=0))
    colnames(ftbl) <- c('pvalue', 'fdr', 'll', 'ul', 'pi0', 'c1', 'S', 'Sp')
    # check only significant p-values
    spvalues <- pvalues[sigpval]
    for (i in 1:nspv) {
	# compute FDR for this p-value
	f <- cit::fdr.od(obs$pvalue,	# USE ALL P-VALUES FOR FDR ESTIMATION
		perml,
		"pvalue",
		nrow(obs),
		spvalues[i]
		)
	# add it to the table
	f<- c(spvalues[i], f)
	ftbl[i,] <- f
    }
    # now sort the table by decreasing p-value (increasing threshold)
    ftbl <- ftbl[order(ftbl$pvalue), ]
    # now sort the table by decreasing p-value (increasing threshold)
    ftbl <- ftbl[order(ftbl$pvalue), ]

    if (verbose) cat("SEARCHING for threshold p-value with FDR >", threshold, '\n')
    # remove NAs
    ftbl = ftbl[!is.na(ftbl$fdr), ]
    # find pvalue that results in an FDR < threshold
    fth <- 0.	# nothing is significant
    for (i in 1:dim(ftbl)[1]) {
	#cat(ftbl[i, 'threshold'], ftbl[i, 'fdr'], '\n')
	if (ftbl[i, 'fdr'] < threshold) {
	    fth <- ftbl$pvalue[i]
	    break
	}
    }
    
    if (verbose) {		   
	print(
	    knitr::kable(ftbl, caption = "FDR: with multiplicity adjustment",
		     digits=3, row.names=FALSE)
	     )
	# and show significant p.values with an FDR < threshold
	# look up in above table at which threshold 
	# does FDR become < our desired rate value, e.g 0.2
	# this will give the first p-value below which FDR is 
	# below the threshold that we want
	print(
	    knitr::kable(obs[obs$pvalue < fth,], 
		    caption = paste(
		    	"P-values with permutation BH FDR <", threshold
			), 
		    digits=4, 
		    row.names=FALSE)
	     )
    }
    
    if (verbose) cat('COMPUTING FDR for p-value threshold', fth, ': ')
    ## FDR results
    f <- cit::fdr.od(obs$pvalue, perml, "pvalue", nrow(obs), fth)
    if (verbose) cat(f['fdr'], '\n')
    
    if (verbose) cat('COMPUTING permutation q-values\n')
    ## q-values
    q <- cit::fdr.q.perm(obs$pvalue, perml, 'pvalue', nrow(obs), cl=0.95)
    
    r <- list()
    r$X <- X
    r$Y <- Y
    r$calc.method <- deparse(substitute(calc.p.by.row))
    r$n.permutations <- n.permutations
    r$p.values <- obs$pvalue
    r$p.threshold <- fth
    r$p.is.significant <- obs$pvalue <= fth
    r$fdr <- f['fdr']
    r$fdr.lower.CI.limit <- f['fdr.ll']
    r$fdr.upper.CI.limit <- f['fdr.ul']
    r$fdr.pi0 <- f['pi.0']
    r$fdr.overdispersion <- f['od']
    r$fdr.observed.positive.tests <- f['s.obs']
    r$fdr.total.positive.tests.over.all.permutations <- f['s.perm']
    r$q.values <- q
    r$q.threshold <- threshold
    r$q.is.significant <- q <= threshold
    return(r)
}


# calculate comparison significance using permutations and the 'fdrci' 
# package
#	NOTE: we assume the same threshold is used for p-value
#	significance as for FDR
# This should be more general as it allows using our own method for
# p-value calculation both initially and in the permutations
#
significance.fdrci <- function(X, Y, 
			       threshold=0.5, 
			       calc.p.by.row=compare.t.test, 
			       n.permutations=1000,
			       verbose=F) {

    if (verbose) cat('COMPUTING p-values\n')
    ## calculate permutation-based FDR, using 'fdrci'
    rslts.obs <- calc.p.by.row(X, Y)

    # add q-values
    if (verbose) cat('COMPUTING q-values\n')
    rslts.obs[,"qvalue"] = p.adjust(rslts.obs[,"pvalue"], method="BH")
    if (verbose) {
        print(
            knitr::kable(rslts.obs[rslts.obs$qvalue < threshold,], 
		     caption = paste(
		     	"Q-values with parametic BH FDR <", threshold
			),
		     digits=4, 
		     row.names=FALSE)
        )
    }
    
    # find out number of statistically significant p-values
    pvalues <- rslts.obs$pvalue
    names(pvalues) <- rslts.obs$ID
    sigpval <- pvalues < threshold
    nspv <- sum(sigpval) 
    if (nspv == 0) {	# nothing is significant
        cat("No significant p-values (pv <", threshold, ")\n")
	r <- list()
        r$X <- X
	r$Y <- Y
        r$calc.method <- deparse(substitute(calc.p.by.row))
        r$n.permutations <- n.permutations
	r$p.values <- pvalues
	r$p.threshold <- threshold
	r$p.is.significant <- sigpval
	r$fdr <- 1
	r$fdr.lower.conf.limit <- NULL
	r$fdr.upper.conf.limit <- NULL
	r$pi0 <- NULL
	r$overdispersion <- NULL
	r$fdr.observed.positive.tests <- 0
	r$fdr.total.positive.tests.over.all.permutations <- 0
	r$q.values <- rslts.obs$qvalue
	r$q.threshold <- threshold
	r$q.is.significant <- rslts.obs$qvalue <= threshold
	return(r)
    }
    
    if (verbose) cat("COMPUTING", n.permutations, "permutations\n")
    ## to calculate a permutation FDR we need to do the permutations
    # make a list vector to hold permutationp-values
    perml <- vector('list', n.permutations)
    
    # repeatedly conduct analysis after permuting X
    for (perm in 1:n.permutations) {
        # generate randomly permutated p-values
	# permutate rows in X
	Xperm = X[order(runif(n.row)),]
	# and compare to Y
	perml[[perm]] = calc.p.by.row(Xperm,Y)
    }
    
    # If there are too many significant p-values, we will estimate
    # a cutoff using 100 increments from min to max
    if (nspv > 100) {
        if (verbose) cat("ESTIMATING threshold in 1/100th intervals\n")
	min.p <- min(pvalues) ; minl <- - log10(min.p)
	max.p <- max(pvalues[pvalues < 0.05]) ; maxl <- - log10(max.p)
	incr <- (minl - maxl) / 100

    	## EXPLORE potential FDR thresholds
	# to compute FDR we must provide a candidate discovery threshold
	# we can explore suitable thresholds by using -log10(p = 1, 2.1, 2.2,...,5
	# to get an output table where each row includes the estimated FDR,
	# 95%CI, pi0, overdispersion estimate, number of discoveries in observed
	# data and total number of discoveries in permuted data (i.e. basically
	# the same as we got with package 'CIT'
	
	if (verbose) {
	    cat("COMPUTING FDR threshold without multiplicity\n")   
	    # estimate FDR and 95% CI without accounting for multiplicity
	    # and see FDR estimate for each threshold (-log10(p))
	    ftbl <- fdrci::fdrTbl(pvalues,
    			   perml,
			   "pvalue",
			   length(pvalues),
			   lowerbound=maxl,	# minimal -log10(p)
			   upperbound=minl,	# maximal -log10(p)
			   incr=incr,	# increment from low to upp
			   cl=.95)		# CI confidence level

	    ftbl = ftbl[!is.na(ftbl$fdr), ]
	    ftbl$threshold <- 10 ^ -ftbl$threshold # show p instead of -log10(p)
	    print(
	        knitr::kable(ftbl, caption = "FDR: no multiplicity adjustment",
			 digits=3, row.names=FALSE)
	         )
	}
        if (verbose) cat("COMPUTING FDR threshold with BH multiplicity\n")
	# use BH for FDR selection and adjustment (see Millstein et al. 2022)
	ftblbh = fdrci::fdrTbl(pvalues,
		       perml,
		       "pvalue",
		       length(pvalues),
		       maxl,	# from -log10(p) = 2 
		       minl,	# to -log10(p) = 5
		       		# in 0.1 intervals
		       incr=incr,
		       correct="BH")
    } else {
        if (verbose) cat("ESTIMATING threshold for all significant p-values\n")
	# check each of the significant p-values one by one
	ftblbh <- data.frame(matrix(ncol=8, nrow=0))
	colnames(ftblbh) <- c('threshold', 'fdr', 'll', 'ul', 'pi0', 'c1', 'S', 'Sp')
	# check only significant p-values
	spvalues <- pvalues[sigpval]
	for (i in 1:nspv) {
	    # compute FDR for this p-value
	    f <- fdrci::fdr_od(pvalues,	# USE ALL P-VALUES
		    perml,
		    "pvalue",
		    nrow(rslts.obs),
		    spvalues[i]
		    )
	    # add it to the table
	    f<- c(-log10(spvalues[i]), f)
	    ftblbh[i,] <- f
	}
	# now sort the table by decreasing p-value (increasing threshold)
	ftblbh <- ftblbh[order(ftblbh$threshold), ]
    }
    
    if (verbose) cat("SEARCHING for threshold p-value with FDR >", threshold, '\n')
    # remove NAs
    ftblbh = ftblbh[!is.na(ftblbh$fdr), ]
    # find pvalue that results in an FDR < threshold
    fth <- 0.	# nothing is significant
    for (i in 1:dim(ftblbh)[1]) {
	#cat(ftblbh[i, 'threshold'], ftblbh[i, 'fdr'], '\n')
	if (ftblbh[i, 'fdr'] < threshold) {
	    fth <- ftblbh[i, 'threshold']
	    break
	}
    }
    # convert to p-value if a valid threshold has been found
    if (fth != 0.) fth <- 10 ^ -fth
    
    # use actual p-values instead of -log10(p)
    ftblbh$threshold <- 10 ^ -ftblbh$threshold # show p instead of -log10(p)
    
    if (verbose) {		   
	print(
	    knitr::kable(ftblbh, caption = "FDR: with multiplicity adjustment",
		     digits=3, row.names=FALSE)
	     )
	# and show significant p.values with an FDR < threshold
	# look up in above table at which threshold 
	# does FDR become < our desired rate value, e.g 0.2
	# this will give the first p-value below which FDR is 
	# below the threshold that we want
	print(
	    knitr::kable(rslts.obs[rslts.obs$pvalue < fth,], 
		    caption = paste(
		    	"P-values with permutation BH FDR <", threshold
			), 
		    digits=4, 
		    row.names=FALSE)
	     )
    }
    
    # do a single FDR calculation for the threshold chosen
    if (verbose) cat('COMPUTING FDR for p-value threshold', fth, ': ')
    f <- fdrci::fdr_od(rslts.obs$pvalue,
                perml,
		"pvalue",
		nrow(rslts.obs),
		fth)
    if (verbose) cat(f['fdr'], '\n')
    
    r <- list()
    r$X <- X
    r$Y <- Y
    r$calc.method <- deparse(substitute(calc.p.by.row))
    r$n.permutations <- n.permutations
    r$p.values <- rslts.obs$pvalue
    r$p.threshold <- fth
    r$p.is.significant <- rslts.obs$pvalue <= fth
    r$fdr <- f['fdr']
    r$fdr.lower.CI.limit <- f['fdr.ll']
    r$fdr.upper.CI.limit <- f['fdr.ul']
    r$fdr.pi0 <- f['pi.0']
    r$fdr.overdispersion <- f['c1']
    r$fdr.observed.positive.tests <- f['S']
    r$fdr.total.positive.tests.over.all.permutations <- f['Sp']
    r$q.values <- rslts.obs$qvalue
    r$q.threshold <- threshold
    r$q.is.significant <- rslts.obs$qvalue <= threshold
    return(r)
}



###### THIS REQUIRES MODIFICATION FOR MORE THAN TWO SAMPLES

### WE NEED A FOR LOOP TO CHECK ALL COMPARISONS ###
samples <- colnames(exprs(syn_data.pep))
samples <- sub(exprColsPrefix, '', samples)	# remove prefix
samples <- sub('.[0-9]+$', '', samples)		# remove replica number
samples <- unique(samples)			# remove duplications


syn_process <- function(syn_data.pep, sample.1, sample.2) {
    ## apply a t-test and extract the p-value
    ###	Here we are comparing the four replicas of sample 1
    ###	to the four replicas of sample 2
    ###	
    #sample.1 <- '260UP'	# select first sample to compare
    #sample.2 <- '260P'	# select second sample to compare
    
    # find each sample columns in the expressions
    group.1 <- grep(sample.1, colnames(exprs(syn_data.pep)))
    group.2 <- grep(sample.2, colnames(exprs(syn_data.pep)))

    # do a t-test between both samples
    pv <- apply(exprs(syn_data.pep), by.row, 
		function(x) { t.test( x[group.1], x[group.2])$p.value } )

    ## calculate q-values
    qv <- qvalue(pv)$qvalues

    ## calculate log2 fold-changes
    ###	exprs are already log2
    ###	so, if we get the mean long2 expression of each group and
    ###	find the difference, we get the log2 fold change
    lfc <- apply(exprs(syn_data.pep), by.row,
        	 function(x) {
	             mean(x[group.1], na.rm=TRUE) - mean(x[group.2], na.rm=TRUE)
		     }
        	)

    ## create a summary table
    syn_res.pep <- data.frame(cbind(exprs(syn_data.pep), pv, qv, lfc))
    ## reorder based on q-values
    syn_res.pep <- syn_res.pep[order(syn_res.pep$qv), ]
    colnames(syn_res.pep) <- c(sub(exprColsPrefix, '', colnames(s)), 
			       "p.value", "q.value", "log2FC")
    
    #knitr::kable(head(round(syn_res.pep, 3)))

    # save
    out=paste("synapter.pep", exprColsPrefix, sample.1, sample.2, "all.tab", sep='.')
    write.table(syn_res.pep, file=out, row.names=T, col.names=T)

    out=paste("synapter.pep", exprColsPrefix, sample.1, sample.2, "p_0.05.tab", sep='.')
    write.table(syn_res.pep[ syn_res.pep$p.value < 0.05, ], 
		file=out, row.names=T, col.names=T)

    # simplify adding protein description
    out=paste("synapter.pep", exprColsPrefix, sample.1, sample.2, "p_0.05.simple.tab", sep='.')
    
    ann <- fData(syn_data.pep)
    syn_res.pep$acc.no <- rownames(syn_res.pep)
    syn_res.pep$accession <- ann$accession[ match(rownames(syn_res.pep), ann$acc.no) ]
    syn_res.pep$Sequence <- ann$Sequence[ match(rownames(syn_res.pep), ann$acc.no) ]
    syn_res.pep$Modifications <- ann$Modifications[ match(rownames(syn_res.pep), ann$acc.no) ]
    
    # annotate the output
    # get protein data
    prot.data <- read.table('../proteins/data.tab', header=T, sep='\t', quote='"')
    
    syn_res.pep$Description <- prot.data$Description[ 
    		match(syn_res.pep$accession, prot.data$Accession) ]
    # save
    write.table(syn_res.pep[ syn_res.pep$p.value <= 0.05, 
                             c("p.value", "q.value", "log2FC",
	                       "Sequence", "Modifications", 
			       "accession", "Description") ],
	       file=out, 
	       row.names=F, 
	       col.names=T)

    # add protein description (repeating multiprotein peptides)
    out=paste("synapter.pep", exprColsPrefix, sample.1, sample.2, "p=0.05.simple.multi.tab", sep='.')
    # find multi-accession rows
    multis <- grepl(";", syn_res.pep$accession)
    # assign single-acc rows
    #syn_res.pep$Description <- prot.data$Description[ match(syn_res.pep$accession, prot.data$accession) ]
    # go over the dataset
    for (i in 1:nrow(syn_res.pep)) {
	# if it is a mult-acc row
	if ( multis[i] ) {
            #print(i)
            # copy it
            row <- syn_res.pep[i, ]
	    #print("original row")
	    #print(row)
	    ## separate the accession ids
            accs <- strsplit(as.character(row$accession), "; ")
	    # for each acc
	    for (j in 1:length(accs[[1]])) {
		#print(paste( "   ", j, accs[[1]][j] ) )
		# set row$accession
		row$accession <- accs[[1]][j]
		# and add it to the dataframe at the end to avoid altering indexing
		syn_res.pep <- rbind(syn_res.pep, row)
	    }
	    # now we could delete the old row but this will alter indexing
	    #syn_res.pep <- syn_res.pep[ -i, ]
	}
    }
    
    # re-identify duplicates in the new data frame
    multis <- grepl(";", syn_res.pep$accession)
    # remove multi-acc rows
    syn_res.pep <- syn_res.pep[ ! multis, ]
    # re-assign all single acc values
    syn_res.pep$Description <- prot.data$Description[ match(syn_res.pep$accession, prot.data$accession) ]
    # and save p-significant results (there are no q-significant results
    write.table(syn_res.pep[ syn_res.pep$p.value <= 0.05, 
    			    c("p.value", "q.value", "log2FC",
	                      "Sequence", "Modifications", 
			      "accession", "Description") ],
	       file=out, 
	       row.names=F, 
	       col.names=T)

    # volcano plot
    plot(syn_res.pep$log2FC, -log10(syn_res.pep$q.value),
	 col = ifelse(grepl(sample.2, colnames(syn_res.pep)),
	   "#4582B3AA",
	   "#A1A1A180"
	   ),
	 pch = 19,
	 xlab = expression(log[2]~fold-change),
	 ylab = expression(-log[10]~q-value))
    grid()
    abline(v = -1, lty = "dotted")
    abline(h = -log10(0.1), lty = "dotted")
    legend("topright", c(sample.2, sample.1),
	   col = c(
	       "#4582B3AA", "#A1A1A1AA"
	           ),
	   pch = 19, bty = "n")

    # heatmap
    msg <- paste("synapter trivial analysis of", 
                 exprColsPrefix, sample.1, sample.2)
    heatmap(as.matrix(syn_res.pep[c(group.1, group.2)]), main=msg)
    heatmap(as.matrix(syn_res.pep[syn_res.pep$p.value <= 0.05,
                                  c(group.1, group.2)]), 
            main=msg)

}


# now we have the sample names, we need a loop
for (s1 in samples) {
    for (s2 in samples){
        if (s1 == s2) next
	# compare both samples, plot, annotate and save the results
    	syn_process(syn_data.pep, s1, s2)
    }
}

# switch level
#
# switch to protein level and repeat analysis
# combine all peptude data into protein data
syn_data.pro <- combineFeatures(syn_data.pep,
                           groupBy = fData(syn_data.pep)$accession,
                           method = sum)


# convert expressions to log2
exprs(syn_data.pro) <- log2(exprs(syn_data.pro))

## apply a t-test and extract the p-value
pv <- apply(exprs(syn_data.pro), 1 ,
            function(x)t.test(x[1:2], x[3:4])$p.value)

## calculate q-values
qv <- qvalue(pv)$qvalues

## calculate log2 fold-changes
lfc <- apply(exprs(syn_data.pro), 1 ,
             function(x) mean(x[1:3], na.rm=TRUE)-mean(x[4:6], na.rm=TRUE))
## create a summary table
syn_res.pro <- data.frame(cbind(exprs(syn_data.pro), pv, qv, lfc))
## reorder based on q-values
syn_res.pro <- syn_res.pro[order(syn_res.pro$q.value), ]
colnames(syn_res.pro) <- c(sub(exprColsPrefix, '', colnames(s)),
			   "p.value", "q.value", "log2FC")
knitr::kable(head(round(syn_res.pro, 3)))

# save
out=paste("synapter.pro", exprColsPrefix, "all.tab", sep='.')
write.table(syn_res.pro, file=out, row.names=T, col.names=T)

out=paste("synapter.pro", exprColsPrefix, "p=0.05.tab", sep='.')
write.table(syn_res.pro[ syn_res.pro$p.value <= 0.05, ], 
	    file=out, row.names=T, col.names=T)

# simplify
out=paste("synapter.pro", exprColsPrefix, "p=0.05.simple.tab", sep='.')
# get protein data
prot.data <- read.table('../proteins/data.tab', header=T, sep='\t', quote='"')
# annotate the output
ann <- fData(syn_data.pro)
syn_res.pro$acc.no <- rownames(syn_res.pro)
syn_res.pro$accession <- ann$accession[ match(rownames(syn_res.pro), ann$acc.no) ]
syn_res.pro$Sequence <- ann$Sequence[ match(rownames(syn_res.pro), ann$acc.no) ]
syn_res.pro$Modifications <- ann$Modifications[ match(rownames(syn_res.pro), ann$acc.no) ]
syn_res.pro$Description <- prot.data$Description[ match(syn_res.pro$accession, prot.data$accession) ]
# save
write.table(syn_res.pro[ syn_res.pro$p.value <= 0.05, c("p.value", "q.value", "log2FC",
	   "Sequence", "Modifications", "accession", "Description") ],
	   file=out, row.names=F, col.names=T)

# add protein description repeating multi-protein peptides
out=paste("synapter.pro", exprColsPrefix, "p=0.05.simple.multi.tab", sep='.')
# find multi-accession rows
multis <- grepl(";", syn_res.pro$accession)
# assign single-acc rows
#syn_res.pro$Description <- prot.data$Description[ match(syn_res.pro$accession, prot.data$accession) ]
# go over the dataset
for (i in 1:nrow(syn_res.pro)) {
    # if it is a mult-acc row
    if ( multis[i] ) {
        #print(i)
        # copy it
        row <- syn_res.pro[i, ]
	#print("original row")
	#print(row)
	## separate the accession ids
        accs <- strsplit(as.character(row$accession), "; ")
	# for each acc
	for (j in 1:length(accs[[1]])) {
	    #print(paste( "   ", j, accs[[1]][j] ) )
	    # set row$accession
	    row$accession <- accs[[1]][j]
	    # and add it to the dataframe at the end to avoid altering indexing
	    syn_res.pro <- rbind(syn_res.pro, row)
	}
	# now we could delete the old row but this will alter indexing
	#syn_res.pro <- syn_res.pro[ -i, ]
    }
}
# re-identify duplicates in the new data frame
multis <- grepl(";", syn_res.pro$accession)
# remove multi-acc rows
syn_res.pro <- syn_res.pro[ ! multis, ]
# re-assign all single acc values
syn_res.pro$Description <- prot.data$Description[ match(syn_res.pro$accession, prot.data$accession) ]
# and save p-significant results (there are no q-significant results
write.table(syn_res.pro[ syn_res.pro$p.value <= 0.05, c("p.value", "q.value", "log2FC",
	   "Sequence", "Modifications", "accession", "Description") ],
	   file=out, row.names=F, col.names=T)

# volcano plot
plot(syn_res.pro$log2FC, -log10(syn_res.pro$q.value),
     col = ifelse(grepl("mut", colnames(syn_res.pro)),
       "#4582B3AA",
       "#A1A1A180"),
     pch = 19,
     xlab = expression(log[2]~fold-change),
     ylab = expression(-log[10]~q-value))
grid()
abline(v = -1, lty = "dotted")
abline(h = -log10(0.1), lty = "dotted")
legend("topright", c("mut", "wt"),
       col = c("#4582B3AA", "#A1A1A1AA"),
       pch = 19, bty = "n")
       
# heatmap
msg <- paste("synapter trivial analysis of", exprColsPrefix)
heatmap(as.matrix(syn_res.pro[1:length(exprCols)]), main=msg)
heatmap(as.matrix(syn_res.pro[syn_res.pro$p.value <= 0.05,1:length(exprCols)]), main=msg)


#
#
#	DEP 
#
#

# We can use any of the protqnt/qntS/qntV/qntV2 datasets 
# calculated above or just process the data from scratch
######

# For use with MSnData instead of a table:
# see https://bioconductor.org/packages/3.18/bioc/html/DEP.html
# convert to SummarizedExperiment
#se <- as(abun, "SummarizedExperiment")
# to convert back, use as(se, "MSnSet")

# default is to start from raw data
from.protein <- F
from.abundances <- F
from.sum.norm <- F
from.sum.vsn.norm <- F
from.vsn.norm <- F

#if (exprColsPrefix == "NormAbun") {
#    from.abundances <- T
#    # data is already Sum normalized, we will do additional VSN below,
#    # so we do not need Sum nor VNS nor Sum+VSN
#}
# we'll prefer the default (start from scratch) so we can keep meta-information

if (from.protein == T) {
    data_se <- as(protqnt, "SummarizedExperiment")
    exprs(data_se) <- log2(exprs(data_se))
} else if (from.abundances == T) {
    data_se <- as(qnt, "SummarizedExperiment")
    exprs(data_se) <- log2(exprs(data_se))
} else if (from.sum.norm == T) {
    data_se <- as(qntS, "SummarizedExperiment")
    exprs(data_se) <- log2(exprs(data_se))
} else if (from.sum.vsn.norm == T) {
    data_se <- as(qntV, "SummarizedExperiment")
    exprs(data_se) <- log2(exprs(data_se))
} else if (from.vsn.norm == T) {
    # preferable for already sum-normalized expression data
    data_se <- as(qntV2, "SummarizedExperiment")
    exprs(data_se) <- log2(exprs(data_se))
} else {		
    # work from scratch
    #
    # from data table
    #
    # https://bioconductor.org/packages/3.18/bioc/vignettes/DEP/inst/doc/DEP.html

    # read data
    data <- read.table('data.tab', header=T, sep='\t', quote='"')
    colnames(data)
    #
    # remove rows with NA
    # data <- na.omit(data)
    #

    # check if there are duplicate proteins/peptides
    data$accession %>% duplicated() %>% any()

    # make names unique and add a name and an ID column
    #data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
    data_unique <- make_unique(data, "accession", "Sequence", delim = ";")

    LFQ_columns <- grep(exprColsPrefix, colnames(data_unique)) # get LFQ column numbers

    # get experimental design (first four rows are Scaled Abundances)
    experimental_design <- read.table('metadata.tab', header=T)
    # fix labels to correct name
    experimental_design$label <- sub("Line", exprColsPrefix, experimental_design$label)
    experimental_design

    # copy "accession" to "ID" and "name" is we didn't use make_unique()
    #data_unique$name <- data_unique$accession
    #data_unique$ID <- data_unique$accession

    # convert to summarized experiment
    data_se <- make_se(data_unique, LFQ_columns, experimental_design)
    # assay data is log2 transformed and rownames correspond to protein names
    # name and ID should have been generated by make_unique
    # colData contains the experimental design in label, condition and replicate
    #	as well as a new ID column
    # IMPORTANT: it applies a log2 transform!!!
    # we do not need to add it later, but if we skip this, we need
    # log2 transformed expressions.

    # we can get the quants from an SE with assay(se) or assay(se)$counts 
    # the experimental design data with colData(se), metadata (if available) 
    # with metadata(se), regions of interest with rowRanges(se)

}

# at this point we have a summarizedExperiment with log2 transformed data
# at the peptide level
#
# if the data comes from previous *qnt* there are no NAs as we have already
# filtered them




########################################
# this should be more or less equivalent to the manual analysis that follows
#	filter based on missing values
#	apply VSN normalization
#	impute eventual remaining missing values
processed <- process(data_se)
pep.dep <- analyze_dep(processed, "control", "L260UP")

plot_all( pep.dep, plots = c("volcano", "heatmap", "single", "freq", "comparison"))
# if any one filas, the pdf() devices remain open and cannot be closed unless
#graphics.off()
# this is harmless unless we are also using pdf() ourselves
# for this reason and to avoid having too many pdf files, 
# we'll do it ourselves in parallel with the other analyses (manual and protein)
########################################
#
# now repeat at the protein level
# 
# combine all peptide data into protein data
data.ms <- as(data_se, "MSnSet")
prot.ms <- combineFeatures(data.ms,
                           groupBy = fData(data.ms)$accession,
                           method = sum)

# we need to sum-normalize first because otherwise process() fails doing VSN
prot.ms <- MSnbase::normalize(prot.ms, 'sum')
prot.dep <- as(prot.ms, "SummarizedExperiment")

prot.proc <- process(prot.dep)
prot.dep <- analyze_dep(processed, "control", "L260P")
plot_all( prot.dep, plots = c("volcano", "heatmap", "single", "freq", "comparison"))
#
# Now we have pep.dep and prot.dep
#
########################################
#
# Now we try to repeat it by hand to get better control
#
# let's have a look at the unprocessed data_se data

# filter missing values (items not identified in all replicas)
# first explore:
# how many proteins are in how many samples
sample.freq <- plot_frequency(data_se)
print(sample.freq)
DT::datatable(sample.freq$data)		# see numbers in browser

# we need to impute missing values, but only for proteins with not too many


# Filter for proteins that are identified in all replicates of at least 
# one condition
data_filt <- filter_missval(data_se, thr = 0)

# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at 
# least one condition (useless here as we only have 2 replicates each)
#data_filt2 <- filter_missval(data_se, thr = 1)

# Plot a barplot of the number of identified items per samples
pltn <- plot_numbers(data_filt)

pltn <- pltn + geom_text(aes(label = sprintf("%.1f", sum), y= sum),  
                                     vjust = 3) 

pltn  + scale_fill_manual("Sample", 
                                               values = c("#fbb4ae", 
                                                          "#b3cde3", 
                                                          "#ccebc5"))
print(pltn)

# Plot a barplot of the protein identification overlap between samples
pltc <- plot_coverage(data_filt)
pltc <- pltc + scale_fill_brewer(palette = "Paired")
pltc <- pltc + geom_text(aes(label=Freq), vjust=1.6, color="white")
#ggsave("./results/Protein coverage.pdf", width = 12, height = 6)
print(pltc)

# combine both plots
pltn / pltc


# The data is background corrected and normalized by variance stabilizing
# transformation (vsn).
# for Sum-normalized data (NormAbun) this will result in Sum+VSN
# Normalize the data
#	we use VSN for it is well behaved in most circumstances
data_norm <- normalize_vsn(data_filt)
pltm <- meanSdPlot(data_norm)

# if you want to customize you can access the plot as shown below
#plt.msd$gg + theme_bw() + scale_fill_distiller(palette = "RdPu")


# Visualize normalization by boxplots for all samples before and after normalization
cond_cols <- c(wt= "#e63946", mut="#2a9d8f")
pltN <- plot_normalization(data_filt, data_norm)
#pltN + scale_fill_manual(values = cond_cols)
print(pltN)
# there are no large differences


# Decide on imputations

# Before replacing NA-values it is important to verify that such values may be
# associated to absent or very low abundances. To do so, we suggest to inspect
# groups of replicate-measurements using matrixNAinspect(). In particular, with
# multiple technical replicates of the same sample it is supposed that any
# variability observed is not linked to the sample itself. So for each NA that
# occurs in the data we suggest to look what was reported for the same protein
# with the other (technical) replicates. This brings us to the term of
# ‘NA-neighbours’ (quantifications for the same protein in replicates).
# When drawing histograms of NA-neighbours one can visually inspect and verify
# that NA-neighbours are typically low abundance values, however, but not
# necessarily the lowest values observed in the entire data-set.
# 
# build a NA inpection plot to detect NA neighbors
matrixNAinspect(assay(data_norm), gr=gl(4,4))

# So only if the hypothesis of NA-neighbours as typically low abundance values
# gets confirmed by visual inspection of the histograms, one may safely proceed
# to replacing them by low random values.


# The remaining missing values in the dataset need to be imputed. The data can
# be missing at random (MAR), for example if proteins are quantified in some
# replicates but not in others. Data can also be missing not at random (MNAR),
# for example if proteins are not quantified in specific conditions (e.g. in
# the control samples). MNAR can indicate that proteins are below the detection
# limit in specific samples, which could be very well the case in proteomics
# experiments. For these different conditions, different imputation methods
# have to be used, as described in the MSnbase vignette and more specifically
# in the impute function descriptions.
# 
# To explore the pattern of missing values in the data, a heatmap is plotted
# indicating whether values are missing (0) or not (1). Only proteins with at
# least one missing value are visualized.

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# To check whether missing values are biased to lower intense proteins, the
# densities and cumulative fractions are plotted for proteins with and without
# missing values.

# Plot intensity distributions and cumulative fraction of proteins with and 
# without missing values
plot_detect(data_filt)

# In this case, missing values  seem to be biased to low-density: 
# it does look MNAR.

# Indeed the proteins with missing values have on average low intensities. This
# data (MNAR and close to the detection limit) should be imputed by a
# left-censored imputation method, such as the quantile regression-based
# left-censored function ("QRILC") or random draws from a
# left-shifted distribution ("MinProb" and "man"). In
# contrast, MAR data should be imputed with methods such as k-nearest neighbor
# ("knn") or maximum likelihood ("MLE") functions. See the
# MSnbase vignette and more specifically the impute function description for
# more information.

# If one uses a unique (very) low value for NA-replacements, this will quickly
# pose a problem at the level of t-tests to look for proteins changing
# abundance between two or more groups of samples. Therefore it is common
# practice to draw random values from a Normal distribution representing this
# lower end of abundance values. Nevertheless, the choice of the parameters of
# this Normal distribution is very delicate.

# All possible imputation methods are printed in an error, if an 
# invalid function name is given. The error will tell us which functions
# are available
#impute(data_norm, fun = "")

# if no imputation, then use filtered data
data_imp_none <- data_filt
# set NAs to zero (if not seen, assume it wasn't present, not that it
# has passed undetected because of its intensity was below equipment 
# sensitivity)
data_imp_zero <- DEP::impute(data_norm, fun = "zero")
# Impute missing data using random draws from a Gaussian distribution 
# centered around a minimal value (for MNAR)
data_imp_mp <- DEP::impute(data_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined 
# left-shifted Gaussian distribution (for MNAR)
data_imp_man <- DEP::impute(data_norm, fun = "man", 
                            shift = 1.8, 	# left shift of the distribution (in SD) from the median of the original distribution
		            scale = 0.3)	# width of the distribution relative to the SD of the original distribution
# note: shift=1 seems to work better

data_imp_qrilc <- DEP::impute(data_norm, fun = "QRILC")

# Impute missing data using the k-nearest neighbour approach (for MAR)
#data_imp_knn <- DEP::impute(data_norm, fun = "knn", rowmax = 0.9)

#data_imp_mle <- DEP::impute(data_norm, fun = "MLE")

# The effect of the imputation on the distributions can be visualized.
# Plot intensity distributions before and after imputation
#plot_imputation(data_norm, data_imp_none)
#plot_imputation(data_norm, data_imp_zero)
plot_imputation(data_norm, data_imp_mp)
plot_imputation(data_norm, data_imp_qrilc)
plot_imputation(data_norm, data_imp_man)
#plot_imputation(data_norm, data_imp_knn)
#plot_imputation(data_norm, data_imp_mle)
#
# can be customized to use previous colors (if done) with
# plti <- plot_imputation(...)
# plt + scale_color_manual(values=cond_cols)

# we will use 'man' imputation (data_imp_man)
data_imp <- data_imp_man


# show correlation matrix
cor_matrix <- plot_cor(data_imp, 
                       significant = F, 
                       lower = 0, 
                       upper = 1, 
                       pal = "GnBu",
                       indicate = c("condition", "replicate"), 
                       plot = F)
print(cor_matrix)
pheatmap(cor_matrix)


### Differential enrichment analysis
# 
# Protein-wise linear models combined with empirical Bayes statistics are used
# for the differential enrichment analysis (or differential expression
# analysis). The test_diff() function introduced here uses limma and
# automatically generates the contrasts to be tested. For the contrasts
# generation, the control sample has to be specified. Additionally, the types
# of contrasts to be produced need to be indicated, allowing the generation of
# all possible comparisons ("all") or the generation of contrasts of
# every sample versus control ("control"). Alternatively, the user
# can manually specify the contrasts to be tested (type = "manual"),
# which need to be specified in the argument test.

# Differential enrichment analysis  based on linear models and empherical 
# Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "L260P")

# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp, type = "all")

# Test manually defined comparisons
data_diff_manual <- test_diff(data_imp, type = "manual", 
                              test = c("L260P_vs_L260UP", "L260UP_vs_L260P"))


# Finally, significant proteins are defined by user-defined cutoffs using 
# add_rejections.

# Denote significant proteins based on user defined cutoffs
man.dep <- add_rejections(data_diff, alpha = p.cutoff, lfc = log2(log2cutoff))



# Visualization of the results
# 
# The results from the previous analysis can be easily visualized by a number
# of functions. These visualizations assist in the determination of the optimal
# cutoffs to be used, highlight the most interesting samples and contrasts, and
# pinpoint differentially enriched/expressed proteins.

# PCA plot
# 
# The PCA plot can be used to get a high-level overview of the data. This can
# be very useful to observe batch effects, such as clear differences between
# replicates.
# Plot the first and second principal components
pltp <- DEP::plot_pca(man.dep, x = 1, y = 2, n = dim(man.dep)[1], point_size = 4)
#pltp <- pltp +   ggtitle("WT/MUT") + 
#  scale_color_manual("Sample type", values = sub_cols)
print(pltp)


DEP::plot_pca(pep.dep, x = 1, y = 2, n=dim(pep.dep)[1], point_size = 4)
# there is another plot_pca in MSnSet.utils which won't work with an SE


# Correlation matrix
# 
# A correlation matrix can be plotted as a heatmap, to visualize the Pearson
# correlations between the different samples.
# Plot the Pearson correlation matrix
plot_cor(man.dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

plot_cor(pep.dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

# Heatmap of all significant proteins
# 
# The heatmap representation gives an overview of all significant proteins
# (rows) in all samples (columns). This allows to see general trends, for
# example if one sample or replicate is really different compared to the
# others. Additionally, the clustering of samples (columns) can indicate closer
# related samples and clustering of proteins (rows) indicates similarly
# behaving proteins. The proteins can be clustered by k-means clustering
# (kmeans argument) and the number of clusters can be defined by argument k.
# Plot a heatmap of all significant proteins with the data centered per 
# protein
plot_heatmap(man.dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"), 
             clustering_distance = "spearman")

plot_heatmap(pep.dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"), 
             clustering_distance = "spearman")

# Alternatively, a heatmap can be plotted using the contrasts, i.e. the direct
# sample comparisons, as columns. Here, this emphasises the enrichment of
# mutant compared to the control sample.

# Plot a heatmap of all significant proteins (rows) and the tested 
# contrasts (columns)
# This is of little use since here we only have one test against the 
# control
plot_heatmap(man.dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = TRUE, 
             clustering_distance = "spearman")
#	     ,
#                       indicate = c("condition", "replicate"), 
#                       show_row_dend= T,
#                       row_dend_side = "right", 
#                       width = 0.5, 
#                       gap = unit(1, "mm"))

plot_heatmap(pep.dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = TRUE, 
             clustering_distance = "spearman")


# Volcano plots of specific contrasts
# 
# Volcano plots can be used to visualize a specific contrast (comparison
# between two samples). This allows to inspect the enrichment of proteins
# between the two samples (x axis) and their corresponding adjusted p value (y
# axis). The add_names argument can be set to FALSE if the protein labels
# should be omitted, for example if there are too many names.
# Plot a volcano plot for the contrast "mut vs wt""

# this is to add if we do not want grids (plot_volcano(...) + remove_grids
# remove_grids <- theme(panel.grid.major = element_blank(), 
#                       panel.grid.minor = element_blank(),
#                       panel.background = element_blank(), 
#                       axis.line = element_line(colour = "black"))



DEP::plot_volcano(man.dep, contrast = "L260UP_vs_L260P", label_size = 4, 
	add_names = TRUE)
## Warning: ggrepel: 58 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps

DEP::plot_volcano(pep.dep, contrast = "L260UP_vs_L260P", label_size = 4, 
	add_names = TRUE)
# gives different (but very similar) results

# Barplots of a protein of interest
# 
# It can also be useful to plot the data of a single protein, for example if
# this protein is of special interest.
# # PSMs
# extreme <- c('Q9S746', 'Q93VG5', 'Q9LJN1',
# 	     'O04311', 'A0A178W0D3')
# # SclAbun
# extreme <- c('Q9LUD4', 'Q93VG5', 'F4KC80', 'P49107', 'O80653', 'Q8GYW0', 'A0A1P8B2M2',
# 	     'Q9SR37', 'O04310', 'O04311')
# # NormAbun
# extreme <- c('F4KC80', 'Q9LUD4', 'P49107', 'Q93VG5', 'O80653', 'Q8GYW0', 'A0A1P8B2M2',
# 	     'Q9SR37', 'O04310', 'O4311')
# general: get significant results IDs
man.extreme <- (get_results(man.dep) %>% filter(significant))$name

pep.extreme <- (get_results(pep.dep) %>% filter(significant))$name

# Plot a barplot for the extreme proteins
plot_single(man.dep, proteins = man.extreme)

plot_single(pep.dep, proteins = pep.extreme)

# Plot a barplot for the proteins with the data centered
plot_single(man.dep, proteins = man.extreme, type = "centered")

plot_single(pep.dep, proteins = pep.extreme, type = "centered")

# Frequency plot of significant proteins and overlap of conditions
# 
# Proteins can be differentially enriched/expressed in multiple comparisons. To
# visualize the distribution of significant conditions per protein and the
# overlap between conditions, the plot_cond function can be used.
#
# Plot a frequency plot of significant proteins for the different conditions
#plot_cond(man.dep)
#
#plot_cond(pep.dep)
# this is of little use as here we only have one condition


# weork in progress!!!
gene.annot <- F
if (gene.annot == T) {
    # no longer works
    # my_protein_ids <- unique(data$accession)
    # results <- POST(url = "https://www.uniprot.org/uploadlists/",
    #                 body = list(from = 'ID',
    #                     	to = 'GENENAME',
    #                     	format = 'tab',
    #                     	query = paste(my_protein_ids, collapse = ' ')))
    #
    # uniprot_results <- content(results, type = 'text/tab-separated-values', 
    #                            col_names = TRUE, 
    #                            col_types = NULL, 
    #                            encoding = "UTF-8")
    # this now should work
    # allToKeys(fromName = "UniProtKB_AC-ID")
    gene.names <- mapUniProt("UniProtKB_AC-ID", "Gene_Name", query = unique(data$accession))
    ### Warning message:
    ### IDs not mapped: P00761 
    gene.ids <- mapUniProt("UniProtKB_AC-ID", "GeneID", query = unique(data$accession))
    ### Warning message:
    ###IDs not mapped: Q8GWI5, A0A2P2CLF9, O77727, Q9LEX8, P02666, P00761 
    #
    # these produced a table with columns "From" and "To"
    #
    # we can get a map from KEGG with KEGGREST
    up2k <- keggConv("ath", "uniprot")
    length(up2k)
    # and make lookups using name 'up:accession'
    kegg_gene_names <- up2k[ paste("up", data$accession, sep=":") ]
    # to get a vector "ath:kegg_gene_name"

    # get mut_wt relevant comparison data
    mut_wt <- data_results[ , c("name","mut_vs_wt_significant", "mut_vs_wt_ratio", "mut_vs_wt_p.adj")]

    # keep only significant differences
    #foldchanges.1 <- subset(mut_wt, mut_vs_wt_p.adj <= 0.05)	# all columns
    #rownames(foldchanges.1) <- foldchanges.1$name

    foldchanges.1 <- mut_wt$mut_vs_wt_ratio		# only ratio
    names(foldchanges.1) = mut_wt$name		# we use names because f.1 is a vector

    foldchanges.1 <- foldchanges.1[ mut_wt$mut_vs_wt_p.adj <= 0.05 ]
    # we use a threshold of -1.2 or + 1.2 (>1.5)
    gene <- names(foldchanges.1)[abs(foldchanges.1) > 1.5]

    # KEGG enrichment
    # Identify KEGG pathways that are enriched.

    mut_wt_kegg <- enrichKEGG(sub('^ath:', '', kegg_gene_names), #gene,
                              organism = 'ath', 
                              pvalueCutoff = 0.05)
			      
    DT::datatable(as.data.frame(mut_wt_kegg))
    
    # Visualize enriched KEGG pathways.
    barplot(mut_wt_kegg, drop = F, showCategory = 12)

    enrichplot::cnetplot(,ut_wt_kegg,categorySize = "pvalue", 
                     foldChange = foldchanges.1, colorEdge= TRUE)
		     
    enrichplot::emapplot(mut_wt_kegg)
    
    heatplot(mut_wt_kegg, 
             foldChange=foldchanges.1,
             showCategory = 10) + ggtitle("Heatplot")

    mkk <- enrichMKEGG(gene = sub('^ath:', '', kegg_gene_names),	#gene,
                   organism = 'ath',
                   pvalueCutoff = 0.25,
                   minGSSize = 5,
                   qvalueCutoff = 0.25)

    # KEGG metabolic map

    # We can also view proteins that were identified in a given route on a
    # KEGG metabolic map. Note: When you run browseKEGG it will open a new window
    # with the default web-browser on your computer.
    #route <- 'xxxxxx'
    #browseKEGG(mut_wt_kegg, route)

}


# Results tables
# 
# To extract a table containing the essential results, the get_results function
# can be used.
# Generate a results table
man_results <- get_results(man.dep)
pep_results <- get_results(pep.dep)

# Number of significant proteins
man_results %>% filter(significant) %>% nrow()
pep_results %>% filter(significant) %>% nrow()
#colnames(man_results)

# Of these columns, the p.val and p.adj columns contain the raw and adjusted p
# values, respectively, for the contrast as depicted in the column name. The
# ratio columns contain the average log2 fold changes. The significant columns
# indicate whether the protein is differentially enriched/expressed, as defined
# by the chosen cutoffs. The centered columns contain the average log2 fold
# changes scaled by protein-wise centering.

# Generate a data.frame from the resulting SummarizedExperiment object
# 
# You might want to obtain an ordinary data.frame of the results. For this purpose, the package provides functions to convert SummarizedExperiment objects to data.frames. get_df_wide will generate a wide table, whereas get_df_long will generate a long table.

# Generate a wide data.frame
man_wide <- get_df_wide(man.dep)
pep_wide <- get_df_wide(pep.dep)
# Generate a long data.frame	(data arranged per sample) (4 x nrows here)
man_long <- get_df_long(man.dep)
pep_long <- get_df_wide(pep.dep)

# Save our results object for reuse
# 
# To facilitate future analysis and/or visualization of our current data,
# saving our analyzed data is highly recommended. We save the final data
# object (dep) as well as intermediates of the analysis, i.e. the initial
# SummarizedExperiment object (data_se), normalized data (data_norm), imputed
# data (data_imp) and differentially expression analyzed data (data_diff). This
# allows us to easily change parameters in future analysis.

# Save analyzed data
out <- paste("DEP", exprColsPrefix, "RData", sep='.')
save(data_se, data_norm, data_imp, data_diff, man.dep, pep.dep, file=out)
# These data can be loaded in future R sessions using this command
#load(out)

# save wide table
out <- paste("DEP.man", exprColsPrefix, "diff.tab", sep='.')
write.table( man_wide[man_wide$significant==T, ], 
	     file=out, row.names=T, col.names=T)

out <- paste("DEP.pep", exprColsPrefix, "diff.tab", sep='.')
write.table( pep_results[pep_results$significant==T, ], 
	     file=out, row.names=T, col.names=T)

# "simplify" manual results
out <- paste("DEP.man", exprColsPrefix, "diff.simple.tab", sep='.')
# get protein data
prot.data <- read.table('../proteins/data.tab', header=T, sep='\t', quote='"')
# add mut_vs_wt_ratio column from man_results
man_wide.a <- man_wide
man_wide.a$FC <- man_results$mut_vs_wt_ratio[ match(man_wide.a$name, man_results$name) ]
man_wide.a$Description <- prot.data$Description[ match(man_wide.a$accession, prot.data$accession) ]
# save
write.table( man_wide.a[man_wide.a$significant==T, 
	c( "mut_vs_wt_p.val", "mut_vs_wt_p.adj", "FC", 
	   "Sequence", "accession", "Description")], 
	file=out, row.names=F, col.names=T)

# expand multiprotein peptides for annotation
out <- paste("DEP.man", exprColsPrefix, "diff.simple.multi.tab", sep='.')
# annotate the dataset (it already has an 'accession' column
# find multi-accession rows
multis <- grepl(";", man_wide.a$accession)
# assign single-acc rows
#man_wide.a$Description <- prot.data$Description[ match(man_wide.a$accession, prot.data$accession) ]
# go over the dataset
for (i in 1:nrow(man_wide.a)) {
    # if it is a mult-acc row
    if ( multis[i] ) {
        #print(i)
        # copy it
        row <- man_wide.a[i, ]
	#print("original row")
	#print(row)
	## separate the accession ids
        accs <- strsplit(as.character(row$accession), "; ")
	# for each acc
	for (j in 1:length(accs[[1]])) {
	    #print(paste( "   ", j, accs[[1]][j] ) )
	    # set row$accession
	    row$accession <- accs[[1]][j]
	    # and add it to the dataframe at the end to avoid altering indexing
	    man_wide.a <- rbind(man_wide.a, row)
	}
	# now we could delete the old row but this will alter indexing
	#man_wide.a <- man_wide.a[ -i, ]
    }
}
# re-identify duplicates in the new data frame to get proper T/F values
multis <- grepl(";", man_wide.a$accession)
# remove multi-acc rows
man_wide.a <- man_wide.a[ ! multis, ]
# re-assign all single acc values
man_wide.a$Description <- prot.data$Description[ match(man_wide.a$accession, prot.data$accession) ]
# save
write.table( man_wide.a[man_wide.a$significant==T, 
	c( "mut_vs_wt_p.val", "mut_vs_wt_p.adj", "Sequence", "accession", "Description")], 
	file=out, row.names=F, col.names=T)

# repeat for pep
#

# "simplify" peptide results
out <- paste("DEP.pep", exprColsPrefix, "diff.simple.tab", sep='.')
pep_wide.a <- pep_wide
# get protein data
prot.data <- read.table('../proteins/data.tab', header=T, sep='\t', quote='"')
# add mut_vs_wt_ratio column from pep_results
pep_wide.a$FC <- pep_results$mut_vs_wt_ratio[ match(pep_wide.a$name, pep_results$name) ]
pep_wide.a$Description <- prot.data$Description[ match(pep_wide.a$accession, prot.data$accession) ]
# save
write.table( pep_wide.a[pep_wide.a$significant==T, 
	c( "mut_vs_wt_p.val", "mut_vs_wt_p.adj", "FC", "Sequence", "accession")], 
	file=out, row.names=F, col.names=T)

# expand multiprotein peptides for annotation
out <- paste("DEP.pep", exprColsPrefix, "diff.simple.multi.tab", sep='.')
# annotate the dataset (it already has an 'accession' column
# find multi-accession rows
multis <- grepl(";", pep_wide.a$accession)
# assign single-acc rows
#pep_wide.a$Description <- prot.data$Description[ match(pep_wide.a$accession, prot.data$accession) ]
# go over the dataset
for (i in 1:nrow(pep_wide.a)) {
    # if it is a mult-acc row
    if ( multis[i] ) {
        #print(i)
        # copy it
        row <- pep_wide.a[i, ]
	#print("original row")
	#print(row)
	## separate the accession ids
        accs <- strsplit(as.character(row$accession), "; ")
	# for each acc
	for (j in 1:length(accs[[1]])) {
	    #print(paste( "   ", j, accs[[1]][j] ) )
	    # set row$accession
	    row$accession <- accs[[1]][j]
	    # and add it to the dataframe at the end to avoid altering indexing
	    pep_wide.a <- rbind(pep_wide.a, row)
	}
	# now we could delete the old row but this will alter indexing
	#pep_wide.a <- pep_wide.a[ -i, ]
    }
}
# re-identify duplicates in the new data frame to get proper T/F values
multis <- grepl(";", pep_wide.a$accession)
# remove multi-acc rows
pep_wide.a <- pep_wide.a[ ! multis, ]
# assign all single acc values
pep_wide.a$Description <- prot.data$Description[ match(pep_wide.a$accession, prot.data$accession) ]
# save
write.table( pep_wide.a[pep_wide.a$significant==T, 
	c( "mut_vs_wt_p.val", "mut_vs_wt_p.adj", "Sequence", "accession", "Description")], 
	file=out, row.names=F, col.names=T)

###########################################################################
#
#	msmsTests using edgeR
#
#

# https://support.bioconductor.org/p/104999/
# we can also use package edgeR through "msmsTests"
# https://bioconductor.org/packages/release/bioc/manuals/msmsTests/man/msmsTests.pdfBiocManager::install("msmsTests")

# we start from the raw data as edgeR will do all the processing for us
pep.msms <- qnt		# base abundance data with NAs removed
pro.msms <- protqnt	# combined by protein

null.f <- "y~replicate"
alt.f <- "y~treatment+replicate"
pep.div <- apply(exprs(pep.msms), by.column, sum)
pro.div <- apply(exprs(pro.msms), by.column, sum)

poisson=F
if (poisson == T) {
    # Differential expression tests on spectral counts
    # 
    # Spectral counts (SpC) is an integer measure which requires of tests suited to
    # compare counts, as the GML methods based in the Poisson distribution, the
    # negative-binomial, or the quasilikelihood [11]
    # 
    # Generally speaking no test is superior to the other. They are just more or
    # less indicated in some cases. The Poisson regression requires the estimation
    # of just one parameter, and is indicated when the number of replicates is low,
    # two or three. A drawback of the Poisson distribution is that its variance
    # equals its mean, and it is not able to explain extra sources of variability a
    # part of the sampling. Quasi-likelihood is a distribution independent GLM, but
    # requires of the estimation of two parameters and hence needs a higher number
    # of replicates, i.e. not less than four. The negative-binomial requires two
    # parameters too, but we may use the implementation of this GLM in the edgeR
    # package [12] which uses an empirical Bayes method to share information across
    # features and may be employed with a restricted number of replicates; in the
    # worst case it limits with the Poisson solution.

    # Poisson GLM regression 
    # 
    # When using the Poisson distribution we implicitly accept a model not
    # sensitive to biological variability [11]. So it is just recommended in
    # cases where we have very few replicates, if any, and we do not expect a
    # signicant biological variability between samples.

    ### Remove all zero rows
    e <- pp.msms.data(pep.msms)
    dim(e)

    ### Null and alternative model
    #null.f <- "y~replica"
    #alt.f <- "y~status+replica"
    null.f <- "y~1"
    alt.f <- "y~condition"
    ### Normalizing condition
    div <- apply(exprs(e),by.column,sum)
    ### Poisson GLM
    pois.res <- msms.glm.pois(e,alt.f,null.f,div=div)
    str(pois.res)

    ### DEPs on unadjusted p-values
    sum(pois.res$p.value<=0.01)

    ### DEPs on multitest adjusted p-values
    adjp <- p.adjust(pois.res$p.value,method="BH")
    sum(adjp<=0.01)

    ### The top features
    o <- order(pois.res$p.value)
    head(pois.res[o,],20)
}

qlglm=F
if (qlglm == T) {
    # Quasi-likelihood GLM regression
    # 
    # The quasi-likelihood is a distribution free model that allows for overdispersion, and could
    # be indicated where an appreciable source of biological variability is expected. In this
    # model, instead of specifying a probability distribution for the data, we just provide a
    # relationship between mean and variance. This relationship takes the form of a function,
    # with a multiplicative factor known as the overdispersion, which has to be estimated from
    # the data [11]. Its use in proteomics has been documented by Li et al. (2010)[13].

    ### Quasi-likelihood GLM
    ql.res <- msms.glm.qlll(pep.msms,alt.f,null.f,div=div)
    str(ql.res)

    ### DEPs on unadjusted p-values
    sum(ql.res$p.value<=0.01)

    ### DEPs on multitest adjusted p-values
    adjp <- p.adjust(ql.res$p.value,method="BH")
    sum(adjp<=0.01)

    ### The top features
    o <- order(ql.res$p.value)
    head(ql.res[o,],20)

}


############################################################################
# edgeR: negative binomial GLM regression example
# 
# The negative-binomial provides another model that allows for overdispersion.
# The im- plementation adopted in this package is entirely based in the
# solution provided by the package edgeR [12] which includes empirical Bayes
# methods to share information among features, and thus may be employed even
# when the number of replicates is as low as two. The negative-binomial is
# downward limited, when no overdispersion is observed, by the Poisson
# distribution.

null.f <- y~replicate
alt.f <- y~treatment.1+replicate

eR.pep <- msms.edgeR(pep.msms,alt.f,null.f,div=pep.div,fnm="treatment.1")
str(eR.pep)
eR.pro <- msms.edgeR(pro.msms,alt.f,null.f,div=pro.div,fnm="treatment.1")
str(eR.pro)


### DEPs on unadjusted p-values
sum(eR.pep$p.value <= 0.01)	# see how many have p <= 0.01
sum(eR.pro$p.value <= 0.01)	# see how many have p <= 0.01

### DEPs on multitest adjusted p-values
adjp <- p.adjust(eR.pep$p.value,method="BH")
sum(adjp <= 0.01)
eR.pep$p.adjust <- adjp

adjp <- p.adjust(eR.pro$p.value,method="BH")
sum(adjp <= 0.01)
eR.pro$p.adjust <- adjp

### The top N features
o <- order(eR.pep$p.value)
head(eR.pep[o,], n.prots)

o <- order(eR.pro$p.value)
head(eR.pro[o,], n.prots)

# Let's now convert the data to counts if needed by multyplying by a factor:
#exprs(eR.pep) <- round(exprs(eR.pep) * 10)	# we got ScaledAbun with one decimal
#x <- msms.edgeR(eR.pep, alt.f, null.f, div = div, fnm = "status")
# not our case


# see fature data for significant differences
head( rownames(eR.pep[ eR.pep$p.value <= p.cutoff, ]) )	# significant diffs
head( rownames(fData(pep.msms)) %in%  rownames(eR.pep[ eR.pep$p.value < p.cutoff, ]) )
head( fData(pep.msms)[rownames(fData(pep.msms)) %in%  rownames(eR.pep[ eR.pep$p.value <= p.cutoff, ]), ] )

head( rownames(eR.pro[ eR.pro$p.value <= p.cutoff, ]) )	# significant diffs
head( rownames(fData(pro.msms)) %in%  rownames(eR.pro[ eR.pro$p.value < p.cutoff, ]) )
head( fData(pro.msms)[rownames(fData(pro.msms)) %in%  rownames(eR.pro[ eR.pro$p.value <= p.cutoff, ]), ] )

# Save 'pep' analysis
# save p <= p.cutoff
sigp <- fData(pep.msms)[rownames(fData(pep.msms)) %in%  rownames(eR.pep[ eR.pep$p.value <= p.cutoff, ]), ] 
# add statistics
sigp <- merge(sigp, eR.pep, by='row.names', all.x=T, all.y=F)[ , -1]
#rownames(sigp) <- sigp$accession
out <- paste("eR.pep.signif", exprColsPrefix, "p", p.cutoff, "tab", sep='.')
write.table( sigp, file=out, row.names=F, col.names=T)


# simplify 'pep' adding protein Description
out <- paste("eR.pep.signif", exprColsPrefix, "p", p.cutoff, "simple.tab", sep='.')
sigp.a <- sigp
# add description to single-protein peptides
# get protein data
prot.data <- read.table('../proteins/data.tab', header=T, sep='\t', quote='"')
sigp.a$Description <- prot.data$Description[ match(sigp.a$accession, prot.data$accession) ]
# save
write.table(sigp.a[ , c("p.value", "p.adjust", "LogFC",
			"Sequence", "Modifications", "accession", "Description")],
	   file=out, row.names=F, col.names=T)

# and now annotate multi-protein peptides
out <- paste("eR.pep.signif", exprColsPrefix, "p", p.cutoff, "simple.multi.tab", sep='.')
# find multi-accession rows
multis <- grepl(";", sigp.a$accession)
# assign single-acc rows
#sigp.a$Description <- prot.data$Description[ match(sigp.a$accession, prot.data$accession) ]
# go over the dataset
for (i in 1:nrow(sigp.a)) {
    # if it is a mult-acc row
    if ( multis[i] ) {
        #print(i)
        # copy it
        row <- sigp.a[i, ]
	#print("original row")
	#print(row)
	## separate the accession ids
        accs <- strsplit(as.character(row$accession), "; ")
	# for each acc
	for (j in 1:length(accs[[1]])) {
	    #print(paste( "   ", j, accs[[1]][j] ) )
	    # set row$accession
	    row$accession <- accs[[1]][j]
	    # and add it to the dataframe at the end to avoid altering indexing
	    sigp.a <- rbind(sigp.a, row)
	}
	# now we could delete the old row but this will alter indexing
	#sigp.a <- sigp.a[ -i, ]
    }
}
# re-identify duplicates in the new data frame
multis <- grepl(";", sigp.a$accession)
# remove multi-acc rows
sigp.a <- sigp.a[ ! multis, ]
# re-assign all single acc values
sigp.a$Description <- prot.data$Description[ match(sigp.a$accession, prot.data$accession) ]
# save
write.table(sigp.a[ , c("p.value", "p.adjust", "LogFC",
			"Sequence", "Modifications", "accession", "Description")],
	   file=out, row.names=F, col.names=T)



# Save 'pro' analysis
# add statistics
sigP <- fData(pro.msms)[rownames(fData(pro.msms)) %in%  rownames(eR.pro[ eR.pro$p.value <= p.cutoff, ]), ] 
sigP <- merge(sigP, eR.pro, by='row.names', all.x=T, all.y=F)[ , -1]
#rownames(sigp) <- sigp$accession
out <- paste("eR.pro.signif", exprColsPrefix, "p", p.cutoff, "tab", sep='.')
write.table( sigP, file=out, row.names=F, col.names=T)

# simplify 'pro' adding protein Description
out <- paste("eR.pro.signif", exprColsPrefix, "p", p.cutoff, "simple.tab", sep='.')
sigpP.a <- sigpP
# add description to single-protein peptides
# get protein data
prot.data <- read.table('../proteins/data.tab', header=T, sep='\t', quote='"')
sigpP.a$Description <- prot.data$Description[ match(sigpP.a$accession, prot.data$accession) ]
# save
write.table(sigpP.a[ , c("p.value", "p.adjust", "LogFC",
			"Sequence", "Modifications", "accession", "Description")],
	   file=out, row.names=F, col.names=T)

# and now annotate multi-protein peptides
out <- paste("eR.pro.signif", exprColsPrefix, "p", p.cutoff, "simple.multi.tab", sep='.')
# find multi-accession rows
multis <- grepl(";", sigpP.a$accession)
# assign single-acc rows
#sigpP.a$Description <- prot.data$Description[ match(sigpP.a$accession, prot.data$accession) ]
# go over the dataset
for (i in 1:nrow(sigpP.a)) {
    # if it is a mult-acc row
    if ( multis[i] ) {
        #print(i)
        # copy it
        row <- sigpP.a[i, ]
	#print("original row")
	#print(row)
	## separate the accession ids
        accs <- strsplit(as.character(row$accession), "; ")
	# for each acc
	for (j in 1:length(accs[[1]])) {
	    #print(paste( "   ", j, accs[[1]][j] ) )
	    # set row$accession
	    row$accession <- accs[[1]][j]
	    # and add it to the dataframe at the end to avoid altering indexing
	    sigpP.a <- rbind(sigpP.a, row)
	}
	# now we could delete the old row but this will alter indexing
	#sigpP.a <- sigpP.a[ -i, ]
    }
}
# re-identify duplicates in the new data frame
multis <- grepl(";", sigpP.a$accession)
# remove multi-acc rows
sigpP.a <- sigpP.a[ ! multis, ]
# re-assign all single acc values
sigpP.a$Description <- prot.data$Description[ match(sigpP.a$accession, prot.data$accession) ]
# save
write.table(sigpP.a[ , c("p.value", "p.adjust", "LogFC",
			"Sequence", "Modifications", "accession", "Description")],
	   file=out, row.names=F, col.names=T)


## PEP
# save padj <= 0.01
sigpa <- fData(pep.msms)[rownames(fData(pep.msms)) %in%  rownames(eR.pep[ eR.pep$p.adjust <= 0.01, ]), ] 
sigpa <- merge(sigpa, eR.pep, by='row.names', all.x=T, all.y=F)[ , -1]
#rownames(sigpa) <- sigpa$accession	# we won't save them

out <- paste("eR.pep.signif", exprColsPrefix, "padj=0.01.tab", sep='.')
write.table( sigpa, file=out, row.names=F, col.names=T)

# simplify
out <- paste("eR.pep.signif", exprColsPrefix, "padj=0.01.simple.tab", sep='.')
sigpa.a <- sigpa
# add description to single-protein peptides
# get protein data
prot.data <- read.table('../proteins/data.tab', header=T, sep='\t', quote='"')
sigpa.a$Description <- prot.data$Description[ match(sigpa.a$accession, prot.data$accession) ]
# save
write.table( sigpa.a[ , c("p.value", "p.adjust", "LogFC",
	   "accession", "Sequence", "Description")], 
	   file=out, row.names=F, col.names=T)

## PRO
sigpaP <- fData(pro.msms)[rownames(fData(pro.msms)) %in%  rownames(eR.pro[ eR.pep$p.adjust <= 0.01, ]), ] 
sigpaP <- merge(sigpaP, eR.pro, by='row.names', all.x=T, all.y=F)[ , -1]
#rownames(sigpaP) <- sigpaP$accession

out <- paste("eR.pro.signif", exprColsPrefix, "padj=0.01.tab", sep='.')
write.table( sigpaP, file=out, row.names=F, col.names=T)

# simplify
out <- paste("eR.pro.signif", exprColsPrefix, "padj=0.01.simple.tab", sep='.')
sigpaP.a <- sigpaP
# add description to single-protein peptides
# get protein data
prot.data <- read.table('../proteins/data.tab', header=T, sep='\t', quote='"')
sigpaP.a$Description <- prot.data$Description[ match(sigpaP.a$accession, prot.data$accession) ]
# save
write.table( sigpaP.a[ , c("p.value", "p.adjust", "LogFC",
	   "accession", "Sequence", "Description")], 
	   file=out, row.names=F, col.names=T)

#
# Reproducibility
# 
# In the omics eld, reproducibility is of biggest concern. Very low p-values
# for a protein in an experiment are not enough to declare that protein as of
# interest as a biomarker. A good biomarker should give as well a reproducible
# signal, and posses a biologically signicant eect size. According to our
# experience, a protein giving less than three counts in the most abundant
# condition results of poor reproducibility, to be declared as statistically
# signicant, experiment after experiment. On the other hand most of the false
# positives in an spiking experiment show log fold changes below 1. These two
# observations [6] allow to improve the results obtained in the previous
# sections. The trick is to flag as relevant those proteins which have low
# p-values, high enough signal, and good effect size. In performing this
# relevance filter we may even accept higher adjusted p-values than usual. This
# flagging is provided by the function test.results.

### Cut-off values for a relevant protein as biomarker
alpha.cut <- 0.05		# now we do not need to use 0.01 for p.adj
SpC.cut <- 2			# min spectral counts considered relevant
lFC.cut <- 1			# min abs log2 fold change
### Relevant peptides according to previous adjustments
eR.pep.tbl <- test.results(eR.pep,			# the dataframe with the stats
			pep.msms,		# the MSnSet with the data
                        pData(pep.msms)$treatment.1,	# the factor used in the tests
                        "UP",			# treatment level name
                        "P",			# contol level name
                        pep.div,		# weights used as divisors
			alpha=alpha.cut,minSpC=SpC.cut,minLFC=lFC.cut,
			method="BH")$tres

(eR.pep.nms <- rownames(eR.pep.tbl)[eR.pep.tbl$DEP])

### Relevant proteins according to previous adjustments
eR.pro.tbl <- test.results(eR.pro,			# the dataframe with the stats
			pro.msms,		# the MSnSet with the data
                        pData(pro.msms)$treatment.1,	# the factor used in the tests
                        "UP",			# treatment level name
                        "P",			# contol level name
                        pro.div,		# weights used as divisors
			alpha=alpha.cut,minSpC=SpC.cut,minLFC=lFC.cut,
			method="BH")$tres

(eR.pro.nms <- rownames(eR.pro.tbl)[eR.pro.tbl$DEP])

# without the post-test filter a relatively low adjusted p-value cut-off is
# required to keep an acceptable number of false positives. The post-test filter
# allows to relax the p-value cut-off improving at the same time both the number
# of true positives and false positives. 

# PEP
# save the results
eR.pep.wide <- merge(fData(pep.msms), eR.pep.tbl, by='row.names', all=F)
rownames(eR.pep.wide) <- eR.pep.wide$accession
eR.pep.wide <- eR.pep.wide[ , -1]

out <- paste("eR.pep.mark", exprColsPrefix, "p", alpha.cut, "tab", sep='.')
write.table(eR.pep.wide, file=out, row.names=F, col.names=T)

# simplify (adding protein names)
out <- paste("eR.pep.mark", exprColsPrefix, "p", alpha.cut, "simple.tab", sep='.')
eR.pep.wide.a <- eR.pep.wide
# get protein data
prot.data <- read.table('../proteins/data.tab', header=T, sep='\t', quote='"')
eR.pep.wide.a$Description <- prot.data$Description[ match(eR.pep.wide.a$accession, prot.data$accession) ]
# save
write.table(eR.pep.wide.a[ , c("p.value", "p.adjust", "LogFC",
	   "accession", "Sequence", "Description")],
	   file=out, row.names=F, col.names=T)

# PRO
eR.pro.wide <- merge(fData(pro.msms), eR.pro.tbl, by='row.names', all=F)
rownames(eR.pro.wide) <- eR.pro.wide$accession
eR.pro.wide <- eR.pro.wide[ , -1]

out <- paste("eR.pro.mark", exprColsPrefix, "p", alpha.cut, "tab", sep='.')
write.table(eR.pro.wide, file=out, row.names=F, col.names=T)

# simplify and add description
out <- paste("eR.pro.mark", exprColsPrefix, "p", alpha.cut, "simple.tab", sep='.')
eR.pro.wide.a <- eR.pro.wide
# get protein data
prot.data <- read.table('../proteins/data.tab', header=T, sep='\t', quote='"')
eR.pro.wide.a$Description <- prot.data$Description[ match(eR.pro.wide.a$accession, prot.data$accession) ]
write.table(eR.pro.wide.a[ , c("p.value", "p.adjust", "LogFC",
	   "accession", "Sequence", "Description")],
	   file=out, row.names=F, col.names=T)


# A useful tool to visualize the global results of dierential expression tests
# is a table of accumulated frequencies of features by p-values in bins of log
# fold changes. It may help in nding the most appropriate post-test lter cut-o
# values in a given experiment.

pval.by.fc(eR.pep.tbl$adjp,eR.pep.tbl$LogFC)
pval.by.fc(eR.pro.tbl$adjp,eR.pro.tbl$LogFC)

### Filtering by minimal signal
flt <- eR.pep.tbl$wt > 2
pval.by.fc(eR.pep.tbl$adjp[flt], eR.pep.tbl$LogFC[flt])

flt <- eR.pro.tbl$wt > 2
pval.by.fc(eR.pro.tbl$adjp[flt], eR.pro.tbl$LogFC[flt])

# Another usual tool is a volcanoplot with the ability to visualize the eect of
# dierent post-test lter cut-o values.
par(mar=c(5,4,0.5,2)+0.1)
res.volcanoplot(eR.pep.tbl,max.pval=0.05,min.LFC=1,maxx=3,maxy=NULL,
		ylbls=1.3)	# 1.3 = -log10(0.05), 2=-log10(0.01)

res.volcanoplot(eR.pro.tbl,max.pval=0.05,min.LFC=1,maxx=3,maxy=NULL,
		ylbls=1.3)	# 1.3 = -log10(0.05), 2=-log10(0.01)


#
#
# wrProteo
#
#

# https://cran.r-project.org/web/packages/wrProteo/vignettes/wrProteoVignette1.html

# PROTEIN LEVEL ANALYSIS

# wrProteo input functions return a list with $raw (initial/raw abundance
# values), $quant with ﬁnal normalized quantitations, $annot (columns ),
# $counts an array with ’PSM’ and ’NoOfRazorPeptides’, $quantNotes,
# $notes and optional setup for meta-data from sdrf; or a data.frame with
# quantitation and annotation if separateAnnot=FALSE

dataSV <- protqntV
# we could convert an MSnSet with
dataWP <- list(raw=exprs(abun),
	       quant=exprs(dataSV),
	       #quant=assay(data_se),
               annot=fData(dataSV),
               notes=list(),  # this should be further detailed
               quantNotes=list()      # ditto
	       )
rownames(dataWP$raw)=fData(abun)$Sequence
rownames(dataWP$quant)=fData(dataSV)$Sequence

# we can get the quants from an SE with assay(se) or assay(se)$counts 
# the experimental design data with colData(se), metadata (if available) 
# with metadata(se), regions of interest with rowRanges(se)

# the protocol is basically the same until normalization, but works for
# LFQ (label-free quantitation) data, not so well for PSM (spectral counting) 
# data, although it can import data from ProteomeDiscovered. But the data
# we got is not complete ProteomeDiscoverer output.

# As an alternative for imputation we can try matrixNAneighbourImpute from
# package wrProteo 

matrixNAinspect(dataWP$quant, gr=gl(2,2), tit="Histogram of Protein Abundances and NA-Neighbours")
matrixNAinspect(assay(data_norm), gr=gl(2,2), tit="Histogram of Protein Abundances and NA-Neighbours")

# This package proposes several related strategies/options for NA-imputation.
# First, the classical imputation of NA-values using Normal distributed random
# data is presented. The mean value for the Normal data can be taken from the
# median or mode of the NA-neighbour values, since (in case of technical
# replicetes) NA-neighbours tell us what these values might have been and thus
# we model a distribution around. Later in this vignette, a more elaborate
# version based on repeated implementations to obtain more robust results will
# be presented.
# 
# The function matrixNAneighbourImpute() proposed in this package offers
# automatic selection of these parameters, which have been tested in a number
# of different projects. However, this choice should be checked by critically
# inspecting the histograms of ‘NA-neighbours’ (ie successful quantitation
# in other replicate samples of the same protein) and the final resulting
# distribution. Initially all NA-neighbours are extracted. It is also worth
# mentioning that in the majority of data-sets encountered, such NA-neighbours
# form skewed distributions.
# 
# The successful quantitation of instances with more than one NA-values per
# group may be considered even more representative, but of course less
# sucessfully quntified values remain. Thus a primary choice is made: If the
# selection of (min) 2 NA-values per group has more than 300 values, this
# distribution will be used as base to model the distribution for drawing
# random values. In this case, by default the 0.18 quantile of the 2
# NA-neighbour distribution will be used as mean for the new Normal
# distribution used for NA-replacements. If the number of 2 NA-neighbours is >=
# 300, (by default) the 0.1 quantile all NA-neighbour values will used. Of
# course, the user has also the possibility to use custom choices for these
# parameters.
# 
# The final replacement is done on all NA values. This also includes proteins
# with are all NA in a given condition as well a instances of mixed successful
# quantitation and NA values.
# 
#data_wP_imp <- matrixNAneighbourImpute(assay(data_norm), gr=gl(2,2), tit="Histogram of Imputed and Final Data")
# (fails)
data_wP_imp <- matrixNAneighbourImpute(dataWP$quant, gr=gl(2,2), tit="Histogram of Imputed and Final Data")
# This function returns a list with $data .. matrix of data where NA are
# replaced by imputed values, $nNA .. number of NA by group, $randParam ..
# parameters used for making random data
data_wP_imp <- matrixNAneighbourImpute(exprs(qntV), gr=c(1,1,2,2), 
	tit="Histogram of Imputed and Final Data",
	plotHist=T)
data_wP_imp <- matrixNAneighbourImpute(assay(data_norm), gr=c(1,1,2,2), 
	tit="Histogram of Imputed and Final Data",
	plotHist=T)
data_wP_imp <- matrixNAneighbourImpute(dataWP$quant, gr=c(1,1,2,2), 
	tit="Histogram of Imputed and Final Data",
	plotHist=T)

#
# However, imputing using normal distributed random data also brings the risk
# of occasional extreme values. In the extreme case it may happen that a given
# protein is all NA in one group, and by chance the random values turn out be
# rather high. Then, the final group mean of imputed values may be even higher
# than the mean of another group with successful quantitations. Of course in
# this case it would be a bad interpretation to consider the protein in
# question upregulated in a sample where all values for this protein were NA.
# To circumvent this problem there are 2 options : 1) one may use special
# filtering schemes to exclude such constellations from final results or 2) one
# could repeat replacement of NA-values numerous times.
# 
# The function testRobustToNAimputation() allows such repeated replacement of
# NA-values. 
# 
# This function replaces NA values based on group neighbours (based on grouping
# of columns in argument gr), following overall assumption of close to Gaussian
# distribution. Furthermore, it is assumed that NA-values originate from
# experimental settings where measurements at or below de- tection limit are
# recoreded as NA. In such cases (eg in proteomics) it is current practice to
# replace NA-values by very low (random) values in order to be able to perform
# t-tests. However, random nor- mal values used for replacing may in rare cases
# deviate from the average (the ’assumed’ value) and in particular, if
# multiple NA replacements are above the average, may look like induced
# biological data and be misinterpreted as so. The statistical testing uses
# eBayes from Bioconductor package limma for robust testing in the context of
# small numbers of replicates. By repeating multiple times the process of
# replacing NA-values and subsequent testing the results can be sumarized
# afterwards by median over all repeated runs to remmove the stochastic effect
# of individual NA-imputation. Thus, one may gain stability towards
# random-character of NA imputations by repeating imputation & test ’nLoop’
# times and summarize p-values by median (results stabilized at 50-100 rounds).
# It is necessary to deﬁne all groups of replicates in gr to obtain all
# possible pair-wise testing (multiple columns in $BH, $lfdr etc). The
# modiﬁed testing-procedure of Bioconductor package ROTS may optionaly be
# included, if desired. This function returns a limma-like S3 list-object
# further enriched by additional ﬁelds/elements
# 
# test_imp <- testRobustToNAimputation(dataWP, 
#                                     gr=gl(2,2),
#				     plotHIst=T
#				     imputMethod="mode2",	# mode2 mode1 datQuant modeAdopt informed non
#				     nLoop=100)
#

# run the comparison tests
testWP <- testRobustToNAimputation(dataWP, 
				   gr=gl(2,2),
        			   plotHist=T,
        			   imputMethod="mode2",       # mode2 mode1 datQuant modeAdopt informed non
        		 	   nLoop=100)
# this will complain of non-unique row-names when used with qntV in dataWP

# head(data_wP_imp$datImp)
# This function returns a limma-type S3 object of class ’MArrayLM’ (which
# can be accessed lika a list); multiple results of testing or multiple testing
# correction types may get included (’p.value’,’FDR’,’BY’,’lfdr’ or ’ROTS.BH’)
# 

# If we had needed imputation, then the imputed data is in
# data_wP_imp$datImp
# then we could
plotPCAw(testWP$datImp, sampleGrp=gl(2,2), tit="PCA on Protein Abundances (MaxQuant,NAs imputed)", rowTyName="proteins", useSymb2=0)
# plot logFC versus average abundance
MAplotW(testWP)
#MAplotW(Mvalue, Avalue=NULL, useComp=1, filtFin=NULL, ProjNA=NULL, FCthrs=NULL,
#        subTxt=NULL, grayIncrem=TRUE, etc.)

# plot the second group of comparisons
#MAplotW(testWP, useComp=2, namesNBest="passFC")

# A Volcano-plot allows to compare the simple fold-change (FC) opposed to the
# outcome of a statistcal test. Frequently we can obsereve, that a some
# proteins show very small FC but enthousiastic p-values and subsequently
# enthousiastic FDR-values. However, generally such proteins with so small FC
# don’t get considered as reliable results, therefore it is common practice
# to add an additional FC-threshold, typically a 1.5 or 2 fold-change.

# FCthrs is the Fold-Change (not log2FC) threshold
# FdrThrs FDR threshold defaults to 0.05
statW <- VolcanoPlotW(testWP, 
        	      useComp=1, 
		      namesNBest="passThr", 	# name N best that past the thresholds
		      FCthrs=3, 		# FC
		      FdrThrs=0.05,		# FDR (FDR, BH, lfdr, BY
		      returnData=T		# df with ID, Mvalue, pValue, FDRvalue, passFilt
		      )

# do another one with FC threshold at 2 and FDR threshold at 0.01
VolcanoPlotW(testWP, useComp=1, namesNBest="passThr", FCthrs=2, FdrThrs=0.01)

# Tables with results can be either directed created using VolcanoPlotW() or,
# as shown below, using the function extractTestingResults().
resWP <- NULL
resWP <- extractTestingResults(testWP, compNo=1, thrsh=0.05, FCthrs=3)
#
# save to file the results for thrsg=0.01 and FCthrs=1.5
out <- paste("wP.pro", exprColsPrefix, "fdr=0.1.fc=1.5.simple.csv", sep='.')
resWP2 <- extractTestingResults(testWP, compNo=1, thrsh=0.01, FCthrs=1.5, 
	annotCol=c("Accession", "EntryName", "GeneName"), addTy=c("AllMeans"),
	filename=out
	)
	
# After FC-filtering for 2-fold (ie change of protein abundance to double or
# half)

# pretty print table using knitr
knitr::kable(resWP[,-1], caption="5%-FDR (BH) Significant results for 1st pairwise set", align="c")

# extend and save results
# rownames in resWP are accession values
# dataWP is our original data containing all relevant info (extras in $annot)
# we'll use it to retrieve the accession by Sequence
ann <- dataWP$annot
prot.data <- read.table('../proteins/data.tab', header=T, sep='\t', quote='"')
# 'match' returns a vector of the positions of (first) matches of
# its first argument in its second.
resWP$accession <- ann$accession[ match(rownames(resWP), ann$accession) ]
resWP$Sequence <- ann$Sequence[ match(rownames(resWP), ann$accession) ]
resWP$Description <- prot.data$Description[ match(resWP$accession, prot.data$accession) ]
out <- paste("wP.pro", exprColsPrefix, "fdr=0.05.fc=3.simple.tab", sep='.')
write.table(resWP[ , c("accession", "Sequence", "FDR", "logFC.1-2")],
	   file=out, row.names=F, col.names=T)

resWP2$accession <- ann$accession[ match(rownames(resWP2), ann$accession) ]
resWP2$Sequence <- ann$Sequence[ match(rownames(resWP2), ann$accession) ]
resWP2$Description <- prot.data$Description[ match(resWP2$accession, prot.data$accession) ]
out <- paste("wP.pro", exprColsPrefix, "fdr=0.1.fc=1.5.simple.tab", sep='.')
write.table(resWP2[ , c("accession", "Sequence", "FDR", "logFC.1-2")],
	   file=out, row.names=F, col.names=T)



# PEPTIDE LEVEL ANALYSIS

# wrProteo input functions return a list with $raw (initial/raw abundance
# values), $quant with ﬁnal normalized quantitations, $annot (columns ),
# $counts an array with ’PSM’ and ’NoOfRazorPeptides’, $quantNotes,
# $notes and optional setup for meta-data from sdrf; or a data.frame with
# quantitation and annotation if separateAnnot=FALSE

dataSV <- qntV2			# we use SVN but not SUM because orig data is
				# already sum-normalized
# we could convert an MSnSet with
dataWp <- list(raw=exprs(abun),
	       quant=exprs(dataSV),
	       #quant=assay(data_se),
               annot=fData(dataSV),
               notes=list(),  # this should be further detailed
               quantNotes=list()      # ditto
	       )
rownames(dataWp$raw)=fData(abun)$Sequence
rownames(dataWp$quant)=fData(dataSV)$Sequence

# we can get the quants from an SE with assay(se) or assay(se)$counts 
# the experimental design data with colData(se), metadata (if available) 
# with metadata(se), regions of interest with rowRanges(se)

# the protocol is basically the same until normalization, but works for
# LFQ (label-free quantitation) data, not so well for PSM (spectral counting) 
# data, although it can import data from ProteomeDiscovered. But the data
# we got is not complete ProteomeDiscoverer output.

# As an alternative for imputation we can try matrixNAneighbourImpute from
# package wrProteo 

matrixNAinspect(dataWp$quant, gr=gl(2,2), tit="Histogram of Peptide Abundances and NA-Neighbours")
# data_norm from DEP analysis
#matrixNAinspect(assay(data_norm), gr=gl(2,2), tit="Histogram of Protein Abundances and NA-Neighbours")

# The function testRobustToNAimputation() allows for repeated replacement of
# NA-values. 
# 
# run the comparison tests
testWp <- testRobustToNAimputation(dataWp, 
				   gr=gl(2,2),
        			   plotHist=T,
        			   imputMethod="mode2",       # mode2 mode1 datQuant modeAdopt informed non
        		 	   nLoop=100)
# this will complain of non-unique row-names when used with qntV in dataWp

# head(data_wP_imp$datImp)
# This function returns a limma-type S3 object of class ’MArrayLM’ (which
# can be accessed lika a list); multiple results of testing or multiple testing
# correction types may get included (’p.value’,’FDR’,’BY’,’lfdr’ or ’ROTS.BH’)
# 

# If we had needed imputation, then the imputed data is in
# data_wP_imp$datImp
# then we could
plotPCAw(testWp$datImp, sampleGrp=gl(2,2), tit="PCA on Protein Abundances (MaxQuant,NAs imputed)", rowTyName="proteins", useSymb2=0)
# plot logFC versus average abundance
MAplotW(testWp)
#MAplotW(Mvalue, Avalue=NULL, useComp=1, filtFin=NULL, ProjNA=NULL, FCthrs=NULL,
#        subTxt=NULL, grayIncrem=TRUE, etc.)

# plot the second group of comparisons
#MAplotW(testWp, useComp=2, namesNBest="passFC")

# A Volcano-plot allows to compare the simple fold-change (FC) opposed to the
# outcome of a statistcal test. Frequently we can obsereve, that a some
# proteins show very small FC but enthousiastic p-values and subsequently
# enthousiastic FDR-values. However, generally such proteins with so small FC
# don’t get considered as reliable results, therefore it is common practice
# to add an additional FC-threshold, typically a 1.5 or 2 fold-change.
layout(1)
# FCthrs is the Fold-Change (not log2FC) threshold
# FdrThrs FDR threshold defaults to 0.05
statWp <- VolcanoPlotW(testWp, 
        	      useComp=1, 
		      namesNBest="passThr", 	# name N best that past the thresholds
		      FCthrs=3, 		# FC
		      FdrThrs=0.05,		# FDR (FDR, BH, lfdr, BY
		      returnData=T		# df with ID, Mvalue, pValue, FDRvalue, passFilt
		      )

# do another one with FC threshold at 2 and FDR threshold at 0.01
VolcanoPlotW(testWp, useComp=1, namesNBest="passThr", FCthrs=2, FdrThrs=0.01)

# Tables with results can be either directed created using VolcanoPlotW() or,
# as shown below, using the function extractTestingResults().
resWp <- NULL
resWp <- extractTestingResults(testWp, compNo=1, thrsh=0.05, FCthrs=3)
#
# save to file the results for thrsg=0.01 and FCthrs=1.5
out <- paste("wP.pep", exprColsPrefix, "fdr=0.1.fc=1.5.simple.csv", sep='.')
resWp2 <- extractTestingResults(testWp, compNo=1, thrsh=0.01, FCthrs=1.5, 
	annotCol=c("Accession", "EntryName", "GeneName"), addTy=c("AllMeans"),
	filename=out
	)
	
# After FC-filtering for 2-fold (ie change of protein abundance to double or
# half)

# pretty print table using knitr
knitr::kable(resWp[,-1], caption="5%-FDR (BH) Significant results for 1st pairwise set", align="c")

# extend and save results
# rownames in resWp are accession values
# dataWp is our original data containing all relevant info (extras in $annot)
# we'll use it to retrieve the accession by Sequence
ann <- dataWp$annot
# 'match' returns a vector of the positions of (first) matches of
# its first argument in its second.
resWp$accession <- ann$accession[ match(rownames(resWp), ann$acc.no) ]
resWp$Sequence <- ann$Sequence[ match(rownames(resWp), ann$acc.no) ]
# get protein data
prot.data <- read.table('../proteins/data.tab', header=T, sep='\t', quote='"')
resWp$Description <- prot.data$Description[ match(resWp$accession, prot.data$accession) ]

out <- paste("wP.pep", exprColsPrefix, "fdr=0.05.fc=3.simple.tab", sep='.')
write.table(resWp[ , c("FDR", "logFC.1-2", "accession", "Sequence", "Description")],
	   file=out, row.names=F, col.names=T)

resWp2$accession <- ann$accession[ match(rownames(resWp2), ann$accession) ]
resWp2$Sequence <- ann$Sequence[ match(rownames(resWp2), ann$accession) ]
resWp2$Description <- prot.data$Description[ match(resWp2$accession, prot.data$accession) ]
out <- paste("wP.pep", exprColsPrefix, "fdr=0.1.fc=1.5.simple.tab", sep='.')
write.table(resWp2[ , c("logFC.1-2", "accession", "Sequence", "FDR", "Description")],
	   file=out, row.names=F, col.names=T)

dev.off()
graphics.off()

