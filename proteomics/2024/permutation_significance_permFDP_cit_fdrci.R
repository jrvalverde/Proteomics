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

####
####
####

cran_packages <- c(
	"coin",
	"dplyr",
	"RVAideMemoire",
	'knitr',
	'fdrci'
	)
	
for (pkg_i in cran_packages) {
  require.cran(pkg_i)
}

library.github('steven-shuken/permFDP')

####
####
####


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
       cat('p.perm.test: returning p=1 to signal invalid comparison\n')
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

# compute p for all rows using RVAideMemoire::perm.t.test
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
    perml = vector('list', n.permutations)
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
	
	    if (verbose)
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
    } else {    # nspv <= 100
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

