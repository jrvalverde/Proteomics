perm.t.test <- function (data,format,B=1000){
  options(scipen=999)
  if (format=="long") {
    unstacked.data <- unstack(data) #requires 'plyr'
    sample1 <- unstacked.data[[1]]
    sample2 <- unstacked.data[[2]]
  } else {
  sample1 <- data[,1]
  sample2 <- data[,2]
  }
  #get some statistics for the two samples
  n1 <- length(sample1)
  n2 <- length(sample2)
  mean1 <- round(mean(sample1), 2)
  mean2 <- round(mean(sample2),2)
  error1 <- qnorm(0.975)*sd(sample1)/sqrt(n1)
  error2 <- qnorm(0.975)*sd(sample2)/sqrt(n2)
  sample1_lci <- round(mean1 - error1,2)
  sample1_uci <- round(mean1 + error1,2)
  sample2_lci <- round(mean2 - error2,2)
  sample2_uci <- round(mean2 + error2,2)
  #get regular t-test results (equal variance)
  p.equal.var <- round(t.test(sample1, sample2, var.equal=TRUE)$p.value, 4)
  #get regular t-test results (unequal variance)
  p.unequal.var <- round(t.test(sample1, sample2, var.equal=FALSE)$p.value, 4)
  #start permutation procedures
  pooledData <- c(sample1, sample2)
  size.sample1 <- length(sample1)
  size.sample2 <- length(sample2)
  size.pooled <- size.sample1+size.sample2
  nIter <- B
  meanDiff <- numeric(nIter+1)
  meanDiff[1] <- round(mean1 - mean2, digits=2)
  for(i in 2:length(meanDiff)){
    index <- sample(1:size.pooled, size=size.sample1, replace=F)
    sample1.perm <- pooledData[index]
    sample2.perm <- pooledData[-index]
    meanDiff[i] <- mean(sample1.perm) - mean(sample2.perm)
  }
  p.value <- round(mean(abs(meanDiff) >= abs(meanDiff[1])), digits=4)
  return (p.value)
  
  plot(density(meanDiff), main="Distribution of permuted mean differences", xlab="", sub=paste("sample 1 (n:", n1,") (95% CI lower bound., mean, 95% CI upper bound.):", sample1_lci, ",", mean1, ",", sample1_uci, "\nsample 2 (n:", n2,") (95% CI lower bound., mean, 95% CI upper bound.):", sample2_lci, ",", mean2, ",", sample2_uci,"\nobserved mean difference (dashed line):", meanDiff[1],"; permuted p.value (2-sided):", p.value, "( number of permutations:", B,")\nregular t-test p-values (2-sided):", p.equal.var, "(equal variance) ;", p.unequal.var, "(unequal variance)"), cex.sub=0.78)
  polygon(density(meanDiff), col="grey")
  rug(meanDiff, col="blue")
  abline(v=meanDiff[1], lty=2, col="red")
}
