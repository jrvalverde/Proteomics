
# Generating data
set.seed(2021)
d <- as.data.frame(cbind(rnorm(1:20, 500, 50), c(rep(0, 10), rep(1, 10))))
treatment <- d$V2
outcome <- d$V1

#Difference in means
original <- diff(tapply(outcome, treatment, mean))
mean(outcome[treatment==1])-mean(outcome[treatment==0])

#Permutation test
permutation.test <- function(treatment, outcome, n){
  distribution=c()
  result=0
  for(i in 1:n){
    distribution[i]=diff(by(outcome, sample(treatment, length(treatment), FALSE), mean))
  }
  result=sum(abs(distribution) >= abs(original))/(n)
  return(list(result, distribution))
}

test1 <- permutation.test(treatment, outcome, 10000)
hist(test1[[2]], breaks=50, col='grey', main="Permutation Distribution", las=1, xlab='')
abline(v=original, lwd=3, col="red")

test1[[1]]


#Compare to t-test
t.test(outcome~treatment)



# Permutation hypothesis test in R
# 
# Exploring a powerful simulation technique with implementation from scratch in
# R
# 
# Serafim Petrov
# Towards Data Science
# Serafim Petrov
# 
# 4 min read
# ·
# Mar 15, 2021
# Photo by Eric Prouzet on Unsplash
# 
# Introduction
# 
# To compare outcomes in experiments, we often use Student’s t-test. It
# assumes that data are randomly selected from the population, arrived in large
# samples (>30), or normally distributed with equal variances between groups.
# If we do not happen to meet these assumptions, we may use one of the
# simulation tests. In this article, we will introduce the Permutation Test.
# Rather than assuming underlying distribution, the permutation test builds its
# distribution, breaking up the associations between or among groups. Often we
# are interested in the difference of means or medians between the groups, and
# the null hypothesis is that there is no difference there. We may ask the
# question: from all the possible permutations, how extreme our data would look
# like? All possible permutations would represent a theoretical distribution.
# In practice, there is no need to perform ALL permutations to build the
# theoretical distribution, but run a reasonable number of simulations to take
# a sample from that distribution. Usually, there are 10k or 100k simulations.
# 
# Example
# 
# Suppose we have a chain of 12 retail stores, and we want to test a new sales
# technique. We may assign or just pick, for instance, 5 stores and try it
# there. Then compare average sales after a certain period with sales in the
# rest 7 stores.
# 
# The data are not randomly collected, presumably not normally distributed; and
# the sample itself has fewer than 30 observations.
# 
# In our dataset, we have two columns: treatment — a binary variable
# indicating whether the store implemented the new approach; outcome — a
# numeric variable recoding a sales amount at the end of a period.
# 
# From the boxplot representation, we may see that the groups are fairly
# different and the variance is approximately equal.
# 
# The histogram of outcome does not appear normally distributed; however, the
# Shapiro-Wilk normality test returns a not significant p-value = 0.3654
# To find the difference in means between groups, we use built-in commands and
# save it to a variable original:
# 
# original <- diff(tapply(outcome, treatment, mean))
# 
# The difference is ~37.8 points. Now we ask the question: if we shuffle the
# data ten thousand times, how often we would observe this or more extreme
# difference?
# 
# We run a simulation drawing 10k shuffles without replacement, and record the
# difference in means for each permutation. Then we perform a two-sided test
# computing the number of records the absolute value of whose exceeded absolute
# values of original difference.
# 
# The function returns the simulated distribution and the p-value.
# If we plot the distribution, we may observe that our original difference is
# not particularly extreme, with an exact p-value of 0.1822818
# 
# In many research projects, this would indicate that there is no statistically
# significant difference between groups: however, in this business case, I
# would pick the option that there is some evidence that the sales approach
# affects the outcome.
# 
# Comparison with other tests
# 
# If we run the Welch t-test, we would get a similar p-value of 0.1815.
# The Shapiro-Wilk normality test indicated that the data were presumably
# obtained from the Normal Distribution, but it is not a general case.
# R also has several libraries to run permutation tests. One of them is
# library(coin), which performs an independence test. It also returns a similar
# p-value of 0.1858.
# 
# Conclusion
# 
# The Permutation test is a powerful tool in measuring effects in experiments.
# It is easy to implement, and it does not rely on many assumptions as other
# tests do. It has not been widely popular until the simulation on computers
# became routinely implemented.
# 
# It is closely related to the other simulation test: the bootstrapping
# hypothesis test, where the samples are drawn with replacement.
# 
# 
