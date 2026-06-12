
# Permutation test - two independent samples
# Notes by Valentin Stefan
# last update: 2021-10-14
# Data
# 
# This is an example about using a permutation test as replacement for a
# t-test. We assume two distributions (samples) that are independent and not
# paired. I used here the sleep data which is used in the examples that come
# with t.test() function. For metadata, see help(sleep).
# 
# It seems that the sleep data are actually paired, that is, each of the 10
# patients was “measured” two times. But we can imagine there were 20
# patients in total, 10 in control (group 1) in 10 in the “treatment” group
# (group 2).

head(sleep)

##   extra group ID
## 1   0.7     1  1
## 2  -1.6     1  2
## 3  -0.2     1  3
## 4  -1.2     1  4
## 5  -0.1     1  5
## 6   3.4     1  6

boxplot(extra ~ group, data = sleep)

# Run a permutation test
# 
# The usual t-test would go like this:

t.test(extra ~ group, data = sleep)

## 
##  Welch Two Sample t-test
## 
## data:  extra by group
## t = -1.8608, df = 17.776, p-value = 0.07939
## alternative hypothesis: true difference in means between group 1 and group 2 is not equal to 0
## 95 percent confidence interval:
##  -3.3654832  0.2054832
## sample estimates:
## mean in group 1 mean in group 2 
##            0.75            2.33

#The results indicate that the difference between the means of the two groups
# is statistically different from 0 (p = 0.079). Mean of group one is 0.75 and
# the mean of group 2 is 2.33. Their difference is 1.58.
# 
# Now we can look at the same question trough the perspective of a permutation
# test.
# 
# We will run n=1000 permutations to get an idea about the distribution of the
# differences of the means between the two groups and also the distribution of
# the t-values. Of course, you can try 10 000 or more, but 1 000 seems a
# generally accepted value.
# 
# By permuting the groups we sample under the null hypothesis. That is, we
# assume that data can come from either group as there should not be any
# difference between the two groups. So, we can too permute (shuffle) the group
# labels for each observation.
# 
# First, we define a vector (allocate memory for a vector) in which we will
# store the differences between the means of the “permuted” groups.
#
dif <- vector(length = 1000)

set.seed(2021-10-14) # for reproducibility when we use the function sample() below
dat <- sleep

for (i in 1:length(dif)){
  dat$group <- sample(dat$group) # shuffle the group labels
  # Compute the means for each group and then take the difference and store it
  dif[i] <- diff(tapply(X = dat$extra, 
                        INDEX = dat$group, 
                        FUN = mean))
}

# Looking at the distribution of the differences of the means between the two
# groups under the Null distribution:

# The observed difference of the original, not-permuted means
obs_dif <- diff(tapply(X = dat$extra, 
                       INDEX = dat$group, 
                       FUN = mean))
hist(dif, col = "grey")
# Add the observed difference at both tails
abline(v = -abs(obs_dif), col = "red", lwd = 2)
abline(v = abs(obs_dif), col = "red", lwd = 2)

# The approximate p-value that comes with the permutation test. Here, we need
# the two-tailed p-value, so we look at the both tails of the “dif”
# distribution.

mean(abs(dif) >= abs(obs_dif))

## [1] 0.061

# Same as above, but more verbose
sum(dif >= abs(obs_dif))/1000 + sum(dif <= -abs(obs_dif))/1000

## [1] 0.061

# The results of the permutation test are similar with the ones of the t-test
# mentioned above.
