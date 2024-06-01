# Text file containing answers for Assignment BIO724P/Statistics for Biologists
# Author? Usama Khan.   Date of Submission? 24/11/2023.
# Given as a single text file as per 0.1.1.


# Problem 1
# Task 1.1

# Defining our sample space as 'G'.
# Defining G as per visual diagram 0.1.2
G <- c("00-00", "00-01", "00-11", "01-00", "01-01", "01-11", "11-00", "11-01", "11-11")

# Alernate allele referenced as "1"
# Count over both sites

countAlternateAlleles <- function(genotype) {
    count <- 0
    for (allele in unlist(strsplit(genotype, ""))) {
        if (allele == "1") {
            count <- count + 1
        }
    }
    return(count)
}

# Applying our function iteratively to each element of G
alleleCounts <- sapply(G, countAlternateAlleles)

# Tabulating the allele frequencies
alleleFrequencies <- table(alleleCounts)

# Total number of genotypes in G?
totalGenotypes <- length(G)


# Calculating the Probability Mass Function (pmf)
# Defining pmf by a calculation of frequencies divided by genotypes
pmf <- alleleFrequencies / totalGenotypes


# Task 1.2. 

# Plotting the pmf

alleleCountValues <- names(pmf)  # x-values (counts of alternate alleles)
probabilityValues <- as.numeric(pmf)  # y-values (probabilities)

plot(
    x = alleleCountValues,
    y = probabilityValues,
    ylab = "Probability",
    xlab = "Count of Alternate Alleles",
    pch = 16,
    cex = 1.5,
    type = "h",
    main = "Probability Mass Function (PMF) of G"
)


# Plotting the cumulative ditribvutive function (CDF)

# Calculating the CDF
cdf_values <- cumsum(probabilityValues)

# Plotting our CDF.
plot(
    x = alleleCountValues,
    y = cdf_values,  # Corrected variable name
    ylab = "Cumulative Probability",
    xlab = "Count of Alternate Alleles",
    pch = 16,
    cex = 1.5,
    type = "s",
    main = "Cumulative Distribution Function (CDF) of G"
)

# End of Problem 1


# Start of Problem 2

# Task 2.1 

# As per diagram in 2.1, writing up genotypes for both sites.
genotypes <- data.frame(
  individual = 1:10,
  first_site = c("00", "00", "11", "01", "11", "00", "01", "00", "11", "11"),
  second_site = c("01", "01", "00", "11", "01", "00", "01", "11", "00", "11")
)

# Counting the alternative alleles
count_alternative_alleles <- function(genotype) {
  sum(unlist(strsplit(genotype, "")) == "1")
}

# Applying the function to each site
genotypes$first_site_count <- sapply(genotypes$first_site, count_alternative_alleles)
genotypes$second_site_count <- sapply(genotypes$second_site, count_alternative_alleles)

# Combining counts for our bar chart
counts <- rbind(genotypes$first_site_count, genotypes$second_site_count)

# Our visualisation of our distribution of G is a
# Barchart 
barplot(counts, beside = TRUE, 
        col = c("orange", "gray"), 
        main = "Alternative Alleles Count for Each Individual at Two Sites",
        xlab = "Individual No.", 
        ylab = "Count of Alternative Alleles", 
        names.arg = genotypes$individual,
        legend.text = c("Count at the First Site", "Count at the Second Site"),
        args.legend = list(title = "Site", cex = 0.8, bty = "n"))

# 2.1 metrics

#Before calculating metrics, a dataframe will be created as a reference point.

# Creating a data frame to store the genotype data for calculating metrics.
genotype_data <- data.frame(
    individual = 1:10,
    first_site = c("00", "00", "11", "01", "11", "00", "01", "00", "11", "11"),
    second_site = c("01", "01", "00", "11", "01", "00", "01", "11", "00", "11")
)

# Function to count alternate alleles in a genotype
count_alleles <- function(genotype) {
    sum(as.numeric(unlist(strsplit(genotype, ""))))
}

# Applying function to each row to calculate the count of alternate alleles
genotype_data$count_alternate_alleles <- apply(genotype_data[, 2:3], 1, function(x) {
    count_alleles(x[1]) + count_alleles(x[2])
})


# Calculating our central tendency as per 2.1
# Calculating the mean of 'G' (the count of alternate alleles)
the_mean_of_G <- mean(genotype_data$count_alternate_alleles)

# Print the_mean_of_G
print(paste("The mean (our measure of central tendency) count of alternate alleles (G) is:", the_mean_of_G))

# Calculating our metric of scale.
# Calculating the range of 'G' (the count of alternate alleles)
range_G <- range(genotype_data$count_alternate_alleles)

# Calculating the sample range (difference between max and min values)
sample_range_G <- diff(range_G)

# Printing our measure of scale.
print(paste("Our measure of scale. Minimum and Maximum of G:", range_G))

# Our metric of skewness
# Calculating mean and standard deviation of 'G'
mean_G <- mean(genotype_data$count_alternate_alleles)
sd_G <- sd(genotype_data$count_alternate_alleles)

# Number of observations
n <- length(genotype_data$count_alternate_alleles)

# Calculate the skewness
skewness_G <- (n / ((n - 1) * (n - 2))) * sum(((genotype_data$count_alternate_alleles - mean_G) / sd_G)^3)

# Print our skewness
print(paste("Skewness of G:", skewness_G))


# Task 2.2
# Comparing the mean of G with the mean of G2

# Hypotheses
# Our null Hypothesis (H0): The mean of G is equal to the mean of G2
# Our alternative Hypothesis (Ha): The mean of G is not equal to the mean of G2

# Given values of G2 as per assignment pdf in 2.2
G2 <- c(2, 4, 2, 3, 3, 3, 2, 3, 3, 4)

# Conducting an independent two-sample t-test
# This two-sample t-test compares the means of G and G2
t_test_result <- t.test(genotype_data$count_alternate_alleles, G2)

# Printing the results of the t-test
# The output includes the t-statistic, degrees of freedom, and p-value
print(t_test_result)

# Report
# Evaluating the difference between the mean count of alternate alleles 
# between G (new population) and G2 (reference population)
# G1 = (1, 1, 2, 3, 3, 0, 2, 2, 2, 4)
# G2 = (2, 4, 2, 3, 3, 3, 2, 3, 3, 4)

# Null Hypothesis (H0): There is no statistical difference between 
# the means of G and G2
# Alternative Hypothesis (Ha): There is a statistical difference 
# between the means of G and G2

# To reject H0, the p-value should be less than 0.05 (p < 0.05)

# Results from the two-sample t-test:
# t-value = -2.0769
# p-value = 0.05505
# Degrees of Freedom (df) = 15.3

# Interpretation:
# Since the p-value (0.05505) is greater than 0.05, we do not have  
# sufficient evidence to reject the null hypothesis.
# Conclusion:
# We continue to uphold the null hypothesis (H0) that there is no 
# statistically significant difference between the mean counts of 
# alternate alleles in G and G2.
# In the context of comparing alleles we may comment on any 
# differences being observed in alleles as being due to chance
# and that for biological relevance, in terms of the alleles,
# being studied- there may be similar genetic characteristics,
# between the two G and G2.

# End of Problem 2

# Problem 3
# Please ensure that the file "assignment.csv"
# is available in the same working directory!
#
if (!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)


df <- read.csv("assignment.csv")

#Making genotypes and ancestry variable into factors
df$genotypes <- as.factor(df$genotypes)
df$ancestry <- as.factor(df$ancestry)

model_genotypes_ancestry_income <- lm(risk ~ genotypes + ancestry + income, data = df)

summary(model_genotypes_ancestry_income)

model_2 <- lm(risk ~ genotypes + ancestry, data = df)
summary(model_2)

ggplot(df, aes(x = genotypes, y = risk)) + 
  geom_point(aes(colour = ancestry)) + 
  geom_smooth(aes(group = genotypes), method = "lm", se = FALSE, colour = "orange") + 
  geom_smooth(aes(group = ancestry, colour = ancestry), method = "lm", se = FALSE) + 
  labs(title = "A Scatter Plot with Regression Lines Showing the Relationship Between Genotype and Ancestry on Disease Risk",
       x = "Genotypes",
       y = "Disease Risk (Arbitrary Units)",
       colour = "Ancestry") + 
  theme_minimal() + 
  scale_colour_brewer(palette = "Set2") 


# Report Statement

# Analysing the data of 200 individuals using an initial General 
# Linear Model (GLM). The data included: Ancestry, genotypes and income.
# Genotypes 2, 3 and 4 had significant p values ((1.33e-05, 9.72e-13, < 2e-16 respectively).
# This indicates a strong association with risk.
# The variables ancestry 2 and ancestry 3 were also significant 
# (p = 0.00779 and 0.00546 respectively).
# The variable income had a p value of 0.90667
# Therefore we did not believe it was a significant predictor of risk
# We therefore had a reduced model which had excluded income to prevent
# confounding. 
# On visualising our reduced model we could see that not only do genotypes
# in general seem to have an association with disease risk but in addition
# as we move from genotype 0 to 4 there seems to be an upward trend in
# disease risk. The graph also seems to indicate that risk varies with 
# ancestry- the regression lines indicate a possible effect of genotype risk
# being influenced by ancestry. 

# Problem 4
 
# Reloading data from csv in case q4 ran as seperately
# Again please enusre assignment file is in the working directory
df <- read.csv('assignment.csv')

# Accessing the variable income
income <- df$income

# The data from income was visually inspected in excel
# the data was split into two distributions 
# rough estimates for the means of these distributions were
# found and these were used as initial parameters.
# To try and find point estimates for the parameters alpha,
# mu1 and mu2 we use Maximum Likelihood Estimate (MLE)

var1 <- 10^2  
var2 <- 10^2  

# Initial parameter estimations: alpha, mu1, and mu2
initial_params <- c(alpha = 0.5, mu1 = 23, mu2 = 43)
mixture_likelihood <- function(params, income) {
  alpha <- params[1]
  mu1 <- params[2]
  mu2 <- params[3]  
  likelihood <- alpha * dnorm(income, mu1, sqrt(var1)) + (1 - alpha) * dnorm(income, mu2, sqrt(var2))  
  -sum(log(likelihood))
}
result <- optim(initial_params, mixture_likelihood, income = income, method = "L-BFGS-B", 
                lower = c(0, -Inf, -Inf), upper = c(1, Inf, Inf))
if (result$convergence == 0) {
  print("Optimisation was successful.")
  print(result$par)  
} else {
  warning("Optimisation failed.")
}
# Now we have point estimates which have been printed off
# We are going to find confidence intervals using a bootstrapping method
# Our number of bootstrap iterations
num_iterations <- 1000
bootstrap_estimates <- matrix(nrow = num_iterations, ncol = 3)  
set.seed(123)  
for (i in 1:num_iterations) {
  sample_indices <- sample(length(income), replace = TRUE)
  bootstrap_sample <- income[sample_indices]
    
  bootstrap_result <- optim(c(alpha = 0.5, mu1 = 23, mu2 = 43), 
                            mixture_likelihood, 
                            income = bootstrap_sample, 
                            method = "L-BFGS-B", 
                            lower = c(0, -Inf, -Inf), 
                            upper = c(1, Inf, Inf))
  
  
  bootstrap_estimates[i, ] <- bootstrap_result$par
}
# Calculating our confidence intervals (CI):
alpha_ci <- quantile(bootstrap_estimates[, 1], probs = c(0.025, 0.975))
mu1_ci <- quantile(bootstrap_estimates[, 2], probs = c(0.025, 0.975))
mu2_ci <- quantile(bootstrap_estimates[, 3], probs = c(0.025, 0.975))
# Printing our CI's: 
print(list(alpha_ci = alpha_ci, mu1_ci = mu1_ci, mu2_ci = mu2_ci))

# End of Problem 4. Both point estimates and confidence intervals
# should have been printed off.
# End of assignment!Thank you for your patience.


