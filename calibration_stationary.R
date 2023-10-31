library(tidyverse)
library(data.table)
library(nloptr)
library(mvtnorm)
library(modeest)

# Let's try this again, but this time in another way. 
#
# More specifically, I'll try to fix my previous issues in the following way:
#   - First, estimate a global measure of error variance across all tags
#   - Then, filter the positions so that you have a fixed measured position for 
#     each tag
#   - Finally, estimate a function between the measured positions and the actual
#     positions

# Read in the data
stationary <- readRDS(file.path("Calibration experiments", 
                                "Data",
                                "preprocessed_stationary.Rds")) %>% 
    ungroup()


# Check again whether the measured positions accurately correspond to a given tag
ggplot(data = stationary, 
       aes(x = x, y = y, color = factor(tag))) +
    geom_point()





################################################################################

######################################
# Task 1: Estimating the error variance

# Pretreat the data to include the mean positions for each of the tags into the 
# data itself as well. The error variance will be computed relative to this mean,
# and thus comes with its undeniable assumptions.
#
# Alternatives can be the mode or median, but needs to be investigated 
# somewhat more. These values are, however, also included in this dataframe.
stationary <- stationary %>% 
    group_by(tag) %>% 
    mutate(mu_x = mean(x),
           mu_y = mean(y),
           median_x = median(x), 
           median_y = median(y), 
           mode_x = mean(mlv(x, method = "mfv")), 
           mode_y = mean(mlv(y, method = "mfv"))) %>% 
    ungroup() %>% 
    mutate(corrected_x = x - mu_x, 
           corrected_y = y - mu_y)

# Check how far off the means are versus the medians and modes. This might 
# give an indication of how justified we are in using the mean instead of the
# other measures
#
# Change x in the ggplot `aes` argument as you please
#
# As a small summary: 
#   - mu_x - median_x: Seems like the mean is slightly bigger than the median 
#                      overall. The difference is, however, very small (most
#                      contained within [-0.001, 0.001], or within 2mm)
#   - mu_y - median_y: This shows the same pattern, but has some outliers that
#                      are quite big (around 0.01, or about 1cm)
#   - mu_x - mode_x: Aside from some outliers, the distribution is quite 
#                    symmetrical. However, the deviation is much bigger 
#                    ([-0.01, 0.01], or within 2cm)
#   - mu_y - mode_y: Same pattern as for x.
stationary %>% 
    mutate(diff_x = mu_x - median_x, 
           diff_y = mu_y - median_y) %>% 
    select(diff_x, diff_y) %>% 
    ggplot(aes(x = diff_y)) +
        geom_histogram(fill = "cornflowerblue", 
                       color = "black", 
                       bins = 20) +
        geom_vline(xintercept = 0, 
                   color = "red",
                   linewidth = 2)

stationary %>% 
    mutate(diff_x = mu_x - mode_x, 
           diff_y = mu_y - mode_y) %>% 
    select(diff_x, diff_y) %>% 
    ggplot(aes(x = diff_y)) +
        geom_histogram(fill = "cornflowerblue", 
                       color = "black", 
                       bins = 20) +
        geom_vline(xintercept = 0, 
                   color = "red",
                   linewidth = 2)

# First attempt: Parametric ######################

# Some assumptions:
#   - The actual measured position is equal to the average of all measured 
#     positions for a given tag
#   - The error variance is equal for each of the tags
#   - This error is normally distributed
#   - Independence of the error between the X and Y dimensions
#   - Independence of error over time

# Create a function that will compute and evaluate the min-log-likelihood
min_log_likelihood <- function(x, mu, sigma){
    # Transpose: `ldmvnorm` assumes that you have one observation in each column
    x <- as.matrix(x) %>% 
        t()
    mu <- as.matrix(mu) %>% 
        t()

    # Transform `sigma` to a `lpMatrix`, as the function `ldmvnorm` requires
    sigma <- ltMatrices(sigma[c(1:2,4)], diag = TRUE)

    # Do the usual: `ldmvnorm` computes the log-likelihood of the given 
    # observations `x` under the means `mu` and Cholesky triangular matrix
    # `sigma`. Under the assumption of independence of the two dimensions
    # (i.e., when the off-diagonal elements in `sigma` are 0), the Cholesky
    # is simply a diagonal matrix with standard deviations as its diagonal 
    # elements
    mvtnorm::ldmvnorm(x, mean = mu, chol = sigma) %>%
        sum() %>% 
        `*` (-1) %>% 
        return()
}

# Create a function that will be minimized in which the variances will be 
# estimated
decomposition_of_variance <- function(x, mu_columns = c("mu_x", "mu_y")){
    # Get the overall standard deviation `sigma`, which will be the measure of 
    # error around the measurements
    sigma <- c(x[1], 0, 0, x[1]) %>% 
        matrix(ncol = 2)

    # Define the min-log-likelihood variable
    MLL <- 0
    
    # Evaluate the min-log-likelihood of the measured positions given the 
    # estimated positions
    min_log_likelihood(stationary[,c("x", "y")],
                       stationary[, mu_columns],
                       sigma) %>% 
        return()
}

# Run the estimation using nloptr
n_tags <- stationary$tag %>% 
    unique() %>% 
    length()
set.seed(35356) # New Storm for Older Lovers - La Dispute
result <- nloptr(x0 = 1,
                 eval_f = function(x) decomposition_of_variance(x),
                 lb = 10^(-15),
                 ub = 1000,
                 opts = list("algorithm" = "NLOPT_LN_BOBYQA", 
                             "ftol_rel" = 10^(-15),
                             "xtol_rel" = 10^(-15),
                             "maxeval" = 10^4))

# Save the results of this estimation
saveRDS(result,
        file.path("Calibration experiments", "Results", "estimate_error_parametric.Rds"))

# Using this result, check how many bins you would need for each dimension to 
# be fairly certain that the position is within this bin.
# 
# If we would bin within the 5 cm, we would be fairly certain the real 
# position would be within this specified bin
binsize = 2 * 1.96 * result$solution

# Time to test some assumptions.
#
# Let's start with normality.
#
# Hump all data of the tags together by mean-correcting them. Then do several 
# tests:
#   - Visual test: Make a histogram and look at the distribution
#   - Visual test: Make a QQ plot and look at its result
#   - Formal test: Shapiro Wilk test
stationary <- stationary %>% 
    mutate(corrected_x = x - mu_x, 
           corrected_y = y - mu_y)

# Let make the histogram: Look quite normal
stationary %>% 
    ggplot(aes(x = corrected_x)) + geom_histogram()
stationary %>% 
    ggplot(aes(x = corrected_y)) + geom_histogram()

# And the QQ-plot: Is quite okay, except at the edges, where it shows large 
# deviations. This is true for both x and y
#
# Solution, use a t-distribution instead? Stays the case apparently.
qqnorm(stationary$corrected_x)
qqline(stationary$corrected_x)

qqnorm(stationary$corrected_y)
qqline(stationary$corrected_y)

car::qqPlot(stationary$corrected_x, distribution = "t", df = length(stationary$x) - 1)
car::qqPlot(stationary$corrected_y, distribution = "t", df = length(stationary$y) - 1)

# And finally, the Shapiro-Wilk test of normality
nortest::ad.test(stationary$corrected_x)
nortest::ad.test(stationary$corrected_y)

ks.test(jitter(stationary$corrected_x), "pnorm")
ks.test(jitter(stationary$corrected_y), "pnorm")

# Now, let's test whether x and y are indeed independent.
#
# This can be done by just computing the correlation between both measured 
# dimensions.
#
# They are significantly correlated
cor.test(stationary$x, stationary$y)
cor.test(stationary$corrected_x, stationary$corrected_y)

# Second attempt: Nonparametric ######################

# Given that the assumptions are not correct, I will need another approach to 
# estimating the error variances. Instead, we will have to compute an error
# covariance matrix, preferably in a nonparametric way. This is what we do in 
# the second attempt.
#
# Some assumptions:
#   - The actual measured position is equal to the average of all measured 
#     positions for a given tag
#   - The error variance is equal for each of the tags
#   - Independence of error over time

# Create a function that will create bootstrapped data in an efficient way. 
# Assumption here is that x is a dataframe that contains, among other, the 
# x and y coordinates that we want to bootstrap.
#
# Importantly, it only bootstraps rows: The relationship between x and y remains
# untouched by this function. This will allow us to not only estimate the 
# variance components for both dimensions, but also the covariance between them.
bootstrap_data <- function(x, iterations) {
    # Get the sample size of the data. Needed to ensure that each of the samples
    # has an equal size to the actual data
    N <- nrow(x)

    # Sample a number of indices for x that is equal to the sample size times 
    # the number of samples one wants to draw
    idx <- sample(1:N, 
                  N * iterations, 
                  replace = TRUE)

    # Extend the dataframe to account for these values and bind them with an 
    # identity number that conveys the sample they are in
    x[idx,] %>% 
        mutate(sample_id = rep(1:iterations, each = N)) %>% 
        return()
}

# Bootstrap the data using this function and immediately compute the necessary
# summary statistics: 2 variances and 1 covariance.
#
# Be careful, this uses a lot of memory!
covariances <- bootstrap_data(stationary, 1000) %>% 
    # Compute the statistics
    group_by(sample_id) %>% 
    mutate(var_x = var(corrected_x), 
           var_y = var(corrected_y), 
           cov_xy = cov(corrected_x, corrected_y)) %>% 
    ungroup() %>% 
    # Delete all other information: Only keep variances and covariance
    group_by(sample_id, var_x, var_y, cov_xy) %>% 
    tidyr::nest() %>% 
    select(-data) %>% 
    ungroup() 

# We can then inspect the histograms of the different statistics under 
# investigation
booted_histogram <- function(col){
    plot_data <- covariances[, col] %>% 
        as.data.frame() %>% 
        setNames(c("X"))

    (ggplot(plot_data, aes(x = X)) +
        geom_histogram(col = "black", 
                       fill = "cornflowerblue")) %>% 
        return()
}

booted_histogram("var_x")
booted_histogram("var_y")
booted_histogram("cov_xy")

# Create several covariance matrices, them being the 99%CI and the mean
covariances <- list("lb" = matrix(c(quantile(covariances$var_x, 0.005), 
                                    quantile(covariances$cov_xy, 0.005),
                                    quantile(covariances$cov_xy, 0.005),
                                    quantile(covariances$var_y, 0.005)),
                                  nrow = 2),
                    "mean" = matrix(c(mean(covariances$var_x), 
                                      mean(covariances$cov_xy),
                                      mean(covariances$cov_xy),
                                      mean(covariances$var_y)),
                                    nrow = 2),
                    "ub" = matrix(c(quantile(covariances$var_x, 0.995), 
                                    quantile(covariances$cov_xy, 0.995),
                                    quantile(covariances$cov_xy, 0.995),
                                    quantile(covariances$var_y, 0.995)),
                                  nrow = 2))

# Save the results of this estimation procedure
saveRDS(covariances,
        file.path("Calibration experiments", "Results", "estimate_error_nonparametric.Rds"))

# And get everything into interpretable quantities.
max_covariance <- covariances[["ub"]]
print(paste0("The standard deviation for the X dimensions is ", 
             max_covariance[1,1] %>% sqrt() %>% round(4), 
             "m, or a total of ",
             max_covariance[1,1] %>% sqrt() %>% round(4) %>% `*` (100), 
             "cm."))
print(paste0("The standard deviation for the Y dimensions is ", 
             max_covariance[2,2] %>% sqrt() %>% round(4), 
             "m, or a total of ",
             max_covariance[2,2] %>% sqrt() %>% round(4) %>% `*` (100), 
             "cm."))
print(paste0("The correlation between the two dimensions is ", 
             max_covariance[1,2] %>% 
                `/` (sqrt(max_covariance[1,1]) * sqrt(max_covariance[2,2])) %>% 
                round(4)))





################################################################################

######################################
# Task 2: Checking the drop-out rate

# Small utility function to extract a piece of a string and convert it to 
# numeric
subnumeric <- function(x, first, last){
    substring(x, first, last) %>% 
        as.numeric() %>% 
        return()
}

# Convert the timestamps of the data into an integer that indicates the number 
# of milliseconds that have passed since the beginning. This is done separately 
# for each of the different tags, as drop-out should also be assessed for each
# one separately.
convert_to_millisecond <- function(x){
    # Only get hours, minutes, seconds, and milliseconds
    times <- x %>% 
        as.POSIXlt(format = "%H:%M:%OS") %>% 
        format("%H:%M:%OS3")

    # Unfortunately, need a dirty hack to be able to work with this further.
    # Assume this format and define hours, minutes, seconds, and milliseconds.
    times %>% 
        (function(x){
            return(3600000 * subnumeric(x, 1, 2) +
                   60000 * subnumeric(x, 4, 5) +
                   1000 * subnumeric(x, 7, 8) +
                   subnumeric(x, 10, 11))
        }) %>% 
        return()
}

# Create a function that will take in the duration, order them according to 
# size, and then compute the mean difference between each of the durations
difference_sampling_times <- function(x){
    x %>% 
        as.numeric() %>% 
        sort() %>% 
        diff() %>% 
        mean() %>% 
        return()
}

# Do this conversion for each of the tags and create a duration column that 
# starts out at 0
stationary <- stationary %>% 
    mutate(duration = convert_to_millisecond(timestamp)) %>% 
    group_by(tag) %>% 
    mutate(duration = duration - min(duration)) %>% 
    ungroup()

# Get the mean difference in the sampled duration between the tags
result <- stationary %>% 
    group_by(tag) %>% 
    mutate(mean_time_between_samples = difference_sampling_times(duration)) %>% 
    ungroup() %>% 
    group_by(tag, mean_time_between_samples) %>% 
    tidyr::nest() %>% 
    select(-data)

# Save these results for future reference
saveRDS(result,
        file.path("Calibration experiments", "Results", "sampling_rate.Rds"))

# And get the mean sampling rate across all tags
mean(result$mean_time_between_samples)





################################################################################

######################################
# Task 3: Defining a function to correct for systematic noise

# Some assumptions:
#   - The actual measured position is equal to the average of all measured 
#     positions for a given tag

########################
# Estimation: 1D

stationary <- readRDS(file.path("Calibration experiments", 
                                "Data", 
                                "preprocessed_stationary.Rds"))

# Create a function that will output a formula for a polynomial of the pwr^th
# degree with y as dependent variable and x as the independent variable
polynomial_formula <- function(x, y, pwr){
    polyn <- paste(y, "~", x)
    for(i in 2:pwr){
        polyn <- paste(polyn, "+", "I(",
                       paste0(x, "^", i), ")")
    }
    polyn %>% 
        formula() %>% 
        return()
}

# Estimate the 10th degree polynomial on the data and find out what the characteristic
# function of the measurements are
result_x <- lm(polynomial_formula("x", "X", 10),
               data = stationary) %>% 
    summary()
result_y <- lm(polynomial_formula("y", "Y", 10),
               data = stationary) %>% 
    summary()

# Save these results for later
saveRDS(result_x,
        file.path("Calibration experiments", 
                  "Results",
                  "polynomial_for_x.Rds"))
saveRDS(result_y,
        file.path("Calibration experiments", 
                  "Results",
                  "polynomial_for_y.Rds"))

# Create a function that will extract all of the coefficients
extract_coefficients <- function(x){
    x$coefficients[,1] %>% 
        as.numeric() %>% 
        return()
}

# Create a function that will correct the measurements that it
# receives
correct_measurements <- function(x, parameters){
    X <- matrix(0, nrow = length(x), ncol = length(parameters))
    for(i in seq_along(parameters)){
        X[,i] <- parameters[i] * x^(i - 1)
    }
    X %>% 
        rowSums() %>% 
        as.numeric() %>% 
        return()
}

# Do the actual correction in the x and y dimensions separately
new_x <- correct_measurements(stationary$mu_x, 
                              extract_coefficients(result_x))
new_y <- correct_measurements(stationary$mu_y, 
                              extract_coefficients(result_y))

# Plot the new coordinates
plot(new_x, new_y)
plot(stationary$mu_x, stationary$mu_y)
plot(stationary$X, stationary$Y)

########################
# Estimation: 1D with Z-scores

# Create a function that will output a formula for a polynomial of the pwr^th
# degree with y as dependent variable and x as the independent variable
polynomial_formula <- function(x, y, pwr){
    polyn <- paste(y, "~", x)
    for(i in 2:pwr){
        polyn <- paste(polyn, "+", "I(",
                       paste0(x, "^", i), ")")
    }
    polyn %>% 
        formula() %>% 
        return()
}

# Create Z-scores
stationary <- stationary %>% 
    mutate(Z_x = scale(mu_x), 
           Z_y = scale(mu_y), 
           Z_X = scale(X), 
           Z_Y = scale(Y))

# Estimate the 10th degree polynomial on the data and find out what the characteristic
# function of the measurements are
result_x <- lm(polynomial_formula("Z_x", "Z_X", 10),
               data = stationary) %>% 
    summary()
result_y <- lm(polynomial_formula("Z_y", "Z_Y", 10),
               data = stationary) %>% 
    summary()

# Save these results for later
saveRDS(result_x,
        file.path("Calibration experiments", 
                  "Results",
                  "polynomial_for_x.Rds"))
saveRDS(result_y,
        file.path("Calibration experiments", 
                  "Results",
                  "polynomial_for_y.Rds"))

# Create a function that will extract all of the coefficients
extract_coefficients <- function(x){
    x$coefficients[,1] %>% 
        as.numeric() %>% 
        return()
}

# Create a function that will correct the measurements that it
# receives
correct_measurements <- function(x, parameters){
    X <- matrix(0, nrow = length(x), ncol = length(parameters))
    for(i in seq_along(parameters)){
        X[,i] <- parameters[i] * x^(i - 1)
    }
    X %>% 
        rowSums() %>% 
        as.numeric() %>% 
        return()
}

# Do the actual correction in the x and y dimensions separately
new_x <- correct_measurements(stationary$Z_x, 
                              extract_coefficients(result_x))
new_y <- correct_measurements(stationary$Z_y, 
                              extract_coefficients(result_y))

# Plot the new coordinates
plot(new_x, new_y)
plot(stationary$Z_x, stationary$Z_y)
plot(stationary$Z_X, stationary$Z_Y)

########################
# Estimation: 2D

degree <- 2

# Create a function that will create the data matrix
data_matrix <- function(x, y, pwr, print_formula = FALSE){
    X <- matrix(0, nrow = length(x), ncol = (pwr + 1)^2)
    f <- 1
    form <- ""
    for(i in 0:pwr){
        for(j in 0:pwr){
            # If the total exponents in the formula exceeds the order of the 
            # polynomial, then you have to skip this one
            if(i + j > pwr){
                next
            }

            # Add the interaction term of x and y to the polynomial
            X[,f] <- x^i * y^j
            f <- f + 1

            # Add the function to the string of the formula
            form <- paste0(form, 
                           ifelse(i + j == 0, "", " + "),
                           "x^", i, " y^", j)
        }
    }
    # Delete the columns in which there are only 0's
    X <- X[,1:(f - 1)]

    # Print formula if requested
    if(print_formula){
        print(form)
    }
    return(X)
}

fit_polynomial <- function(x, y, pwr){
    # Unlist x and y
    x <- as.matrix(x)
    y <- as.matrix(y)

    # Create the data matrix based on x and y
    xy <- data_matrix(x[,1], x[,2], pwr)

    # Get beta through the formula of the least-squares method
    beta <- solve( t(xy) %*% xy ) %*% t(xy) %*% y
    
    return(beta)
}

# Fit the measured coordinates to the real coordinates with a 4th degree 
# polynomial (after 5 becomes singular)
result <- fit_polynomial(stationary[,c("mu_x", "mu_y")], 
                         stationary[,c("X", "Y")], 
                         degree)

# Save these results for later
saveRDS(result,
        file.path("Calibration experiments", 
                  "Results",
                  "polynomial_for_2D.Rds"))

# Create a function that will correct the measurements that it
# receives
correct_measurements <- function(x, parameters, pwr){
    # Unlist x: Otherwise error in data_matrix
    x <- as.matrix(x)

    # Create the data matrix based on x and y. Unfortunately, difficult to extract
    # the order of the polynomial from only its parameters, which is why it has
    # to be provided as an argument
    XY <- data_matrix(x[,1], x[,2], pwr)

    # Do the actual conversion of the measured positions to the newly estimated
    # ones
    ( XY %*% parameters ) %>% 
        return()
}

# Do the actual correction in the x and y dimensions separately
new_xy <- correct_measurements(stationary[,c("mu_x", "mu_y")], 
                               result,
                               degree)

# Plot the new coordinates
plot(new_xy[,1], new_xy[,2])
plot(stationary$mu_x, stationary$mu_y)
plot(stationary$X, stationary$Y)
