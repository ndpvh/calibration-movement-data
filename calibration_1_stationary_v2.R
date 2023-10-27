library(tidyverse)
library(data.table)
library(nloptr)
library(mvtnorm)

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
                                "preprocessed_stationary.Rds"))

# 1) Estimating the error variance
#
# Some assumptions:
#   - The actual measured position is equal to the average of all measured 
#     positions for a given tag
#   - The error variance is equal for each of the tags
#   - This error is normally distributed
#   - Independence of the error between the X and Y dimensions
#   - Independence of error over time
#
# The model below is completely free. It can also be useful to take rows
# and columns of (X,Y) into account, already putting them on a line.
#
# I am going to have to change this: It is far too difficult to be able to 
# pull this off

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
decomposition_of_variance <- function(x){
    # Get the overall standard deviation `sigma` and the means associated to 
    # each of the tags `mu`
    sigma <- c(x[1], 0, 0, x[1]) %>% 
        matrix(ncol = 2)
    mu <- x[-1] %>% 
        matrix(ncol = 2)

    # Define the min-log-likelihood variable
    MLL <- 0

    # Put the estimated values for the positions into the dataframe
    stationary <- stationary %>% 
        ungroup() %>% 
        # Add the estimated means to the dataframe
        tidyr::nest(data = -tag) %>% 
        mutate(mu_x = mu[,1],
               mu_y = mu[,2]) %>%
        tidyr::unnest(data)
    
    # Evaluate the min-log-likelihood of the measured positions given the 
    # estimated positions
    min_log_likelihood(stationary[,c("x", "y")],
                       stationary[,c("mu_x", "mu_y")],
                       sigma) %>% 
        return()


    # Get each of the tags and loop over each of them
    tags <- unique(stationary$tag)
    for(i in seq_along(tags)){
        # Select the data of one specific tag
        y <- stationary %>% 
            ungroup() %>% 
            filter(tag == tags[i]) %>% 
            select(x, y) %>% 
            as.matrix()
            
        # Evaluate the likelihood of the data for a tag under the normal model
        # with an estimated mean value and an estimated error variance that is 
        # the same for each of the tags
        MLL <- dmvnorm(y, mu[i,], sigma) %>% 
            log() %>% 
            sum() %>% 
            `*` (-1) %>% 
            `+` (MLL)
    }

    return(MLL)
}

# Run the estimation using nloptr
n_tags <- stationary$tag %>% 
    unique() %>% 
    length()
set.seed(35356) # New Storm for Older Lovers - La Dispute
result <- nloptr(x0 = c(1, rnorm(n_tags * 2)),
                 eval_f = decomposition_of_variance,
                 lb = c(10^(-5), rep(-1000, each = n_tags * 2)),
                 ub = rep(1000, each = n_tags * 2 + 1),
                 opts = list("algorithm" = "NLOPT_GN_DIRECT", 
                             "ftol_rel" = 10^(-15),
                             "xtol_rel" = 10^(-15),
                             "maxeval" = 10^5))

