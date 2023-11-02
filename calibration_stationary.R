################################################################################
# Purpose: Analyze the preprocessed stationary data gathered on 14/10/2023     #
#          and 21/10/2023. The main goal of this analysis is to inform us on   #
#          ways to deal with systematic and unsystematic noise in our data.    #
#          Consequently, the results from this analysis were used to create    #
#          the `filter_measurements` function in the "utility" folder.         #
################################################################################

source(file.path("utility", "utility.R"))

# Get the dates of both stationary dataframes, which will be used both for 
# reading in the data and for saving the results of these data.
stationary <- c("14-10-2023",
                "21-10-2023")

#-------------------------------------------------------------------------------
# Some utilities
#-------------------------------------------------------------------------------

# Function that will read in a given stationary dataset. Here, `x` is one of the 
# two dates contained within the `stationary` vector.
load_stationary <- function(x){
    readRDS(file.path("data", 
                      paste0("preprocessed_stationary_", x, ".Rds"))) %>% 
        return()
}

# Function that will save the results for a given stationary dataset. Here, `x` 
# is again one of the dates contained within `stationary`. `name` is the name 
# you want to give this specific result. Finally, `result` is the thing you want 
# to save under that name.
save_result <- function(result, name, x){
    saveRDS(result, 
            file.path("results", 
                      paste0(name, "_", x, ".Rds")))
}





#-------------------------------------------------------------------------------
# Estimating the error variance
#-------------------------------------------------------------------------------

# Here, we estimate the error covariances of the tags in the stationary data. 
# 
# While deleted in this version of the scripts, we did notice some issues with 
# the data. More specifically, we found that several parametric assumptions
# (e.g., normality) do not hold for our data, making it difficult to estimate
# the error variance in a parametric way. We therefore use a nonparametric way 
# for estimating the error covariance matrix, which has the following 
# assumptions:
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

# Loop over the different dates and do all your estimation 
set.seed(39) # The Messenger - Thrice
for(i in seq_along(stationary)){
    # Load data from specific date and change it somewhat
    data <- load_stationary(stationary[i]) %>% 
        # Compute several center statistics, namely mean, median, and mode per
        # tag
        group_by(tag) %>%         
        mutate(mu_x = mean(x),
               mu_y = mean(y),
               median_x = median(x), 
               median_y = median(y), 
               mode_x = mean(mlv(x, method = "mfv")), 
               mode_y = mean(mlv(y, method = "mfv"))) %>% 
        ungroup() %>% 
        # Create a "corrected" version of the x- and y-positions using the mean.
        # This one will be used in the estimation of the error covariances.
        #
        # Can be replaced with median or mode too, if desired.
        mutate(corrected_x = x - mu_x, 
               corrected_y = y - mu_y)

    # Bootstrap the data using this function and immediately compute the necessary
    # summary statistics: 2 variances and 1 covariance.
    #
    # Be careful, this uses a lot of memory!
    covariances <- bootstrap_data(data, 2500) %>% 
        # Compute the covariances based on the corrected x- and y-positions. 
        # Importantly, this is done for each separate bootstrapped sample.
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

    # Create summary statistics for each of the covariances, and more specifically 
    # given quantiles of the bootstrapped distribution. This should give us an 
    # idea of how badly off we are.
    #
    # Here, we use a 99% CI, just to be sure
    result <- covariances %>% 
        # Summarize the different variables into CI and mean
        summarize(lb_var_x = quantile(var_x, 0.005), 
                  lb_var_y = quantile(var_y, 0.005),
                  lb_cov_xy = quantile(cov_xy, 0.005), 
                  m_var_x = mean(var_x), 
                  m_var_y = mean(var_y), 
                  m_cov_xy = mean(cov_xy), 
                  ub_var_x = quantile(var_x, 0.995), 
                  ub_var_y = quantile(var_y, 0.995),
                  ub_cov_xy = quantile(cov_xy, 0.995)) %>% 
        # Restructure the dataframe to be more useful: Put each of the 
        # covariances as a row, and the lower bounds, mean, and upper bounds 
        # as the columns 
        unlist() %>% 
        matrix(nrow = 3, ncol = 3) %>% 
        as.data.frame() %>% 
        setNames(c("lb", "mean", "ub")) %>% 
        cbind(covariance = c("var_x", "var_y", "cov_xy"))  

    # Save these results
    save_result(result, 
                "error_covariance_", 
                stationary[i])

    # Create histograms of the bootstrapped samples for the covariances and 
    # save these in the specified location.
    #
    # First create a function that will make the plot
    booted_histogram <- function(col){
        plot_data <- covariances[, col] %>% 
            as.data.frame() %>% 
            setNames(c("X"))

        (ggplot(plot_data, aes(x = X)) +
            geom_histogram(col = "black", 
                        fill = "cornflowerblue")) %>% 
            return()
    }

    # Make three plots and bind them together
    plt <- ggarrange(booted_histogram("var_x") +
                        labs(title = "VAR(X)", 
                             x = "Bootstrapped value"),
                     booted_histogram("var_y") +
                        labs(title = "VAR(Y)", 
                             x = "Bootstrapped value"),
                     booted_histogram("cov_xy") +
                        labs(title = "COV(X,Y)", 
                             x = "Bootstrapped value"),
                     nrow = 1)

    # Finally save them
    ggsave(file.path("figures", "calibration_stationary", 
                     paste0("error_covariance_", stationary[i], ".png")),
           plot = plt,
           units = "px",
           width = 3000, 
           height = 1300)

    # Release the memory that is held up by the bootstrapped data and the data 
    # itself.
    rm(data, covariances, result, plt)
    gc()
}

# Interpretation of the results: 
#   





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
