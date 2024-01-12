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
                "21-10-2023",
                "4_anchors_22-12-2023",
                "6_anchors_22-12-2023")

#-------------------------------------------------------------------------------
# Some utilities
#-------------------------------------------------------------------------------

# Function that will read in a given stationary dataset. Here, `x` is one of the 
# two dates contained within the `stationary` vector.
load_stationary <- function(x){
    readRDS(file.path("data", 
                      "calibration",
                      "stationary",
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
                      "calibration",
                      "stationary",
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
    # Be careful, this uses a lot of memory! Currently, 1000 samples is all my 
    # computer can manage (vectorized). Could consider going for unvectorized at 
    # a later stage
    covariances <- bootstrap_data(data, 1000) %>% 
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
                "error_covariance", 
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
    ggsave(file.path("figures", 
                     "calibration", 
                     "stationary", 
                     paste0("error_covariance_", stationary[i], ".png")),
           plot = plt,
           units = "px",
           width = 3000, 
           height = 1300)

    # Release the memory that is held up by the bootstrapped data and the data 
    # itself.
    rm(data, covariances, result, plt)
}

# Interpretation of the results: 
#   - Covariances between x and y are as good as 0, so not relationship in the 
#     error between both dimensions
#   - Variances are twice as high for the second calibration period, which might 
#     suggest that there is more error on the sides than in the center 
#     (in calibration period 2, only 3 tags were put in the center)
#   - Worst case scenario -- which is calibration period 2, upper bound of the
#     99%CI -- there is about 2.25cm of standard error in the x-direction, and
#     1.80cm of standard error in the y-direction 
#
# Additional comments after calibration of 22-12-2023
#   - Error on this day is a lot bigger. More specifically, the error in standard
#     deviations is: 
#       - 4 anchors: 8.7cm (x) and 3.9cm (y)
#       - 6 anchors: 6.9cm (x) and 8.3cm (y)





#-------------------------------------------------------------------------------
# Sampling rate
#-------------------------------------------------------------------------------

# Here, we estimate the attained sampling rate for our stationary data.

# Create a function that will take in the duration, order them according to 
# size, and then compute the mean difference between each of the durations
sampling_rate <- function(x){
    x %>% 
        as.numeric() %>% 
        sort() %>% 
        diff() %>% 
        mean() %>% 
        return()
}

# Again loop over each of the stationary datasets
for(i in seq_along(stationary)){
    # Load the stationary data for a given date and convert the timestamps to 
    # milliseconds
    data <- load_stationary(stationary[i]) %>% 
        mutate(duration = convert_to_millisecond(timestamp))

    # Compute the mean sampling rate for each tag separately and make the 
    # dataframe somewhat easier to interpret
    result <- data %>% 
        # Get the sampling rate in msec
        group_by(tag) %>% 
        mutate(msec = sampling_rate(duration)) %>% 
        ungroup() %>% 
        # Convert the sampling rate to Hz: 1 / s -> 1000 / msec
        mutate(Hz = 1000 / msec) %>% 
        # Only retain those variables that you want to interpret
        group_by(tag, msec, Hz) %>% 
        tidyr::nest() %>% 
        select(-data)

    # Save this result
    save_result(result, 
                "sampling_rate", 
                stationary[i])

    # Make another histogram to visualize the result across tags 
    plt <- ggplot(result, 
                  aes(x = Hz)) +
        geom_histogram(color = "black", 
                       fill = "cornflowerblue") +
        labs(title = "Sampling rate per tag", 
             x = "Hz") +
        geom_vline(xintercept = 5, 
                   color = "red")

    ggsave(file.path("figures", 
                     "calibration", 
                     "stationary",
                     paste0("sampling_rate_", stationary[i], ".png")),
           plot = plt, 
           units = "px", 
           width = 1000, 
           height = 1100)
}

# Interpretation of the result:
#   - On the 14th of October, sampling rates varied substantially between 2 and 
#     4Hz instead of the expected 5Hz. Sampling rates should thus be increased 
#     for our future experiments.
#   - On the 21st of October, we found a higher than expected sampling rate. This 
#     might suggest that our initial attempts of changing the sampling rate from 
#     5Hz to 8Hz on that day might have been successful, leading to an observed 
#     sampling rate of 6-7Hz 
#   - Some tags seem to perform worse and only send out responses every once in 
#     a while.
#
# Additional comments after calibration on 22-12-2023
#   - On this day, sampling with the 6 anchors remained relatively similar to 
#     the sampling rate on 21-10-2023 (although with a slight loss of frequency).
#   - Sampling with 4 anchors had an effect on the sampling frequency. Need to 
#     find out how to increase it again.





#-------------------------------------------------------------------------------
# Systematic distortions
#-------------------------------------------------------------------------------

# Here, we estimate a function to filter out systematic distortions in the data. 
# We do this by estimating a polynomial on the stationary data and checking
# the differences between the functions estimated on both days. We have two 
# major assumption in this analysis, which is: 
#   - The actual measured position is equal to the average of all measured 
#     positions for a given tag
#   - The systematic error in the x- and y-dimensions is independent


# Create a function that will output a formula for a polynomial of the degree^th
# degree with `dependent` as dependent variable and `independent` as the 
# independent variable.
polynomial_formula <- function(independent, dependent, degree){
    # Start with putting the dependent variable against the independent variable
    polyn <- paste(dependent, "~", independent)

    # Add each term of the polynomial to this string. Should be of the format 
    # I(x^degree)
    for(i in 2:degree){
        polyn <- paste(polyn, "+", "I(",
                       paste0(independent, "^", i), ")")
    }

    # Return this string as a formula to be used in `lm`
    polyn %>% 
        formula() %>% 
        return()
}

# Create a formula that will do the actual estimation, based on data, a dependent 
# variable, an independent variable, and a degree
fit_polynomial <- function(data, independent, dependent, degree){
    polynomial_formula(independent, dependent, 10) %>% 
        lm(data = data) %>% 
        summary() %>% 
        return()
}

# Create a function that will extract all of the coefficients
extract_coefficients <- function(x){
    x$coefficients[,1] %>% 
        as.numeric() %>% 
        return()
}

# As one can predict already: Loop over all stationary data again.
for(i in seq_along(stationary)){
    # Load the data for a given date and standardize all data: Measured positions
    # and supposed real positions. Standardization is performed to make the
    # function scalable to other data in which the absolute positions are not 
    # the same as these of the calibration phase.
    #
    # Importantly, the kind of standardization that is performed is not the 
    # standard way of doing this. Rather, we want to transform all measured positions 
    # that fall within the measurement space (i.e., the space bounded by the 
    # anchors) to fall within -1 and 1. This will allow us to more readily 
    # translate the calibration on one day to the measurements on another, as 
    # you are not dependent on the actual measurements for the standardization
    # (in contrast to the scaled equivalent of the Z-score). Importantly, the 
    # `minmax_standardize` function does this transformation, as defined in the 
    # utility functions.
    #
    # Unfortunately, we don't have the anchors' positions, so we will assume that 
    # they lie at the sides of the "real positions" grid. In reality, this will 
    # probably be off by a few tens of centimeters, but we will have to do 
    # another, more precisely done calibration to combat the issues that this 
    # brings.
    data <- load_stationary(stationary[i]) %>% 
        mutate(standardized_x = minmax_standardize(x, 
                                                   min_x = min(X), 
                                                   max_x = max(X)), 
               standardized_y = minmax_standardize(y, 
                                                   min_x = min(Y), 
                                                   max_x = max(Y)), 
               standardized_X = minmax_standardize(X), 
               standardized_Y = minmax_standardize(Y))

    # Estimate the 10th degree polynomial on the standardized data
    pars <- list("x" = fit_polynomial(data, 
                                      "standardized_x", 
                                      "standardized_X", 
                                      10),
                 "y" = fit_polynomial(data, 
                                      "standardized_y", 
                                      "standardized_Y", 
                                      10))

    # Save these results
    save_result(pars, 
                "polynomial", 
                stationary[i])
                 
    # Create a dataframe to be used in `correct_distortion` which will 
    # contain all of the parameters for the separate dimensions
    result <- cbind("x" = extract_coefficients(pars[["x"]]), 
          "y" = extract_coefficients(pars[["y"]])) %>% 
        as.data.frame()

    # Save these parameters as well 
    save_result(result, 
                "parameters_polynomial", 
                stationary[i])

    # Check how well this polynomial is able to correct for the distortion 
    # in the data.
    coordinates <- data %>% 
        ungroup() %>% 
        select(standardized_x, standardized_y) %>% 
        correct_distortion(result)

    # Create a plot of the locations for the distorted and the corrected data.
    #
    # First, create a function for this
    location_plot <- function(data, x, y, title){
        # Create plot data with the string column names `x` and `y`
        plot_data <- data[, c(x, y)] %>% 
            as.data.frame() %>% 
            setNames(c("X", "Y"))

        # Create the plot itself
        plt <- ggplot(plot_data, 
                      aes(x = X, y = Y)) +
            geom_point(size = 2, 
                       color = "black") +
            labs(x = "Standardized X", 
                 y = "Standardized Y", 
                 title = title)

        return(plt)
    }
    
    # Create and bind the plots of the distorted and corrected data
    plt <- ggarrange(location_plot(data, 
                                   "standardized_x", 
                                   "standardized_y", 
                                   "Original"), 
                     location_plot(coordinates, 
                                   "standardized_x", 
                                   "standardized_y", 
                                   "Corrected"),
                     nrow = 1)

    # Save this plot 
    ggsave(file.path("figures", 
                     "calibration", 
                     "stationary",
                     paste0("correct_distortion_", stationary[i], ".png")), 
           plot = plt, 
           units = "px", 
           width = 2000, 
           height = 1200)
}

# Interpretation of the results: 
#   - For both datasets, the distortions seems to disappear somewhat
#   - However, they do still show some distortion of the distances between each 
#     of the measured tags
#   - Question: is this a problem in our measurements, or a problem in where we 
#     put the tags (i.e., is our assumption of the supposed coordinates X and Y 
#     correct, or did we make a human error here)
#   - Parameters between the two datasets are similar, though do not match 
#     exactly. Given that their are more measured positions in the data from 
#     the 14th of October, we will assume these estimates have a higher power, 
#     and will use those for the correction of our distortion.
#
# Additional comments after calibration on 22-12-2023
#   - The 4-anchor data is corrected better than the 6-anchor data, but both 
#     are worse than the 14-10-2023 counterpart.