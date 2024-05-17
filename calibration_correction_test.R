################################################################################
# Purpose: Test whether the `preprocess_measurements` function enhances the    #
#          quality of the data that we have. For this, we compare several      #
#          features of the raw data against the same features for the          #
#          preprocessed/corrected data. The features to be checked are the     #
#          following:                                                          #
#              - Measurement error                                             #
#              - Distance from real location                                   #
#                                                                              #
#          To estimate these features, we will make use of stationary data     #
#          from the 21-10-2023 and from 22-12-2023, as these data are well-    #
#          behaved (no movement) and not from the same day as the data         #
#          `preprocess_measurements` is based on (14-10-2023)                  #
################################################################################

source(file.path("utility", "utility.R"))

# Get the dates of both stationary dataframes, which will be used both for 
# reading in the data and for saving the results of these data.
stationary <- c("21-10-2023",
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
                      "preprocess_measurements",
                      paste0(name, "_", x, ".Rds")))
}





#-------------------------------------------------------------------------------
# Estimating the error variance
#-------------------------------------------------------------------------------

# Here, we use a same approach as in `calibration_stationary.R` to estimate the 
# error variance of the tags. We do this for each dataset separately, furthermore
# distinguishing between "corrected" and "uncorrected" data.

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
set.seed(3347456) # All the Best Cowboys Have Daddy Issues - Senses Fail
correction <- c("corrected", "uncorrected")
for(i in seq_along(stationary)){
    covariances <- list()

    for(j in seq_along(correction)) {
        # Load the data and use the `preprocess_measurements` function or not, 
        # depending on the value of `correction`
        data <- load_stationary(stationary[i])
        if(correction[j] == "corrected") {
            data <- preprocess_measurements(data) %>% 
                filter(!(is.na(x)), 
                       !(is.na(y)))
        }

        # Load data from specific date and change it somewhat
        data <- data %>% 
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
        covariances[[j]] <- bootstrap_data(data, 1000) %>% 
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
        result <- covariances[[j]] %>% 
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
                    paste0("error_covariance_",
                           correction[j]), 
                    stationary[i])

        rm(data, result)
        gc()
    }

    # Create histograms of the bootstrapped samples for the covariances and 
    # save these in the specified location.
    #
    # First create a function that will make the plot
    booted_histogram <- function(idx, col){
        plot_data <- covariances[[idx]][, col] %>% 
            as.data.frame() %>% 
            setNames(c("X"))

        (ggplot(plot_data, aes(x = X)) +
            geom_histogram(col = "black", 
                           fill = "cornflowerblue")) %>% 
            return()
    }

    multiple_histograms <- function(idx) {
        ggarrange(booted_histogram(idx, "var_x") + labs(title = "VAR(X)", x = ""), 
                  booted_histogram(idx, "var_y") + labs(title = "VAR(Y)", x = ""),
                  booted_histogram(idx, "cov_xy") + labs(title = "COV(X,Y)", x = ""),
                  nrow = 1) %>% 
            return()
    }

    # Make six plots and bind them together
    plt <- ggarrange(plotlist = list(NULL, 
                                     multiple_histograms(1),
                                     NULL,
                                     multiple_histograms(2)),
                     nrow = 2,
                     ncol = 2,
                     widths = c(0.25, 1),
                     labels = c("Corrected", "", "Uncorrected", ""))

    # Finally save them
    ggsave(file.path("figures", 
                     "calibration", 
                     "preprocess_measurements", 
                     paste0("error_covariance_", stationary[i], ".png")),
           plot = plt,
           units = "px",
           width = 3000, 
           height = 2600)

    # Release the memory that is held up by the bootstrapped data and the data 
    # itself.
    rm(covariances, plt)
    gc()
}


# Conclusions:
#   - 21-10-2023: Error is reduced significantly in the x-direction, but also 
#                 increased significantly in the y-direction. In standard 
#                 deviations, this increase and decrease is equal to (approx.):
#                     - x: 2.24cm (uncorrected) vs. 2.00cm (corrected)
#                     - y: 1.73cm (uncorrected) vs. 2.24cm (corrected)
#
#   - 22-12-2023: Error is reduced significantly in both the x- and y-
#                 direction, showing that our preprocessing accomplishes its 
#                 goal for these datasets. This decreases is equal to (approx.):
#                 
#                 ~ 4 anchors      
#                     - x: 8.60cm (uncorrected) vs. 7.42cm (corrected)
#                     - y: 3.87cm (uncorrected) vs. 3.32cm (corrected)
#
#                 ~ 6 anchors      
#                     - x: 6.71cm (uncorrected) vs. 4.58cm (corrected)
#                     - y: 8.19cm (uncorrected) vs. 7.28cm (corrected)
#
# Reductions in error variance are in general small, suggesting we should find 
# other ways to handle this kind of error. It is a start though.





#-------------------------------------------------------------------------------
# Estimating absolute error
#-------------------------------------------------------------------------------

# Here, we will estimate the absolute positional error in measuring the x and y 
# positions. Given that the stationary tags always rest at a same location, we 
# can estimate how far off the measured positions are from the real ones. 
#
# Importantly, this will give us an idea of bias in the positions, not to be 
# confused with the errors we computed in the previous section.