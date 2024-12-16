#' Compute summary statistics
#' 
#' Overhead function that takes in the data and then computes all summary 
#' statistics contained in \code{fx}.
#' 
#' @param data Data.frame containing the columns in \code{.vars} and for which 
#' to compute the summary statistics.
#' @param fx List of functions that compute the summary statistics of interest.
#' Should take in two numeric vectors.
#' @param .vars Character vector or matrix containing the columns for which to 
#' compute the summary statistics. If matrix, each row should contain the columns
#' to compare.
#' @param .by Character vector denoting the columns to group by. Defaults to
#' \code{"id"}.
#' 
#' @return Data.frame containing the summary statistics of interest.
#' 
#' @export 
summary_statistics <- function(data, 
                               fx, 
                               .vars,
                               .by = "id") {

    # Check and/or transform the .vars argument
    if(!is.matrix(.vars)) {
        .vars <- matrix(.vars, ncol = 2)
    }

    # Check whether the .by arguments are all available in the data. If not, 
    # delete those that cannot be found
    idx <- sapply(.by, \(x) !is.null(data[, x]))
    .by <- .by[idx]
    
    # Group the data and create a new data.frame that will contain all results
    grouped_data <- data %>% 
        dplyr::group_by(tidyselect::all_of(.by))

    results <- grouped_data %>% 
        dplyr::summarize() %>% 
        list()

    # Loop over all functions and add their values to the dataframe.
    for(i in fx) {
        for(j in seq_len(nrow(.vars))) {
            # Get the name of the function. Will be used as column name for its 
            # output
            column <- i %>% 
                substitute() %>% 
                as.character() %>% 
                paste0("_", .vars[j, 1], vars[j, 2])
    
            # Add the summary statistics to the results list
            my_result <- grouped_data %>% 
                dplyr::summarize(result = i(tidyselect::all_of(.vars[j, 1]), 
                                            tidyselect::all_of(.vars[j, 2]))) %>% 
                dplyr::select(result) %>% 
                setNames(column)
    
            results <- append(results, my_result)
        }
    }

    # Combine and return
    return(do.call("cbind", results))
}

# Create a function in which we will compute the summary statistics of interest.
# Key here is to compare the positions that result from preprocessing with the 
# real positions that are contained in the `original` data set.
#
# This function takes in the arguments `data` -- the preprocessed data set -- 
# and `kind` -- the key that defines the original data. It then computes the 
# difference between preprocessed and actual positions and provides some 
# summary statistics.
# compute_summary_statistics <- function(data, kind) {
#     original[[kind]] %>% 
#         # Add the original dataset into the one it is being compared to and 
#         # give `X` and `Y` as labels for the real positions (similar to 
#         # stationary preprocessing)
#         dplyr::rename(X = x, 
#                       Y = y) %>% 
#         dplyr::full_join(data, by = c("nsim", "time", "id")) %>% 

#         # Compute the difference between filtered and expected positions. Used 
#         # to measure the extent to which systematic error is present in the data
#         dplyr::mutate(difference_x = X - x, 
#                       difference_y = Y - y, 
#                       distance = sqrt((X - x)^2 + (Y - y)^2)) %>% 

#         # Compute the statistics of interest per simulation and id. This 
#         # will allow for a more broad view on where it still might go awry
#         dplyr::group_by(nsim, id) %>% 
#         dplyr::arrange(time) %>% 
#         dplyr::summarize(# Statistics about how close we are to the actual 
#                          # positions
#                          mean_diff_x = mean(difference_x, na.rm = TRUE), 
#                          mean_diff_y = mean(difference_y, na.rm = TRUE), 
#                          mean_dist = mean(distance, na.rm = TRUE),
#                          q025_diff_x = quantile(difference_x, probs = 0.025, na.rm = TRUE),
#                          q025_diff_y = quantile(difference_y, probs = 0.025, na.rm = TRUE),
#                          q025_diff_y = quantile(distance, probs = 0.025, na.rm = TRUE),
#                          q975_diff_x = quantile(difference_x, probs = 0.975, na.rm = TRUE),
#                          q975_diff_y = quantile(difference_y, probs = 0.975, na.rm = TRUE), 
#                          q975_diff_y = quantile(distance, probs = 0.975, na.rm = TRUE), 

#                          # Statistics about the size of the measurement error
#                          # (compared to the actual positions)
#                          rmse_diff_x = sd(difference_x, na.rm = TRUE), 
#                          rmse_diff_y = sd(difference_y, na.rm = TRUE), 
#                          rmse_dist = sd(distance, na.rm = TRUE),
#                          mae_diff_x = mean(abs(difference_x), na.rm = TRUE), 
#                          mae_diff_y = mean(abs(difference_y), na.rm = TRUE), 
#                          mae_dist = mean(abs(distance), na.rm = TRUE), 
                         
#                          # Autocorrelation in the residuals
#                          auto_x = cor(difference_x[2:length(difference_x)], 
#                                       difference_x[2:length(difference_x) - 1],
#                                       use = "pairwise.complete.obs"), 
#                          auto_y = cor(difference_y[2:length(difference_y)], 
#                                       difference_y[2:length(difference_y) - 1], 
#                                       use = "pairwise.complete.obs"), 
#                          auto_dist = cor(distance[2:length(distance)], 
#                                          distance[2:length(distance) - 1], 
#                                          use = "pairwise.complete.obs")) %>% 
#         dplyr::ungroup() %>% 
#         suppressMessages() %>% 
#         return()
# }





################################################################################
# SUMMARY STATISTICS

#' Compute Bias
#' 
#' @param x Numeric vector of values to be compared to \code{y}.
#' @param y Numeric vector of values to be compared to \code{x}.
#' 
#' @return Numeric denoting the mean difference between the two vectors
#' 
#' @export
bias <- function(x, y) {
    (x - y) %>% 
        mean() %>% 
        return()
}

#' Compute RMSE
#' 
#' @param x Numeric vector of values to be compared to \code{y}.
#' @param y Numeric vector of values to be compared to \code{x}.
#' 
#' @return Numeric denoting the RMSE of the difference between the two vectors
#' 
#' @export
rmse <- function(x, y) {
    (x - y) %>% 
        sd() %>% 
        return()
}

#' Compute MAE
#' 
#' @param x Numeric vector of values to be compared to \code{y}.
#' @param y Numeric vector of values to be compared to \code{x}.
#' 
#' @return Numeric denoting the MAE of the difference between the two vectors
#' 
#' @export
mae <- function(x, y) {
    (x - y) %>% 
        abs() %>% 
        mean() %>% 
        return()
}