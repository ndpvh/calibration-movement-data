#' Execute function on a moving window
#' 
#' Create a moving window of a given span and execute the specified function on 
#' the data within that window. In this project, a moving window is used to 
#' smooth the data.
#' 
#' Importantly, the span is a relative variable, that is it represents the 
#' number of neighbours to take in to account in each direction relative to the 
#' central data point to be smoothed. In other words, if you take `span = 2`, 
#' then you will smooth over the indices `c(-2, -1, 0, 1, 2)` where index `0` 
#' represents the datapoint to be smoothed. These kinds of indices are also 
#' provided to the researcher in the function `fx`, allowing them to assign 
#' different weights to each of the neighbours around a central observation. 
#' The relative indices in a window are found under the column `index`.
#' 
#' @param data Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' @param span Integer denoting the number of neighbours in each direction to 
#' take into account within each moving window.
#' @param fx Function to execute on the data within the moving window. Should 
#' take in a similar data structure to `data` and output a single-row dataframe
#' with columns `x` and `y` (representing the smoothed values of both variables).
#' @param .by String denoting whether the moving window should be taken with 
#' respect to a given grouping variable. Defaults to `NULL`
#' 
#' @return Smoothed dataframe with a similar structure as `data`
#' 
#' @export 
moving_window <- function(data, 
                          span, 
                          fx,
                          .by = NULL) {

    # Check whether there is a grouping variable to account for. If not, then 
    # we just use a dummy
    if(!is.null(.by)) {
        group <- unique(data[,.by])
    } else {
        group <- 1
    }

    # Go over each of the data points and smooth the data using the moving 
    # window. We dispatch/loop over all the different possibilities in the 
    # grouping variable, for each of which we will create an individual 
    # dataframe to be smoothed.
    #
    # Within this loop, we do the following:
    #   - Create a group-specific dataframe
    #   - Determine which time points to contain within each moving window based
    #     on the specified span
    #   - Apply the specified function to the data
    #   - Put these data in a list
    all_data <- list()
    for(i in seq_along(group)) {
        # Get the data specific to the group (if .by is not NULL)
        if(!is.null(.by)) {
            individual_data <- data[data[, .by] == group[i],]
        } else {
            individual_data <- data
        }

        # Add an integer index to the data. This will allow researchers some 
        # flexibility with regard to how to assign weights to the different
        # observations
        individual_data <- individual_data %>% 
            dplyr::mutate(index = dplyr::row_number())

        # Define the length of the data and create the time points to be taken 
        # into account using the `span` argument
        N <- nrow(individual_data)
        idx <- cbind(from = 1:N - span, 
                     to = 1:N + span,
                     obs = 1:N) %>% 
            as.data.frame()

        # Delete nonexisting indices from this list
        idx[idx < 1] <- 1
        idx[idx > N] <- N

        # Add nested tibbles that contain the data for the specified spans to
        # the indices. Make the indices in `index` relative to the primary 
        # observation index in 'obs'
        individual_data <- idx %>% 
            dplyr::rowwise() %>% 
            dplyr::mutate(data = individual_data[from:to,] %>% 
                                     dplyr::mutate(index = index - obs) %>% 
                                     dplyr::arrange(time) %>% 
                                     tidyr::nest()) %>% 
            dplyr::ungroup() %>% 
            cbind(individual_data) %>% 
            dplyr::select(-x, -y)

        # Apply the specified function to the selected data
        all_data[[i]] <- individual_data %>% 
            dplyr::rowwise() %>% 
            dplyr::mutate(data = data$data %>% 
                                     as.data.frame() %>% 
                                     fx()) %>% 
            tidyr::unnest(data) %>% 
            dplyr::ungroup()
    }

    # Bind all data together and sort based on time
    data <- do.call("rbind", all_data) %>% 
        dplyr::select(-from, -to, -obs, -index) %>% 
        dplyr::select(time:id, x:y) %>% 
        dplyr::arrange(time) %>% 
        as.data.frame()    

    return(data)
}





################################################################################
# STATISTICS TO USE WITH THE MOVING WINDOW

#' Average observations within a moving window
#' 
#' @param data Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' 
#' @return Single-row dataframe with transformed columns `x` and `y`
#' 
#' @export
average <- function(data) {
    data.frame(x = mean(data$x), 
               y = mean(data$y)) %>% 
        return()
}

#' Average observations within a moving window
#' 
#' @param data Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' @param weights Either a numeric vector of weights or a function that creates 
#' the weights to be assigned to each of the data points. Defaults to 
#' `\(x) dnorm(x, mean = 0, sd = 1)`
#' @param .by String denoting the column by which the weights should be determined
#' (i.e., of which the values are passed on to `weights`). Ignored if a numeric 
#' is provided for `weights`. Defaults to `"index"`.
#' 
#' @return Single-row dataframe with transformed columns `x` and `y`
#' 
#' @export
weighted_average <- function(data, 
                             weights = \(x) dnorm(x, mean = 0, sd = 1), 
                             .by = "index") {
    
    # Create the weights and normalize them (should sum to 1). Depends on whether
    # the weights are a numeric vector or a function
    if(class(weights) == "function") {
        weights <- weights(data[, .by])
    } else if(class(weights) != "numeric") {
        stop("Weights should either by numeric or a function.")
    }
    weights <- weights / sum(weights)

    # Throw error if the weights are not of equal length to the data
    if(length(weights) != nrow(data)) {
        stop("Weights are not of same size as the data.")
    }

    # Compute the weighted means
    data.frame(x = sum(data$x * weights), 
               y = sum(data$y * weights)) %>% 
        return()
}