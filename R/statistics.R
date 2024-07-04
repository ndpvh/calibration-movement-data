#' Average observations within a moving window
#' 
#' @param data Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' 
#' @return Single-row dataframe with transformed columns `x` and `y`
#' 
#' @export
average <- function(data, 
                    cols = NULL) {
    
    # Compute the average for x and y (minimal requirement)
    result <- data.frame(x = mean(data$x), 
                         y = mean(data$y))

    # Do the same for any other variable that is asked for
    if(!is.null(cols)) {
        for(i in cols) {
            result[,i] <- mean(data[,i])
        }
    }

    return(result)
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

    # Compute the weighted means for x and y
    result <- data.frame(x = sum(data$x * weights), 
                         y = sum(data$y * weights))

    # If the weighted means should be computed for other columns, do so
    if(!is.null(cols)) {
        for(i in cols) {
            result[,i] <- sum(data[,i] * weights)
        }
    }

    return(result)
}