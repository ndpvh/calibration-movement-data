#' Filter data through LOESS
#' 
#' Uses the \code{loess} function in R to filter the data. Is a generalization 
#' and improvement of \code{linear} and \code{parabola} used with the 
#' \code{moving_window}.
#' 
#' @param data Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' @param degree Integer denoting the degree of the polynomial. 
#' 
#' @return Single-row dataframe with transformed columns `x` and `y`
#' 
#' @export
local_regression <- function(data, 
                             ...) {
    
    # Perform a loess of a given degree
    result <- loess(formula = y ~ x, 
                    data = data, 
                    ...)

    # Once done, we can replace the results of the data with these results
    result <- data %>%
        dplyr::mutate(x = result$x, 
                      y = result$y)

    return(result)
}