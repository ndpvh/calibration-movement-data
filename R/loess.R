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
                             .by = "id",
                             span_obs = NULL,
                             span = 0.75,
                             ...) {

    # If `span_obs` is defined, we want to change the `span` argument to account
    # for the number of observations within each window
    if(!is.null(span_obs)) {
        span <- round(span_obs / nrow(data), digits = 4)
    }

    # Do a regression for each value of .by
    data %>% 
        dplyr::group_by_at(.by) %>% 
        dplyr::mutate(x = loess(formula = x ~ time, 
                                span = span,
                                ...) %>% 
                          predict(), 
                      y = loess(formula = y ~ time, 
                                span = span,
                                ...) %>% 
                          predict()) %>% 
        dplyr::ungroup() %>% 
        return()
}