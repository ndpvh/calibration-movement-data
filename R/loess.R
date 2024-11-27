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
                             .by = NULL,
                             degree = 1,
                             ...) {

    # Create an internal function that will take in a single dataset (grouped by 
    # the argument provided here), do a cross-validation and then do the loess 
    # regression. Builds upon the locfit package for this.
    loess <- function(grouped_data) {
        # Use a cross-validation method to find out which value of the span you would
        # like to use for the loess. Try a whole range between 10% to 90%. Given that 
        # locfit is quite fast, we can do this quite exhaustively. (9,000,000 rows, 
        # 5% between each span; 2min 40sec).
        #
        # Here, we assume that the span-value will be the same for both x and y 
        # direction. To ensure this is the case, we bind both types of in one big 
        # data.frame under another name z.
        spans <- seq(0.1, 0.9, 0.05)
        xy_data <- data.frame(z = c(grouped_data$x, grouped_data$y), 
                              time = c(grouped_data$time, max(grouped_data$time) + grouped_data$time))
    
        fits <- sapply(spans, 
                       \(x) locfit::gcv(z ~ locfit::lp(time, 
                                                       deg = degree, 
                                                       nn = x), 
                                        data = xy_data,
                                        ...)) 
    
        # Select that span that maximizes the likelihood, as this one is most closely 
        # related to the RMSE within this package 
        span <- spans[fits[1,] == max(fits[1,])]
        if(length(span) > 1) {
            span <- mean(span)
        }

        # Now that the span has been selected for, we can finally fit the loess
        # and extract the data that we want. 
        grouped_data$x <- locfit::locfit(x ~ locfit::lp(time, 
                                                        deg = degree,
                                                        nn = span), 
                                         data = grouped_data, 
                                         ...) %>% 
            predict(locfit::lp(grouped_data$time, 
                               deg = degree,
                               nn = span))
        grouped_data$y <- locfit::locfit(y ~ locfit::lp(time, 
                                                        deg = degree,
                                                        nn = span), 
                                         data = grouped_data, 
                                         ...) %>% 
            predict(locfit::lp(grouped_data$time, 
                               deg = degree, 
                               nn = span))

        return(grouped_data)
    }

    if(is.null(.by)) {
        return(loess(data))
    } else {
        groups <- unique(data[, .by])
        filtered <- lapply(groups, 
                           \(x) loess(data[data[, .by] == x, ]))
                           
        return(do.call("rbind", filtered))
    }
}