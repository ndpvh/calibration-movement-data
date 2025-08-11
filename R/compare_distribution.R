#' Compare two distributions
#' 
#' This function will check whether two distributions differ significantly. Done
#' by first taking the difference in the distribution of the values provided in
#' \code{x} and \code{y}, then defining the (1 - \code{alpha}) * 100% confidence
#' interval of this distributions, and finally checking whether 0 is a part of 
#' that interval. 
#' 
#' If 0 is part of this interval, the two distributions are said to not differ 
#' significantly, and vice-versa if 0 is not a part of this interval. 
#'
#' This function therefore allows an analysis analogous to a posterior check, or 
#' otherwise absolute goodness-of-fit with a parametric bootstrap procedure.
#' 
#' @param x Numeric vector of test values.
#' @param y Numeric vector of reference values.
#' @param alpha Numeric denoting the signficance level.
#' @param bootstrapped Integer denoting the number of samples to have in a 
#' bootstrap. If not defined, confidencen intervals will be selected based on raw
#' difference scores. Defaults to \code{NA}
#' 
#' @return Single-rowed data.frame containing the separate confidence intervals 
#' for \code{x} and \code{y}, as well as the confidence interval for the difference 
#' \code{x - y} and a logical denoting whether 0 is contained within the latter
#' interval.
#' 
#' @export
compare_distribution <- function(x, 
                                 y,
                                 alpha = 0.05,
                                 bootstrapped = NA) {
    
    bounds <- c(alpha/2, 0.5, 1 - alpha/2)
    diff <- x - y

    if(!is.na(bootstrapped)) {
        N <- length(x)

        idx <- sample(1:N, N * bootstrapped, replace = TRUE) 

        x <- x[idx] %>% 
            matrix(nrow = N, ncol = bootstrapped) %>% 
            colMeans()
        y <- y[idx] %>% 
            matrix(nrow = N, ncol = bootstrapped) %>% 
            colMeans()
        diff <- diff[idx] %>% 
            matrix(nrow = N, ncol = bootstrapped) %>% 
            colMeans()
    }

    ci_x <- quantile(x, probs = bounds)
    ci_y <- quantile(y, probs = bounds)
    ci <- quantile(diff, probs = bounds)

    return(data.frame(mean_x = mean(x, na.rm = TRUE),
                      median_x = ci_x[2],
                      ci_x_lower = ci_x[1], 
                      ci_x_upper = ci_x[3],
                      mean_y = mean(y, na.rm = TRUE),
                      median_y = ci_y[2],
                      ci_y_lower = ci_y[1], 
                      ci_y_upper = ci_y[3],
                      mean_diff = mean(x - y, na.rm = TRUE),
                      median_diff = ci[2],
                      ci_diff_lower = ci[1],
                      ci_diff_upper = ci[3],
                      significant = (0 < ci[1]) | 0 > ci[3]))
}