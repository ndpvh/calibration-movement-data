#' Transform data to -1 and 1
#' 
#' Transform the delivered data to the -1 and 1 range based on the provided 
#' minimal and maximal value that the data can take.
#' 
#' @param x Vector of numerical data to be transformed
#' @param min_x Numeric denoting the minimal value that the data `x` can take.
#' Defaults to the minimal observed value in `x`.
#' @param max_x Numeric denoting the maximal value that the data `x` can take.
#' Defaults to the maximal observed value in `x`.
#' 
#' @export 
minmax_standardize <- function(x, min_x = NULL, max_x = NULL){
    # If `min_x` and/or `max_x` are not provided, default to the minimal and/or
    # maximal values of `x`
    if(is.null(min_x)){
        min_x <- min(x)
    }

    if(is.null(max_x)){
        max_x <- max(x)
    }

    # Once done, you can start standardizing the data
    y <- 2 * (x - min_x)/(max_x - min_x) - 1
    return(y)
}

#' Transform data from -1 and 1 to normal range
#' 
#' Transform the delivered data from the -1 and 1 range back to their original 
#' range, based on the provided minimal and maximal value that the data can take.
#' 
#' @param x Vector of numerical data to be transformed back to their normal range
#' @param min_x Numeric denoting the minimal value that the data originally 
#' could take. Has no default and should therefore be provided in order for the 
#' transformation to happen.
#' @param max_x Numeric denoting the maximal value that the data originally 
#' could take. Has no default and should therefore be provided in order for the 
#' transformation to happen.
#' 
#' @export
minmax_destandardize <- function(x, min_x = NULL, max_x = NULL){
    # If `min_x` or `max_x` is not provided, just return `x` with the warning 
    # that the data were not destandardized.
    if(is.null(min_x) || is.null(max_x)){
        paste0("`min_x` or `max_x` is not provided to the `minmax_destandardize` ", 
               "function. Without these arguments, the requested transformation ", 
               "is not possible. Provided data is returned as is.") %>% 
            warnings()
        return(x)
    }

    # If `min_x` and `max_x` are known, we can do the opposite transformation 
    # than the one defined in `minmax_standardize`
    y <- 0.5 * (x + 1) * (max_x - min_x) + min_x 
    return(y)
}
