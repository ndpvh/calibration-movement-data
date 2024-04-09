#' Create a perfect rectangle
#' 
#' A function that creates locations on within a rectangle. Can be used to 
#' create a perfect location-wise rectangle on which the stationary data 
#' can be estimated.
#' 
#' @param size A numeric vector containing the width and height of the 
#' rectangle.
#' @param delta_x A numeric denoting the space inbetween the points of the 
#' rectangle. Defaults to `1`.
#' @param number_tags A boolean denoting whether each position should receive
#' a unique number. Can later be used to match measured tags to them. Defaults
#' to `TRUE`.
#' 
#' @export
rectangle <- function(size,
                      delta_x = 1,
                      number_tags = TRUE){
    # Create X and Y series for the given size and given space inbetween each
    # point
    X <- seq(0, size[1], delta_x)
    Y <- seq(0, size[2], delta_x)

    # Now bind them together in a single tibble
    XY <- cbind(X = rep(X, each = length(Y)),
                Y = rep(Y, times = length(X))) %>% 
        as_tibble()

    # If needed, create tag numbers that can be used if needed.
    if(number_tags){
        XY <- XY %>% 
            mutate(tag = row_number())
    }
    
    return(XY)
}
