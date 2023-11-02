#' Correct distorted measurements
#' 
#' This function will take in the measured positions for a given dimension and 
#' correct these measurements based on the parameters of the parameters of 
#' a polynomial. Returns the corrected measurements.
#' 
#' @param x A vector or matrix of measurements. If a matrix, the columns should 
#' represent the measured dimensions, while rows contain the actual measurements.
#' @param parameters A vector or matrix of parameters of the polynomial. The
#' format of the parameters depends on the format of `x`, as both should be in 
#' the same format. If a matrix, the columns should therefore again represent 
#' the measured dimensions. Parameters should be ordered according to the 
#' exponent of the polynomial that they correspond to: you start with the
#' intercept (x^0), then with x (x^1), then x^2, and so on.
#' 
#' @export
correct_distortion <- function(x, parameters){
    # Dispatch based on the format of `x`
    if(ncol(x) > 1){
        # Convert `x` to a matrix 
        x <- as.matrix(x)

        # If a matrix, then we should loop over the columns to correct each 
        # dimension separately
        for(i in seq_len(ncol(x))){
            x[,i] <- correct_single_dimension(x[,i], parameters[,i])
        }
    } else {
        # If a vector, we can just correct this one dimension
        x <- correct_single_dimension(x, parameters)
    }
    return(x)
}

# This lower-level function will do the actual correcting
correct_single_dimension <- function(x, parameters){
    # Reserve memory for the different terms of the polynomial
    X <- matrix(0, nrow = length(x), ncol = length(parameters))

    # Loop over all terms in the polynomial and compute the value for `x` for 
    # each of these
    for(i in seq_along(parameters)){
        X[,i] <- parameters[i] * x^(i - 1)
    }

    # Sum all terms and convert to a numeric
    X %>% 
        rowSums() %>% 
        as.numeric() %>% 
        return()
}
