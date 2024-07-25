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
                             .by = "index", 
                             cols = NULL) {
    
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

#' Take the middle observation within a moving window
#' 
#' @param data Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' 
#' @return Single-row dataframe with transformed columns `x` and `y`
#' 
#' @export
middle <- function(data, 
                   cols = NULL) {

    # Get the middle of the time-bin
    data <- data %>% 
        dplyr::arrange(time) %>% 
        dplyr::mutate(index = dplyr::row_number())

    idx <- median(data$index)

    # Dispatch on whether the index is a whole number or not. If it is a whole 
    # number, we simple take the middle x and y. If it is not a whole number, 
    # then we take the mean of the two x's and y's that make up the middle of 
    # the interval
    if(idx %in% data$index) {
        result <- data.frame(x = data$x[idx], 
                             y = data$y[idx])
    } else {
        result <- data.frame(x = mean(data$x[c(floor(idx), ceiling(idx))]), 
                             y = mean(data$y[c(floor(idx), ceiling(idx))]))
    }

    # Do the same for any other variable that it is asked for
    if(!is.null(cols)) {
        for(i in cols) {
            if(idx %in% data$index) {
                result[,i] <- data[idx,i]
            } else {
                result[,i] <- mean(data[c(floor(idx), ceiling(idx)), i])
            }
        }
    }

    return(result)
}

#' Approximate positions with linear regression
#' 
#' Estimate a linear regression in both the x- and the y-direction (i.e., 
#' with y ~ x and x ~ y) and take the best fitting one. We then approximate the 
#' real position by filling out the observed x and y in the equation and 
#' averaging observed (x, y) coordinates with the predicted ones. This should 
#' (hopefully) decrease the error in both dimensions instead of keeping the 
#' error in one dimension to predict the location in the other.
#' 
#' @param data Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' Additionally, the variable `"index"` should be provided (automatically the 
#' case when using `moving_window`), as this will tell us which of the 
#' observations to use for our predictions.
#' 
#' @return Single-row dataframe with transformed columns `x` and `y`
#' 
#' @export
linear <- function(data, 
                   cols = NULL) {

    rotated <- FALSE

    # Do two regressions and check which one has the highest r-squared
    lm_y <- lm(y ~ x, data = data)
    lm_x <- lm(x ~ y, data = data)

    # Get the index of the data for which the `index == 0`: These coordinates 
    # needed
    idx <- data$index == 0
    co <- as.numeric(data[idx, c("x", "y")])

    # Special case: vertical or horizontal predictivity. If we find such a case
    # We will first try to rotate the data and try again
    if(is.na(lm_y$coefficients[2]) | is.na(lm_x$coefficients[2])) {
        rotated <- TRUE 

        # Rotate the data with 45 degrees
        R <- matrix(c(cos(pi / 4), sin(pi / 4), -sin(pi / 4), cos(pi / 4)), 
                    nrow = 2, 
                    ncol = 2)
        rotated_data <- (R %*% t(as.matrix(data[, c("x", "y")]))) %>% 
            t() %>% 
            as.data.frame() %>% 
            setNames(c("x", "y"))

        # Now that the data have been rotated, do the regressions again
        lm_y <- lm(y ~ x, data = rotated_data)
        lm_x <- lm(x ~ y, data = rotated_data)

        # Replace the coordinate with the rotated one
        co <- as.numeric(rotated_data[idx, c("x", "y")])
    }

    # Generate predictions for the (x, y) coordinates
    if(summary(lm_y)$r.squared >= summary(lm_x)$r.squared) {
        # Extract coefficients
        b0 <- lm_y$coefficients[1]
        b1 <- lm_y$coefficients[2]

        # Create predicted (x, y) coordinates based on the observed coordinates
        y_hat <- as.numeric(b0 + b1 * co[1])
        x_hat <- as.numeric((co[2] - b0) / b1)
    } else {
        # Extract coefficients
        b0 <- lm_x$coefficients[1]
        b1 <- lm_x$coefficients[2]

        # Create predicted (x, y) coordinates based on the observed coordinates
        y_hat <- as.numeric((co[1] - b0) / b1)
        x_hat <- as.numeric(b0 + b1 * co[2])
    }

    # If the data were rotated, rotate them back now and replace the predicted 
    # values with these rotated ones
    if(rotated) {
        # Rotate back
        R <- matrix(c(cos(-pi / 4), sin(-pi / 4), -sin(-pi / 4), cos(-pi / 4)), 
                    nrow = 2, 
                    ncol = 2)
        rotated_co <- R %*% matrix(c(x_hat, y_hat), nrow = 2, ncol = 1)

        # Replace the values
        x_hat <- rotated_co[1]
        y_hat <- rotated_co[2]

        # Get the real coordinates back
        co <- as.numeric(data[idx, c("x", "y")])
    }
    
    # Compute the average of the observed and predicted coordinates. If the 
    # data are not linearly related (even after rotation), then just return 
    # the raw coordinates
    if(is.na(lm_y$coefficients[2]) | is.na(lm_x$coefficients[2]) | 
       is.na(x_hat) | is.na(y_hat) |
       round(lm_y$coefficients[2], 10) == 0 | round(lm_x$coefficients[2], 10) == 0) {
        result <- data.frame(x = co[1], 
                             y = co[2])  
    } else {
        result <- data.frame(x = mean(c(co[1], x_hat)), 
                             y = mean(c(co[2], y_hat)))
    }

    # Use the `index == 0` values of the data for each of the columns that 
    # should be exported together with these coordinates
    if(!is.null(cols)) {
        for(i in cols) {
            result[,i] <- data[data$index == 0, i]
        }
    }

    return(result)
}

#' Approximate positions with a parabola
#' 
#' Estimate a linear regression in both the x- and the y-direction (i.e., 
#' with y ~ x and x ~ y) and take the best fitting one. We then approximate the 
#' real position by filling out the observed x and y in the equation and 
#' averaging observed (x, y) coordinates with the predicted ones. This should 
#' (hopefully) decrease the error in both dimensions instead of keeping the 
#' error in one dimension to predict the location in the other.
#' 
#' @param data Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' Additionally, the variable `"index"` should be provided (automatically the 
#' case when using `moving_window`), as this will tell us which of the 
#' observations to use for our predictions.
#' 
#' @return Single-row dataframe with transformed columns `x` and `y`
#' 
#' @export
parabola <- function(data, 
                     cols = NULL) {

    rotated <- FALSE

    # Do two regressions and check which one has the highest r-squared
    lm_y <- lm(y ~ x + I(x^2), data = data)
    lm_x <- lm(x ~ y + I(y^2), data = data)

    # Select the (x, y) coordinate at index 0
    idx <- data$index == 0
    co <- as.numeric(data[idx, c("x", "y")])

    # Special case: vertical or horizontal predictivity. If we find such a case
    # We will first try to rotate the data and try again
    if((is.na(lm_y$coefficients[3]) | is.na(lm_x$coefficients[3])) &
       (is.na(lm_y$coefficients[2]) | is.na(lm_x$coefficients[2]))) {
        rotated <- TRUE 

        # Rotate the data with 45 degrees
        R <- matrix(c(cos(pi / 4), sin(pi / 4), -sin(pi / 4), cos(pi / 4)), 
                    nrow = 2, 
                    ncol = 2)
        rotated_data <- (R %*% t(as.matrix(data[, c("x", "y")]))) %>% 
            t() %>% 
            as.data.frame() %>% 
            setNames(c("x", "y"))

        # Now that the data have been rotated, do the regressions again
        lm_y <- lm(y ~ x + I(x^2), data = rotated_data)
        lm_x <- lm(x ~ y + I(y^2), data = rotated_data)

        # Replace the coordinate with the rotated one
        co <- as.numeric(rotated_data[idx, c("x", "y")])
    }

    # Generate predictions for the (x, y) coordinates
    if(summary(lm_y)$r.squared >= summary(lm_x)$r.squared) {
        # Extract coefficients
        b0 <- lm_y$coefficients[1]
        b1 <- lm_y$coefficients[2]
        b2 <- lm_y$coefficients[3]

        # Special case: If the approximation yields no parabola, but rather a linear 
        # line, pass it on to `linear`. This is done so that the functions are a bit
        # more sparse, instead of building all potential cases into this one
        if((round(b2, 10) == 0) | is.na(b2)) {
            return(linear(data, cols = cols))
        }

        # Create predicted (x, y) coordinates based on the observed coordinates.
        # For the inverse one, a small trick is needed: You subtract y from the 
        # complete equation and then find the roots of this equation. This 
        # ensures at least one solution for x. If there are two solutions, pick 
        # the one that is closest to the actual observation of x
        y_hat <- b0 + b1 * co[1] + b2 * co[1]^2

        D <- b1^2 - 4 * b2 * (b0 - co[2])
        x_hat <- c((-b1 + sqrt(D)) / (2 * b2), 
                   (-b1 - sqrt(D)) / (2 * b2))
        x_hat <- x_hat[which.min((x_hat - co[1])^2)]
    } else {
        # Extract coefficients
        b0 <- lm_x$coefficients[1]
        b1 <- lm_x$coefficients[2]
        b2 <- lm_y$coefficients[3]

        # Special case: If the approximation yields no parabola, but rather a linear 
        # line, pass it on to `linear`. This is done so that the functions are a bit
        # more sparse, instead of building all potential cases into this one
        if((round(b2, 10) == 0) | is.na(b2)) {
            return(linear(data, cols = cols))
        }

        # Create predicted (x, y) coordinates based on the observed coordinates, 
        # using a similar approach as above
        D <- b1^2 - 4 * b2 * (b0 - co[1])
        y_hat <- c((-b1 + sqrt(D)) / (2 * b2), 
                   (-b1 - sqrt(D)) / (2 * b2))
        y_hat <- y_hat[which.min((y_hat - co[2])^2)]

        x_hat <- b0 + b1 * co[2] + b2 * co[2]^2
    }

    # If the data were rotated, rotate them back now and replace the predicted 
    # values with these rotated ones
    if(rotated) {
        # Rotate back
        R <- matrix(c(cos(-pi / 4), sin(-pi / 4), -sin(-pi / 4), cos(-pi / 4)), 
                    nrow = 2, 
                    ncol = 2)
        rotated_co <- R %*% matrix(c(x_hat, y_hat), nrow = 2, ncol = 1)

        # Replace the values
        x_hat <- rotated_co[1]
        y_hat <- rotated_co[2]

        # Get the real coordinates back
        co <- as.numeric(data[idx, c("x", "y")])
    }
    
    # Compute the average of the observed and predicted coordinates. If the 
    # data are not linearly related (even after rotation), then just return 
    # the raw coordinates
    if(is.na(x_hat) | is.na(y_hat)) {
        result <- data.frame(x = co[1], 
                             y = co[2])  
    } else {
        result <- data.frame(x = mean(c(co[1], x_hat)), 
                             y = mean(c(co[2], y_hat)))
    }

    # Use the `index == 0` values of the data for each of the columns that 
    # should be exported together with these coordinates
    if(!is.null(cols)) {
        for(i in cols) {
            result[,i] <- data[data$index == 0, i]
        }
    }

    return(result)
}