library(tidyverse)
library(data.table)
library(nloptr)

# Load in both the data about the experiments
exp <- fread(file.path("Calibration experiments", "Data", "experiments.csv"),
             data.table = FALSE) %>% 
    # Select only a few columns and rename `id` so that you can later join
    # this dataframe with the other one
    select(id, name) %>% 
    rename(experiment_id = id)

# Load and preprocess the data about the positions
pos <- fread(file.path("Calibration experiments", "Data", "datapoints.csv"),
             data.table = FALSE) %>%
    # Column names from NUC
    setNames(c("id", "created", "modified", "tag_id", "experiment_id",
               "timestamp", "x", "y", "person_id")) %>% 
    # Join the position data with the experimental data
    plyr::join(exp, by = "experiment_id") %>% 
    # Delete and rename columns
    select(timestamp, tag_id, x:y, name) %>% 
    rename(tag = tag_id,
           experiment = name) %>% 
    # Delete data in which no experiment was running
    filter(!(is.na(experiment))) 

# Experiments that have been done are the following:
#   - test experiment: walking around
#   - my_experiment:
#   - stationary 2, 3, 5, 6: laying down the tags on the net and letting them
#                                    sit there
#   - movement L X V: moving with X Volts across the Lth line of the Y-axis
#                     (movement itself is along X-axis)  
#   - movement Lperp X V: same principle, but X and Y are inversed
#   - movement diag X V: same principle, but this time along the diagonals

###################################################
# Stationary data
###################################################

# Just visualizing the data a bit
visualize_positions <- function(name){
    plot_data <- pos %>% 
        filter(experiment == name) %>% 
        select(x, y)

    plot(plot_data$x, plot_data$y)
}

visualize_positions("stationary 2")
visualize_positions("stationary 3")
visualize_positions("stationary 5")
visualize_positions("stationary 6")

# Lets make a sum of them all and plot these as well
to_sum <- c("stationary 2", 
            "stationary 3",
            "stationary 5",
            "stationary 6")
stationary <- pos %>% 
    group_by(experiment) %>% 
    tidyr::nest() %>% 
    mutate(stationary = experiment %in% to_sum) %>% 
    tidyr::unnest(data) %>% 
    ungroup()

stationary %>% 
    filter(stationary) %>% 
    ggplot(aes(x = x, y = y)) + 
        geom_point()

########################
# Preparing the data

# Create a function that will create a dotted rectangle of a given length and 
# width. Importantly, both should be natural numbers, as there is one meter 
# apart between each dot
create_rectangle <- function(w, h){
    w <- w + 1 ; h <- h + 1
    cbind(X = rep((1:w) - 1, each = h),
          Y = rep((1:h) - 1, times = w)) %>% 
        as_tibble() %>% 
        mutate(tag = row_number()) %>% 
        return()
}

# Create a function that will classify the measurements into several discrete
# groups, based on their continuous value.
# Assumptions are that each of the tags are enough apart so that the measurement 
# error is not severe enough to not measure their locations accurately
classify_dimension <- function(x, numcat){
    (x - min(x)) %>% 
        `/` (max(x) - min(x)) %>%
        `*` (numcat) %>% 
        round() %>% 
        return()
}

# Delete the original tags and replace them with the ones created by the grid
# Furthermore assign the "real" grid coordinates to each of the tags as to be 
# able to estimate the bias and variation of the measurement
stationary <- stationary %>% 
    # Only get data from the stationary measurements and delete tags
    filter(stationary) %>% 
    select(-tag) %>% 
    # Correct the measurements to begin at the origin for the best of your
    # abilities
    mutate(x = x,
           y = y) %>%
    # Assign each of the measurements to discrete categories as defined in 
    # classify_dimension
    mutate(X = classify_dimension(x, 10),
           Y = classify_dimension(y, 8))

# Quick test of whether this works: Apparently it does! (Change X to Y for the 
# other dimension)
ggplot(data = stationary, aes(x = x, y = y, color = factor(X))) +
    geom_point()

# Join the stationary dataframe with a rectangle that is created with the same 
# specifications
stationary <- stationary %>% 
    plyr::join(create_rectangle(10, 8), 
               by = c("X", "Y"))

# Finally, correct the "real positions" by adding the minimum of x and y to it.
# While this is not going to be a perfect fit, it does give us a close approximation
# to the actual positions of X and Y, which were not measured in our first 
# calibration experiment.
stationary <- stationary %>% 
    mutate(X = X + min(x),
           Y = Y + min(y))

# Check whether it worked with some plots. It actually seems to work, and the
# relationship is quite okay
ggplot(data = stationary, aes(x = x, y = X)) +
    geom_point() +
    geom_smooth()
ggplot(data = stationary, aes(x = y, y = Y)) +
    geom_point() +
    geom_smooth()

# Also visualize the relationship between x and y to make sure there are no
# interactions happening: Otherwise 2D estimation
ggplot(data = stationary, aes(x = x, y = y)) +
    geom_point() + 
    geom_smooth()

# Save these data as being preprocessed
saveRDS(stationary,
        file.path("Calibration experiments", 
                  "Data", 
                  "preprocessed_stationary.Rds"))

########################
# Estimation: 1D

stationary <- readRDS(file.path("Calibration experiments", 
                                "Data", 
                                "preprocessed_stationary.Rds"))

# Create a function that will output a formula for a polynomial of the pwr^th
# degree with y as dependent variable and x as the independent variable
polynomial_formula <- function(x, y, pwr){
    polyn <- paste(y, "~", x)
    for(i in 2:pwr){
        polyn <- paste(polyn, "+", "I(",
                       paste0(x, "^", i), ")")
    }
    polyn %>% 
        formula() %>% 
        return()
}

# Estimate the 10th degree polynomial on the data and find out what the characteristic
# function of the measurements are
result_x <- lm(polynomial_formula("x", "X", 10),
               data = stationary) %>% 
    summary()
result_y <- lm(polynomial_formula("y", "Y", 10),
               data = stationary) %>% 
    summary()

# Save these results for later
saveRDS(result_x,
        file.path("Calibration experiments", 
                  "Results",
                  "polynomial_for_x.Rds"))
saveRDS(result_y,
        file.path("Calibration experiments", 
                  "Results",
                  "polynomial_for_y.Rds"))

########################
# Estimation: 2D

stationary <- readRDS(file.path("Calibration experiments", 
                                "Data", 
                                "preprocessed_stationary.Rds"))

# Create a function that will create the data matrix
data_matrix <- function(x, y, pwr, print_formula = FALSE){
    X <- matrix(0, nrow = length(x), ncol = (pwr + 1)^2)
    f <- 1
    form <- ""
    for(i in 0:pwr){
        for(j in 0:pwr){
            # If the total exponents in the formula exceeds the order of the 
            # polynomial, then you have to skip this one
            if(i + j > pwr){
                next
            }

            # Add the interaction term of x and y to the polynomial
            X[,f] <- x^i * y^j
            f <- f + 1

            # Add the function to the string of the formula
            form <- paste0(form, 
                           ifelse(i + j == 0, "", " + "),
                           "x^", i, " y^", j)
        }
    }
    # Delete the columns in which there are only 0's
    X <- X[,1:(f - 1)]

    # Print formula if requested
    if(print_formula){
        print(form)
    }
    return(X)
}

fit_polynomial <- function(x, y, pwr, 
                           deterministic = TRUE,        # Whether analytic formula of least-squares should be used
                           analytic_inverse = TRUE){    # Whether inverse should be computed analytically or whether a pseudoinverse can be used
    # Unlist x and y
    x <- as.matrix(x)
    y <- as.matrix(y)

    # Do the least-squares
    if(deterministic){
        # Create the data matrix based on x and y
        xy <- data_matrix(x[,1], x[,2], pwr)

        # Get beta through the formula of the least-squares method
        if(analytic_inverse){
            beta <- solve( t(xy) %*% xy ) %*% t(xy) %*% y
        } else {
            beta <- pseudoinverse( t(xy) %*% xy ) %*% t(xy) %*% y
        }
    } else {
        # THIS IS A DEAD END: YOU CANNOT TURN BACK TO YOUR ORIGINAL MEASUREMENTS
        # IF YOU WERE TO ESTIMATE THE POLYNOMIAL AND EVALUATE IT ON ANOTHER 
        # DATA MATRIX. MAYBE THIS CAN BE DONE IN ANOTHER WAY, BUT IT IS 
        # NOT THE WAY IN WHICH I WAS THINKING ABOUT IT.

        # Get the data matrix of the real positions
        # XY <- data_matrix(y[,1], y[,2], pwr)

        # # Define the function to be minimized
        # to_minimize <- function(x){
        #     Y_hat <- (xy %*% x)

        #     print(dim(Y_hat))
        #     print(dim(y))
        #     (y - Y_hat)^2 %>% 
        #         sum() %>% 
        #         return()
        # }

        # # Define the nloptr algorithm options and run the minimization with
        # # only one initial condition (everything is 0)
        # results <- nloptr(x0 = numeric(length = ncol(XY)),
        #                   eval_f = to_minimize,
        #                   lb = rep(-Inf, each = ncol(XY)),
        #                   ub = rep(Inf, each = ncol(XY)),
        #                   opts = list("algorithm" = "NLOPT_GN_DIRECT", 
        #                               "ftol_rel" = 10^(-15),
        #                               "xtol_rel" = 10^(-15)))

        # beta <- result$solution %>% 
        #     matrix(ncol = 2)
    }
    return(beta)
}

# Fit the measured coordinates to the real coordinates with a 4th degree 
# polynomial (after 5 becomes singular)
result <- fit_polynomial(stationary[,c("x", "y")], 
                         stationary[,c("X", "Y")], 
                         10)

# Save these results for later
saveRDS(result,
        file.path("Calibration experiments", 
                  "Results",
                  "polynomial_for_2D.Rds"))

########################
# Checking how far you get: 1D

stationary <- readRDS(file.path("Calibration experiments", 
                                "Data", 
                                "preprocessed_stationary.Rds"))
result_x <- readRDS(file.path("Calibration experiments", 
                              "Results",
                              "polynomial_for_x.Rds"))
result_y <- readRDS(file.path("Calibration experiments", 
                              "Results",
                              "polynomial_for_y.Rds"))

# Create a function that will extract all of the coefficients
extract_coefficients <- function(x){
    x$coefficients[,1] %>% 
        as.numeric() %>% 
        return()
}

# Create a function that will correct the measurements that it
# receives
correct_measurements <- function(x, parameters){
    X <- matrix(0, nrow = length(x), ncol = length(parameters))
    for(i in seq_along(parameters)){
        X[,i] <- parameters[i] * x^(i - 1)
    }
    X %>% 
        rowSums() %>% 
        as.numeric() %>% 
        return()
}

# Do the actual correction in the x and y dimensions separately
new_x <- correct_measurements(stationary$x, 
                              extract_coefficients(result_x))
new_y <- correct_measurements(stationary$y, 
                              extract_coefficients(result_y))

# Plot the new coordinates
plot(new_x, new_y)
plot(stationary$x, stationary$y)
plot(stationary$X, stationary$Y)

########################
# Checking how far you get: 2D

stationary <- readRDS(file.path("Calibration experiments", 
                                "Data", 
                                "preprocessed_stationary.Rds"))
result <- readRDS(file.path("Calibration experiments", 
                            "Results",
                            "polynomial_for_2D.Rds"))

# Create a function that will correct the measurements that it
# receives
correct_measurements <- function(x, parameters, pwr){
    # Unlist x: Otherwise error in data_matrix
    x <- as.matrix(x)

    # Create the data matrix based on x and y. Unfortunately, difficult to extract
    # the order of the polynomial from only its parameters, which is why it has
    # to be provided as an argument
    XY <- data_matrix(x[,1], x[,2], pwr)

    # Do the actual conversion of the measured positions to the newly estimated
    # ones
    ( XY %*% parameters ) %>% 
        return()
}

# Do the actual correction in the x and y dimensions separately
new_xy <- correct_measurements(stationary[,c("x", "y")], 
                               result,
                               10)

# Plot the new coordinates
plot(new_xy[,1], new_xy[,2])
plot(stationary$x, stationary$y)
plot(stationary$X, stationary$Y)

########################
# Estimation: Splines

stationary <- readRDS(file.path("Calibration experiments", 
                                "Data", 
                                "preprocessed_stationary.Rds"))

# Do the splines estimations separately for x and y
# This is a cubic spline, as defined by the function `ns`, and the knots being 
# the typical interquartal quantiles
nbins <- 2
quants <- seq(0 + 0.5/nbins, 1 - 0.5/nbins, 1/nbins)
result_x <- lm(data = stationary, 
               X ~ ns(x, knots = quantile(x, quants)))
result_y <- lm(data = stationary, 
               Y ~ ns(y, knots = quantile(y, quants)))

# Save these results for later
saveRDS(result_x,
        file.path("Calibration experiments", 
                  "Results",
                  "splines_for_x.Rds"))
saveRDS(result_y,
        file.path("Calibration experiments", 
                  "Results",
                  "splines_for_y.Rds"))

# Read them back in (just a formality: Keeps you from having to run everything)
result_x <- readRDS(file.path("Calibration experiments", 
                              "Results",
                              "splines_for_x.Rds"))
result_y <- readRDS(file.path("Calibration experiments", 
                              "Results",
                              "splines_for_x.Rds"))

# Use the spline fits to create new values for the measured x and y, which 
# should hopefully correspond to the actual X and Y
X_recovered <- predict(result_x, newdata = list(x = stationary$x)) %>% 
    as.numeric()
Y_recovered <- predict(result_y, newdata = list(x = stationary$y)) %>% 
    as.numeric()

# Check the correspondence between the actual and the recovered positions
plot(X_recovered, stationary$X)
plot(Y_recovered, stationary$Y) 

plot(X_recovered, stationary$x)
plot(Y_recovered, stationary$y)

plot(stationary$x, stationary$X)
plot(stationary$y, stationary$Y)

# And check the grid itself
stationary %>% 
    ungroup() %>% 
    mutate(X_recovered = X_recovered,
           Y_recovered = Y_recovered) %>% 
    ggplot(aes(x = X_recovered, y = Y_recovered)) +
        geom_point()

########################
# Estimation: Linear regression

stationary <- readRDS(file.path("Calibration experiments", 
                                "Data", 
                                "preprocessed_stationary.Rds"))

# Do the splines estimations separately for x and y
# This is a cubic spline, as defined by the function `ns`, and the knots being 
# the typical interquartal quantiles
result_x <- lm(data = stationary, 
               X ~ x)
result_y <- lm(data = stationary, 
               Y ~ y)

# Save these results for later
saveRDS(result_x,
        file.path("Calibration experiments", 
                  "Results",
                  "linear_for_x.Rds"))
saveRDS(result_y,
        file.path("Calibration experiments", 
                  "Results",
                  "linear_for_y.Rds"))

# Read them back in (just a formality: Keeps you from having to run everything)
result_x <- readRDS(file.path("Calibration experiments", 
                              "Results",
                              "linear_for_x.Rds"))
result_y <- readRDS(file.path("Calibration experiments", 
                              "Results",
                              "linear_for_x.Rds"))

# Use the spline fits to create new values for the measured x and y, which 
# should hopefully correspond to the actual X and Y
X_recovered <- predict(result_x, newdata = list(x = stationary$x)) %>% 
    as.numeric()
Y_recovered <- predict(result_y, newdata = list(x = stationary$y)) %>% 
    as.numeric()

# Check the correspondence between the actual and the recovered positions
plot(X_recovered, stationary$X)
plot(Y_recovered, stationary$Y) 

plot(X_recovered, stationary$x)
plot(Y_recovered, stationary$y)

plot(stationary$x, stationary$X)
plot(stationary$y, stationary$Y)

# And check the grid itself
stationary %>% 
    ungroup() %>% 
    mutate(X_recovered = X_recovered,
           Y_recovered = Y_recovered) %>% 
    ggplot(aes(x = X_recovered, y = Y_recovered)) +
        geom_point()

