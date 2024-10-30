#' Smooth using a Kalman filter
#' 
#' This is a higher-level function that will first look at whether the data needs
#' to be processed by a given variable. Then it will run the Kalman filter on 
#' these grouped data.
#' 
#' @param data Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' @param model String denoting which model to use.
#' @param .by String denoting whether the moving window should be taken with 
#' respect to a given grouping variable. Defaults to `NULL`.
#' 
#' @return Smoothed dataframe with a similar structure as `data`
#' 
#' @export 
kalman_filter <- function(data,     
                          reverse = TRUE,
                          internal = FALSE,
                          model = "constant_velocity", 
                          .by = NULL,
                          assumed_variance = 0.031^2,
                          check = FALSE) {
    
    # Dispatch on whether to group the data by a given variable or not
    if(is.null(.by)) {
        return(kalman_filter_individual(data, 
                                        reverse = reverse, 
                                        model = model, 
                                        check = check, 
                                        assumed_variance = assumed_variance,
                                        internal = internal))
    } else {
        data %>% 
            dplyr::group_by(.dots = .by) %>% 
            tidyr::nest() %>% 
            dplyr::mutate(data = data %>% 
                as.data.frame() %>% 
                kalman_filter_individual(reverse = reverse, 
                                         model = model, 
                                         check = check,
                                         assumed_variance = assumed_variance, 
                                         internal = internal) %>% 
                list()) %>% 
            tidyr::unnest(data) %>% 
            dplyr::ungroup() %>% 
            as.data.frame() %>% 
            return()
    }
}

#' Use Kalman filter on the data 
#' 
#' @param data Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' @param model String denoting which model to use.
#' 
#' @return Smoothed dataframe with a similar structure as `data`
#' 
#' @export 
kalman_filter_individual <- function(data, 
                                     reverse = TRUE,
                                     internal = FALSE,
                                     model = "constant_velocity", 
                                     assumed_variance = 0.031^2, 
                                     check = FALSE) {

    # Robustness against too little data. When there was only 1 row, errors arose
    if(nrow(data) <= 5) {
        return(data)
    }

    # Get the model parameters and initial conditions
    parameters <- kalman_models[[model]](data, 
                                         reverse = reverse, 
                                         internal = internal,
                                         assumed_variance = assumed_variance)

    # If you want to use another package for the estimation. Ideally, both methods 
    # would converge on the same results, but if not, then I should figure out
    # why they diverge
    if(!internal) {
        # Extract the data to be explained and the columns to be changed
        y <- parameters[["y"]] %>% 
            dplyr::select(x, y) %>% 
            as.matrix() %>% 
            t()

        cols <- parameters[["cols_of_interest"]]

        # Delete them from the list
        idx <- !(names(parameters) %in% c("y", "cols_of_interest"))
        parameters <- parameters[idx]

        # Do the estimation
        smoothed_y <- kalmanfilter::kalman_filter(parameters, 
                                                  y, 
                                                  smooth = reverse)

        # Adjust the data
        data[, cols] <- t(smoothed_y$y_tt)


    # If you want to use the internal mechanisms, then we should proceed in the 
    # normal way.
    } else {
        # Extract some of the more useful parameters, namely data and initial 
        # conditions
        y <- parameters[["y"]]
        x0 <- parameters[["x0"]]
        F0 <- parameters[["F0"]]
    
        cols <- parameters[["cols_of_interest"]]
    
        # Make a copy of y that will contain the smoothed data. Futhermore create 
        # lists that will hold the estimation covariance and the Kalman gain at 
        # each iteration
        smoothed_y <- y
        P <- list()
        K <- list()
        z <- list()
        
        # Iterate over the data to smooth it
        for(i in seq_len(nrow(y))) {
            # Perform the three steps of the Kalman filter: 
            #   (a) kf_predict: Predict the next time step t + 1
            #   (b) kf_innovation: Compute the Kalman gain
            #   (c) kf_update: Update the initial prediction with the measurement
            #
            # Check whether this is the first iteration. If this is the case, we 
            # should use the initial (prior) guess as prediction. If not, then we 
            # can use the previously predicted values to create new predictions
            if(i == 1) {
                prediction <- list("x" = x0, 
                                   "F" = F0)
            } else {
                # Extract and create the parameters that change each iteration
                # (i.e., those that depend on the time that has passed)
                A <- parameters[["A"]](y$Delta_t[i])
                W <- parameters[["W"]](y$Delta_t[i])
    
                # Do the prediction
                prediction <- kf_predict(x0, 
                                         A, 
                                         parameters[["u"]][i], 
                                         parameters[["B"]], 
                                         W, 
                                         F0)
            }
    
            innovation <- kf_innovation(matrix(as.numeric(y[i, cols]), ncol = 1),
                                        prediction[["x"]],
                                        parameters[["H"]],
                                        parameters[["V"]],
                                        prediction[["F"]])
    
            result <- kf_update(prediction[["x"]], 
                                innovation[["z"]], 
                                parameters[["H"]],
                                parameters[["V"]],
                                prediction[["F"]],
                                innovation[["K"]])
    
            # Save the results in the smoothed dataset and in the list of 
            # estimation uncertainties
            smoothed_y[i, cols] <- result[["x"]][which(cols %in% c("x", "y"))]
            P[[i]] <- t(result[["F"]]) %*% result[["F"]]
            K[[i]] <- innovation[["K"]]
            z[[i]] <- innovation[["z"]]
    
            # Overwrite the initial conditions with the newly acquired values
            x0 <- result[["x"]]
            F0 <- result[["F"]]
        }
    
        # If you want to check the autocorrelation assumptions in the innovations, 
        # print out the results
        if(check) {
            idx <- y$original
            z <- do.call("cbind", z[idx])
    
            message(paste0("Correlations between the innovations are ", 
                           cor(z[1, 2:ncol(z) - 1], z[1, 2:ncol(z)]), 
                           " and ", 
                           cor(z[2, 2:ncol(z) - 1], z[2, 2:ncol(z)])))
    
            dist_1 <- data %>% 
                dplyr::mutate(dist = sqrt((x_original - x)^2 + (y_original - y)^2)) %>% 
                dplyr::select(dist) %>% 
                unlist() %>% 
                as.numeric()
    
            idx <- smoothed_y$original
            dist_2 <- data %>% 
                dplyr::mutate(dist = sqrt((x_original - smoothed_y$x[idx])^2 + (y_original - smoothed_y$y[idx])^2)) %>% 
                dplyr::select(dist) %>% 
                unlist() %>% 
                as.numeric()
    
            print(sd(dist_1))
            print(sd(dist_2))
            browser()
        }
    
        # If you reversed the data, delete the reversed data and only keep the new 
        # (smoothed) values for the original ones
        if(reverse) {
            smoothed_y <- smoothed_y %>% 
                dplyr::filter(original) %>% 
                dplyr::select(-original)
        }

        # Replace the original dataset with the smoothed ones
        data <- data %>% 
            dplyr::select(-x, -y) %>% 
            dplyr::full_join(smoothed_y, by = "time") %>% 
            dplyr::select(-Delta_t, -index) %>% 
            dplyr::relocate(time, x, y)
    }        
        
    return(data)
}

#' Predict step in the Kalman filter
#' 
#' Assumptions in the model:
#'  a) Mean of the process noise is equal to 0
#' 
#' @param x Vector of values of the movement equation at time t.
#' @param A Transition matrix, relating values of X at time t to those at t.
#' @param u Vector of values for the external variables at time t.
#' @param B Matrix that connects the external variables in u to the values in x.
#' @param W Cholesky decomposition of the process noise covariance matrix.
#' @param F Cholesky decomposition of the estimation covariance matrix at time t.
#' 
#' @return List of predicted values for x and F at time t + 1
#' 
#' @export 
kf_predict <- function(x, 
                       A,
                       u, 
                       B,
                       W,
                       F) {
    
    # Predict values of the mean and estimation covariance. Use the square root 
    # of the covariances to ensure that they will lead to positive definite 
    # matrices. For this, use and predict values of the Cholesky decomposition 
    # and use the R matrix of a QR decomposition to update this Cholesky.
    x <- A %*% x + B %*% u
    # F <- rbind(F %*% t(A), W) %>% 
    #     qr() %>% 
    #     qr.R()
    F <- (A %*% t(F) %*% F %*% t(A) + W) %>% 
        chol()

    return(list("x" = matrix(x, ncol = 1), 
                "F" = F))
}

#' Innovation step in the Kalman filter
#' 
#' Assumptions in the model:
#'  a) Mean of the process noise is equal to 0
#' 
#' @param y Vector of measured values at time t + 1
#' @param x Vector of predicted values of the movement equation at time t + 1.
#' @param H Matrix relating the values in x to the values in y.
#' @param V Cholesky decomposition of the measurement noise covariance matrix.
#' @param F Cholesky decomposition of the estimation covariance matrix at time t + 1.
#' 
#' @return Smoothed dataframe with a similar structure as `data`
#' 
#' @export 
kf_innovation <- function(y, 
                          x, 
                          H,
                          V,
                          F) {

    # Compute the innovation z and its covariance matrix's decomposition
    z <- y - H %*% x
    G <- rbind(F %*% t(H), V) %>% 
        qr() %>% 
        qr.R()

    # Compute the Kalman gain
    K <- t(solve(G) %*% t(solve(G)) %*% H %*% t(F) %*% F)

    return(list("z" = z,
                "K" = K))
}

#' Update step in the Kalman filter
#' 
#' @param x Vector of predicted values of the movement equation at time t + 1.
#' @param z Vector of innovations at time t + 1
#' @param H Matrix relating the values in x to the values in y.
#' @param V Cholesky decomposition of the measurement noise covariance matrix.
#' @param F Cholesky decomposition of the estimation covariance matrix at time t + 1.
#' @param K Kalman gain at this updating step
#' 
#' @return Smoothed dataframe with a similar structure as `data`
#' 
#' @export 
kf_update <- function(x, 
                      z, 
                      H, 
                      V, 
                      F, 
                      K) {
    x <- x + K %*% z
    F <- rbind(F %*% t(diag(nrow(F)) - K %*% H), V %*% t(K)) %>% 
        qr() %>% 
        qr.R()
    
    return(list("x" = x,
                "F" = F))
}

# Constant velocity model: Transform data to and create the parameters
constant_velocity <- function(data,
                              reverse = TRUE, 
                              internal = FALSE,
                              assumed_variance = 0.031^2) {
    # Measurements
    y <- data %>% 
        dplyr::select(time, x, y) %>% 
        dplyr::arrange(time) %>% 
        dplyr::mutate(index = dplyr::row_number(), 
                      Delta_t = c(0, time[2:length(time)] - time[2:length(time) - 1]), 
                      original = TRUE)

    # If you want to smooth the data forwards and backwards, add the reversed 
    # data to `y`. In these data, \Delta t should still be positive, as time 
    # cannot be negative in the constant velocity model
    if(reverse & internal) {
        reversed_y <- data %>% 
            dplyr::select(time, x, y) %>% 
            dplyr::arrange(time) %>% 
            dplyr::mutate(index = dplyr::row_number()) %>% 
            dplyr::arrange(desc(time)) %>% 
            dplyr::mutate(Delta_t = abs(c(0, time[2:length(time)] - time[2:length(time) - 1])),
                          original = FALSE)

        y <- rbind(reversed_y, dplyr::filter(y, Delta_t != 0))
        y$original <- c(rep(FALSE, nrow(data) - 1),
                        rep(TRUE, nrow(data)))
    }

    # Define the columns in `y` that should be taken into account when using the 
    # filter
    cols_of_interest <- c("x", "y")

    # Create the transition matrix A, which depends on the data
    A <- function(Delta_t) {
        c(1, 0, Delta_t, 0,
          0, 1, 0, Delta_t,
          0, 0, 1, 0,
          0, 0, 0, 1) %>% 
            matrix(nrow = 4, ncol = 4, byrow = TRUE) %>% 
            return()
    }

    # Create B, which in this case is empty
    B <- matrix(0, nrow = 4, ncol = 1)

    # Define the process noise as the Random Velocity process noise, for which 
    # we will use empirically informed values using two steps. 
    #
    # Step 1 consists of creating the empirically informed values. First, we 
    # assume a given variance based on the measurement error observed in the 
    # calibration experiments. Additionally, we compute observed velocities and 
    # accelerations from the data and compute the variance at each level. For 
    # the constant velocity model, it is sufficient to use the observed 
    # variance of the acceleration for the creation of the W matrix. Note that
    # we subtract the assumed variance at the position level from this 
    # acceleration variance to get a more accurate estimate of this variance 
    # (allowing for a more accurate Kalman filter).
    #
    # Step 2 consists of using these values to create the W matrix, consisting 
    # of the movement variances at each level at each time step. Note that in 
    # most sources, this matrix is called Q instead of W, while the allowed 
    # variation at the acceleration level itself is often denoted w. I chose to 
    # be consistent with the lowercase notation of this variation. Furthermore 
    # note that we delete the derived measurement error observed at the position 
    # level from the derived variance in the acceleration, as in our derivation 
    # the measurement error seeps through.
    observed_data <- y[y$original, ]
    velocity <- data.frame(x = diff(observed_data$x) / abs(diff(observed_data$time)), 
                           y = diff(observed_data$y) / abs(diff(observed_data$time)), 
                           Delta_t = abs(diff(observed_data$time)))
    
    var_w <- c(var(velocity$x), var(velocity$y)) - 2 * mean(velocity$Delta_t)^(-2) * assumed_variance
    var_w <- ifelse(var_w <= 1e-10, 1e-10, var_w)

    cov_w <- cov(velocity$x, velocity$y)

    # In some sources called Q, while W is reserved for just the errors 
    # themselves
    W <- function(Delta_t) {
        c(Delta_t^2 * var_w[1], 0, Delta_t * var_w[1], 0, 
          0, Delta_t^2 * var_w[2], 0, Delta_t * var_w[2], 
          Delta_t * var_w[1], 0, var_w[1], 0,
          0, Delta_t * var_w[2], 0, var_w[2]) %>% 
            matrix(nrow = 4, ncol = 4, byrow = TRUE) %>% 
            return()
    }

    # Create the measurement matrix H. Only positions x and y are measured
    H <- c(1, 0, 0, 0,
           0, 1, 0, 0) %>% 
        matrix(nrow = 2, byrow = TRUE)

    # Define the measurement error covariances
    V <- matrix(c(assumed_variance, 0, 0, assumed_variance), nrow = 2, ncol = 2) %>% 
        chol()

    # Define the initial conditions. Very vague but data-informed priors
    x0 <- c(mean(observed_data$x, na.rm = TRUE), 
            mean(observed_data$y, na.rm = TRUE),
            mean(velocity$x, na.rm = TRUE),
            mean(velocity$y, na.rm = TRUE)) %>% 
        matrix(ncol = 1)
    
    F0 <- cov(cbind(observed_data$x, 
                    observed_data$y, 
                    c(NA, velocity$x), 
                    c(NA, velocity$y)), 
              use = "pairwise.complete.obs") %>% 
        diag() %>% 
        diag() %>% 
        chol()

    # Put everything in a list and return. This list looks different for the 
    # internal functions than for the kalman_filter function of the 
    # package kalmanfilter.
    if(internal) {
        return(list("y" = y,                    # Data to smooth
                    "u" = numeric(nrow(y)),     # External variables
                    "x0" = x0,                  # Prior mean
                    "F0" = F0,                  # Prior variance
                    "A" = A,                    # Transition matrix movement equation
                    "B" = B,                    # Slope for external variables
                    "W" = W,                    # Covariance matrix movement equation
                    "H" = H,                    # Measurement matrix
                    "V" = V,                    # Covariance matrix measurement equation
                    "cols_of_interest" = cols_of_interest))
    } else {
        return(list("y" = y, 
                    "B0" = x0, 
                    "P0" = t(F0) %*% F0,
                    "Dm" = lapply(1:nrow(y), \(x) matrix(0, nrow = 4, ncol = 1)) %>% 
                        make_array(),
                    "Am" = lapply(1:nrow(y), \(x) matrix(0, nrow = 2, ncol = 1)) %>% 
                        make_array(),
                    "Fm" = lapply(1:nrow(y), \(i) A(y$Delta_t[i])) %>% 
                        make_array(),
                    "Qm" = lapply(1:nrow(y), \(i) W(y$Delta_t[i])) %>% 
                        make_array(),
                    "Rm" = lapply(1:nrow(y), \(x) t(V) %*% V) %>% 
                        make_array(),
                    "Hm" = lapply(1:nrow(y), \(x) H) %>% 
                        make_array(),
                    "cols_of_interest" = cols_of_interest))
    }
}

# Constant acceleration model: Transform data to and create the parameters
constant_acceleration <- function(data,
                                  reverse = TRUE, 
                                  internal = FALSE,
                                  assumed_variance = 0.031^2) {
    # Measurements
    y <- data %>% 
        dplyr::select(time, x, y) %>% 
        dplyr::arrange(time) %>% 
        dplyr::mutate(index = dplyr::row_number(), 
                      Delta_t = c(0, time[2:length(time)] - time[2:length(time) - 1]), 
                      original = TRUE)

    # If you want to smooth the data forwards and backwards, add the reversed 
    # data to `y`. In these data, \Delta t should still be positive, as time 
    # cannot be negative in the constant velocity model
    if(reverse) {
        reversed_y <- data %>% 
            dplyr::select(time, x, y) %>% 
            dplyr::arrange(time) %>% 
            dplyr::mutate(index = dplyr::row_number()) %>% 
            dplyr::arrange(desc(time)) %>% 
            dplyr::mutate(Delta_t = abs(c(0, time[2:length(time)] - time[2:length(time) - 1])),
                          original = FALSE)

        y <- rbind(reversed_y, dplyr::filter(y, Delta_t != 0))
        y$original <- c(rep(FALSE, nrow(data) - 1),
                        rep(TRUE, nrow(data)))
    }

    # Define the columns in `y` that should be taken into account when using the 
    # filter
    cols_of_interest <- c("x", "y")

    # Create the transition matrix A, which depends on the data
    A <- function(Delta_t) {
        c(1, 0, Delta_t, 0,
          0, 1, 0, Delta_t,
          0, 0, 1, 0,
          0, 0, 0, 1) %>% 
            matrix(nrow = 4, ncol = 4, byrow = TRUE) %>% 
            return()
    }

    # Create B, which in this case is empty
    B <- matrix(0, nrow = 4, ncol = 1)

    # Define the process noise as the Random Velocity process noise, for which 
    # we will use empirically informed values using two steps. 
    #
    # Step 1 consists of creating the empirically informed values. First, we 
    # assume a given variance based on the measurement error observed in the 
    # calibration experiments. Additionally, we compute observed velocities and 
    # accelerations from the data and compute the variance at each level. For 
    # the constant velocity model, it is sufficient to use the observed 
    # variance of the acceleration for the creation of the W matrix. Note that
    # we subtract the assumed variance at the position level from this 
    # acceleration variance to get a more accurate estimate of this variance 
    # (allowing for a more accurate Kalman filter).
    #
    # Step 2 consists of using these values to create the W matrix, consisting 
    # of the movement variances at each level at each time step. Note that in 
    # most sources, this matrix is called Q instead of W, while the allowed 
    # variation at the acceleration level itself is often denoted w. I chose to 
    # be consistent with the lowercase notation of this variation. Furthermore 
    # note that we delete the derived measurement error observed at the position 
    # level from the derived variance in the acceleration, as in our derivation 
    # the measurement error seeps through.
    observed_data <- y[y$original, ]
    velocity <- data.frame(x = diff(observed_data$x) / abs(diff(observed_data$time)), 
                           y = diff(observed_data$y) / abs(diff(observed_data$time)), 
                           Delta_t = abs(diff(observed_data$time)))
    acceleration <- data.frame(x = diff(velocity$x) / abs(diff(observed_data$time, lag = 2)), 
                               y = diff(velocity$y) / abs(diff(observed_data$time, lag = 2)), 
                               Delta_t = abs(diff(observed_data$time, lag = 2)))
    
    var_w <- c(var(acceleration$x), var(acceleration$y)) - 4 * mean(acceleration$Delta_t)^(-4) * assumed_variance
    var_w <- ifelse(var_w <= 1e-10, 1e-10, var_w)

    # In some sources called Q, while W is reserved for just the errors 
    # themselves
    W <- function(Delta_t) {
        c(Delta_t^4 / 4 * var_w[1], 0, Delta_t^3 / 2 * var_w[1], 0, 
          0, Delta_t^4 / 4 * var_w[2], 0, Delta_t^3 / 2 * var_w[2], 
          Delta_t^3 / 2 * var_w[1], 0, Delta_t^2 * var_w[1], 0,
          0, Delta_t^3 / 2 * var_w[2], 0, Delta_t^2 * var_w[2]) %>% 
            matrix(nrow = 4, ncol = 4) %>% 
            return()
    }

    # Create the measurement matrix H. Only positions x and y are measured
    H <- c(1, 0, 0, 0,
           0, 1, 0, 0) %>% 
        matrix(nrow = 2, byrow = TRUE)

    # Define the measurement error covariances
    V <- matrix(c(assumed_variance, 0, 0, assumed_variance), nrow = 2, ncol = 2) %>% 
        chol()

    # Define the initial conditions. Very vague but data-informed priors
    x0 <- c(mean(data$x, na.rm = TRUE), 
            mean(data$y, na.rm = TRUE),
            mean(velocity$x, na.rm = TRUE),
            mean(velocity$y, na.rm = TRUE)) %>% 
        matrix(ncol = 1)
    
    F0 <- cov(cbind(observed_data$x, 
                    observed_data$y, 
                    c(NA, velocity$x), 
                    c(NA, velocity$y)), 
              use = "pairwise.complete.obs") %>% 
        diag() %>% 
        diag() %>% 
        chol()

    # Put everything in a list and return
    # Put everything in a list and return. This list looks different for the 
    # internal functions than for the kalman_filter function of the 
    # package kalmanfilter.
    if(internal) {
        return(list("y" = y,                    # Data to smooth
                    "u" = numeric(nrow(y)),     # External variables
                    "x0" = x0,                  # Prior mean
                    "F0" = F0,                  # Prior variance
                    "A" = A,                    # Transition matrix movement equation
                    "B" = B,                    # Slope for external variables
                    "W" = W,                    # Covariance matrix movement equation
                    "H" = H,                    # Measurement matrix
                    "V" = V,                    # Covariance matrix measurement equation
                    "cols_of_interest" = cols_of_interest))
    } else {
        return(list("y" = y, 
                    "B0" = x0, 
                    "P0" = t(F0) %*% F0,
                    "Dm" = lapply(1:nrow(y), \(x) matrix(0, nrow = 4, ncol = 1)) %>% 
                        make_array(),
                    "Am" = lapply(1:nrow(y), \(x) matrix(0, nrow = 2, ncol = 1)) %>% 
                        make_array(),
                    "Fm" = lapply(1:nrow(y), \(i) A(y$Delta_t[i])) %>% 
                        make_array(),
                    "Qm" = lapply(1:nrow(y), \(i) W(y$Delta_t[i])) %>% 
                        make_array(),
                    "Rm" = lapply(1:nrow(y), \(x) t(V) %*% V) %>% 
                        make_array(),
                    "Hm" = lapply(1:nrow(y), \(x) H) %>% 
                        make_array(),
                    "cols_of_interest" = cols_of_interest))
    }
}

# List of all models that exist
kalman_models <- list("constant_velocity" = constant_velocity, 
                      "constant_acceleration" = constant_acceleration)

# Utility function that will bind a list of matrices together in a 3D array. 
# Allows for changes in the transition matrix etc based on the Delta t of 
# those trials.
make_array <- function(x) {
    return(do.call(abind::abind, 
                   c(x, along = 3)))
}
