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
                          model = "constant_velocity", 
                          .by = NULL) {
    
    # Dispatch on whether to group the data by a given variable or not
    if(is.null(.by)) {
        return(kalman_filter_individual(data, reverse = reverse, model = model))
    } else {
        data %>% 
            dplyr::group_by(.dots = .by) %>% 
            tidyr::nest() %>% 
            dplyr::mutate(data = data %>% 
                as.data.frame() %>% 
                kalman_filter_individual(reverse = reverse, model = model) %>% 
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
                                     model = "constant_velocity") {

    # Get the model parameters and initial conditions
    parameters <- kalman_models[[model]](data, reverse = reverse)

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

        # Overwrite the initial conditions with the newly acquired values
        x0 <- result[["x"]]
        F0 <- result[["F"]]
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
                              reverse = TRUE) {
    # Measurements
    y <- data %>% 
        dplyr::select(time, x, y) %>% 
        dplyr::arrange(time) %>% 
        dplyr::mutate(index = dplyr::row_number(), 
                      Delta_t = c(0, time[2:length(time)] - time[2:length(time) - 1]), 
                      original = TRUE)

    # If you want to smooth the data forwards and backwards, add the reversed 
    # data to `y`
    if(reverse) {
        reversed_y <- data %>% 
            dplyr::select(time, x, y) %>% 
            dplyr::arrange(time) %>% 
            dplyr::mutate(index = dplyr::row_number()) %>% 
            dplyr::arrange(desc(time)) %>% 
            dplyr::mutate(Delta_t = c(0, time[2:length(time)] - time[2:length(time) - 1]),
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
        c(1, 0, Delta_t, 0, Delta_t^2 / 2, 0,
          0, 1, 0, Delta_t, 0, Delta_t^2 / 2,
          0, 0, 1, 0, Delta_t, 0,
          0, 0, 0, 1, 0, Delta_t, 
          0, 0, 0, 0, 1, 0, 
          0, 0, 0, 0, 0, 1) %>% 
            matrix(nrow = 6, ncol = 6, byrow = TRUE) %>% 
            return()
    }

    # Create B, which in this case is empty
    B <- matrix(0, nrow = 6, ncol = 1)

    # Define the process noise as the Random Acceleration process noise
    # (see Saho (2018)). Also depends on the data and on an arbitrary value 
    # to be given for the variance that is expected. Here, we use the empirical
    # value of the acceleration variation and correct for the measurement error
    # that we have on the positions.
    var_q <- c(y$x[y$original] %>% diff() %>% diff() %>% var(),
               y$y[y$original] %>% diff() %>% diff() %>% var()) - 3 * 0.031^2
    var_q <- ifelse(var_q <= 1e-6, 1e-6, var_q)

    W <- function(Delta_t) {
        c(Delta_t^4 * var_q[1] / 4, 0, Delta_t^3 * var_q[1] / 2, 0, 0, 0,
          0, Delta_t^4 * var_q[2] / 4, 0, Delta_t^3 * var_q[2] / 2, 0, 0,
          Delta_t^3 * var_q[1] / 2, 0, Delta_t^2 * var_q[1], 0, 0, 0,
          0, Delta_t^3 * var_q[2] / 2, 0, Delta_t^2 * var_q[2], 0, 0,
          0, 0, 0, 0, var_q[1], 0,
          0, 0, 0, 0, 0, var_q[2]) %>% 
            matrix(nrow = 6, ncol = 6) %>% 
            return()
    }

    # Create the measurement matrix H. Only positions x and y are measured
    H <- c(1, 0, 0, 0, 0, 0,
           0, 1, 0, 0, 0, 0) %>% 
        matrix(nrow = 2, byrow = TRUE)

    # Define the measurement error covariances
    V <- matrix(c(0.031^2, 0, 0, 0.031^2), nrow = 2, ncol = 2) %>% 
        chol()

    # Define the initial conditions. Very vague but data-informed priors
    x0 <- c(mean(data$x, na.rm = TRUE), 
            mean(data$y, na.rm = TRUE),
            mean(diff(data$x), na.rm = TRUE),
            mean(diff(data$y), na.rm = TRUE),
            mean(diff(diff(data$x)), na.rm = TRUE), 
            mean(diff(diff(data$y)), na.rm = TRUE)) %>% 
        matrix(ncol = 1)
    F0 <- cov(cbind(data$x, data$y, 
                    c(0, diff(data$x)), c(0, diff(data$y)),
                    c(0, 0, diff(diff(data$x))), c(0, 0, diff(diff(data$y)))), 
              use = "pairwise.complete.obs") %>% 
        chol()

    # Put everything in a list and return
    return(list("y" = y, 
                "u" = numeric(nrow(y)),
                "x0" = x0, 
                "F0" = F0,
                "A" = A, 
                "B" = B,
                "W" = W, 
                "H" = H, 
                "V" = V,
                "cols_of_interest" = cols_of_interest))
}

# List of all models that exist
kalman_models <- list("constant_velocity" = \(x, reverse) constant_velocity(x, reverse = reverse))
