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
            dplyr::group_by(dplyr::across(tidyselect::all_of(.by))) %>% 
            tidyr::nest() %>% 
            dplyr::mutate(data = data %>% 
                as.data.frame() %>% 
                kalman_filter_individual(reverse = reverse, 
                                         model = model, 
                                        #  check = check,
                                         assumed_variance = assumed_variance#, 
                                        #  internal = internal
                                        ) %>% 
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
#' @param model String denoting which of the models to use.
#' 
#' @return Smoothed dataframe with the same structure as `data`
#' 
#' @export 
kalman_filter_individual <- function(data, 
                                     reverse = TRUE, 
                                     model = "constant_velocity", 
                                     assumed_variance = 0.031^2) {

    # Robustness against too little data. When there was only 1 row, errors arose
    if(nrow(data) <= 5) {
        return(data)
    }

    # Get the model parameters and initial conditions
    parameters <- kalman_models[[model]](data, 
                                         reverse = reverse, 
                                         internal = TRUE,
                                         assumed_variance = assumed_variance)

    # Extract some of the more useful parameters, namely data and initial 
    # conditions.
    z <- parameters[["z"]]
    x0 <- parameters[["x"]]
    P0 <- parameters[["P"]]

    cols <- parameters[["cols_of_interest"]]
    
    # Make a copy of y that will contain the smoothed data. Futhermore create 
    # lists that will hold the estimation covariance and the Kalman gain at 
    # each iteration
    smoothed <- z
    P <- list()
    K <- list()
    y <- list()
    
    # Iterate over the data to smooth it. In this smoothing, you take the three
    # steps of the Kalman filter, namely:
    #   1) Prediction: Using the movement equation to predict the next state of 
    #                  x based on the initial conditions x0 and P0
    #   2) Innovation: Use the prediction to defined how "wrong" the measurement 
    #                  may be and compute the Kalman gain
    #   3) Update: Make a guess about the latent state of x based on the 
    #              prediction and measurement as weighted by the Kalman gain
    for(i in seq_len(nrow(z))) {
        # Perform the prediction step of the Kalman filter. Note that in the 
        # first iteration, we use the initial (prior) guess as prediction.
        if(i == 1) {
            prediction <- list("x" = x0, 
                               "P" = P0)
        } else {
            # Extract and create the parameters that change each iteration
            # (i.e., those that depend on the time that has passed)
            F <- parameters[["F"]](z$Delta_t[i])
            W <- parameters[["W"]](z$Delta_t[i])

            # Do the prediction. Uses the initial conditions x0 and P0, the 
            # movement transition matrix F, and the movement error covariance 
            # matrix W.
            prediction <- kf_predict(x0, 
                                     P0, 
                                     F, 
                                     W,
                                     u = parameters[["u"]][i, , drop = FALSE], 
                                     B = parameters[["B"]])
        }

        # Do the inovation step. Uses the measured state z, the predicted state 
        # x and its covariance P, and finally the measurement matrix H and the 
        # assumed measurement covariance matrix R.
        innovation <- kf_innovation(matrix(as.numeric(z[i, cols]), ncol = 1),
                                    prediction[["x"]],
                                    prediction[["P"]],
                                    parameters[["H"]],
                                    parameters[["R"]])

        # Finally, update the value of x
        result <- kf_update(prediction[["x"]],
                            prediction[["P"]], 
                            innovation[["y"]], 
                            parameters[["H"]],
                            innovation[["K"]])

        # Save the results in the smoothed dataset.
        smoothed[i, cols] <- result[["x"]][seq_along(cols)]

        # Also save some of the intermediate results in a separate list, allowing
        # us to check them in case we want to
        # P[[i]] <- t(result[["P"]]) %*% result[["P"]]
        P[[i]] <- result[["P"]]
        K[[i]] <- innovation[["K"]]
        y[[i]] <- innovation[["y"]]

         # Overwrite the initial conditions with the newly acquired values
        x0 <- result[["x"]]
        P0 <- result[["P"]]
    }

    # Replace the original dataset with the smoothed ones
    smoothed <- dplyr::filter(smoothed, original)
    data <- data %>% 
        dplyr::select(-x, -y) %>% 
        dplyr::full_join(smoothed, by = "time") %>% 
        dplyr::select(-Delta_t, -index, -original) %>% 
        dplyr::relocate(time, x, y)
        
    return(data)
}





################################################################################
# STEPS IN KALMAN FILTER

#' Predict step in the Kalman filter
#' 
#' In this step, we use the movement equation to predict the new values of x 
#' based on the initial condition \code{x0}, namely through the equation
#' 
#' \begin{equation}
#'     x = F x0 + G u
#' \end{equation}
#' where F is the movement transition matrix, G is the error matrix, and u is 
#' the error at that time. For the constant velocity model, this error is 
#' equal to random changes in acceleration.
#' 
#' We also compute the certainty around this prediction based on the initial 
#' estimate of this covariance and updated with the error variance associated to 
#' u, so that:
#' 
#' \begin{equation}
#'     P = F P0 F^T + G U G^T
#' \end{equation}
#' where U is the covariance matrix of u. We additionally contract $G u$ and 
#' $G U G^T$ into the singular error $w$ and its covariance matrix $W$ for 
#' simplicity purposes.
#' 
#' Note that in this implementation, we always use Cholesky decomposition 
#' instead of actual covariances.
#' 
#' @param x0 Vector of values of the movement equation at time t.
#' @param P0 Cholesky decomposition of the covariance of x0 at time t.
#' @param F Transition matrix, relating values of X at time t to those at t.
#' @param W Cholesky decomposition of the process noise covariance matrix.
#' @param u Vector of values for the external variables at time t. By default an 
#' empty vector.
#' @param B Matrix that connects the external variables in u to the values in x.
#' By default an empty matrix.
#' 
#' @return List of predicted values for x and P at time t + 1.
#' 
#' @export 
kf_predict <- function(x0,
                       P0, 
                       F, 
                       W, 
                       u = matrix(0, nrow = length(x0), ncol = 1), 
                       B = matrix(0, nrow = length(u), ncol = length(u))) {

    # Predict values of the mean and estimation covariance. Use the square root 
    # of the covariances to ensure that they will lead to positive definite 
    # matrices. For this, use and predict values of the Cholesky decomposition 
    # and use the R matrix of a QR decomposition to update this Cholesky.
    x <- F %*% x0 + B %*% u
    P <- F %*% P0 %*% t(F) + W

    # Return the predicted value of x as well as the Cholesky decomposition of 
    # the covariance.
    return(list("x" = matrix(x, ncol = 1), 
                "P" = P)) 
}

#' Innovation step in the Kalman filter
#' 
#' In the innovation step, we will define the innovation (or error) of the 
#' measurement if the movement equation is correct in its prediction. This 
#' expected innovation is defined as:
#' 
#' \begin{equation}
#'     y = z - Hx
#' \end{equation}
#' where z is the measurement, H is the measurement matrix, and x is the 
#' prediction.
#' 
#' Additionally, we compute the covariance of this innovation, which is defined
#' as:
#' 
#' \begin{equation}
#'     S = H P H^T + R
#' \end{equation}
#' where R is the (assumed) measurement covariance. It is interesting to note 
#' that the value of S is always symmetric for each value of H. Whether its 
#' Cholesky decomposition exists, however, will depend on H as well as P.
#' 
#' Finally, we use the covariance of the prediction and the measurement 
#' covariance to define the Kalman gain, here computed as:
#' 
#' \begin{equation}
#'     K = P H^T S^{-1}
#' \end{equation}
#' making the Kalman gain a measure of how much we can trust the prediction 
#' (pitting prediction covariance over the total variance, i.e. prediction 
#' covariance + measurement covariance).
#' 
#' @param z Vector of measured values at time t + 1
#' @param x Vector of predicted values of the movement equation at time t + 1.
#' @param P Cholesky decomposition of the estimation covariance matrix at time 
#' t + 1.
#' @param H Matrix relating the values in x to the values in y.
#' @param R Cholesky decomposition of the measurement noise covariance matrix.#' 
#' 
#' @return List of the innovation ("y"), the Cholesky decomposition of its 
#' covariance matrix ("S"), and the Kalman gain ("K")
#' 
#' @export 
kf_innovation <- function(z, 
                          x, 
                          P,
                          H,
                          R) {

    # Transform the relevant matrices to a covariance matrix using their 
    # Cholesky decomposition.
    R <- t(R) %*% R

    # Compute the innovation y and its covariance matrix decomposition
    y <- z - H %*% x
    S <- H %*% P %*% t(H) + R

    # Compute the Kalman gain
    K <- P %*% t(H) %*% solve(S)

    return(list("y" = y,
                "S" = chol(S),
                "K" = K))
}

#' Update step in the Kalman filter
#' 
#' In the update step, we use the measurement and the predicted value of the 
#' variable x and use both to make a good guess of what the "real" value of 
#' x might be (as well as what the uncertainty around this guess is). 
#' 
#' We do this using the following equation:
#' 
#' \begin{equation}
#'     x_t = x + K y
#' \end{equation}
#' where x is the predicted value for the variable (coming from the prediction 
#' step), K is the Kalman gain, and y is the innovation (both coming from the 
#' innovation step).
#' 
#' For the covariance, we use the following equation:
#' 
#' \begin{equation}
#'     P_t = (I - K H) P
#' \end{equation}
#' where I is the identity matrix and P is the prediction covariance coming from 
#' the prediction step.
#' 
#' @param x Vector of predicted values of the movement equation at time t + 1.
#' @param P Cholesky decomposition of the estimation covariance matrix at time 
#' t + 1.
#' @param y Vector of innovations at time t + 1
#' @param H Measurement matrix connecting measurement to movement.
#' @param K Kalman gain
#' 
#' @return List containing the smoothed value of x ("x") together with the  
#' Cholesky decomposition of its covariance ("P")
#' 
#' @export 
kf_update <- function(x, 
                      P, 
                      y, 
                      H, 
                      K) {

    # Compute x and P         
    x <- x + K %*% y
    P <- (diag(nrow(P)) - K %*% H) %*% P
    
    return(list("x" = x,
                "P" = P))
}





################################################################################
# MODELS

# Constant velocity model
#
# GIVE SOME ADDITIONAL INFORMATION ABOUT THIS MODEL, WHAT IT IMPLIES ETC
constant_velocity <- function(data,
                              reverse = TRUE, 
                              internal = TRUE,
                              assumed_variance = 0.031^2) {

    # PREPARING DATA


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





    # MOVEMENT EQUATION PARAMETERS

    # F: Transition matrix of the movement equation, which depends on the values
    #    of the time interval
    F <- function(Delta_t) {
        c(1, 0, Delta_t, 0,
          0, 1, 0, Delta_t,
          0, 0, 1, 0,
          0, 0, 0, 1) %>% 
            matrix(nrow = 4, ncol = 4, byrow = TRUE) %>% 
            return()
    }

    # W: The process noise here is chosen to be come in through random 
    #    random accelerations that may occur, which are assumed to follow a
    #    zero-mean Gaussian distribution. For its actual values, we use 
    #    empirical estimations.
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
        tmp <- c(Delta_t^2 * var_w[1], 0, Delta_t * var_w[1], 0, 
          0, Delta_t^2 * var_w[2], 0, Delta_t * var_w[2], 
          Delta_t * var_w[1], 0, var_w[1], 0,
          0, Delta_t * var_w[2], 0, var_w[2]) %>% 
            matrix(nrow = 4, ncol = 4, byrow = TRUE)

        # tryCatch(chol(tmp), error = function(e) browser())

        tmp %>% 
            # chol() %>% 
            return()
    }

    # B and u: In the case of the constant velocity model, we assume that these
    #          are empty
    B <- matrix(0, nrow = 4, ncol = 1)
    u <- matrix(0, nrow = nrow(y), ncol = 1)





    # MEASUREMENT EQUATION PARAMETERS

    # H: Measurement matrix that connects variables of the movement equation to 
    #    the ones of the measurement equation. Here, only positions x and y are 
    #    measured
    H <- c(1, 0, 0, 0,
           0, 1, 0, 0) %>% 
        matrix(nrow = 2, byrow = TRUE)

    # R: The measurement error covariance matrix. Here, we use the value
    #    provided to the function as the assumed value of the measurement error.
    R <- c(assumed_variance, 0, 0, assumed_variance) %>% 
        matrix(nrow = 2, ncol = 2) %>% 
        chol()





    # INITIAL CONDITIONS

    # Define the initial conditions. Very vague but data-informed priors
    x0 <- c(mean(observed_data$x, na.rm = TRUE), 
            mean(observed_data$y, na.rm = TRUE),
            mean(velocity$x, na.rm = TRUE),
            mean(velocity$y, na.rm = TRUE)) %>% 
        matrix(ncol = 1)
    
    P0 <- cov(cbind(observed_data$x, 
                    observed_data$y, 
                    c(NA, velocity$x), 
                    c(NA, velocity$y)), 
              use = "pairwise.complete.obs") %>% 
        diag() %>% 
        diag()

    # Put everything in a list and return. This list looks different for the 
    # internal functions than for the kalman_filter function of the 
    # package kalmanfilter.
    if(internal) {
        return(list("z" = y,          # Data to smooth
                    "x" = x0,         # Current value of x (prior mean)
                    "P" = P0,         # Current covariance of x (prior covariance)
                    "F" = F,          # Movement transition matrix
                    "W" = W,          # Movement covariance matrix
                    "B" = B,          # External variable transition matrix
                    "u" = u,          # External variables themselves
                    "H" = H,          # Measurement matrix
                    "R" = R,          # Measurement covariance matrix
                    "cols_of_interest" = cols_of_interest))
                    
    } else {
        return(list("y" = y,                    # Data to smooth
                    "B0" = x0,                  # Prior mean
                    "P0" = t(F0) %*% F0,        # Prior covariance
                    "Dm" = lapply(1:nrow(y), \(x) matrix(0, nrow = 4, ncol = 1)) %>%    # Intercept movement equation
                        make_array(),
                    "Am" = lapply(1:nrow(y), \(x) matrix(0, nrow = 2, ncol = 1)) %>%    # Intercept measurement equation
                        make_array(),
                    "Fm" = lapply(1:nrow(y), \(i) A(y$Delta_t[i])) %>%                  # Transition matrix movement equation
                        make_array(),
                    "Qm" = lapply(1:nrow(y), \(i) W(y$Delta_t[i])) %>%                  # Covariance matrix movement equation
                        make_array(),
                    "Rm" = lapply(1:nrow(y), \(x) t(V) %*% V) %>%                       # Covariance matrix measurement equation
                        make_array(),
                    "Hm" = lapply(1:nrow(y), \(x) H) %>%                                # Measurement matrix
                        make_array(),
                    "cols_of_interest" = cols_of_interest))
    }
}

# List of all models that exist
kalman_models <- list("constant_velocity" = constant_velocity)

# Utility function that will bind a list of matrices together in a 3D array. 
# Allows for changes in the transition matrix etc based on the Delta t of 
# those trials.
make_array <- function(x) {
    return(do.call(abind::abind, 
                   c(x, along = 3)))
}
