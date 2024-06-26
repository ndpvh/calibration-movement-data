#' Use a dynamical model to filter data
#' 
#' This function uses VAR(1) parameters to estimate the true positions of the 
#' data. It does so through means of min-log-likelihood estimation where: 
#' 
#' \bm{\mu}_t = \bm{y}_t - B \bm{y}_{t - 1} - \bm{\epsilon}_t
#' \bm{epsilon}_t \sim N(\bm{0}, \Sigma)
#' 
#' In this model, the values of $\bm{\mu}$ are estimates using the data and the 
#' "known" values of $\bm{y}$, $B$ and $\Sigma$. The values of the parameters 
#' $B$ and $\Sigma$ are estimated on the stationary data. 
#' 
#' Another approach to this kind of filtering might be more individualistic, 
#' requiring us to do a stationary calibration before each experiment and, 
#' ultimately, allowing us to estimate parameters for each tag for each session
#' separately. Something to be discussed, depending on the efficacy of this 
#' method in the synthetic simulation study.
#' 
#' @param data Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' @param parameters List containing user-specified parameters. Defaults to 
#' the parameters estimated on the stationary calibration data.
#' @param ... Arguments passed on to DEoptim.control
#' 
#' @return Smoothed dataframe with a similar structure as `data`
#' 
#' @export 
dynamic_filter <- function(data, 
                           parameters = NULL,
                           ...) {
    # Preliminaries that are used in the objecive function, but can be 
    # initialized beforehand. Consists of:
    #   - The parameters to be used, consisting of the autoregressive effect and 
    #     the lower-triangular Cholesky decomposition of the covariance matrix,
    #     as required by ldmvnorm
    #   - The predicted values of the movement y_hat according to these 
    #     parameters
    if(is.null(parameters)) {
        B <- c(7.77e-1, -1.37e-3, 4.98e-3, 7.85e-1) %>% 
            matrix(nrow = 2, ncol = 2)
        S <- c(9.66e-5, 7.25e-7, 7.25e-7, 9.32e-5) %>% 
            matrix(nrow = 2, ncol = 2)
    } else {
        B <- parameters[["B"]]
        S <- parameters[["S"]]
    }
    C <- S %>% 
        chol() %>% 
        t()
    C <- mvtnorm::ltMatrices(C[lower.tri(C, diag = TRUE)], 
                             diag = TRUE)

    # Create a y_hat that will contain the means for the normal distribution
    y_hat <- t(y[2:nrow(y), ]) - B %*% t(y[2:nrow(y) - 1,])

    # Create a local minimization functions in which the values of $\mu$ are 
    # estimated    
    objective_function <- function(y, x) {
        # Compute the min-log-likelihood
        x %>% 
            matrix(nrow = 2) %>% 
            mvtnorm::ldmvnorm(mean = y_hat, 
                              chol = C, 
                              logLik = FALSE) %>% 
            `*` (-1) %>% 
            sum(na.rm = TRUE) %>% 
            return()
    }

    # Create another function that estimates and extracts the parameters
    estimate <- function(y) {
        # First arrange the values in y according to time and select only the 
        # values of x and y
        y <- y %>% 
            dplyr::arrange(time) %>% 
            dplyr::select(x, y)

        # Do the estimation and extract the results
        params <- DEoptim::DEoptim(\(x) objective_function(y, x), 
                                   lower = rep(-1e5, (nrow(y) - 1) * 2),
                                   upper = rep(1e5, (nrow(y) - 1) * 2),
                                   control = DEoptim::DEoptim.control(...))

        params$optim$bestmem %>% 
            matrix(nrow = 2) %>% 
            t() %>% 
            return()
    }

    # Now use these two functions to estimate the values of $\bm{\mu}$. Dispatch
    # based on id
    ids <- unique(data$id) %>% 
        sort()
    results <- list()
    for(i in seq_along(ids)) {
        cat("\rEstimation for ID ", ids[i])

        tmp <- data %>% 
            dplyr::filter(id == ids[i]) %>% 
            dplyr::arrange(time)

        results[[i]] <- tmp %>% 
            estimate() %>% 
            as.data.frame() %>% 
            setNames(c("x_filtered", "y_filtered")) %>% 
            cbind(tmp[2:nrow(tmp),]) %>% 
            dplyr::select(time:y, x_filtered, y_filtered)
    }
    cat("\n")

    results <- do.call("rbind", results)
    return(results)
}