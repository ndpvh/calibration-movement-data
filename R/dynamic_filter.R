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

    # Create a local minimization functions in which the values of $\mu$ are 
    # estimated    
    objective_function <- function(y, x) {
        # Get the residuals from subtracting x (the estimated means) from y (the
        # data). It is on the residuals that the carry-over effect is defined, 
        # not on the raw data itself.
        x <- matrix(x, ncol = 2)
        residuals <- (y - x) %>% 
            as.matrix() %>% 
            t()

        # Compute the min-log-likelihood
        mvtnorm::ldmvnorm(residuals[, 2:nrow(residuals)], 
                          mean = B %*% residuals[, 2:nrow(residuals) - 1], 
                          chol = C, 
                          logLik = TRUE) %>% 
            `*` (-1) %>%
            return()
    }

    # Create another function that estimates and extracts the parameters
    estimate <- function(y) {
        # First arrange the values in y according to time and select only the 
        # values of x and y
        y <- y %>% 
            dplyr::arrange(time) %>% 
            dplyr::select(x, y) %>% 
            as.matrix()

        # Do the estimation and extract the results
        # params <- DEoptim::DEoptim(\(x) objective_function(y, y_hat, x), 
        #                            lower = rep(-1e2, (nrow(y) - 1) * 2),
        #                            upper = rep(1e2, (nrow(y) - 1) * 2),
        #                            control = DEoptim::DEoptim.control(...))
        params <- nloptr::nloptr(rnorm(length(y), mean = as.numeric(y), sd = 0.05), 
                                 \(x) objective_function(y, x), 
                                 lb = rep(-1e2, nrow(y) * 2),
                                 ub = rep(1e2, nrow(y) * 2),
                                 opts = list("algorithm" = "NLOPT_LN_BOBYQA", 
                                             ...))

        # params$optim$bestmem %>% 
        params$solution %>% 
            matrix(ncol = 2) %>% 
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
            cbind(tmp)
    }
    cat("\n")

    results <- do.call("rbind", results)
    return(results)
}

#' Use a equilibrium model to filter data
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
equilibrium_filter <- function(data, 
                               parameters = NULL,
                               ...) {
    # Preliminaries that are used in the objecive function, but can be 
    # initialized beforehand. Consists of:
    #   - The parameters to be used, consisting of the autoregressive effect and 
    #     the lower-triangular Cholesky decomposition of the covariance matrix,
    #     as required by ldmvnorm
    if(is.null(parameters)) {
        S <- readRDS(file.path("results", "stationary", "unsystematic_error.Rds"))
        S <- lapply(S, 
                    \(x) matrix(x$mean[c(1, 3, 3, 2)], nrow = 2, ncol = 2))
        S <- Reduce("+", S) / length(S)
    } else {
        S <- parameters[["S"]]
    }
    C <- S %>% 
        chol() %>% 
        t()
    C <- mvtnorm::ltMatrices(C[lower.tri(C, diag = TRUE)], 
                             diag = TRUE)

    # Create a local minimization functions in which the values of $\mu$ are 
    # estimated    
    objective_function <- function(y, x) {
        # Compute the min-log-likelihood
        mvtnorm::ldmvnorm(t(y), 
                          mean = t(matrix(x, ncol = 2)), 
                          chol = C, 
                          logLik = TRUE) %>% 
            `*` (-1) %>%
            return()
    }

    # Create another function that estimates and extracts the parameters
    estimate <- function(y) {
        # First arrange the values in y according to time and select only the 
        # values of x and y
        y <- y %>% 
            dplyr::arrange(time) %>% 
            dplyr::select(x, y) %>% 
            as.matrix()

        # Do the estimation and extract the results
        # params <- DEoptim::DEoptim(\(x) objective_function(y, y_hat, x), 
        #                            lower = rep(-1e2, (nrow(y) - 1) * 2),
        #                            upper = rep(1e2, (nrow(y) - 1) * 2),
        #                            control = DEoptim::DEoptim.control(...))
        params <- nloptr::nloptr(rnorm(length(y), mean = as.numeric(y), sd = 0.05), 
                                 \(x) objective_function(y, x), 
                                 lb = rep(-1e2, nrow(y) * 2),
                                 ub = rep(1e2, nrow(y) * 2),
                                 opts = list("algorithm" = "NLOPT_LN_BOBYQA", 
                                             ...))

        # params$optim$bestmem %>% 
        params$solution %>% 
            matrix(ncol = 2) %>% 
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
            cbind(tmp)
    }
    cat("\n")

    results <- do.call("rbind", results)
    return(results)
}