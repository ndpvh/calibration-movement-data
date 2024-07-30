testthat::test_that("Kalman filter: Invididual, Original, CV", {
    # Create some test data
    data <- list(cbind(time = 1:10, 
                       x = 1:10, 
                       y = 1:10) %>% 
                     as.data.frame(), 
                 cbind(time = 1:10, 
                       x = rep(1, 10) + 5 * c(seq(-1, 1, length.out = 5), seq(1, -1, length.out = 5)), 
                       y = rep(1, 10)) %>% 
                     as.data.frame())

    # Add some noise to these data
    set.seed(1)
    S <- matrix(c(0.031^2, 0, 0, 0.031^2), nrow = 2, ncol = 2)
    noisy_data <- lapply(data, 
                         \(x) x + cbind(rep(0, 10), MASS::mvrnorm(nrow(x), mu = c(0, 0), Sigma = S)))

    # Use the Kalman filter on these data
    tst <- lapply(noisy_data, 
                  \(x) nameless::kalman_filter_individual(x, reverse = FALSE, model = "constant_velocity"))

    # Check how big the sd is compared to the original 0.031
    std <- sapply(seq_along(tst), 
                  \(i) sd(c(tst[[i]]$x - data[[i]]$x, 
                            tst[[i]]$y - data[[i]]$y)))

    testthat::expect_equal(std[1] / 0.031, 0.70, tolerance = 1e-2)
    testthat::expect_equal(std[2] / 0.031, 56.31, tolerance = 1e-2)
})

testthat::test_that("Kalman filter: Individual, Reverse, CV", {
    # Create some test data
    data <- list(cbind(time = 1:10, 
                       x = 1:10, 
                       y = 1:10) %>% 
                     as.data.frame(), 
                 cbind(time = 1:10, 
                       x = rep(1, 10) + 5 * c(seq(-1, 1, length.out = 5), seq(1, -1, length.out = 5)), 
                       y = rep(1, 10)) %>% 
                     as.data.frame())

    # Add some noise to these data
    set.seed(1)
    S <- matrix(c(0.031^2, 0, 0, 0.031^2), nrow = 2, ncol = 2)
    noisy_data <- lapply(data, 
                         \(x) x + cbind(rep(0, 10), MASS::mvrnorm(nrow(x), mu = c(0, 0), Sigma = S)))

    # Use the Kalman filter on these data
    tst <- lapply(noisy_data, 
                  \(x) nameless::kalman_filter_individual(x, reverse = TRUE, model = "constant_velocity"))

    # Check how big the sd is compared to the original 0.031
    std <- sapply(seq_along(tst), 
                  \(i) sd(c(tst[[i]]$x - data[[i]]$x, 
                            tst[[i]]$y - data[[i]]$y)))

    testthat::expect_equal(std[1] / 0.031, 0.30, tolerance = 1e-2)
    testthat::expect_equal(std[2] / 0.031, 73.13, tolerance = 1e-2)
})