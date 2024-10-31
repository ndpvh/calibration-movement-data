# Create some test data
data <- list(# Linear and circular movement at same velocity
             data.frame(time = 1:10, 
                        x = 1:10, 
                        y = 1:10,
                        id = rep("test 1", 10)), 
             data.frame(time = 1:10,
                        x = cos(seq(0, 2 * pi, length.out = 10)), 
                        y = sin(seq(0, 2 * pi, length.out = 10)),
                        id = rep("test 2", 10)),
             # Linear and circular movement at increasing velocity
             data.frame(time = 1:10, 
                        x = seq(0, 2, length.out = 10)^2, 
                        y = seq(0, 2, length.out = 10)^2,
                        id = rep("test 3", 10)),
             data.frame(time = 1:10,
                        x = cos(seq(0, 1, length.out = 10)^2 * 2 * pi), 
                        y = sin(seq(0, 1, length.out = 10)^2 * 2 * pi),
                        id = rep("test 4", 10)),
             # Linear and circular movement at increasing and decreasing velocity
             data.frame(time = 1:10, 
                        x = rep(c(-1, 1), each = 5) * seq(-1, 1, length.out = 10)^2, 
                        y = rep(c(-1, 1), each = 5) * seq(-1, 1, length.out = 10)^2,
                        id = rep("test 5", 10)),
             data.frame(time = 1:10, 
                        x = cos(c(seq(0, 1, length.out = 5)^2, seq(1, 0, length.out = 5)^2) * 2 * pi), 
                        y = sin(c(seq(0, 1, length.out = 5)^2, seq(1, 0, length.out = 5)^2) * 2 * pi),
                        id = rep("test 6", 10)),
             # Linear movement in one direction
             data.frame(time = 1:10, 
                        x = rep(1, each = 10), 
                        y = 1:10,
                        id = rep("test 7", 10)), 
             data.frame(time = 1:10, 
                        x = 1:10, 
                        y = rep(1, each = 10),
                        id = rep("test 8", 10)), 
             data.frame(time = 1:10, 
                        x = rep(1, each = 10), 
                        y = seq(0, 1, length.out = 10)^2,
                        id = rep("test 9", 10)), 
             data.frame(time = 1:10, 
                        x = seq(0, 1, length.out = 10)^2, 
                        y = rep(1, each = 10),
                        id = rep("test 10", 10)), 
             # No movement 
             data.frame(time = 1:10, 
                        x = rep(1, each = 10), 
                        y = rep(2, each = 10),
                        id = rep("test 11", 10)),
             # A lot of data
             data.frame(time = 1:100, 
                        x = rep(1, 100) + 5 * c(seq(-1, 1, length.out = 50), seq(1, -1, length.out = 50)), 
                        y = rep(1, 100),
                        id = rep("test 12", 100)))

# Add some noise to these data
set.seed(1)
S <- matrix(c(0.031^2, 0, 0, 0.031^2), nrow = 2, ncol = 2)

noisy_data <- lapply(data, 
                     function(x) {
                        x[,c("x", "y")] <- x[,c("x", "y")] + MASS::mvrnorm(nrow(x), mu = c(0, 0), Sigma = S)
                        return(x)
                     })

testthat::test_that("Kalman filter: Invididual, Original, CV", {
    # Use the Kalman filter on these data
    tst <- lapply(noisy_data, 
                  \(x) nameless::kalman_filter_individual(x, reverse = FALSE, model = "constant_velocity"))

    # Check how big the sd is compared to the original 0.031
    std <- sapply(seq_along(tst), 
                  \(i) sd(c(tst[[i]]$x - data[[i]]$x, 
                            tst[[i]]$y - data[[i]]$y))) / 0.031

    testthat::expect_true(all(std < 1.15 & std > 0.59))
})

testthat::test_that("Kalman filter: Individual, Reverse, CV", {
    # Use the Kalman filter on these data
    tst <- lapply(noisy_data, 
                  \(x) nameless::kalman_filter_individual(x, reverse = TRUE, model = "constant_velocity"))

    # Check how big the sd is compared to the original 0.031
    std <- sapply(seq_along(tst), 
                  \(i) sd(c(tst[[i]]$x - data[[i]]$x, 
                            tst[[i]]$y - data[[i]]$y))) / 0.031

    testthat::expect_true(all(std < 1.19 & std > 0.48))
})

testthat::test_that("Kalman filter: Grouped, Original, CV", {
    original <- do.call("rbind", data)
    bound_data <- do.call("rbind", noisy_data)

    # Use the Kalman filter on these data
    tst <- kalman_filter(bound_data, reverse = FALSE, model = "constant_velocity", .by = "id")

    # Check how big the sd is compared to the original 0.031
    tst$x <- tst$x - original$x 
    tst$y <- tst$y - original$y
    std <- tst %>% 
        dplyr::group_by(id) %>% 
        dplyr::summarize(std = sd(c(x, y)) / 0.031) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(std) %>% 
        unlist() %>% 
        as.numeric() 

    testthat::expect_true(all(std < 1.15 & std > 0.59))
})

testthat::test_that("Kalman filter: Grouped, Reverse, CV", {
    original <- do.call("rbind", data)
    bound_data <- do.call("rbind", noisy_data)

    # Use the Kalman filter on these data
    tst <- kalman_filter(bound_data, reverse = TRUE, model = "constant_velocity", .by = "id")

    # Check how big the sd is compared to the original 0.031
    tst$x <- tst$x - original$x 
    tst$y <- tst$y - original$y
    std <- tst %>% 
        dplyr::group_by(id) %>% 
        dplyr::summarize(std = sd(c(x, y)) / 0.031) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(std) %>% 
        unlist() %>% 
        as.numeric() 

    testthat::expect_true(all(std < 1.19 & std > 0.48))
})

testthat::test_that("Kalman filter: Preserves other variables", {
    bound_data <- do.call("rbind", noisy_data) %>% 
        dplyr::mutate(another_variable = id)

    # Use the Kalman filter on these data
    tst <- kalman_filter(bound_data, reverse = TRUE, model = "constant_velocity", .by = "id")

    testthat::expect_true("another_variable" %in% colnames(tst))
})

testthat::test_that("Kalman filter: Prediction step", {
    # Create a function that will compute the manual x value. Only works for 
    # 2D arrays, assuming that if Kalman prediction works for 2D it will work 
    # for higher dimensions as well.
    manual_x <- function(x, F) {
        return(matrix(c(F[1, 1] * x[1] + F[1, 2] * x[2], 
                        F[2, 1] * x[1] + F[2, 2] * x[2]), 
                      ncol = 1))
    }

    # Create a function that will compute the manual P value. Only works for 
    # 2D arrays, assuming that if Kalman prediction works for 2D it will work 
    # for higher dimensions as well.
    manual_P <- function(P, F, W) {
        c(F[1, 1]^2 * P[1, 1] + 2 * F[1, 2] * F[1, 1] * P[1, 2] + F[1, 2] * P[2, 2], 
          F[1, 1] * (F[2, 1] * P[1, 1] + F[2, 2] * P[1, 2]) + F[1, 2] * (F[2, 1] * P[2, 1] + F[2, 2] * P[2, 2]), 
          F[1, 1] * (F[2, 1] * P[1, 1] + F[2, 2] * P[1, 2]) + F[1, 2] * (F[2, 1] * P[2, 1] + F[2, 2] * P[2, 2]), 
          F[2, 1]^2 * P[1, 1] + 2 * F[2, 1] * F[2, 2] * P[1, 2] + F[2, 2]^2 * P[2, 2]) %>% 
            matrix(nrow = 2, ncol = 2) %>% 
            `+` (W) %>% 
            return()
    }

    # Create lists of potential possibilities
    x0 <- list(c(1, 1), c(1, 0), c(0, 1), c(0, 0))
    P0 <- list(matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
               matrix(c(2, 0, 0, 2), nrow = 2, ncol = 2), 
               matrix(c(2, 0.5, 0.5, 2), nrow = 2, ncol = 2))
    F <- list(matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
              matrix(c(1, 0.5, 0, 1), nrow = 2, ncol = 2), 
              matrix(c(1, 0, 0, 0.5), nrow = 2, ncol = 2))
    W <- P0

    # Do the manual and nonmanual translations and put the results in gigantic 
    # lists
    ref <- list()
    tst <- list()
    f <- 1
    for(i in seq_along(x0)) {
        for(j in seq_along(P0)) {
            for(k in seq_along(F)) {
                for(l in seq_along(W)) {
                    # Manual translations
                    ref[[f]] <- list("x" = manual_x(x0[[i]], F[[k]]), 
                                     "P" = manual_P(P0[[j]], F[[k]], W[[l]]) %>% 
                                        chol())

                    # Through kf_predict
                    tst[[f]] <- nameless::kf_predict(x0[[i]],
                                                     P0[[j]] %>% chol(), 
                                                     F[[k]], 
                                                     W[[l]] %>% chol())

                    # Update the index f
                    f <- f + 1
                }
            }
        }
    }

    # Do the test
    testthat::expect_equal(tst, ref)
})
