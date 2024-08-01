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