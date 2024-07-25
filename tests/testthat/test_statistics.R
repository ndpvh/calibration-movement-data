# Moving average
testthat::test_that("Moving window: Average", {
    # Create different predictable data sets
    data <- list(cbind(rep(1, each = 5),
                       rep(1, each = 5), 
                       1:5,
                       6:10), 
                 cbind(rep(1, each = 5), 
                       rep(1, each = 5), 
                       rep(1, each = 5), 
                       rep(2, each = 5)),
                 cbind(rep(1, each = 5), 
                       rep(1, each = 5), 
                       c(1, 1, NA, 1, 1), 
                       c(1, 2, 3, 4, NA)))
                       
    # Apply the average function on each of these data
    tst <- lapply(data, 
                  \(x) x %>% 
                      as.data.frame() %>% 
                      setNames(c("time", "id", "x", "y")) %>% 
                      nameless::average())

    # Define the references
    ref <- list(data.frame(x = 3, y = 8),
                data.frame(x = 1, y = 2), 
                data.frame(x = as.numeric(NA), y = as.numeric(NA)))

    # Tests
    testthat::expect_equal(tst, ref)
})

# Moving average
testthat::test_that("Moving window: Weighted average", {
    
    # Create different predictable data sets
    data <- list(cbind(rep(1, each = 5),
                       rep(1, each = 5), 
                       1:5,
                       6:10,
                       1:5), 
                 cbind(rep(1, each = 5), 
                       rep(1, each = 5), 
                       rep(1, each = 5), 
                       rep(2, each = 5),
                       1:5),
                 cbind(rep(1, each = 5), 
                       rep(1, each = 5), 
                       c(1, 1, NA, 1, 1), 
                       c(1, 2, 3, 4, NA),
                       1:5))

    ########################
    # Test error

    testthat::expect_error(nameless::weighted_average(data, weights = rep("test", 5)))
    testthat::expect_error(nameless::weighted_average(data, weights = rep("test", 3)))
    testthat::expect_error(nameless::weighted_average(data, weights = 1))

    ########################
    # Test with numeric

    # Symmetric weights
    weights <- c(0.10, 0.20, 0.40, 0.20, 0.10)

    tst <- lapply(1:3, 
                  \(i) data[[i]] %>%  
                      as.data.frame() %>% 
                      setNames(c("time", "id", "x", "y", "index")) %>% 
                      nameless::weighted_average(weights = weights))

    ref <- list(data.frame(x = 3, y = 8),
                data.frame(x = 1, y = 2), 
                data.frame(x = as.numeric(NA), y = as.numeric(NA)))

    testthat::expect_equal(tst, ref)

    # Asymmetric weights
    weights <- c(0.10, 0.20, 0.40, 0.60, 0.70)

    tst <- lapply(1:3, 
                  \(i) data[[i]] %>%  
                      as.data.frame() %>% 
                      setNames(c("time", "id", "x", "y", "index")) %>% 
                      nameless::weighted_average(weights = weights))

    ref <- list(data.frame(x = 3.8, y = 8.8),
                data.frame(x = 1, y = 2), 
                data.frame(x = as.numeric(NA), y = as.numeric(NA)))

    testthat::expect_equal(tst, ref)

    ########################
    # Test with function

    # By index
    weights <- \(x) x

    tst <- lapply(1:3, 
                  \(i) data[[i]] %>%  
                      as.data.frame() %>% 
                      setNames(c("time", "id", "x", "y", "index")) %>% 
                      nameless::weighted_average(weights = weights, .by = "index"))

    ref <- list(data.frame(x = 3.67, y = 8.67),
                data.frame(x = 1, y = 2), 
                data.frame(x = as.numeric(NA), y = as.numeric(NA)))

    testthat::expect_equal(tst, ref, tolerance = 1e-2)            
    
    # By time 
    weights <- \(x) x

    tst <- lapply(1:3, 
                  \(i) data[[i]] %>%  
                      as.data.frame() %>% 
                      setNames(c("time", "id", "x", "y", "index")) %>% 
                      nameless::weighted_average(weights = weights, .by = "time"))

    ref <- list(data.frame(x = 3, y = 8),
                data.frame(x = 1, y = 2), 
                data.frame(x = as.numeric(NA), y = as.numeric(NA)))

    testthat::expect_equal(tst, ref)    
})

# Moving middle of bin
testthat::test_that("Moving window: Middle", {
    # Create different predictable data sets
    data <- list(cbind(1:5,
                       rep(1, each = 5), 
                       1:5,
                       6:10), 
                 cbind(1:5, 
                       rep(1, each = 5), 
                       rep(1, each = 5), 
                       rep(2, each = 5)),
                 cbind(1:5, 
                       rep(1, each = 5), 
                       c(1, 1, NA, 1, 1), 
                       c(1, 2, 3, 4, NA)))
                       
    # Apply the average function on each of these data
    tst <- lapply(data, 
                  \(x) x %>% 
                      as.data.frame() %>% 
                      setNames(c("time", "id", "x", "y")) %>% 
                      nameless::middle())

    # Define the references
    ref <- list(data.frame(x = 3, y = 8),
                data.frame(x = 1, y = 2), 
                data.frame(x = as.numeric(NA), y = 3))

    # Tests
    testthat::expect_equal(tst, ref)
})

# Moving middle of bin
testthat::test_that("Moving window: Linear regression", {
    # Create different predictable data sets
    data <- list(cbind(1:5,
                       rep(1, each = 5), 
                       1:5,
                       1:5, 
                       seq(-2, 2, 1)), 
                 cbind(1:5, 
                       rep(1, each = 5), 
                       c(1, 2, 3, 4, 5), 
                       c(1, 2, 2, 3, 4),
                       seq(-2, 2, 1)),
                 cbind(1:5, 
                       rep(1, each = 5), 
                       rep(1, each = 5), 
                       rep(2, each = 5), 
                       seq(-2, 2, 1)),
                 cbind(1:5, 
                       rep(1, each = 5), 
                       rep(1, each = 5), 
                       1:5, 
                       seq(-2, 2, 1)),
                 cbind(1:5, 
                       rep(1, each = 5), 
                       1:5, 
                       rep(1, each = 5), 
                       seq(-2, 2, 1)),
                 cbind(1:5, 
                       rep(1, each = 5), 
                       c(1, 2, 2, 2, 3), 
                       c(2, 1, 2, 3, 2), 
                       seq(-2, 2, 1)),
                 cbind(1:5, 
                       rep(1, each = 5), 
                       c(1, 1, NA, 1, 1), 
                       c(1, 2, 3, 4, NA), 
                       seq(-2, 2, 1)))
                       
    # Apply the average function on each of these data
    tst <- lapply(data, 
                  \(x) x %>% 
                      as.data.frame() %>% 
                      setNames(c("time", "id", "x", "y", "index")) %>% 
                      nameless::linear() %>% 
                      suppressWarnings())

    # Define the references
    ref <- list(data.frame(x = 3, y = 3),
                data.frame(x = 2.73, y = 2.2),   # Computed by hand
                data.frame(x = 1, y = 2), 
                data.frame(x = 1, y = 3),
                data.frame(x = 3, y = 1),
                data.frame(x = 2, y = 2),
                data.frame(x = as.numeric(NA), y = 3))

    # Tests
    testthat::expect_equal(tst, ref, tolerance = 1e-2)
})

# Moving middle of bin
testthat::test_that("Moving window: Quadratic regression", {
    # Create different predictable data sets
    data <- list(cbind(1:5,
                       rep(1, each = 5), 
                       1:5,
                       1:5, 
                       seq(-2, 2, 1)), 
                 cbind(1:5, 
                       rep(1, each = 5), 
                       c(1, 2, 3, 4, 5), 
                       c(1, 2, 2, 3, 4),
                       seq(-2, 2, 1)),
                 cbind(1:5, 
                       rep(1, each = 5), 
                       rep(1, each = 5), 
                       rep(2, each = 5), 
                       seq(-2, 2, 1)),
                 cbind(1:5, 
                       rep(1, each = 5), 
                       rep(1, each = 5), 
                       1:5, 
                       seq(-2, 2, 1)),
                 cbind(1:5, 
                       rep(1, each = 5), 
                       1:5, 
                       rep(1, each = 5), 
                       seq(-2, 2, 1)),
                 cbind(1:5, 
                       rep(1, each = 5), 
                       c(1, 2, 2, 2, 3), 
                       c(2, 1, 2, 3, 2), 
                       seq(-2, 2, 1)),
                 cbind(1:5, 
                       rep(1, each = 5), 
                       c(1, 1, NA, 1, 1), 
                       c(1, 2, 3, 4, NA), 
                       seq(-2, 2, 1)), 
                 cbind(1:5, 
                       rep(1, each = 5),
                       c(-2, -1, 0, 1, 2),
                       c(4, 1, 0, 1, 4), 
                       seq(-2, 2, 1)),
                 cbind(1:5, 
                       rep(1, each = 5),
                       c(4, 1, 0, 1, 4),
                       c(-2, -1, 0, 1, 2),                        
                       seq(-2, 2, 1)))
                       
    # Apply the average function on each of these data
    tst <- lapply(data, 
                  \(x) x %>% 
                      as.data.frame() %>% 
                      setNames(c("time", "id", "x", "y", "index")) %>% 
                      nameless::parabola() %>% 
                      suppressWarnings())

    # Define the references
    ref <- list(data.frame(x = 3, y = 3),
                data.frame(x = 2.8, y = 2.12),
                data.frame(x = 1, y = 2), 
                data.frame(x = 1, y = 3),
                data.frame(x = 3, y = 1),
                data.frame(x = 2, y = 2),
                data.frame(x = as.numeric(NA), y = 3),
                data.frame(x = 0, y = 0), 
                data.frame(x = 0, y = 0))

    # Tests
    testthat::expect_equal(tst, ref, tolerance = 1e-2)
})
