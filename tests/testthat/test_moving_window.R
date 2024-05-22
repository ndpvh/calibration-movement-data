# Test without grouping
testthat::test_that("Moving window: No grouping", {
    # Create predictable data. No grouping should only be used when you have 
    # data at multiple time points for only a single individual. Hence the 
    # mock-data will reflect this
    data <- cbind(1:10,
                  rep(1, each = 10), 
                  1:10,
                  1:10) %>% 
        as.data.frame() %>% 
        setNames(c("time", "id", "x", "y"))

    # Create function that will take the mean of a given bin
    average <- function(x) {
        matrix(c(mean(x$x), mean(x$y)), 
               nrow = 1) %>% 
            as.data.frame() %>% 
            setNames(c("x", "y")) %>% 
            return()
    }

    # Apply the moving average for several cases
    tst_0 <- nameless::moving_window(data, span = 0, fx = average)
    tst_1 <- nameless::moving_window(data, span = 1, fx = average)
    tst_2 <- nameless::moving_window(data, span = 2, fx = average)

    # Create the references
    ref_0 <- cbind(x = 1:10, y = 1:10) %>% 
        as.data.frame() %>% 
        cbind(data[, c("time", "id")]) %>% 
        dplyr::select(time, id, x, y) %>% 
        dplyr::arrange(time)
    ref_1 <- cbind(x = c(1.5, 2:9, 9.5), y = c(1.5, 2:9, 9.5)) %>% 
        as.data.frame() %>% 
        cbind(data[, c("time", "id")]) %>% 
        dplyr::select(time, id, x, y) %>% 
        dplyr::arrange(time)
    ref_2 <- cbind(x = c(2, 2.5, 3:8, 8.5, 9), y = c(2, 2.5, 3:8, 8.5, 9)) %>% 
        as.data.frame() %>% 
        cbind(data[, c("time", "id")]) %>% 
        dplyr::select(time, id, x, y) %>% 
        dplyr::arrange(time)

    # Tests
    testthat::expect_equal(tst_0, ref_0)
    testthat::expect_equal(tst_1, ref_1)
    testthat::expect_equal(tst_2, ref_2)
})

# Test with grouping
testthat::test_that("Moving window: Grouping", {
    # Create predictable data. When there is grouping, you will have multiple 
    # points for each individual, preferably at the same time. Hence the 
    # mock-data will reflect this
    data <- cbind(rep(1:5, times = 2),
                  rep(1:2, each = 5), 
                  1:10,
                  1:10) %>% 
        as.data.frame() %>% 
        setNames(c("time", "id", "x", "y"))

    # Create function that will take the mean of a given bin
    average <- function(x) {
        matrix(c(mean(x$x), mean(x$y)), 
               nrow = 1) %>% 
            as.data.frame() %>% 
            setNames(c("x", "y")) %>% 
            return()
    }

    # Apply the moving average for several cases
    tst_0 <- nameless::moving_window(data, span = 0, fx = average, .by = "id")
    tst_1 <- nameless::moving_window(data, span = 1, fx = average, .by = "id")
    tst_2 <- nameless::moving_window(data, span = 2, fx = average, .by = "id")

    # Create the references
    ref_0 <- cbind(x = 1:10, y = 1:10) %>% 
        as.data.frame() %>% 
        cbind(data[, c("time", "id")]) %>% 
        dplyr::select(time, id, x, y) %>% 
        dplyr::arrange(time)
    ref_1 <- cbind(x = c(1.5, 2:4, 4.5, 6.5, 7:9, 9.5), 
                   y = c(1.5, 2:4, 4.5, 6.5, 7:9, 9.5)) %>% 
        as.data.frame() %>% 
        cbind(data[, c("time", "id")]) %>% 
        dplyr::select(time, id, x, y) %>% 
        dplyr::arrange(time)
    ref_2 <- cbind(x = c(2, 2.5, 3, 3.5, 4, 7, 7.5, 8, 8.5, 9), 
                   y = c(2, 2.5, 3, 3.5, 4, 7, 7.5, 8, 8.5, 9)) %>% 
        as.data.frame() %>% 
        cbind(data[, c("time", "id")]) %>% 
        dplyr::select(time, id, x, y) %>% 
        dplyr::arrange(time)

    # Tests
    testthat::expect_equal(tst_0, ref_0)
    testthat::expect_equal(tst_1, ref_1)
    testthat::expect_equal(tst_2, ref_2)
})





################################################################################
# STATISTICS TO USE WITH THE MOVING WINDOW

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
