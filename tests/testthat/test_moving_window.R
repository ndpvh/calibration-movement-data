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
    local_average <- function(x) {
        matrix(c(mean(x$x), mean(x$y)), 
               nrow = 1) %>% 
            as.data.frame() %>% 
            setNames(c("x", "y")) %>% 
            return()
    }

    # Apply the moving average for several cases
    tst_0 <- nameless::moving_window(data, span = 0, fx = local_average)
    tst_1 <- nameless::moving_window(data, span = 1, fx = local_average)
    tst_2 <- nameless::moving_window(data, span = 2, fx = local_average)

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
    local_average <- function(x) {
        matrix(c(mean(x$x), mean(x$y)), 
               nrow = 1) %>% 
            as.data.frame() %>% 
            setNames(c("x", "y")) %>% 
            return()
    }

    # Apply the moving average for several cases
    tst_0 <- nameless::moving_window(data, span = 0, fx = local_average, .by = "id")
    tst_1 <- nameless::moving_window(data, span = 1, fx = local_average, .by = "id")
    tst_2 <- nameless::moving_window(data, span = 2, fx = local_average, .by = "id")

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