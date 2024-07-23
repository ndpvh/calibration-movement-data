# Test without grouping
testthat::test_that("Bin: No grouping", {
    # Create predictable data. No grouping should only be used when you have 
    # data at multiple time points for only a single individual. Hence the 
    # mock-data will reflect this
    data <- cbind(1:10,
                  rep(1, each = 10), 
                  1:10,
                  1:10) %>% 
        as.data.frame() %>% 
        setNames(c("time", "id", "x", "y"))

    # Apply the moving average for several cases
    tst_2 <- nameless::bin(data, span = 2, fx = nameless::average)
    tst_5 <- nameless::bin(data, span = 5, fx = nameless::average)
    tst_10 <- nameless::bin(data, span = 10, fx = nameless::average)

    # Create the references
    ref_2 <- data.frame(time = seq(1, 9, 2), 
                        id = rep(1, each = 5), 
                        x = seq(1.5, 9.5, 2), 
                        y = seq(1.5, 9.5, 2)) %>% 
        dplyr::arrange(time)
    ref_5 <- data.frame(time = c(2.5, 7.5), 
                        id = rep(1, each = 2), 
                        x = c(3, 8), 
                        y = c(3, 8)) %>% 
        dplyr::arrange(time)
    ref_10 <- data.frame(time = 5, 
                         id = 1, 
                         x = 5.5, 
                         y = 5.5) %>% 
        dplyr::arrange(time)

    # Tests
    testthat::expect_equal(tst_2, ref_2)
    testthat::expect_equal(tst_5, ref_5)
    testthat::expect_equal(tst_10, ref_10)
})

# Test with grouping
testthat::test_that("Bin: Grouping", {
    # Create predictable data. When there is grouping, you will have multiple 
    # points for each individual, preferably at the same time. Hence the 
    # mock-data will reflect this
    data <- cbind(rep(1:5, times = 2),
                  rep(1:2, each = 5), 
                  1:10,
                  1:10) %>% 
        as.data.frame() %>% 
        setNames(c("time", "id", "x", "y"))

    # Apply the moving average for several cases
    tst_2 <- nameless::bin(data, span = 2, fx = nameless::average, .by = "id")
    tst_5 <- nameless::bin(data, span = 5, fx = nameless::average, .by = "id")

    # Create the references
    ref_2 <- data.frame(time = rep(c(1, 3, 5), times = 2), 
                        id = rep(1:2, each = 3), 
                        x = c(1.5, 3.5, 5, 6.5, 8.5, 10), 
                        y = c(1.5, 3.5, 5, 6.5, 8.5, 10)) %>% 
        dplyr::arrange(time)
    ref_5 <- data.frame(time = c(2.5, 2.5), 
                        id = 1:2, 
                        x = c(3, 8), 
                        y = c(3, 8)) %>% 
        dplyr::arrange(time)

    # Tests
    testthat::expect_equal(tst_2, ref_2)
    testthat::expect_equal(tst_5, ref_5)
})
