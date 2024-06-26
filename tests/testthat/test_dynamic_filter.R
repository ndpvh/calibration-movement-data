# Dynamic filtering without grouping
testthat::test_that("Dynamic filter: No grouping", {
    # Create predictable data. No grouping should only be used when you have 
    # data at multiple time points for only a single individual. Hence the 
    # mock-data will reflect this
    data <- cbind(1:10,
                  rep(1, each = 10), 
                  1:10,
                  1:10) %>% 
        as.data.frame() %>% 
        setNames(c("time", "id", "x", "y"))

    # Create the reference 
    ref <- data %>% 
        dplyr::filter(time != 1) %>%
        dplyr::select(x, y) %>% 
        as.matrix()
    rownames(ref) <- NULL
    colnames(ref) <- NULL

    # Add error to the x and y data that corresponds to the autoregressive model
    # that was being created (this way, we can do a kind of recovery study)
    B <- c(7.77e-1, -1.37e-3, 4.98e-3, 7.85e-1) %>% 
            matrix(nrow = 2, ncol = 2)
    S <- c(9.66e-5, 7.25e-7, 7.25e-7, 9.32e-5) %>% 
        matrix(nrow = 2, ncol = 2)

    cols <- c("x", "y")

    set.seed(1)
    data[1, cols] <- 
        data[1, cols] + 
        MASS::mvrnorm(1, c(0, 0), S)

    for(i in 2:nrow(data)) {
        data[i, cols] <- 
            data[i, cols] + 
            t(B %*% t(data[i - 1, cols])) + 
            MASS::mvrnorm(1, c(0, 0), S)
    }

    # Create the test, filtered data
    set.seed(1)
    tst <- nameless::dynamic_filter(data, 
                                    itermax = 5e2, 
                                    NP = 100, 
                                    trace = FALSE) %>% 
        suppressWarnings() %>% 
        dplyr::select(x_filtered, y_filtered) %>% 
        as.matrix()
    rownames(tst) <- NULL
    colnames(tst) <- NULL  

    # Tests
    testthat::expect_equal(tst, ref, tolerance = 2e-1)
})

# Dynamic filtering with grouping
testthat::test_that("Dynamic filter: Grouping", {
    # Create predictable data. No grouping should only be used when you have 
    # data at multiple time points for only a single individual. Hence the 
    # mock-data will reflect this
    data <- cbind(rep(1:10, times = 2),
                  rep(1:2, each = 10), 
                  rep(1:10, times = 2),
                  rep(1:10, times = 2)) %>% 
        as.data.frame() %>% 
        setNames(c("time", "id", "x", "y"))

    # Create the reference 
    ref <- data %>% 
        dplyr::filter(time != 1) %>%
        dplyr::select(x, y) %>% 
        as.matrix()
    rownames(ref) <- NULL
    colnames(ref) <- NULL

    # Add error to the x and y data that corresponds to the autoregressive model
    # that was being created (this way, we can do a kind of recovery study). The
    # error should be added for each of the id's separately
    B <- c(7.77e-1, -1.37e-3, 4.98e-3, 7.85e-1) %>% 
            matrix(nrow = 2, ncol = 2)
    S <- c(9.66e-5, 7.25e-7, 7.25e-7, 9.32e-5) %>% 
        matrix(nrow = 2, ncol = 2)

    cols <- c("x", "y")

    ids <- unique(data$id)

    set.seed(1)
    data_list <- list()
    for(i in ids) {
        ind_data <- data[data$id == i,]

        ind_data[1, cols] <- 
            ind_data[1, cols] + 
            MASS::mvrnorm(1, c(0, 0), S)

        for(j in 2:nrow(ind_data)) {
            ind_data[j, cols] <- 
                ind_data[j, cols] + 
                t(B %*% t(ind_data[j - 1, cols])) + 
                MASS::mvrnorm(1, c(0, 0), S)
        }
        data_list[[i]] <- ind_data
    }
    data <- do.call("rbind", data_list)

    # Create the test, filtered data
    set.seed(1)
    tst <- nameless::dynamic_filter(data, 
                                    itermax = 5e2, 
                                    NP = 100, 
                                    trace = FALSE) %>% 
        suppressWarnings() %>% 
        dplyr::select(x_filtered, y_filtered) %>% 
        as.matrix()
    rownames(tst) <- NULL
    colnames(tst) <- NULL  

    # Tests
    testthat::expect_equal(tst, ref, tolerance = 4e-1)
})
