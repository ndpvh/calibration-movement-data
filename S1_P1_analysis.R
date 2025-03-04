################################################################################
# Purpose: Analyze the raw data gathered on 14/10/2023, 21/10/2023,            #
#          22/12/2023, and 02/12/2024. In these experiments, we measured       #
#          the positions of tags that were fixed in space, allowing us to      #
#          take a closer look at the measurement error.                        #
#                                                                              #
#          The main goal of this analysis is to inform us on the structure of  #
#          the error, allowing us to simulate data with a similar type of      #
#          error as those of the real data.                                    #
#                                                                              # 
#          The structure of this file is as follows:                           # 
#              Lx-y: Loading in datafiles                                      # 
#              Lx-y: Examining bias                                            # 
#              Lx-y: Examining variance                                        #     
#              Lx-y: Examining sampling frequency                              #
#                                                                              #
#          Note that it is possible to skip the first part and continue from   #
#          second part onwards (Bias)                                          #
################################################################################

devtools::load_all()

################################################################################
# DATAFILES

# Load in the datafiles in which all results are saved. 
data <- data.table::fread(
    file.path("data", "raw_data", "datapoints.csv"),
    data.table = FALSE,
    fill = TRUE
)

# Load the experiment names and join together with the bigger datafile. This 
# will make it easier for us to select on specific dates on which experiments 
# were ran, as well as on different conditions for each of these experiments.
data <- data.table::fread(
    file.path("data", "raw_data", "experiments.csv"),
    data.table = FALSE
) %>% 
    dplyr::rename(experiment_id = id) %>% 
    dplyr::select(experiment_id, name) %>% 
    dplyr::full_join(data, by = "experiment_id") %>% 
    dplyr::select(-experiment_id) %>% 
    dplyr::rename(experiment = name)

# Define the datasets that contain the measurements of the fixed locations on 
# the different days. Add those all in a list with the date as its identifier 
# and save the list in a separate file. This is done for several reasons:
#   - Will make it easier to perform the analyses for the separate days
#   - Makes it possible for people who look at the code to reproduce our results
#     without them having access to the complete datafile (and without having to 
#     interpret the names that we provided to the experiments directly).
#   - For some datasets, some additional processing is necessary because of tags 
#     that were emitting a signal while not lying on the grid.
experiments <- list(
    "14-10-2023" = paste("stationary", c(2:3, 5:6)),
    "21-10-2023" = paste("STATIONARY", 1:2),
    "22-12-2023" = paste("stationarity", 1:4, "- 22-12-2023"),
    "16-11-2024" = paste("rec e", 1:9)
)
anchors <- list(
    "14-10-2023" = readRDS(file.path("data", "anchor_position_14-10-2023.Rds")),
    "21-10-2023" = readRDS(file.path("data", "anchor_position_21-10-2023.Rds")),
    "22-12-2023" = readRDS(file.path("data", "anchor_position_22-12-2023.Rds")),
    "16-11-2024" = readRDS(file.path("data", "anchor_position_16-11-2024.Rds"))
)

data_list <- lapply(
    names(experiments),
    function(name) {
        idx <- data$experiment %in% experiments[[name]]
        idy <- !is.na(data$x) & !is.na(data$y)

        data[idx & idy, ] %>% 
            dplyr::mutate(
                day = name,
                anchor_xmin = min(anchors[[name]][, 2]),
                anchor_xmax = max(anchors[[name]][, 2]),
                anchor_ymin = min(anchors[[name]][, 3]),
                anchor_ymax = max(anchors[[name]][, 3]),
            ) %>% 
            return()
    }
) %>% 
    `names<-` (names(experiments))
    
data_list[[3]] <- data_list[[3]][data_list[[3]]$x < 13.5, ]
    
saveRDS(
    data_list,
    file.path("data", "study 1", "data_list.Rds")
)





################################################################################
# SYSTEMATIC ERROR

# Adding real positions ########################################################

data_list <- readRDS(file.path("data", "study 1", "data_list.Rds"))

# Link the real positions of the tags to the data, allowing us to create a 
# metric of the bias around tags in a given position. This is done in a few 
# steps by creating a grid of the required size and finding out to which of the 
# grids intersections each tag can be assigned. A function is created for this.
add_locations <- function(data, 
                          size){
    # Create X and Y series for the given size and given space inbetween each
    # point
    XY <- data.frame(
        X = rep(
            seq(0, size[1], 1), 
            each = size[2] + 1
        ),
        Y = rep(
            seq(0, size[2], 1),
            times = size[1] + 1
        )
    ) %>% 
        dplyr::mutate(tag = dplyr::row_number())

    # Assign rows and columns to the tag_id's in the data. This is done through
    # a standardization and then a guess of where the tag might be. Importantly, 
    # we center the "real" positions on the measured positions.
    #
    # Note that the filter for NA values is imposed to account for when not all 
    # positions of the grid (XY) have been measured, as is the case in the 
    # experiment of the 21-10-2023.
    assign_row <- function(x, n_rows){
        (x - min(x)) %>% 
            `/` (max(x) - min(x)) %>%
            `*` (n_rows) %>% 
            round() %>% 
            return()
    }

    data <- data %>% 
        dplyr::mutate(
            X = assign_row(x, size[1]),
            Y = assign_row(y, size[2])
        ) %>% 
        dplyr::select(-tag_id) %>% 
        dplyr::full_join(
            XY,
            by = c("X", "Y")
        ) %>% 
        dplyr::filter(!is.na(x) & !is.na(y)) %>% 
        dplyr::mutate(
            X = X - mean(X) + mean(x),
            Y = Y - mean(Y) + mean(y)
        )

    return(data)    
}

sizes <- list(
    "14-10-2023" = c(10, 8),
    "21-10-2023" = c(10, 8),
    "22-12-2023" = c(10, 7),
    "16-11-2024" = c(7, 10)
)

for(i in names(data_list)) {
    data_list[[i]] <- add_locations(
        data_list[[i]],
        sizes[[i]]
    )
}


# Visualize error ##############################################################

# Now that this is all done, we can visualize how far off the measurements are
# of the real positions. First step: Doing this on average. 
dist <- lapply(
    seq_along(data_list),
    \(i) data_list[[i]] %>% 
        dplyr::mutate(dist = sqrt((x - X)^2 + (y - Y)^2)) %>% 
        dplyr::summarize(
            study = names(data_list)[i],
            mean = mean(dist),
            sd = sd(dist),
            q025 = quantile(dist, probs = c(0.025)),
            q975 = quantile(dist, probs = c(0.975))
        )
)
dist <- do.call("rbind", dist)
View(dist)

# Now do this based on the distance from the center. In the initial plots, one 
# could see that the measurements were somewhat pushed to the center.
dist <- lapply(
    seq_along(data_list),
    \(i) data_list[[i]] %>% 
        dplyr::mutate(
            dist = sqrt((x - X)^2 + (y - Y)^2),
            x_group = X - mean(x),
            y_group = Y - mean(y)
        ) %>% 
        dplyr::group_by(tag) %>% 
        dplyr::summarize(
            study = names(data_list)[i],
            mean = mean(dist),
            sd = sd(dist),
            q025 = quantile(dist, probs = c(0.025)),
            q975 = quantile(dist, probs = c(0.975)),
            x_group = x_group[1],
            y_group = y_group[1]            
        )
)
dist <- do.call("rbind", dist)
View(dist)

# Visualize this last step.
bias_plot <- function(data, 
                      title) {
    # Get distances from the mean inside of the data
    data <- data %>% 
        dplyr::mutate(
            x_group = X - mean(x),
            y_group = Y - mean(y)
        )

    # Make a plot for x and y, summarizing over the other dimension.
    plots <- lapply(
        list(
            c("x", "x_group", "X", title, ""), 
            c("y", "y_group", "Y", "", "Position")
        ),
        function(x) {
            plot_data <- data[, x[1:3]] %>% 
                setNames(c("X", "M", "actual"))

            plt <- ggplot2::ggplot(plot_data, 
                                   ggplot2::aes(x = X, 
                                                color = factor(M), 
                                                fill = factor(M))) +
                ggplot2::geom_density(alpha = 0.5) +
                ggplot2::geom_vline(ggplot2::aes(xintercept = actual,
                                                 color = factor(M)),
                                    linewidth = 0.5) +
                ggplot2::labs(x = x[5],
                              y = "Density",
                              title = x[4]) +
                ggplot2::theme_minimal() +
                ggplot2::theme(plot.title = ggplot2::element_text(size = 40,
                                                                  hjust = 0.5),
                               axis.title = ggplot2::element_text(size = 30),
                               axis.text = ggplot2::element_text(size = 20), 
                               panel.background = ggplot2::element_rect(fill = NA, 
                                                                        linewidth = 1.5),
                               legend.position = "none")
            
            return(plt)
        }
    )

    ggpubr::ggarrange(
        plotlist = plots, 
        ncol = 1
    ) %>% 
        return()
}

plots <- lapply(
    seq_along(data_list), 
    \(i) bias_plot(
        data_list[[i]],
        names(data_list)[i]
    )
)
plots <- append(
    list(
        ggpubr::ggarrange(
            nameless::name_plot("X", size = 17),
            nameless::name_plot("Y", size = 17),
            ncol = 1
        )
    ),
    plots
)

plt <- ggpubr::ggarrange(
    plotlist = plots,
    nrow = 1,
    widths = c(0.1, rep(0.2, 4))
)
ggplot2::ggsave(
    file.path("figures", "study 1", "systematic error, uncorrected.png"),
    plt,
    width = 1500 * 4,
    height = 750 * 4,
    units = "px"
)


# Polynomial ###################################################################

# Now as to how to handle this. Bind all datasets together, normalize all data 
# so that they fall between -1 and 1 depending on the locations of the anchors, 
# and then fit a 4th degree multilvel polynomial on the result. 
#
# The fixed effects of this polynomial will be taken as the parameters of the 
# model.
all_data <- do.call(
    "rbind",
    data_list
) %>% 
    dplyr::mutate(
        x = 2 * (x - anchor_xmin) / (anchor_xmax - anchor_xmin) - 1,
        y = 2 * (y - anchor_ymin) / (anchor_ymax - anchor_ymin) - 1,
        X = 2 * (X - anchor_xmin) / (anchor_xmax - anchor_xmin) - 1,
        Y = 2 * (Y - anchor_ymin) / (anchor_ymax - anchor_ymin) - 1,
    )

# Create a function that will take in the data and compute the parameters of 
# interest through least-squares. We allow the degree of the polynomial to be 
# specified by the user, allowing us to compare different polynomial fits to 
# each other through AIC and cross-validation.
#
# Separate functions for the creation of the primary variables in the least-
# squares and for the computation. This allows us to use the same functions 
# in the cross-validation, reducing the redundancy in the code.
create_xy <- function(x, degree) {
    Y <- cbind(x$X, x$Y)

    X <- matrix(
        1, 
        nrow = nrow(x),
        ncol = 1
    )
    for(i in 1:degree) {
        X <- cbind(X, x$x^i, x$y^i)
    }

    return(
        list(
            "Y" = Y,
            "X" = X
        )
    )
}

polynomial <- function(x, degree) {
    # Create the Y and X variables to be used in the polynomial computation
    variables <- create_xy(x, degree)
    Y <- variables$Y 
    X <- variables$X

    # Compute the parameters through the analytic solution of the least-squares
    B <- tryCatch(
        solve(t(X) %*% X) %*% (t(X) %*% Y),
        error = function(e) {
            return(NA)
        }
    )

    # Get the SSE
    if(any(is.na(B))) {
        SSE <- NA
    } else {
        SSE <- sum((Y - X %*% B)^2)
    }
    
    return(
        list(
            "Y" = Y,
            "X" = X,
            "B" = B, 
            "SSE" = SSE
        )
    )
}

# Define the AIC function. Takes in the data and degree, and will return the 
# resulting AIC based on the SSE.
aic <- function(x, degree) {
    # Compute the polynomial
    result <- polynomial(x, degree)

    # Compute and return the AIC
    k <- length(result$B)
    n <- nrow(x)

    return(n * log(result$SSE / n) + 2 * (k + 1))
}

# Define the cross-validation function. Takes in data and degree, and will return
# the MSE for the test data.
#
# Used a leave-one-out cross-validation on the tags, so that all data of a single
# tag (location) will be used as test data.
cross_validation <- function(x, degree) {
    # Define the tags to leave out of the training procedure at each iteration
    tags <- unique(x$tag)
    iter <- length(tags)

    # Loop over each of the iterations
    MSE <- numeric(iter)
    for(i in 1:iter) {
        # Get all non-NA indices to be used in this iteration
        idy <- which(x$tag == tags[i])

        # Compute the polynomial on the training set
        data_i <- x[-idy, ]
        params <- polynomial(data_i, degree)

        if(is.na(params$SSE)) {
            MSE[i] <- NA
            next
        }

        # Create the variables of interest based on the test set
        variables <- create_xy(x[idy, ], degree)
        Y <- variables$Y 
        X <- variables$X 

        # Compute the MSE for this iteration
        B <- params$B 
        MSE[i] <- sum((Y - X %*% B)^2) / length(idy)
    }

    # Return the mean MSE across all iterations as the metric for fit
    return(mean(MSE, na.rm = TRUE))
}

# With all necessary functions defined, loop over all data-files and all 
# polynomial degrees of interest and compute AICs and MSEs from the 
# cross-validation
studies <- unique(dist$day)

set.seed(10) # Retrograde - Silverstein
result <- lapply(
    studies, 
    function(x) {
        # Select the data for only that day
        data_x <- all_data[all_data$day == x, ]

        # Perform aic and cross-validation for each of the degrees of interest
        AIC <- MSE <- numeric(15)
        for(i in 1:15) {
            print(i)
            AIC[i] <- aic(data_x, i)
            MSE[i] <- cross_validation(data_x, i)
        }

        # Create a dataframe combining these results for this specific day
        result <- cbind(1:15, AIC, MSE) %>% 
            as.data.frame() %>% 
            setNames(c("degree", "aic", "mse"))
        result$day <- x 

        return(result)
    }
)

# Bind the results together and inspect them
result <- do.call("rbind", result)
result <- rbind(
    result,
    result %>% 
        dplyr::group_by(degree) %>% 
        dplyr::summarize(
            day = "average",
            degree = degree[1],
            aic = mean(aic, na.rm = TRUE),
            mse = mean(mse, na.rm = TRUE)
        )
)
data.table::fwrite(
    result,
    file.path("results", "study 1", "polynomial, degree comparison.csv")
)

plt_1 <- ggplot2::ggplot(result,
                         ggplot2::aes(x = degree, 
                                      y = aic, 
                                      color = factor(day))) +
    ggplot2::geom_line(linewidth = 2) +
    ggplot2::geom_point(size = 4) +
    ggplot2::scale_color_manual(values = c("14-10-2023" = "cornflowerblue", 
                                           "21-10-2023" = "salmon",
                                           "22-12-2023" = "goldenrod",
                                           "16-11-2024" = "darkolivegreen4",
                                           "average" = "black")) +
    ggplot2::scale_y_continuous(labels = scales::scientific) +
    ggplot2::labs(x = "Degree of the polynomial",
                  y = "AIC",
                  title = "Fit",
                  color = "Day") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 40,
                                                      hjust = 0.5),
                   axis.title = ggplot2::element_text(size = 30),
                   axis.text = ggplot2::element_text(size = 20), 
                   legend.title = ggplot2::element_text(size = 20),
                   legend.text = ggplot2::element_text(size = 17),
                   panel.background = ggplot2::element_rect(fill = NA, 
                                                            linewidth = 1.5))

plt_2 <- ggplot2::ggplot(result,
                         ggplot2::aes(x = degree, 
                                      y = mse, 
                                      color = factor(day))) +
    ggplot2::geom_line(linewidth = 2) +
    ggplot2::geom_point(size = 4) +
    ggplot2::scale_color_manual(values = c("14-10-2023" = "cornflowerblue", 
                                           "21-10-2023" = "salmon",
                                           "22-12-2023" = "goldenrod",
                                           "16-11-2024" = "darkolivegreen4",
                                           "average" = "black")) +
    ggplot2::scale_y_continuous(labels = scales::scientific) +
    ggplot2::labs(x = "Degree of the polynomial",
                  y = "MSE",
                  title = "Cross-validation",
                  color = "Day") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 40,
                                                      hjust = 0.5),
                   axis.title = ggplot2::element_text(size = 30),
                   axis.text = ggplot2::element_text(size = 20), 
                   legend.title = ggplot2::element_text(size = 20),
                   legend.text = ggplot2::element_text(size = 17),
                   panel.background = ggplot2::element_rect(fill = NA, 
                                                            linewidth = 1.5))

plt <- ggpubr::ggarrange(
    plt_1, 
    plt_2,
    nrow = 1,
    common.legend = TRUE,
    legend = "right"
)

ggplot2::ggsave(
    file.path("figures", "study 1", "systematic error, polynomial degree comparison.png"),
    plt,
    width = 1500 * 3,
    height = 700 * 3,
    unit = "px"
)

# Result:
#   - AIC does not show sufficient punishing for more complex models and is 
#     therefore discarded
#   - Cross-validation shows the best fit across days for the polynomial of 
#     degree 9. This one is selected for the next steps.



# Error reduction ##############################################################

# Let's check whether this works. 
correct <- function(x) {
    # Estimate a polynomial of the 9th degree
    params <- polynomial(x, 9)

    # Retrieve the independent variables X and the parameters B
    X <- params$X 
    B <- params$B

    # Compute the result Y and add it to the dataframe
    Y <- X %*% B
    x$x_corrected <- Y[, 1]
    x$y_corrected <- Y[, 2]

    return(x)
}

# Correct the data and transform the positions back to their original locations.
# Additionally compute the distance of the raw and the corrected measurements
# to the real positions.
corrected_data <- lapply(
    data_list, 
    \(x) x %>% 
        # Transform to -1, 1 range based on anchor positions
        dplyr::mutate(
            x = 2 * (x - anchor_xmin) / (anchor_xmax - anchor_xmin) - 1,
            y = 2 * (y - anchor_ymin) / (anchor_ymax - anchor_ymin) - 1,
            X = 2 * (X - anchor_xmin) / (anchor_xmax - anchor_xmin) - 1,
            Y = 2 * (Y - anchor_ymin) / (anchor_ymax - anchor_ymin) - 1
        ) %>% 
        # Use the polynomial to correct the distortion
        correct() %>% 
        # Transform back to original scale
        dplyr::mutate(
            x = (anchor_xmax - anchor_xmin) * (x + 1) / 2 + anchor_xmin,
            y = (anchor_ymax - anchor_ymin) * (y + 1) / 2 + anchor_ymin,
            X = (anchor_xmax - anchor_xmin) * (X + 1) / 2 + anchor_xmin,
            Y = (anchor_ymax - anchor_ymin) * (Y + 1) / 2 + anchor_ymin,
            x_corrected = (anchor_xmax - anchor_xmin) * (x_corrected + 1) / 2 + anchor_xmin,
            y_corrected = (anchor_ymax - anchor_ymin) * (y_corrected + 1) / 2 + anchor_ymin
        ) %>% 
        # Compute distance as a metric of success of the transformation
        dplyr::mutate(
            dist = sqrt((x - X)^2 + (y - Y)^2),
            dist_cor = sqrt((x_corrected - X)^2 + (y_corrected - Y)^2),
            x_group = X - mean(x),
            y_group = Y - mean(y) 
        )
)
names(corrected_data) <- names(data_list)

# Let's compare the efficacy of getting rid of the distortion with the current
# method. For this, we use a nonparametric bootstrapping procedure, bootstrapping
# the mean distance of the measured/corrected position vs the real position.
# We do this for each of the datasets separately.
studies <- unique(dist$day)
x_group <- unique(dist$x_group)
y_group <- unique(dist$y_group)
results <- matrix(
    0, 
    nrow = length(studies), 
    ncol = 9
)
N <- 1000

set.seed(11) # Hagiophobia - Trophy Scars
for(i in seq_along(studies)) {
    # Select the data of interest
    data_i <- corrected_data[[studies[i]]]
    tags <- unique(data_i$tag)

    # Instantiate two vectors that will contain the bootstrapped result. Despite
    # me wanting to, I didn't vectorize this to spare the memory of my pc, as 
    # vectorized nonparametric boostraps are quite heavy on memory with these 
    # data.
    #
    # Specifications of the bootstrap:
    #   - 1000 samples
    #   - Relationship between uncorrected and corrected remains
    #   - Sample sizes equal for bootstrap and real data within tags
    #
    # First determine which indices to take for each of the bootstraps. Then 
    # loop over each of the bootstraps and compute the statistic of interest.
    bootstrapped <- matrix(
        0, 
        nrow = nrow(data_i), 
        ncol = N
    )
    idx <- 1
    for(j in seq_along(tags)) {
        # Select tag data
        sample_idx <- which(data_i$tag == tags[j])

        # Bootstrap indices and save them in the general index list
        bootstrapped[idx:(idx + length(sample_idx) - 1), ] <- sample(
            sample_idx,
            length(sample_idx) * N,
            replace = TRUE
        )

        # Update index idx
        idx <- idx + length(sample_idx)
    }

    # Now that we have the indices, bootstrap the data itself and compute
    # the mean distance from the real positions
    uncorrected <- corrected <- difference <- numeric(N)
    for(j in seq_len(N)) {
        idx <- bootstrapped[, j]
        uncorrected[j] <- dist$dist[idx] %>% 
            mean()
        corrected[j] <- dist$dist_cor[idx] %>% 
            mean()
        difference[j] <- (dist$dist[idx] - dist$dist_cor[idx]) %>% 
            mean()
    }

    # Save the statistics of interest
    results[i, ] <- c(
        mean(uncorrected),
        quantile(uncorrected, 0.025),
        quantile(uncorrected, 0.975),
        mean(corrected), 
        quantile(corrected, 0.025),
        quantile(corrected, 0.975),
        mean(difference), 
        quantile(difference, 0.025),
        quantile(difference, 0.975)
    )

    # # Transform to a dataframe, add information on the x- and y-groups and 
    # # compute the statistics of interest
    # tmp <- data.frame(
    #     uncorrected = uncorrected,
    #     corrected = corrected,
    #     difference = difference,
    #     x_group = x_group,
    #     y_group = y_group
    # ) 
    
    # tmp_x <- tmp %>% 
    #     dplyr::group_by(x_group) %>% 
    #     dplyr::summarize(
    #         x_group = x_group[1],
    #         y_group = NA,
    #         m_uncorrected = mean(uncorrected),
    #         q025_uncorrected = quantile(uncorrected, 0.005),
    #         q975_uncorrected = quantile(uncorrected, 0.995),
    #         m_corrected = mean(corrected), 
    #         q025_corrected = quantile(corrected, 0.005),
    #         q975_corrected = quantile(corrected, 0.995),
    #         m_difference = mean(difference), 
    #         q025_difference = quantile(difference, 0.005),
    #         q975_difference = quantile(difference, 0.995)
    #     )

    # tmp_y <- tmp %>% 
    #     dplyr::group_by(y_group) %>% 
    #     dplyr::summarize(
    #         x_group = NA,
    #         y_group = y_group[1],
    #         m_uncorrected = mean(uncorrected),
    #         q025_uncorrected = quantile(uncorrected, 0.005),
    #         q975_uncorrected = quantile(uncorrected, 0.995),
    #         m_corrected = mean(corrected), 
    #         q025_corrected = quantile(corrected, 0.005),
    #         q975_corrected = quantile(corrected, 0.995),
    #         m_difference = mean(difference), 
    #         q025_difference = quantile(difference, 0.005),
    #         q975_difference = quantile(difference, 0.995)
    #     )

    # tmp_xy <- tmp %>% 
    #     dplyr::group_by(x_group, y_group) %>% 
    #     dplyr::summarize(
    #         x_group = x_group[1],
    #         y_group = y_group[1],
    #         m_uncorrected = mean(uncorrected),
    #         q025_uncorrected = quantile(uncorrected, 0.005),
    #         q975_uncorrected = quantile(uncorrected, 0.995),
    #         m_corrected = mean(corrected), 
    #         q025_corrected = quantile(corrected, 0.005),
    #         q975_corrected = quantile(corrected, 0.995),
    #         m_difference = mean(difference), 
    #         q025_difference = quantile(difference, 0.005),
    #         q975_difference = quantile(difference, 0.995)
    #     )

    # # Bind the results together and put them in the list
    # tmp <- rbind(
    #     tmp_x, 
    #     tmp_y,
    #     tmp_xy
    # )
    # tmp$significance <- (tmp$q025_difference > 0) | (tmp$q975_difference < 0)
    # tmp$full_ridance <- (tmp$q025_corrected < 0) & (tmp$q975_corrected > 0)

    # results[[studies[i]]] <- tmp
}

results <- as.data.frame(results) %>% 
    setNames(
        c(
            "m_uncorrected",
            "q025_uncorrected",
            "q975_uncorrected",
            "m_corrected",
            "q025_corrected",
            "q975_corrected",
            "m_difference",
            "q025_difference",
            "q975_difference"
        )
    )
results$significance <- (results$q025_difference > 0) | (results$q975_difference < 0)
results$full_ridance <- (results$q025_difference < 0) & (results$q975_difference > 0)

# Significant reduction in the systematic error, but no riddance yet. Instead of
# an (unsigned) bias of around 24cm, we reduce the (unsigned) bias to around 
# 11cm. Unclear, unhowever, how much bias reduction there is on the level of the 
# tags.
#
# TO DO: Check whether relationship to x_group and y_group. A bit more difficult
#        to achieve with the bootstrap, so maybe real ANOVA?



# Visualization ################################################################

# Visualize the bias when correction with the polynomial equation is done.
plots <- lapply(
    seq_along(studies), 
    function(i) {
        tmp <- corrected_data[[studies[i]]]
        tmp$x <- tmp$x_corrected
        tmp$y <- tmp$y_corrected
        
        return(
            bias_plot(
                tmp,
                studies[i]
            )
        )
    }
)
plots <- append(
    list(
        ggpubr::ggarrange(
            nameless::name_plot("X", size = 17),
            nameless::name_plot("Y", size = 17),
            ncol = 1
        )
    ),
    plots
)

plt <- ggpubr::ggarrange(
    plotlist = plots,
    nrow = 1,
    widths = c(0.1, rep(0.2, 4))
)
ggplot2::ggsave(
    file.path("figures", "study 1", "systematic error, corrected.png"),
    plt,
    width = 1500 * 4,
    height = 750 * 4,
    units = "px"
)

# Some additional visualization of the real and measured positions before and 
# after the correction.
#
# Create the function that will handle the plotting.
grid <- function(data, title = TRUE) {
    # Uncorrected
    plt_1 <- ggplot2::ggplot(data, 
                             ggplot2::aes(x = x,
                                          y = y)) +
        ggplot2::geom_point(size = 3, 
                            color = "black") +
        ggplot2::annotate("point", 
                          x = data$X, 
                          y = data$Y,
                          size = 3,
                          color = "cornflowerblue") +
        ggplot2::annotate("point", 
                          x = rep(c(data$anchor_xmin, data$anchor_xmax), each = 2),
                          y = rep(c(data$anchor_ymin, data$anchor_ymax), times = 2),
                          shape = 17,
                          size = 5, 
                          color = "salmon") +
        ggplot2::lims(x = range(c(data$anchor_xmin, data$anchor_xmax)), 
                      y = range(c(data$anchor_ymin, data$anchor_ymax))) +
        ggplot2::labs(x = "x",
                      y = "y",
                      title = ifelse(title, "Before", " ")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 40,
                                                          hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 30),
                       axis.text = ggplot2::element_text(size = 20), 
                       panel.background = ggplot2::element_rect(fill = NA, 
                                                                linewidth = 1.5),
                       legend.position = "none") +
        ggplot2::coord_equal()
            
    # Corrected
    plt_2 <- ggplot2::ggplot(data, 
                             ggplot2::aes(x = x_corrected,
                                          y = y_corrected)) +
        ggplot2::geom_point(size = 3, 
                            color = "black") +
        ggplot2::annotate("point", 
                          x = data$X, 
                          y = data$Y,
                          size = 3,
                          color = "cornflowerblue") +
        ggplot2::annotate("point", 
                          x = rep(c(data$anchor_xmin, data$anchor_xmax), each = 2),
                          y = rep(c(data$anchor_ymin, data$anchor_ymax), times = 2),
                          shape = 17,
                          size = 5, 
                          color = "salmon") +
        ggplot2::lims(x = range(c(data$anchor_xmin, data$anchor_xmax)), 
                      y = range(c(data$anchor_ymin, data$anchor_ymax))) +
        ggplot2::labs(x = "x",
                      y = "y",
                      title = ifelse(title, "After", " ")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 40,
                                                          hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 30),
                       axis.text = ggplot2::element_text(size = 20), 
                       panel.background = ggplot2::element_rect(fill = NA, 
                                                                linewidth = 1.5),
                       legend.position = "none") +
        ggplot2::coord_equal()

    return(
        ggpubr::ggarrange(
            plt_1,
            plt_2,
            nrow = 1
        )
    )
}

plots <- lapply(
    seq_along(studies), 
    function(i) {
        tmp <- corrected_data[[studies[i]]]
        tmp$x <- tmp$x_corrected
        tmp$y <- tmp$y_corrected
        
        return(grid(tmp, i == 1))
    }
)

plt <- ggpubr::ggarrange(
    # Name plots
    ggpubr::ggarrange(
        plotlist = lapply(
            studies,
            \(x) nameless::name_plot(x, size = 17)
        ),
        ncol = 1
    ),
    # Actual grids
    ggpubr::ggarrange(
        plotlist = plots,
        nrow = length(studies)
    ),
    widths = c(0.2, 0.8)
)
ggplot2::ggsave(
    file.path("figures", "study 1", "systematic error, grid.png"),
    plt,
    width = 1500 * 4,
    height = 3000 * 4,
    units = "px"
)





################################################################################
# UNSYSTEMATIC ERROR

# PURPOSE: Estimate the error covariance of the stationary data. For this, we 
#          use a nonparametric bootstrap of the standard deviation around the 
#          central measured position (mean, mode, or median).
#
#          We have two variations: One for the overall variation across tags, 
#          and one that is tag-specific. For the second one, we will check whether
#          the error relates to the position of the tag in space.
#
#          Importantly, the analysis makes the assumption of independence of 
#          error across time. This assumption is tested beforehand.





#------------------------------------------------------------------------------#
# Assumption of time-independence

# Create a function that will compute the autocorrelation for a given variable
autocorr <- function(x) {
    # If too few data points, skip
    if(nrow(x) < 3) {
        return(NA)
    }

    # Arrange data according to time
    x <- x %>% 
        dplyr::arrange(time) %>% 
        dplyr::select(x) %>% 
        unlist() %>% 
        as.numeric()

    return(cor(x[2:length(x)], x[2:length(x) - 1]))
}

# Test the time-independence of the measured positions per experiment per tag
corrs <- lapply(data_list, 
                \(x) x %>% 
                    dplyr::arrange(experiment, id, time) %>% 
                    dplyr::group_by(experiment, id) %>% 
                    tidyr::nest() %>% 
                    dplyr::mutate(auto_x = purrr::map(data, 
                                                      \(X) X %>% 
                                                          dplyr::select(time, x) %>% 
                                                          setNames(c("time", "x")) %>% 
                                                          autocorr()), 
                                  auto_y = purrr::map(data, 
                                                      \(X) X %>% 
                                                          dplyr::select(time, y) %>% 
                                                          setNames(c("time", "x")) %>% 
                                                          autocorr())) %>% 
                    dplyr::select(-data))
names(correlations) <- data_files

# Results: 
#   There is quite a high autocorrelation, meaning that errors are not 
#   independent over time. 
#
# Consequence: 
#   In our estimation of the error covariances, we need to account for this 
#   time dependence
#
# How: 
#   In our approach, we will use the analytic least-squares method to estimate 
#   the mean, autoregressive component, and the covariance matrix of a simple 
#   VAR(1) on the data, defined as:
#
#       \bm{y}_t = \bm{\alpha} + B \bm{y}_{t - 1} + \bm{\epsilon}_t
#       \bm{\epsilon}_t \sim N(\bm{0}, \Sigma).
#
#   We will then take the mean value for the transition and 
#   covariance matrix as an approximation for the actual time-dependence and 
#   random error. This info can then be used to estimate the actual positions 
#   of the tags at time t as: 
#
#       \bm{\mu}_t = \bm{y}_t - B \bm{y}_{t - 1} - \bm{\epsilon}_t
#       \bm{\epsilon}_t \sim N(\bm{0}, \Sigma),
#
#   where the values of \bm{\mu}_t are estimated using the min-log-likelihood.
#   Note that we use the data to estimate an equal amount of parameters. This 
#   may not be optimal, but that's left to practice to figure out. Furthermore 
#   note that while we use the same data sets to estimate the parameters of the 
#   VAR(1) and the values of \bm{\mu}_t, in the actual data we won't do this, 
#   but rather use the values that we get from these stationary data to estimate
#   \bm{\mu}_t

# Let's create a function to estimate the parameters of the VAR(1)
autoregression <- function(x) {
    # Arrange the variables according to time
    x <- x %>% 
        dplyr::arrange(time)

    # Prepare the variables
    Y <- x %>% 
        dplyr::mutate(x = ifelse(time == max(time), NA, x), 
                      y = ifelse(time == max(time), NA, y)) %>% 
        dplyr::filter(!is.na(x)) %>% 
        dplyr::select(x, y) %>% 
        as.matrix() %>% 
        t()

    X <- x %>% 
        dplyr::mutate(x = ifelse(time == min(time), NA, x), 
                      y = ifelse(time == min(time), NA, y),
                      intercept = 1) %>% 
        dplyr::filter(!is.na(x)) %>% 
        dplyr::select(intercept, x, y) %>% 
        as.matrix() %>% 
        t()

    # Do least-squares
    B <-  Y %*% t(X) %*% solve(X %*% t(X))

    # Compute the residuals of the model and use them to estimate the covariance
    # matrix.
    e <- Y - B %*% X
    S <- cov(t(e))

    # Extract the other parameters and return in a vector
    return(c(as.vector(B[,2:3]), as.vector(S)))
}

# Apply the function to each of the data sets
results <- lapply(data_list, 
                  \(x) x %>% 
                      dplyr::group_by(experiment, id) %>% 
                      dplyr::summarize(experiment = experiment[1], 
                                       id = id[1],
                                       data = cbind(x, y, time) %>% 
                                           as.data.frame() %>% 
                                           setNames(c("x", "y", "time")) %>% 
                                           autoregression() %>% 
                                           t() %>% 
                                           as.data.frame()) %>% 
                      tidyr::unnest(data) %>% 
                      dplyr::ungroup() %>% 
                      setNames(c("experiment", "id", 
                                 "auto_x", "cross_yx", "cross_xy", "auto_y", 
                                 "var_x", "cov_yx", "cov_xy", "var_y")))

# Get the overall quantiles 0.025, 0.50, and 0.975 for each of the parameters
results <- do.call("rbind", results) %>% 
    dplyr::select(auto_x:var_y) %>% 
    as.matrix() %>% 
    matrixStats::colQuantiles(probs = c(0.025, 0.5, 0.975))

saveRDS(results, 
        file.path("results", "stationary", "var_params.Rds"))

# As a small (quick) test, check whether the dynamic filtering approach would 
# work on these stationary data (recovery done in tests, but gives some
# weird values in the simulated data: This thus serves as a sanity check)
n_cores <- parallel::detectCores()
# results <- parallel::mclapply(data_list, 
results <- lapply(list(data_list[[1]]),
                              \(x) equilibrium_filter(x, 
                                                      maxeval = 1e3, 
                                                      print_level = 1) %>% 
                                  suppressWarnings()))

results[[1]] %>% 
    dplyr::mutate(x = x_filtered, 
                  y = y_filtered) %>% 
    plot(per_iteration = FALSE)

# Currently seems to be infeasible due to time and computational constraints, 
# unfortunately. Way to go seems to be through simple reduction of measurement
# error. Maybe you could try something else than DEoptim?
#
# -> At this moment, filtering done on the measurement level, but maybe I should 
#    only filter the residuals. Might give different results. Flagged this 
#    because values of y_hat are weird in the optimizer





#------------------------------------------------------------------------------#
# Overall measurement error

# Create a function that will create bootstrapped data in an efficient way. 
# Assumption here is that x is a dataframe that contains, among other, the 
# x and y coordinates that we want to bootstrap.
#
# Importantly, it only bootstraps rows: The relationship between x and y remains
# untouched by this function. This will allow us to not only estimate the 
# variance components for both dimensions, but also the covariance between them.
bootstrapped_covariance <- function(x, 
                                    iterations,
                                    vectorized_iterations = 100) {

    # Get the sample size of the data. Needed to ensure that each of the samples
    # has an equal size to the actual data
    N <- nrow(x)

    # Determine how many times you will have to run the `vectorized_iterations` 
    # to attain the `iterations`
    whole_number <- floor(iterations / vectorized_iterations) 
    iters <- c(rep(vectorized_iterations, each = whole_number), 
               iterations %% vectorized_iterations)

    # Remove iterations that are equal to 0 (only the case if 
    # vectorized_iterations) is a diviser of iterations
    iters <- iters[iters != 0]

    # Do a mix of vectorized and unvectorized bootstrapping to spare your system's 
    # memory.
    results <- list() ; f <- 1
    for(i in seq_along(iters)) {
        # Sample a number of indices for x that is equal to the sample size times 
        # the number of samples one wants to draw
        idx <- sample(1:N, 
                      N * iters[i], 
                      replace = TRUE)

        # Extend the dataframe to account for these values and bind them with an 
        # identity number that conveys the sample they are in
        results[[i]] <- x[idx,] %>% 
            dplyr::mutate(sample_id = rep(f:(f + iters[i] - 1), each = N)) %>% 
            # Compute the covariances based on the corrected x- and y-positions. 
            # Importantly, this is done for each separate bootstrapped sample.
            dplyr::group_by(sample_id) %>% 
            dplyr::mutate(var_x = var(x), 
                          var_y = var(y), 
                          cov_xy = cov(x, y)) %>% 
            dplyr::ungroup() %>% 
            # Delete all other information: Only keep variances and covariance
            dplyr::group_by(sample_id, var_x, var_y, cov_xy) %>% 
            tidyr::nest() %>% 
            dplyr::select(-data) %>% 
            dplyr::ungroup()

        f <- f + iters[i]
    } 
    return(do.call("rbind", results))
}

# Loop over the different dates and do all your estimation 
set.seed(39) # The Messenger - Thrice
results <- list()
for(i in seq_along(data_list)){
    print(data_files[i])

    # Load data from specific date and change it somewhat
    data <- data_list[[i]] %>% 
        # Compute several center statistics, namely mean, median, and mode per
        # tag
        dplyr::group_by(id) %>%         
        dplyr::mutate(mu_x = mean(x),
                      mu_y = mean(y),
                      median_x = median(x), 
                      median_y = median(y), 
                      mode_x = mean(modeest::mlv(x, method = "mfv")), 
                      mode_y = mean(modeest::mlv(y, method = "mfv"))) %>% 
        dplyr::ungroup() %>% 
        # Create a "corrected" version of the x- and y-positions using the mean.
        # This one will be used in the estimation of the error covariances.
        #
        # Can be replaced with median or mode too, if desired.
        dplyr::mutate(x = x - mu_x, 
                      y = y - mu_y)

    # Bootstrap the data using this function and immediately compute the necessary
    # summary statistics: 2 variances and 1 covariance.
    covariances <- bootstrapped_covariance(data, 10000, vectorized_iterations = 250) 

    # Create summary statistics for each of the covariances, and more specifically 
    # given quantiles of the bootstrapped distribution. This should give us an 
    # idea of how badly off we are.
    #
    # Here, we use a 99% CI, just to be sure
    results[[data_files[i]]] <- covariances %>% 
        # Summarize the different variables into CI and mean
        dplyr::summarize(lb_var_x = quantile(var_x, 0.005), 
                         lb_var_y = quantile(var_y, 0.005),
                         lb_cov_xy = quantile(cov_xy, 0.005), 
                         m_var_x = mean(var_x), 
                         m_var_y = mean(var_y), 
                         m_cov_xy = mean(cov_xy), 
                         ub_var_x = quantile(var_x, 0.995), 
                         ub_var_y = quantile(var_y, 0.995),
                         ub_cov_xy = quantile(cov_xy, 0.995)) %>% 
        # Restructure the dataframe to be more useful: Put each of the 
        # covariances as a row, and the lower bounds, mean, and upper bounds 
        # as the columns 
        unlist() %>% 
        matrix(nrow = 3, ncol = 3) %>% 
        as.data.frame() %>% 
        setNames(c("lb", "mean", "ub")) %>% 
        cbind(covariance = c("var_x", "var_y", "cov_xy")) 

    # Release the memory that is held up by the bootstrapped data and the data 
    # itself.
    rm(data, covariances)
}

# Save the results
saveRDS(results, 
        file.path("results", "stationary", "unsystematic_error.Rds"))

# Interpretation of the results: 
#   - Covariances between x and y are as good as 0, so not relationship in the 
#     error between both dimensions
#   - Variances are about: 
#         - 14-10-2023: x: 99%CI = [0.00020, 0.00022], mean = 0.00021
#                       y: 99%CI = [0.00018, 0.00020], mean = 0.00019
#         - 21-10-2023: x: 99%CI = [0.00046, 0.00050], mean = 0.00048
#                       y: 99%CI = [0.00031, 0.00032], mean = 0.00031
#         - 22-12-2023: x: 99%CI = [0.00734, 0.00762], mean = 0.00748
#                       y: 99%CI = [0.00804, 0.00824], mean = 0.00814
#   - Variances are lowest for the first calibration period (14-10-2023), twice 
#     as high for the second calibration period (21-10-2023), and highest in the 
#     third calibration period (22-12-2023). This has some consequences:
#         - There may be more error on the sides than in the center, as the 
#           second calibration period only measured the sides
#         - Either the precision with which you measure depends on the room, or 
#           something went wrong in the third calibration session
#   - In the worst case scenario -- which is calibration period 3, upper bound 
#     of the 99%CI -- there is about 8.73cm of standard error in the x-direction, 
#     and 9.08cm of standard error in the y-direction 





#------------------------------------------------------------------------------#
# Measurement error per tag

# Loop over the different dates and do all your estimation, but this time 
# dispatching on the id
set.seed(7244) # Falling on Deaf Ears - Hail the Sun
results <- list()
for(i in seq_along(data_list)){
    print(data_files[i])

    # Load data from specific date and change it somewhat
    data <- data_list[[i]] %>% 
        # Compute several center statistics, namely mean, median, and mode per
        # tag
        dplyr::group_by(id) %>%         
        dplyr::mutate(mu_x = mean(x),
                      mu_y = mean(y),
                      median_x = median(x), 
                      median_y = median(y), 
                      mode_x = mean(modeest::mlv(x, method = "mfv")), 
                      mode_y = mean(modeest::mlv(y, method = "mfv"))) %>% 
        dplyr::ungroup() %>% 
        # Create a "corrected" version of the x- and y-positions using the mean.
        # This one will be used in the estimation of the error covariances.
        #
        # Can be replaced with median or mode too, if desired.
        dplyr::mutate(x = x - mu_x, 
                      y = y - mu_y)

    # Bootstrap the data using this function and immediately compute the necessary
    # summary statistics: 2 variances and 1 covariance.
    covariances <- data %>% 
        dplyr::group_by(experiment, id) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(covariance = bootstrapped_covariance(as.data.frame(data), 
                                                           10000, 
                                                           vectorized_iterations = 10000) %>% 
                          list()) %>% 
        dplyr::select(-data) %>% 
        tidyr::unnest(covariance) %>% 
        dplyr::ungroup() 

    # Create summary statistics for each of the covariances, and more specifically 
    # given quantiles of the bootstrapped distribution. This should give us an 
    # idea of how badly off we are.
    #
    # Here, we use a 99% CI, just to be sure
    results[[data_files[i]]] <- covariances %>% 
        # Summarize the different variables into CI and mean
        dplyr::group_by(experiment, id) %>% 
        dplyr::summarize(lb_var_x = quantile(var_x, 0.005), 
                         lb_var_y = quantile(var_y, 0.005),
                         lb_cov_xy = quantile(cov_xy, 0.005), 
                         m_var_x = mean(var_x), 
                         m_var_y = mean(var_y), 
                         m_cov_xy = mean(cov_xy), 
                         ub_var_x = quantile(var_x, 0.995), 
                         ub_var_y = quantile(var_y, 0.995),
                         ub_cov_xy = quantile(cov_xy, 0.995)) %>% 
        dplyr::ungroup()

    # Release the memory that is held up by the bootstrapped data and the data 
    # itself.
    rm(data, covariances)
}

# Add the X and Y locations to the dataframes. Needed to be able to infer when 
# standard errors are the worst
results <- lapply(seq_along(data_list), 
                  \(x) data_list[[x]] %>% 
                      dplyr::group_by(experiment, id) %>% 
                      dplyr::summarize(X = mean(X), 
                                       Y = mean(Y), 
                                       x = mean(x), 
                                       y = mean(y)) %>% 
                      dplyr::inner_join(results[[x]]))
names(results) <- data_files

# Save the results
saveRDS(results, 
        file.path("results", "stationary", "unsystematic_error_per_tag.Rds"))

# Visualize the location of the coordinates together with their 99%CI based on 
# the mean error variance
error_plot <- function(x, 
                       title = "",
                       linewidth = 1,
                       ...) {
    # Create polygons that define a circle with center equal to the measured 
    # position and a radius equal to the distance from the mean to the bounds 
    # of a 99%CI (computed with the mean error standard deviation)
    compute_normal <- function(vx, vy, cv, mx, my) {
        co <- cbind(rep(seq(mx - 5 * sqrt(vx), mx + 5 * sqrt(vx), length.out = 250), 
                        each = 250),
                    rep(seq(my - 5 * sqrt(vy), my + 5 * sqrt(vy), length.out = 250), 
                        times = 250))
        
        S <- matrix(c(vx, cv, cv, vy), 
                    nrow = 2, 
                    ncol = 2)

        return(cbind(co, 
                     mvtnorm::dmvnorm(co,
                                      mean = c(mx, my),
                                      sigma = S)))
    }

    plot_data <- x %>% 
        dplyr::group_by(experiment, id) %>% 
        dplyr::mutate(data = compute_normal(m_var_x, 
                                            m_var_y, 
                                            m_cov_xy, 
                                            x, 
                                            y) %>% 
                          as.data.frame() %>% 
                          setNames(c("x", "y", "z")) %>% 
                          dplyr::mutate(z = z / sum(z)) %>% 
                          list()) %>% 
        dplyr::select(experiment, id, data) %>% 
        tidyr::unnest(data) %>% 
        dplyr::mutate(experiment_id = paste0(experiment, "_", id)) %>% 
        dplyr::ungroup()

    # Create the plot itself. For some reason, it doesn't want to make the plot
    # by using geom_contour (most probably because I don't have continuous 
    # values for z across the whole range of x and y). Therefore all contours 
    # added in an annotate
    expid <- unique(plot_data$experiment_id)

    plt <- ggplot2::ggplot()
    for(i in expid) {
        tmp <- dplyr::filter(plot_data, experiment_id == i)
        plt <- plt + 
            ggplot2::stat_contour(data = tmp, 
                                  ggplot2::aes(x = x, y = y, z = z),
                                  color = "black",
                                  linewidth = linewidth,
                                  alpha = 1)
    }
    plt <- plt +
        ggplot2::coord_equal() +
        ggplot2::labs(x = "x", 
                      y = "y", 
                      title = title) +
        ggplot2::theme_minimal() +
        ggplot2::theme(...)

    return(plt)
}

plots <- lapply(data_files, 
                \(x) error_plot(results[[x]], 
                                title = x, 
                                linewidth = 0.01,
                                axis.title = ggplot2::element_text(size = 20),
                                plot.title = ggplot2::element_text(size = 30,
                                                                   hjust = 0.5),
                                panel.border = ggplot2::element_rect(color = "black", 
                                                                     fill = NA,
                                                                     linewidth = 2)))

ggplot2::ggsave(file.path("figures", "stationary", "unsystematic_error_per_tag.jpg"), 
                ggpubr::ggarrange(plotlist = plots,
                                  ncol = 1),
                width = 5000, 
                height = 20000, 
                unit = "px",
                limitsize = FALSE)

# Let's do an additionaly analysis: Examine the relationship between the distance
# of a tag to the center and the error that we observe
results <- lapply(results, 
                  \(x) x %>% 
                      dplyr::mutate(distance_idealized = sqrt((X - mean(X))^2 + (Y - mean(Y))^2), 
                                    distance_idealized_x = abs(X - mean(X)),
                                    distance_idealized_y = abs(Y - mean(Y)),
                                    distance_measured = sqrt((x - mean(x))^2 + (y - mean(y))^2), 
                                    distance_measured_x = abs(x - mean(x)),
                                    distance_measured_y = abs(y - mean(y))) %>% 
                      dplyr::rowwise() %>% 
                      dplyr::mutate(lb_var = mean(c(lb_var_x, lb_var_y)), 
                                    m_var = mean(c(m_var_x, m_var_y)),
                                    ub_var = mean(c(ub_var_x, ub_var_y))) %>% 
                      dplyr::ungroup())

# Save the results
saveRDS(results, 
        file.path("results", "stationary", "unsystematic_error_per_tag.Rds"))

# Create plots of the relationship between each
dists <- c("distance_idealized", "distance_idealized_x", "distance_idealized_y", 
           "distance_measured", "distance_measured_x", "distance_measured_y")
spread <- c("m_var", "m_var_x", "m_var_y", 
            "m_var", "m_var_x", "m_var_y")

plots <- list() ; f <- 1
for(i in data_files) {
    for(j in seq_along(dists)) {
        plot_data <- results[[i]] %>% 
            dplyr::select(dists[j], spread[j]) %>% 
            setNames(c("x", "y"))
        plots[[f]] <- ggplot2::ggplot(data = plot_data, 
                                      ggplot2::aes(x = x, 
                                                   y = y)) +
            ggplot2::geom_point(size = 4, 
                                shape = 21, 
                                color = "cornflowerblue") +
            ggplot2::labs(title = i,
                          x = dists[j], 
                          y = spread[j])
        f <- f + 1
    }
}

ggplot2::ggsave(file.path("figures", "stationary", "unsystematic_error_distance.jpg"), 
                ggpubr::ggarrange(plotlist = plots, 
                                  ncol = length(dists), 
                                  nrow = length(data_files)), 
                width = 6000, 
                height = 5500, 
                unit = "px")

# Interpretation of the results: 
#   - There is some variation in the unsystematic error per tag
#   - There is evidence that this variation is related to the position of the 
#     tag in the grid: The closer to the edges, the greater the measurement 
#     error (in all datasets), and in the middle there seems to be a problem for
#     the 22-12-2023 data
#       - This only seems to be the case slightly, and primarily in the 
#         y-direction





################################################################################
# OTHER STUFF OF INTEREST

# PURPOSE: Check other things that may be of interest in the stationary 
#          calibration data. This includes:
#              - Sampling rate: Will give us an idea of how many datapoints will
#                               be averaged over when binning the data





#------------------------------------------------------------------------------#
# Sampling rate

# Create a function that will take in the duration, order them according to 
# size, and then compute the mean difference between each of the durations
sampling_rate <- function(x){
    x %>% 
        as.numeric() %>% 
        sort() %>% 
        diff() %>% 
        mean() %>% 
        return()
}

# Again loop over each of the stationary datasets
for(i in seq_along(stationary)){
    # Load the stationary data for a given date and convert the timestamps to 
    # milliseconds
    data <- load_stationary(stationary[i]) %>% 
        mutate(duration = convert_to_millisecond(timestamp))

    # Compute the mean sampling rate for each tag separately and make the 
    # dataframe somewhat easier to interpret
    result <- data %>% 
        # Get the sampling rate in msec
        group_by(tag) %>% 
        mutate(msec = sampling_rate(duration)) %>% 
        ungroup() %>% 
        # Convert the sampling rate to Hz: 1 / s -> 1000 / msec
        mutate(Hz = 1000 / msec) %>% 
        # Only retain those variables that you want to interpret
        group_by(tag, msec, Hz) %>% 
        tidyr::nest() %>% 
        select(-data)

    # Save this result
    save_result(result, 
                "sampling_rate", 
                stationary[i])

    # Make another histogram to visualize the result across tags 
    plt <- ggplot(result, 
                  aes(x = Hz)) +
        geom_histogram(color = "black", 
                       fill = "cornflowerblue") +
        labs(title = "Sampling rate per tag", 
             x = "Hz") +
        geom_vline(xintercept = 5, 
                   color = "red")

    ggsave(file.path("figures", 
                     "calibration", 
                     "stationary",
                     paste0("sampling_rate_", stationary[i], ".png")),
           plot = plt, 
           units = "px", 
           width = 1000, 
           height = 1100)
}

# Interpretation of the result:
#   - On the 14th of October, sampling rates varied substantially between 2 and 
#     4Hz instead of the expected 5Hz. Sampling rates should thus be increased 
#     for our future experiments.
#   - On the 21st of October, we found a higher than expected sampling rate. This 
#     might suggest that our initial attempts of changing the sampling rate from 
#     5Hz to 8Hz on that day might have been successful, leading to an observed 
#     sampling rate of 6-7Hz 
#   - Some tags seem to perform worse and only send out responses every once in 
#     a while.
#
# Additional comments after calibration on 22-12-2023
#   - On this day, sampling with the 6 anchors remained relatively similar to 
#     the sampling rate on 21-10-2023 (although with a slight loss of frequency).
#   - Sampling with 4 anchors had an effect on the sampling frequency. Need to 
#     find out how to increase it again.

