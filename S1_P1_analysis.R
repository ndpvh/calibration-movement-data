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
#              L23-171:     Loading in datafiles                               # 
#              L177-1020:   Examining bias                                     # 
#              L1026-1641:  Examining variance                                 #     
#              L1647-1707:  Examining sampling frequency                       #
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

# Add a time-variable to the data_list. Is a transformation of the timestamp 
# provided in the actual datapoints, but now as a numeric in msec
for(i in seq_along(data_list)) {
    data_list[[i]] <- data_list[[i]] %>% 
        dplyr::group_by(experiment) %>% 
        dplyr::mutate(
            time = withr::with_options(
                list(digits = 16),
                as.numeric(timestamp)
            ), 
            time = time - min(time)
        ) %>% 
        dplyr::ungroup()
}

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
    # a standardization and then a guess of where the tag might be.
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

# Finally, save the data_list to be used in all next steps   
saveRDS(
    data_list,
    file.path("data", "study 1", "data_list.Rds")
)





################################################################################
# SYSTEMATIC ERROR

# PURPOSE: Uncover the structure of the systematic bias that exists in the 
#          measurements. The main strategy is to compare different polynomials 
#          of varying degrees to each other through model comparison tools. 
#          Once a model is selected, we test its efficacy on the data at hand.

data_list <- readRDS(file.path("data", "study 1", "data_list.Rds"))

# Visualize error ##############################################################

# Now that this is all done, we can visualize how far off the measurements are
# of the real positions. First step: Doing this on average. 
dist <- lapply(
    seq_along(data_list),
    \(i) data_list[[i]] %>% 
        dplyr::mutate(dist = sqrt((x - X)^2 + (y - Y)^2)) %>% 
        dplyr::summarize(
            study = names(data_list)[i],
            mean = mean(dist, na.rm = TRUE),
            sd = sd(dist, na.rm = TRUE),
            q025 = quantile(dist, probs = c(0.025), na.rm = TRUE),
            q975 = quantile(dist, probs = c(0.975), na.rm = TRUE)
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
        paste0("Day ", i)
    )
)
plots <- append(
    list(
        ggpubr::ggarrange(
            nameless::name_plot("x", size = 17),
            nameless::name_plot("y", size = 17),
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
create_xy <- function(x, 
                      degree, 
                      interactions) {

    # Create Y
    Y <- cbind(x$X, x$Y)

    # Create the start of the matrix X where you already append the intercepts
    X <- matrix(
        1, 
        nrow = nrow(x),
        ncol = 1
    )

    # Check whether you want interactions to happen. If so, then we need a 
    # slightly different approach from the case where we only want main effects
    # in the polynomial.
    if(interactions) {
        # Create a degree matrix that combines all main and interaction effects
        # This makes sure that we are not defining identical terms at any point 
        # in the equation.
        #
        # Add a check that the total sum of the rows does not exceed the degree:
        # in a polynomial, this would be impossible. Also delete the intercept 
        # condition where the degree for both x and y is equal to 0.
        degrees <- 0:degree
        degrees <- cbind(
            rep(degrees, each = length(degrees)),
            rep(degrees, times = length(degrees))
        )

        degrees <- degrees[degrees[, 1] + degrees[, 2] <= degree, , drop = FALSE]
        degrees <- degrees[degrees[, 1] + degrees[, 2] != 0, , drop = FALSE]

        # Loop over each row and add the effects to the X matrix
        for(i in seq_len(nrow(degrees))) {
            X <- cbind(X, x$x^degrees[i, 1] * x$y^degrees[i, 2])
        }

    # Approach when having only main effects
    } else {
        for(i in 1:degree) {
            X <- cbind(X, x$x^i, x$y^i)
        }
    }

    return(
        list(
            "Y" = Y,
            "X" = X
        )
    )
}

polynomial <- function(x, 
                       degree, 
                       interactions) {

    # Create the Y and X variables to be used in the polynomial computation
    variables <- create_xy(
        x, 
        degree, 
        interactions
    )

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
aic <- function(x, ...) {
    # Compute the polynomial
    result <- polynomial(x, ...)

    # Compute and return the AIC
    k <- length(result$B)
    n <- nrow(x)

    return(n * log(result$SSE / n) + 2 * (k + 1))
}

# Define the BIC function. Takes in the data and degree, and will return the 
# resulting BIC based on the SSE.
bic <- function(x, ...) {
    # Compute the polynomial
    result <- polynomial(x, ...)

    # Compute and return the BIC
    k <- length(result$B)
    n <- nrow(x)

    return(log(result$SSE / n) + (k / n) * log(n))
}

# Define the cross-validation function. Takes in data and degree, and will return
# the MSE for the test data.
#
# Used a leave-one-out cross-validation on the tags, so that all data of a single
# tag (location) will be used as test data.
cross_validation <- function(x, ...) {
    # Define the tags to leave out of the training procedure at each iteration
    tags <- unique(x$tag)
    iter <- length(tags)

    # Loop over each of the iterations
    d <- numeric(iter)
    for(i in 1:iter) {
        # Get all non-NA indices to be used in this iteration
        idy <- which(x$tag == tags[i])

        # Compute the polynomial on the training set
        data_i <- x[-idy, ]
        params <- polynomial(data_i, ...)

        if(is.na(params$SSE)) {
            d[i] <- NA
            next
        }

        # Create the variables of interest based on the test set
        variables <- create_xy(x[idy, ], ...)
        Y <- variables$Y 
        X <- variables$X 

        # Compute the MSE for this iteration
        B <- params$B 
        Y_hat <- X %*% B 
        d[i] <- mean(sqrt((Y[,1] - Y_hat[,1])^2 + (Y[,2] - Y_hat[,2])^2))
    }

    # Return the mean MSE across all iterations as the metric for fit
    return(mean(d, na.rm = TRUE))
}

# With all necessary functions defined, loop over all data-files and all 
# polynomial degrees of interest and compute AICs and MSEs from the 
# cross-validation
studies <- unique(dist$study)
degrees <- 1:15

set.seed(10) # Retrograde - Silverstein
result <- lapply(
    studies, 
    function(x) {
        print(x)

        # Select the data for only that day
        data_x <- all_data[all_data$day == x, ]

        # Perform aic and cross-validation for each of the degrees of interest
        AIC <- BIC <- d <- numeric(length(degrees) * 2)
        f <- 1
        for(i in degrees) {
            print(i)
            for(j in c(TRUE, FALSE)) {
                AIC[f] <- aic(data_x, i, j)
                BIC[f] <- bic(data_x, i, j)
                d[f] <- cross_validation(data_x, i, j)
                f <- f + 1
            }
        }

        # Create a dataframe combining these results for this specific day
        result <- data.frame(
            degree = rep(degrees, each = 2),
            interaction = rep(c(TRUE, FALSE), times = length(degrees)), 
            aic = AIC, 
            bic = BIC,
            d = d
        ) 
        result$day <- x 

        return(result)
    }
)

# Bind the results together and inspect them
result <- do.call("rbind", result)
result <- rbind(
    result,
    result %>% 
        dplyr::group_by(degree, interaction) %>% 
        dplyr::summarize(
            day = "average",
            degree = degree[1],
            interaction = interaction[1],
            aic = mean(aic, na.rm = TRUE),
            bic = mean(bic, na.rm = TRUE),
            d = mean(d, na.rm = TRUE)
        ) %>% 
        dplyr::ungroup()
)
data.table::fwrite(
    result,
    file.path("results", "study 1", "polynomial, degree comparison.csv")
)

# Visualize the results through a line plot for AIC and CV
combos <- data.frame(
    rep(c("aic", "bic", "d"), times = 2),
    rep(c(FALSE, TRUE), each = 3)
)

plots <- lapply(
    seq_len(nrow(combos)),
    \(i) ggplot2::ggplot(result[result$interaction == combos[i, 2], ],
                         ggplot2::aes(x = degree, 
                                      y = .data[[combos[i, 1]]], 
                                      color = factor(day))) +
    ggplot2::geom_line(linewidth = 2) +
    ggplot2::geom_point(size = 4) +
    ggplot2::scale_color_manual(values = c("14-10-2023" = "cornflowerblue", 
                                           "21-10-2023" = "salmon",
                                           "22-12-2023" = "goldenrod",
                                           "16-11-2024" = "darkolivegreen4",
                                           "average" = "black"),
                                labels = c("14-10-2023" = 1, 
                                           "21-10-2023" = 2,
                                           "22-12-2023" = 3,
                                           "16-11-2024" = 4,
                                           "average" = "average")) +
    ggplot2::scale_y_continuous(labels = \(x) scales::scientific(x, digits = 1)) +
    ggplot2::labs(x = ifelse(combos[i, 2], "Degree", " "),
                  y = "",
                  title = ifelse(
                    i %in% 4:6, 
                    "",
                    ifelse(
                      combos[i, 1] == "aic", 
                      "AIC", 
                      ifelse(
                          combos[i, 1] == "bic", 
                          "BIC", 
                          "d"
                      )
                    )
                  ),
                  color = "Day") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 55,
                                                      hjust = 0.5),
                   axis.title = ggplot2::element_text(size = 35),
                   axis.text = ggplot2::element_text(size = 25), 
                   legend.title = ggplot2::element_text(size = 35),
                   legend.text = ggplot2::element_text(size = 30),
                   panel.background = ggplot2::element_rect(fill = NA, 
                                                            linewidth = 1.5))
)

ggplot2::ggsave(
    file.path("figures", "study 1", "systematic error, polynomial degree comparison.png"),
    ggpubr::ggarrange(
        ggpubr::ggarrange(
            nameless::name_plot("No interaction", size = 20),
            nameless::name_plot("Interaction", size = 20),
            ncol = 1
        ),
        ggpubr::ggarrange(
            plotlist = plots,
            nrow = 2,
            ncol = 3,
            common.legend = TRUE,
            legend = "right",
            widths = c(0.505, 0.495, 0.495)
        ),
        nrow = 1,
        widths = c(1/6, 5/6)
    ),
    width = 3100 * 3,
    height = 1500 * 3,
    unit = "px"
)

# Result:
#   - AIC/BIC do not show sufficient punishing for more complex models and are 
#     therefore discarded
#   - Cross-validation shows the best fit across days for the polynomial of 
#     degree 9 and without interaction effects. This one is selected for the 
#     next steps.



# Error reduction ##############################################################

# Let's check whether this works. 
correct <- function(x) {
    # Estimate a polynomial of the 9th degree
    params <- polynomial(x, 9, FALSE)

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
studies <- unique(dist$study)
x_group <- unique(dist$x_group)
y_group <- unique(dist$y_group)
results <- matrix(
    0, 
    nrow = length(studies), 
    ncol = 9
)
N <- 10000

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
    #   - 10000 samples
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
        uncorrected[j] <- data_i$dist[idx] %>% 
            mean()
        corrected[j] <- data_i$dist_cor[idx] %>% 
            mean()
        difference[j] <- (data_i$dist[idx] - data_i$dist_cor[idx]) %>% 
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

results <- cbind(studies, results) %>% 
    as.data.frame() %>% 
    setNames(
        c(
            "study",
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

data.table::fwrite(
    results,
    file.path("results", "study 1", "systematic error, reduction.csv")
)

# Significant reduction in the systematic error, but no riddance yet. Instead of
# an (unsigned) bias of around 24cm, we reduce the (unsigned) bias to around 
# 11cm. Unclear, unhowever, how much bias reduction there is on the level of the 
# tags.
#
# TO DO: Check whether relationship to x_group and y_group. A bit more difficult
#        to achieve with the bootstrap, so maybe real ANOVA?



# Finally, get some descriptives for the corrected and uncorrected data with 
# regard to their distance from the actual positions.
systematic <- lapply(
    seq_along(studies), 
    \(i) corrected_data[[studies[i]]] %>% 
        dplyr::rename(study = day) %>% 
        dplyr::group_by(study) %>% 
        dplyr::summarize(
            mean_raw = mean(dist, na.rm = TRUE),
            sd_raw = sd(dist, na.rm = TRUE),
            q025_raw = quantile(dist, prob = 0.025, na.rm = TRUE),
            q975_raw = quantile(dist, prob = 0.975, na.rm = TRUE),
            mean_cor = mean(dist_cor, na.rm = TRUE),
            sd_cor = sd(dist_cor, na.rm = TRUE),
            q025_cor = quantile(dist_cor, prob = 0.025, na.rm = TRUE),
            q975_cor = quantile(dist_cor, prob = 0.975, na.rm = TRUE),
            x_group = x_group[1], 
            y_group = y_group[1]
        ) %>% 
        dplyr::ungroup() %>% 
        return()
)
systematic <- do.call("rbind", systematic) 

data.table::fwrite(
    systematic,
    file.path("results", "study 1", "systematic error, descriptives.csv")
)


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
                paste0("Day", i)
            )
        )
    }
)
plots <- append(
    list(
        ggpubr::ggarrange(
            nameless::name_plot("x", size = 17),
            nameless::name_plot("y", size = 17),
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
#          error across time. This assumption is tested at the end of this 
#          section in the code.

data_list <- readRDS(file.path("data", "study 1", "data_list.Rds"))

# Overall measurement error ####################################################

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
    iters <- c(
        rep(vectorized_iterations, each = whole_number), 
        iterations %% vectorized_iterations
    )

    # Remove iterations that are equal to 0 (only the case if 
    # vectorized_iterations) is a diviser of iterations
    iters <- iters[iters != 0]

    # Do a mix of vectorized and unvectorized bootstrapping to spare your system's 
    # memory.
    results <- list() ; f <- 1
    for(i in seq_along(iters)) {
        # Sample a number of indices for x that is equal to the sample size times 
        # the number of samples one wants to draw
        idx <- sample(
            1:N, 
            N * iters[i], 
            replace = TRUE
        )

        # Extend the dataframe to account for these values and bind them with an 
        # identity number that conveys the sample they are in
        results[[i]] <- x[idx,] %>% 
            dplyr::mutate(sample_id = rep(f:(f + iters[i] - 1), each = N)) %>% 
            # Compute the covariances based on the corrected x- and y-positions. 
            # Importantly, this is done for each separate bootstrapped sample.
            dplyr::group_by(sample_id) %>% 
            dplyr::mutate(
                var_x = var(x), 
                var_y = var(y), 
                cov_xy = cov(x, y)
            ) %>% 
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
    print(names(data_list)[i])

    # Load data from specific date and change it somewhat
    data <- data_list[[i]] %>% 
        # Compute several center statistics, namely mean, median, and mode per
        # tag
        dplyr::group_by(tag) %>%         
        dplyr::mutate(
            mu_x = mean(x),
            mu_y = mean(y),
            median_x = median(x), 
            median_y = median(y), 
            mode_x = mean(modeest::mlv(x, method = "mfv")), 
            mode_y = mean(modeest::mlv(y, method = "mfv"))
        ) %>% 
        dplyr::ungroup() %>% 
        # Create a "corrected" version of the x- and y-positions using the mean.
        # This one will be used in the estimation of the error covariances.
        #
        # Can be replaced with median or mode too, if desired.
        dplyr::mutate(
            x = x - mu_x, 
            y = y - mu_y
        )

    # Bootstrap the data using this function and immediately compute the necessary
    # summary statistics: 2 variances and 1 covariance.
    covariances <- bootstrapped_covariance(
        data, 
        10000, 
        vectorized_iterations = 250
    ) 

    # Create summary statistics for each of the covariances, and more specifically 
    # given quantiles of the bootstrapped distribution. This should give us an 
    # idea of how badly off we are.
    #
    # Here, we use a 99% CI, just to be sure
    results[[names(data_list)[i]]] <- covariances %>% 
        # Summarize the different variables into CI and mean
        dplyr::summarize(
            lb_var_x = quantile(var_x, 0.005), 
            lb_var_y = quantile(var_y, 0.005),
            lb_cov_xy = quantile(cov_xy, 0.005), 
            m_var_x = mean(var_x), 
            m_var_y = mean(var_y), 
            m_cov_xy = mean(cov_xy), 
            ub_var_x = quantile(var_x, 0.995), 
            ub_var_y = quantile(var_y, 0.995),
            ub_cov_xy = quantile(cov_xy, 0.995)
        ) %>% 
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
saveRDS(
    results, 
    file.path("results", "study 1", "unsystematic error, overall covariance.Rds")
)

# Interpretation of the results: 
#   - Covariances between x and y are as good as 0, so not relationship in the 
#     error between both dimensions
#   - Variances are about: 
#         - 14-10-2023: x: 99%CI = [0.00020, 0.00022], mean = 0.00021
#                       y: 99%CI = [0.00018, 0.00020], mean = 0.00019
#         - 21-10-2023: x: 99%CI = [0.00046, 0.00050], mean = 0.00048
#                       y: 99%CI = [0.00031, 0.00032], mean = 0.00031
#         - 22-12-2023: x: 99%CI = [0.00418, 0.00447], mean = 0.00432
#                       y: 99%CI = [0.00043, 0.00045], mean = 0.00044
#         - 16-11-2024: x: 99%CI = [0.00049, 0.00051], mean = 0.00050
#                       y: 99%CI = [0.00292, 0.00314], mean = 0.00303
#   - Variances are lowest for the first calibration period (14-10-2023) and 
#     highest in the third calibration period (22-12-2023). 
#   - In the worst case scenario -- which is calibration period 3, upper bound 
#     of the 99%CI -- there is about 6.69cm of standard error in the x-direction, 
#     and 2.1cm of standard error in the y-direction 



# Measurement error per tag ####################################################

# Loop over the different dates and do all your estimation, but this time 
# dispatching on the id
set.seed(7244) # Falling on Deaf Ears - Hail the Sun
results <- list()
for(i in seq_along(data_list)){
    print(names(data_list)[i])

    # Load data from specific date and change it somewhat
    data <- data_list[[i]] %>% 
        # Compute several center statistics, namely mean, median, and mode per
        # tag
        dplyr::group_by(tag) %>%         
        dplyr::mutate(
            mu_x = mean(x),
            mu_y = mean(y),
            median_x = median(x), 
            median_y = median(y), 
            mode_x = mean(modeest::mlv(x, method = "mfv")), 
            mode_y = mean(modeest::mlv(y, method = "mfv"))
        ) %>% 
        dplyr::ungroup() %>% 
        # Create a "corrected" version of the x- and y-positions using the mean.
        # This one will be used in the estimation of the error covariances.
        #
        # Can be replaced with median or mode too, if desired.
        dplyr::mutate(
            x = x - mu_x, 
            y = y - mu_y
        )

    # Bootstrap the data using this function and immediately compute the necessary
    # summary statistics: 2 variances and 1 covariance.
    covariances <- data %>% 
        dplyr::group_by(tag) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(
            covariance = bootstrapped_covariance(
                as.data.frame(data), 
                10000, 
                vectorized_iterations = 250
            ) %>% 
                list()
        ) %>% 
        dplyr::select(-data) %>% 
        tidyr::unnest(covariance) %>% 
        dplyr::ungroup() 

    # Create summary statistics for each of the covariances, and more specifically 
    # given quantiles of the bootstrapped distribution. This should give us an 
    # idea of how badly off we are.
    #
    # Here, we use a 99% CI, just to be sure
    results[[names(data_list)[i]]] <- covariances %>% 
        # Summarize the different variables into CI and mean
        dplyr::group_by(tag) %>% 
        dplyr::summarize(
            lb_var_x = quantile(var_x, 0.005), 
            lb_var_y = quantile(var_y, 0.005),
            lb_cov_xy = quantile(cov_xy, 0.005), 
            m_var_x = mean(var_x), 
            m_var_y = mean(var_y), 
            m_cov_xy = mean(cov_xy), 
            ub_var_x = quantile(var_x, 0.995), 
            ub_var_y = quantile(var_y, 0.995),
            ub_cov_xy = quantile(cov_xy, 0.995)
        ) %>% 
        dplyr::ungroup()

    # Release the memory that is held up by the bootstrapped data and the data 
    # itself.
    rm(data, covariances)
}

# Add information about the mean of the tags to the results. This will be used
# in the analyses to come.
for(i in names(data_list)) {
    # Means per tag
    means <- data_list[[i]] %>% 
        dplyr::group_by(tag) %>% 
        dplyr::summarize(
            mu_x = mean(x), 
            mu_y = mean(y),
            mu_X = mean(X),
            mu_Y = mean(Y)
        ) %>% 
        dplyr::ungroup()

    # Bind them to the results
    results[[i]] <- results[[i]] %>% 
        dplyr::full_join(means, by = "tag")
}

# Save the results
saveRDS(
    results, 
    file.path("results", "study 1", "unsystematic error, per tag.Rds")
)

# Investigate whether the relationship between the observed error and the 
# distance of a tag to the center of the grid. First compute the "real" and 
# the measured distance of a tag to the center of the grid.
results <- lapply(
    results, 
    \(x) x %>% 
        dplyr::mutate(
            distance_idealized = sqrt((mu_X - mean(mu_X))^2 + (mu_Y - mean(mu_Y))^2), 
            distance_idealized_x = abs(mu_X - mean(mu_X)),
            distance_idealized_y = abs(mu_Y - mean(mu_Y)),
            distance_measured = sqrt((mu_x - mean(mu_x))^2 + (mu_y - mean(mu_y))^2), 
            distance_measured_x = abs(mu_x - mean(mu_x)),
            distance_measured_y = abs(mu_y - mean(mu_y))
        ) %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(
            lb_var = mean(c(lb_var_x, lb_var_y)), 
            m_var = mean(c(m_var_x, m_var_y)),
            ub_var = mean(c(ub_var_x, ub_var_y))
        ) %>% 
        dplyr::ungroup()
)

# Perform a linear regression between the distance from a tag to the center and 
# the variance within each of the dimensions
results_lm <- lapply(
    results,
    function(x) {
        return(
            list(
                "x" = lm(sqrt(x$m_var_x) ~ x$distance_idealized_x),
                "y" = lm(sqrt(x$m_var_y) ~ x$distance_idealized_y),
                "xy" = lm(sqrt(x$m_var) ~ x$distance_idealized)
            )
        )
    }
)
names(results_lm) <- names(results)

saveRDS(
    results_lm,
    file.path("results", "study 1", "unsystematic error, regressions per tag.Rds")
)

# Interpretation of the results:
#   - Overall, there doesn't seem to be an effect of distance from center on the 
#     variance of the error. Seems like we can consider each tag as more or less
#     having the same variance.
#   - In the y-direction, 2 days do display some significance, namely 14-10-2023
#     and 22-12-2023. Let's plot to make sense of the results.

# Create plots of the relationship between each
dists <- paste0("distance_idealized", c("", "_x", "_y"))
locations <- paste0("m_var", c("", "_x", "_y")) 

plots <- lapply(
    names(results),
    function(x) {
        data <- results[[x]]

        plt <- lapply(
            seq_along(dists), 
            function(i) {
                plot_data <- data[, c(dists[i], locations[i])] %>% 
                    setNames(c("X", "Y"))

                return(
                    ggplot2::ggplot(plot_data, 
                                    ggplot2::aes(x = X, 
                                                 y = Y)) +
                        ggplot2::geom_point(size = 4, 
                                            shape = 21, 
                                            color = "cornflowerblue") +
                        ggplot2::labs(title = ifelse(i == 1, x, ""),
                                      x = dists[i], 
                                      y = locations[i])
                )
            }
        )

        return(
            ggpubr::ggarrange(
                plotlist = plt,
                ncol = 1
            )
        )
    }
)

ggplot2::ggsave(
    file.path("figures", "study 1", "unsystematic error, per tag distance.png"), 
    ggpubr::ggarrange(
        plotlist = plots, 
        nrow = 1
    ), 
    width = 6000, 
    height = 5500, 
    unit = "px"
)

# Interpretation of the results: 
#   - If there was significance, it seemed to be mostly driven by outliers.



# Assumption of time-independence ##############################################

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
corrs <- lapply(
    data_list, 
    \(x) x %>% 
        dplyr::arrange(experiment, tag, time) %>% 
        dplyr::group_by(experiment, tag) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(
            auto_x = purrr::map(
                data, 
                \(z) z %>% 
                    dplyr::select(time, x) %>% 
                    setNames(c("time", "x")) %>% 
                    autocorr()
            ) %>% 
                unlist(), 
            auto_y = purrr::map(
                data, 
                \(z) z %>% 
                    dplyr::select(time, y) %>% 
                    setNames(c("time", "x")) %>% 
                    autocorr()
            ) %>% 
                unlist()
        ) %>% 
        dplyr::select(-data)
)
names(corrs) <- names(data_list)

# Visualize this autocorrelation
plots <- lapply(
    seq_along(corrs), 
    function(i) {
        data_i <- corrs[[i]]
        plot_data <- data.frame(
            X = c(data_i$auto_x, data_i$auto_y),
            M = rep(c("x", "y"), each = nrow(data_i))
        )

        plt <- ggplot2::ggplot(plot_data,
                               ggplot2::aes(x = X, 
                                            fill = factor(M), 
                                            color = factor(M))) +
            ggplot2::geom_density(alpha = 0.5) +
            ggplot2::scale_color_manual(values = c("x" = "cornflowerblue",
                                                   "y" = "salmon"),
                                        guide = "none") +
            ggplot2::scale_fill_manual(values = c("x" = "cornflowerblue",
                                                  "y" = "salmon")) +
            ggplot2::labs(x = "Autocorrelation",
                          y = "Density",
                          title = names(corrs)[i],
                          fill = "Dimension") +
            ggplot2::lims(x = c(0, 1)) +
            ggplot2::theme_minimal() +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 40,
                                                              hjust = 0.5),
                           axis.title = ggplot2::element_text(size = 30),
                           axis.text = ggplot2::element_text(size = 20), 
                           legend.title = ggplot2::element_text(size = 30),
                           legend.text = ggplot2::element_text(size = 20),
                           panel.background = ggplot2::element_rect(fill = NA, 
                                                                    linewidth = 1.5))

        return(plt)
    } 
)

plt <- ggpubr::ggarrange(
    plotlist = plots,
    nrow = 1,
    legend = "bottom",
    common.legend = TRUE
)

ggplot2::ggsave(
    file.path("figures", "study 1", "unsystematic error, autocorrelation.png"),
    plt,
    width = 4 * 500 * 3, 
    height = 650 * 3,
    unit = "px"
)

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
#
# This dynamic approach will be explored next

# Let's create a function to estimate the parameters of the VAR(1)
autoregression <- function(x) {
    # Return all NAs if not enough data is provided
    if(nrow(x) < 11) {
        return(rep(NA, 10))
    }

    # Arrange the variables according to time
    x <- x %>% 
        dplyr::arrange(time)

    # Prepare the variables. Delete the first observation from the Y matrix and 
    # the last observation from the X matrix
    Y <- x %>% 
        dplyr::mutate(
            x = ifelse(time == min(time), NA, x), 
            y = ifelse(time == min(time), NA, y)
        ) %>% 
        dplyr::filter(!is.na(x)) %>% 
        dplyr::select(x, y) %>% 
        as.matrix()

    X <- x %>% 
        dplyr::mutate(
            x = ifelse(time == max(time), NA, x), 
            y = ifelse(time == max(time), NA, y),
            intercept = 1
        ) %>% 
        dplyr::filter(!is.na(x)) %>% 
        dplyr::select(intercept, x, y) %>% 
        as.matrix()

    # Do least-squares
    B <- tryCatch(solve(t(X) %*% X) %*% t(X) %*% Y,
                  error = function(e) browser())

    # Compute the residuals of the model and use them to estimate the covariance
    # matrix.
    S <- cov(Y - X %*% B)

    return(
        c(
            as.numeric(t(B)),
            as.numeric(S)
        )
    )
}

# Apply the function to each of the data sets
results <- lapply(
    data_list, 
    \(x) x %>% 
        dplyr::group_by(experiment, tag) %>% 
        dplyr::mutate(
            x = x - mean(x),
            y = y - mean(y)
        ) %>% 
        dplyr::summarize(
            experiment = experiment[1], 
            tag = tag[1],
            data = cbind(x, y, time) %>% 
                as.data.frame() %>% 
                setNames(c("x", "y", "time")) %>% 
                autoregression() %>% 
                matrix(nrow = 1) %>% 
                as.data.frame()
        ) %>% 
        tidyr::unnest(data) %>% 
        dplyr::ungroup() %>% 
        setNames(
            c(
                "experiment", 
                "tag", 
                paste0("intercept", c("_x", "_y")),
                paste0("auto", c("_x", "_yx", "_xy", "_y")),
                paste0("sigma", c("_x", "_yx", "_xy", "_y"))
            )
        )
)

# Get the 99%CI for each of the parameters, providing us with some type of idea
# of what to account for
results <- do.call("rbind", results) %>% 
    dplyr::select(intercept_x:sigma_y) %>% 
    as.matrix() %>% 
    matrixStats::colQuantiles(
        probs = c(0.005, 0.5, 0.995),
        na.rm = TRUE
    )

saveRDS(
    results, 
    file.path("results", "study 1", "unsystematic error, autoregression parameters.Rds")
)





################################################################################
# OTHER STUFF OF INTEREST

# PURPOSE: Check other things that may be of interest in the stationary 
#          calibration data. This includes:
#              - Sampling rate: Will give us an idea of how many datapoints will
#                               be averaged over when binning the data

data_list <- readRDS(file.path("data", "study 1", "data_list.Rds"))

# Sampling rate ################################################################

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

# Again loop over each of the datasets. Note that in the results, you will find 
# an NA value for tag 66 in experiment "stationarity 3 - 22-12-2023". This is 
# because that tag only has a single value attached to itself, meaning that we 
# can safely discard this NA from the results
results <- lapply(
    names(data_list),
    \(x) data_list[[x]] %>% 
        dplyr::group_by(experiment, tag) %>% 
        dplyr::summarize(bin_size = sampling_rate(time)) %>% 
        dplyr::ungroup() %>% 
        # Convert the sampling rate to Hz
        dplyr::mutate(Hz = 1 / bin_size) %>% 
        dplyr::select(experiment, tag, bin_size, Hz) %>% 
        dplyr::mutate(day = x)
)
results <- do.call("rbind", results)

data.table::fwrite(
    results,
    file.path("results", "study 1", "sampling rate.csv")
)

results %>% 
    dplyr::group_by(day) %>% 
    dplyr::summarize(
        mean = mean(Hz, na.rm = TRUE),
        var = var(Hz, na.rm = TRUE), 
        q025 = quantile(Hz, 0.025, na.rm = TRUE),
        q975 = quantile(Hz, 0.975, na.rm = TRUE),
        min = min(Hz, na.rm = TRUE), 
        max = max(Hz, na.rm = TRUE)
    ) %>% 
    View()
