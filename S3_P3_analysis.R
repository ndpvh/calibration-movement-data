################################################################################
# Purpose: Perform the formal analyses for study 3 and describe how we got to  #
#          our conclusions.                                                    #
################################################################################

devtools::load_all()


# Define the names of all the summary files, which contains the values of the 
# summary statistics for each simulation and condition. This will be used as the 
# basis for the main analyses.
files <- paste0(
    "data_",
    c("R10", "T10", "R6R", "T6R", "R6N", "T6N")
)

data <- lapply(
    files, 
    \(x) data.table::fread(
        file.path("results", "study 3", paste0("summary_", x, ".csv")), 
        data.table = FALSE
    )
) %>% 
    `names<-` (files)

# Define the columns in which we are interested
columns <- c("bias_dist", "rmse_dist")

# Loop over the different datafiles and compare the values of the summary 
# statistics for each of the conditions to the values of these same statistics
# before any preprocessing was done. 
distribution <- lapply(
    names(data), 
    function(x) {
        print(x)

        # Separate the data from before and after the preprocessing
        before <- dplyr::filter(data[[x]], preprocessed == "before")
        after <- dplyr::filter(data[[x]], preprocessed == "after")

        # Nest the data based on the unique conditions in the 
        # preprocessed data
        after <- after %>%
            dplyr::group_by(preprocessing_function) %>% 
            tidyr::nest()

        # Loop over the different columns and do the comparisons 
        # we want to do.
        result <- lapply(
            columns, 
            function(y) {
                # Create reference
                reference <- before[, y]
 
                # Use mutate on the nested data
                after <- after %>% 
                    dplyr::rowwise() %>% 
                    dplyr::mutate(
                        data_2 = data %>%
                            as.data.frame() %>% 
                            dplyr::select(tidyselect::matches(y)) %>% 
                            unlist() %>% 
                            as.numeric() %>% 
                            nameless::compare_distribution(
                                reference, 
                                bootstrapped = 10000
                            ) %>% 
                            list()) %>% 
                    dplyr::ungroup() %>% 
                    dplyr::select(-data) %>% 
                    tidyr::unnest(data_2) 

                nameless::compare_distribution(
                    reference,
                    reference, 
                    bootstrapped = 10000
                ) %>% 
                    dplyr::mutate(preprocessing_function = "before") %>% 
                    rbind(after) %>% 
                    dplyr::relocate(preprocessing_function) %>% 
                    return()
            }
        )

        # Bind together into one big dataframe with an additional 
        # column specifying the statistic that it used
        statistics <- sapply(
            seq_along(result), 
            \(i) rep(columns[i], nrow(result[[i]]))
        ) %>% 
            as.character() 
        result <- do.call("rbind", result)
        result$statistic <- statistics

        return(result)
    }
) %>% 
    `names<-` (names(data))

saveRDS(
    distribution, 
    file.path("results", "study 3", "difference_distribution.Rds")
)





#-------------------------------------------------------------------------------
# Selecting pipelines
#-------------------------------------------------------------------------------

# Let's check which conditions were significant for the "fixed" and "movement"
# datasets
significant <- data.frame(
    preprocessing_function = distribution[[1]]$preprocessing_function,
    statistics = distribution[[1]]$statistic
)
for(i in names(data)) {
    significant[, i] <- distribution[[i]]$significant
    significant[, paste0(i, "_effect")] <- distribution[[i]]$ci_diff_upper
}

# Inspect the results visually and try to come up with some pipelines to use in
# the second simulation. 
#
# To select the pipelines, we inspect the statistics of interest and use the
# following criteria for the inclusion of a singular pipeline:
# on the following criteria:
#   - Negative value of the statistic (i.e., measurement became better)
#   - Significance of the effect 
#   - Success in at least 2/3 of the datasets. This is an arbitrary criterion, 
#     but ensures some pipelines get picked even when all pipelines performed 
#     bad in some datasets (e.g., time-related error has no significant pipelines
#     that reduced the mean distance)
statistics <- c("bias_dist", "rmse_dist")
results <- data.frame(
    preprocessing_function = significant %>% 
        dplyr::filter(statistics == "bias_dist") %>% 
        dplyr::select(preprocessing_function) %>% 
        unlist() %>% 
        as.vector()
)

for(i in statistics) {
    tmp <- significant %>% 
        dplyr::filter(statistics == i) %>% 
        dplyr::mutate(
            R10 = data_R10 & data_R10_effect < 0,
            R6R = data_R6R & data_R6R_effect < 0,
            R6N = data_R6N & data_R6N_effect < 0,
            T10 = data_T10 & data_T10_effect < 0,
            T6R = data_T6R & data_T6R_effect < 0,
            T6N = data_T6N & data_T6N_effect < 0
        ) %>%
        dplyr::rowwise() %>% 
        dplyr::mutate(
            selected = sum(dplyr::across(R10:T6N)) / 6,
            selected = selected >= 2/3
        ) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(-statistics) %>% 
        dplyr::relocate(
            preprocessing_function,
            selected
        )

    tmp <- setNames(
        tmp, 
        c(
            "preprocessing_function", 
            i,
            paste0(i, "_", colnames(tmp[, -c(1:2)]))
        )
    )

    results <- dplyr::full_join(
        results, 
        tmp, 
        by = "preprocessing_function"
    )    
}

results <- dplyr::relocate(
    results,
    preprocessing_function,
    bias_dist,
    rmse_dist
)

View(results)
data.table::fwrite(
    results, 
    file.path(".", "results", "study 3", "selected_filters.csv")
)

# For some reason, LOESS/LOWESS performs really bad now!




#-------------------------------------------------------------------------------
# Some descriptives
#-------------------------------------------------------------------------------

# Create one big file of distributional differences
results <- lapply(
    names(distribution),
    \(x) distribution[[x]] %>% 
        dplyr::mutate(
            file = x,
            type = stringr::str_split(x, pattern = "_")[[1]][1],
            error = stringr::str_split(x, pattern = "_")[[1]][2]
        )
)
results <- do.call("rbind", results)

# Look at mean reductions of the good pipelines (Kalman, time, and combo)
results[results$preprocessing_function == "kalm", ] %>% 
    dplyr::group_by(type, statistic) %>% 
    dplyr::summarize(
        mean = mean(mean_diff, na.rm = TRUE),
        ci_lower = mean(ci_diff_lower, na.rm = TRUE), 
        ci_upper = mean(ci_diff_upper, na.rm = TRUE)
    )

results[results$preprocessing_function == "time", ] %>% 
    dplyr::group_by(type, statistic) %>% 
    dplyr::summarize(
        mean = mean(mean_diff, na.rm = TRUE),
        ci_lower = mean(ci_diff_lower, na.rm = TRUE), 
        ci_upper = mean(ci_diff_upper, na.rm = TRUE)
    )

results[results$preprocessing_function %in% c("time_kalm", "kalm_time"), ] %>% 
    dplyr::group_by(type, statistic) %>% 
    dplyr::summarize(
        mean = mean(mean_diff, na.rm = TRUE),
        ci_lower = mean(ci_diff_lower, na.rm = TRUE), 
        ci_upper = mean(ci_diff_upper, na.rm = TRUE)
    )

results[!c(results$preprocessing_function %in% c("before", "time", "kalm", "kalm_time", "time_kalm")), ] %>% 
    dplyr::group_by(type, statistic) %>% 
    dplyr::summarize(
        mean = mean(mean_diff, na.rm = TRUE),
        ci_lower = mean(ci_diff_lower, na.rm = TRUE), 
        ci_upper = mean(ci_diff_upper, na.rm = TRUE)
    )

results[results$preprocessing_function == "before", ] %>% 
    dplyr::group_by(type, statistic) %>% 
    dplyr::summarize(
        mean = mean(mean_x, na.rm = TRUE),
        ci_lower = mean(ci_x_lower, na.rm = TRUE), 
        ci_upper = mean(ci_x_upper, na.rm = TRUE)
    )






#-------------------------------------------------------------------------------
# ANOVA
#-------------------------------------------------------------------------------

# Define the names of all the summary files, which contains the values of the 
# summary statistics for each simulation and condition. This will be used as the 
# basis for the main analyses.
files <- paste(
    "data",
    c("R10", "T10", "R6R", "T6R", "R6N", "T6N"),
    sep = "_"
)

# Create the data with differences in it already
data <- lapply(
    files, 
    function(x) {
        # Read in the datafile and add information on the type of file it is
        data <- data.table::fread(
            file.path("results", "study 3", paste0("summary_", x, ".csv")), 
            data.table = FALSE
        ) %>% 
            dplyr::mutate(
                file = x,
                type = stringr::str_split(x, pattern = "_")[[1]][1],
                error = stringr::str_split(x, pattern = "_")[[1]][2]
            )

        # Create a reference and test to be used to create differences
        reference <- data[data$preprocessed == "before", ]
        test <- data[data$preprocessed != "before", ]

        # Once done, we can make differences as follows
        test %>% 
            dplyr::group_by(preprocessing_function) %>% 
            tidyr::nest() %>% 
            dplyr::mutate(
                data_2 = data[[1]] %>% 
                    dplyr::select(
                        nsim, 
                        id,
                        bias_dist, 
                        rmse_dist, 
                        mae_dist
                    ) %>% 
                    dplyr::rename(
                        bias_test = bias_dist,
                        rmse_test = rmse_dist,
                        mae_test = mae_dist
                    ) %>% 
                    dplyr::full_join(
                        reference,
                        by = c("nsim", "id")
                    ) %>% 
                    dplyr::mutate(
                        bias_diff = bias_test - bias_dist, 
                        rmse_diff = rmse_test - rmse_dist,
                        mae_diff = mae_test - mae_dist
                    ) %>% 
                    dplyr::select(
                        nsim, 
                        id,
                        bias_diff, 
                        rmse_diff, 
                        mae_diff,
                        file, 
                        type, 
                        error
                    ) %>% 
                    list()
            ) %>% 
            dplyr::select(-data) %>% 
            tidyr::unnest(data_2) %>% 
            dplyr::ungroup() %>% 
            return()
    }
)  
data <- do.call("rbind", data)

# Add all types of dummy variables to the mix, allowing us to create an ANOVA-type
# of analysis
components <- c(
    "av", 
    "idx", 
    "time", 
    "kalm", 
    "loess-1", 
    "loess-2",
    "loess-3"
)
for(i in components) {
    data[, i] <- stringr::str_detect(data$preprocessing_function, i)
}

data$moving_window <- rowSums(data[, c("av", "idx", "time")]) > 0
data$loess <- rowSums(data[, c("loess-1", "loess-2", "loess-3")]) > 0

# Perform ANOVAs of interest for the data: Not a lot of interesting stuff that
# comes out, other than Kalman filter having a non-significant effect
model <- aov(
    bias_diff ~ kalm, 
    data = data
) 
lsr::etaSquared(model)

model <- aov(
    rmse_diff ~ moving_window * loess * kalm, 
    data = data
) 
lsr::etaSquared(model)







#-------------------------------------------------------------------------------
# Error covariance
#-------------------------------------------------------------------------------

# Read in all trajectory data
types <- paste(
    "data",
    c("R10", "R6R", "R6N", "T10", "T6R", "T6N"),
    sep = "_"
)
data_list <- lapply(
    types,
    \(x) data.table::fread(
        file = file.path("results", "study 3", paste0("trajectory_", x, ".csv")),
        data.table = FALSE
    )
) %>% 
    `names<-` (types)

# Loop over the different dates and do all your estimation 
results <- list()
for(i in seq_along(data_list)){
    print(names(data_list)[i])

    # Load data from specific date and change it somewhat
    data <- data_list[[i]] %>% 
        dplyr::mutate(
            x = x - x_actual, 
            y = y - y_actual
        )
    data <- data[data$preprocessing_function == "kalm", ]

    # Define which pipelines you're looking at
    covariances <- lapply(
        unique(data$preprocessing_function),
        function(name) {
            # Bootstrap the data using this function and immediately compute the necessary
            # summary statistics: 2 variances and 1 covariance.
            data %>% 
                dplyr::group_by(nsim, preprocessing_function, id) %>% 
                dplyr::summarize(
                    var_x = var(x, na.rm = TRUE),
                    var_y = var(y, na.rm = TRUE), 
                    cov_xy = cov(x, y, use = "pairwise.complete.obs")
                ) %>% 
                dplyr::ungroup() %>% 
                dplyr::group_by(preprocessing_function) %>% 
                tidyr::nest() %>% 
                dplyr::mutate(
                    data_2 = data[[1]] %>% 
                        dplyr::ungroup() %>% 
                        dplyr::summarize(
                            lb_var_x = quantile(var_x, 0.025), 
                            lb_var_y = quantile(var_y, 0.025),
                            lb_cov_xy = quantile(cov_xy, 0.025), 
                            m_var_x = mean(var_x), 
                            m_var_y = mean(var_y), 
                            m_cov_xy = mean(cov_xy), 
                            ub_var_x = quantile(var_x, 0.975), 
                            ub_var_y = quantile(var_y, 0.975),
                            ub_cov_xy = quantile(cov_xy, 0.975)
                        ) %>% 
                        unlist() %>% 
                        matrix(nrow = 3, ncol = 3) %>% 
                        as.data.frame() %>% 
                        setNames(c("lb", "mean", "ub")) %>% 
                        cbind(covariance = c("var_x", "var_y", "cov_xy")) %>% 
                        list()
                ) %>% 
                dplyr::select(-data) %>% 
                tidyr::unnest(data_2) %>% 
                dplyr::ungroup() %>% 
                suppressWarnings() %>% 
                suppressMessages() %>% 
                return() 
        }
    )

    # Add to the results list
    results[[names(data_list)[i]]] <- do.call("rbind", covariances)

    # Release the memory that is held up by the bootstrapped data and the data 
    # itself.
    rm(data, covariances)
}

# Save the results
saveRDS(
    results, 
    file.path("results", "study 3", "unsystematic error, overall covariance.Rds")
)

# Bind all together
results <- lapply(
    names(results),
    \(x) results[[x]] %>% 
        dplyr::mutate(
            file = x,
            type = stringr::str_split(x, pattern = "_")[[1]][1],
            error = stringr::str_split(x, pattern = "_")[[1]][2]
        )
)
results <- do.call("rbind", results)

# Check mean values for the error. Do separately for the different types, so 
# to be able to filter out the problematic pipelines
results %>% 
    dplyr::group_by(covariance) %>% 
    dplyr::summarize( 
        lb = mean(lb),
        mean = format(mean(mean), nsmall = 6), 
        ub = mean(ub)
    ) %>% 
    View()




#-------------------------------------------------------------------------------
# Temporal error
#-------------------------------------------------------------------------------

# Read in all trajectory data
types <- paste(
    "data",
    c("R10", "R6R", "R6N", "T10", "T6R", "T6N"),
    sep = "_"
)
data_list <- lapply(
    types,
    \(x) data.table::fread(
        file = file.path("results", "study 3", paste0("trajectory_", x, ".csv")),
        data.table = FALSE
    )
) %>% 
    `names<-` (types)

# Let's create a function to estimate the parameters of the VAR(1)
autoregression <- function(x) {
    # Return all NAs if not enough data is provided
    if(nrow(x) < 9) {
        return(rep(NA, 8))
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
            y = ifelse(time == max(time), NA, y)
        ) %>% 
        dplyr::filter(!is.na(x)) %>% 
        dplyr::select(x, y) %>% 
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

# Loop over the different dates and do all your estimation 
dims <- c("_x", "_yx", "_xy", "_y")
results <- list()
for(i in seq_along(data_list)){
    print(names(data_list)[i])

    # Load data from specific date and change it somewhat
    data <- data_list[[i]] %>% 
        dplyr::mutate(
            x = x - x_actual, 
            y = y - y_actual
        )
    data <- data[data$preprocessing_function == "kalm", ]

    # Define which pipelines you're looking at
    parameters <- parallel::mclapply(
        unique(data$preprocessing_function),
        function(name) {
            print(name)

            # Bootstrap the data using this function and immediately compute the necessary
            # summary statistics: 2 variances and 1 covariance.
            data %>% 
                dplyr::group_by(nsim, preprocessing_function, id) %>% 
                tidyr::nest() %>% 
                dplyr::mutate(
                    data_2 = data[[1]] %>% 
                        autoregression() %>% 
                        matrix(nrow = 1) %>% 
                        as.data.frame() %>% 
                        setNames(
                            c(
                                paste0("auto", dims),
                                paste0("sigma", dims)
                            )
                        ) %>% 
                        list()
                ) %>% 
                dplyr::select(-data) %>% 
                tidyr::unnest(data_2) %>% 
                dplyr::ungroup() %>% 
                dplyr::group_by(preprocessing_function) %>% 
                tidyr::nest() %>% 
                dplyr::mutate(
                    data_2 = data[[1]] %>% 
                        dplyr::ungroup() %>% 
                        dplyr::summarize(
                            lb_auto_x = quantile(auto_x, 0.025), 
                            lb_auto_yx = quantile(auto_yx, 0.025), 
                            lb_auto_xy = quantile(auto_xy, 0.025), 
                            lb_auto_y = quantile(auto_y, 0.025),
                            lb_sigma_x = quantile(sigma_x, 0.025), 
                            lb_sigma_yx = quantile(sigma_yx, 0.025), 
                            lb_sigma_xy = quantile(sigma_xy, 0.025), 
                            lb_sigma_y = quantile(sigma_y, 0.025),

                            m_auto_x = mean(auto_x), 
                            m_auto_yx = mean(auto_yx), 
                            m_auto_xy = mean(auto_xy), 
                            m_auto_y = mean(auto_y),
                            m_sigma_x = mean(sigma_x), 
                            m_sigma_yx = mean(sigma_yx),
                            m_sigma_xy = mean(sigma_xy),
                            m_sigma_y = mean(sigma_y),  

                            ub_auto_x = quantile(auto_x, 0.975), 
                            ub_auto_yx = quantile(auto_yx, 0.975), 
                            ub_auto_xy = quantile(auto_xy, 0.975), 
                            ub_auto_y = quantile(auto_y, 0.975),
                            ub_sigma_x = quantile(sigma_x, 0.975), 
                            ub_sigma_yx = quantile(sigma_yx, 0.975),
                            ub_sigma_xy = quantile(sigma_xy, 0.975),
                            ub_sigma_y = quantile(sigma_y, 0.975)
                        ) %>% 
                        unlist() %>% 
                        matrix(nrow = 8, ncol = 3) %>% 
                        as.data.frame() %>% 
                        setNames(c("lb", "mean", "ub")) %>% 
                        cbind(
                            covariance = c(
                                paste0("auto", dims),
                                paste0("sigma", dims)
                            )
                        ) %>% 
                        list()
                ) %>% 
                dplyr::select(-data) %>% 
                tidyr::unnest(data_2) %>% 
                dplyr::ungroup() %>% 
                suppressWarnings() %>% 
                suppressMessages() %>% 
                return() 
        },
        mc.cores = 3
    )

    # Add to the results list
    results[[names(data_list)[i]]] <- do.call("rbind", parameters)

    # Release the memory that is held up by the bootstrapped data and the data 
    # itself.
    rm(data, parameters)
}

# Save the results
saveRDS(
    results, 
    file.path("results", "study 3", "unsystematic error, autoregression.Rds")
)

# Bind all together
results <- lapply(
    names(results),
    \(x) results[[x]] %>% 
        dplyr::mutate(
            file = x,
            type = stringr::str_split(x, pattern = "_")[[1]][1],
            error = stringr::str_split(x, pattern = "_")[[1]][2]
        )
)
results <- do.call("rbind", results)

# Check mean values for the error. Do separately for the different types, so 
# to be able to filter out the problematic pipelines
results %>% 
    dplyr::filter(stringr::str_detect(error, "T")) %>% 
    dplyr::group_by(covariance) %>% 
    dplyr::summarize( 
        lb = mean(lb),
        mean = mean(mean), 
        ub = mean(ub)
    ) %>% 
    View()





#-------------------------------------------------------------------------------
# Visualization of difference distributions
#-------------------------------------------------------------------------------

# Select four cases, namely movement/fixed vs R10/T6N
cases <- paste(
    "data",
    c("R10", "T6N"),
    sep = "_"
)

# Define titles for each case
statistic <- list(
    "rmse_dist" = "RMSE",
    "bias_dist" = "Bias"
)

# Loop over each of these cases and create a plot
for(i in c("rmse_dist", "bias_dist")) {
    plt <- lapply(
        cases,
        function(x) {
            # Select the data
            data <- distribution[[x]]

            # Adjust the data 
            data <- data[data$statistic == i, ] %>% 
                dplyr::mutate(
                    M = preprocessing_function == "before",
                    M = factor(M)
                ) %>% 
                dplyr::select(
                    preprocessing_function,
                    M,
                    mean_x, 
                    ci_x_lower,
                    ci_x_upper
                ) %>% 
                dplyr::rename(
                    y = mean_x,
                    ymin = ci_x_lower,
                    ymax = ci_x_upper
                ) %>% 
                dplyr::arrange(y)

            levels <- rev(data$preprocessing_function)
            
            # With the data defined, we can create a plot
            plt <- ggplot2::ggplot(
                data = data,
                ggplot2::aes(
                    x = factor(preprocessing_function, levels = levels),
                    y = y,
                    ymin = y,
                    ymax = ymax,
                    fill = M
                )
            ) +
                ggplot2::geom_errorbar(
                    linewidth = 0.5,
                    width = 0.5
                ) +
                ggplot2::geom_bar(
                    stat = "identity"
                ) +
                ggplot2::coord_flip() +
                ggplot2::scale_fill_manual(
                    values = c(
                        "TRUE" = "salmon",
                        "FALSE" = "cornflowerblue"
                    )
                ) +
                ggplot2::labs(
                    title = "",
                    y = statistic[[i]],
                    x = ""
                ) +
                ggplot2::scale_y_continuous(expand = c(0, 0)) +
                ggplot2::theme(
                    panel.background = ggplot2::element_rect(
                        fill = "white"
                    ),
                    panel.border = ggplot2::element_rect(
                        fill = NA,
                        color = "black",
                        linewidth = 1.5
                    ),
                    panel.grid.major.x = ggplot2::element_line(
                        color = "gray75"
                    ),
                    panel.grid.major.y = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank(),
                    plot.title = ggplot2::element_text(
                        hjust = 0.5, 
                        size = 30
                    ),
                    axis.text.y = ggplot2::element_text(size = 8),
                    axis.title = ggplot2::element_text(size = 25),
                    legend.position = "none"
                )

            return(plt)
        }
    )

    # Bind all plots together
    plt <- ggpubr::ggarrange(
        plotlist = plt,
        nrow = 1,
        ncol = 2
    )

    ggplot2::ggsave(
        file.path("figures", "study 3", paste0("result_", i, ".png")),
        plt,
        width = 5500, 
        height = 2500,
        unit = "px"
    )
}






#-------------------------------------------------------------------------------
# Visualization of movement in x- and y-direction
#-------------------------------------------------------------------------------

# Read in data from all types of error
files <- paste(
    "data", 
    c("R10", "R6R", "R6N", "T10", "T6R", "T6N"),
    sep = "_"
)

data <- lapply(
    files, 
    \(x) data.table::fread(
        file.path("data", "study 3", paste0(x, ".csv")), 
        data.table = FALSE
    ) %>% 
        dplyr::filter(nsim == 1)
) %>%  
    `names<-` (files)

# With the data read in, we can start visualizing the trajectories in the x and 
# y direction for a given participant (just some type of visualization is okay)
id <- "csvbe"
for(i in seq_along(data)) {
    # Select the relevant data
    plot_data <- data[[i]][data[[i]]$id == id, ] %>% 
        dplyr::select(x, y, time)
    
    # Loop over the columns of interest
    cols <- c("x", "y")
    plt <- lapply(
        cols, 
        function(x) {
            tmp <- plot_data[, c("time", x)] %>% 
                setNames(c("X", "Y"))

            plt <- ggplot2::ggplot(
                data = tmp, 
                ggplot2::aes(
                    x = X, 
                    y = Y
                )
            ) +
                ggplot2::geom_line(
                    col = "black", 
                    linewidth = 2
                ) + 
                ggplot2::labs(
                    x = "Time", 
                    y = x
                ) +
                ggplot2::theme(
                    panel.background = ggplot2::element_rect(
                        fill = "white"
                    ),
                    panel.border = ggplot2::element_rect(
                        fill = NA,
                        color = "black",
                        linewidth = 1.5
                    ),
                    panel.grid.major = ggplot2::element_line(
                        color = "gray75"
                    ),
                    plot.title = ggplot2::element_text(
                        hjust = 0.5, 
                        size = 30
                    ),
                    axis.text = ggplot2::element_text(size = 8),
                    axis.title = ggplot2::element_text(size = 25)
                )
            
            return(plt)
        }
    )

    # Bind together and save
    plt <- ggpubr::ggarrange(
        plotlist = plt, 
        nrow = 2
    )

    ggplot2::ggsave(
        file.path("figures", "study 3", paste0(files[i], "_timeseries.png")),
        plt,
        width = 5500, 
        height = 3500, 
        unit = "px"
    )

    # Let's also look at the actual movements
    plot_data <- cbind(
        plot_data[2:nrow(plot_data), c("x", "y")], 
        plot_data[2:nrow(plot_data) - 1, c("x", "y")]
    ) %>% 
        setNames(c("x", "y", "xend", "yend"))

    plt_2 <- ggplot2::ggplot(
        data = plot_data, 
        ggplot2::aes(
            x = x, 
            y = y,
            xend = xend, 
            yend = yend
        )
    ) +
        ggplot2::geom_segment(
            col = "black", 
            linewidth = 2
        ) + 
        ggplot2::labs(
            x = "x", 
            y = "y"
        ) +
        ggplot2::coord_equal() +
        ggplot2::theme(
            panel.background = ggplot2::element_rect(
                fill = "white"
            ),
            panel.border = ggplot2::element_rect(
                fill = NA,
                color = "black",
                linewidth = 1.5
            ),
            panel.grid.major = ggplot2::element_line(
                color = "gray75"
            ),
            plot.title = ggplot2::element_text(
                hjust = 0.5, 
                size = 30
            ),
            axis.text = ggplot2::element_text(size = 8),
            axis.title = ggplot2::element_text(size = 25)
        )

    ggplot2::ggsave(
        file.path("figures", "study 3", paste0(files[i], "_space.png")),
        plt_2
    )

    # Bind them together, because why not
    plt <- ggpubr::ggarrange(
        plotlist = list(
            plt_2 + ggplot2::theme(aspect.ratio = 1), 
            plt
        ),
        ncol = 2, 
        widths = c(1/3, 2/3),
        labels = c("A", "B"),
        font.label = list(
            size = 30
        )
    )

    ggplot2::ggsave(
        file.path("figures", "study 3", paste0(files[i], ".png")),
        plt,
        width = 8000 * 0.7, 
        height = 3500 * 0.7,
        units = "px"
    )
}
