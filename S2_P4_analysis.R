################################################################################
# Purpose: Perform the formal analyses for study 2 and describe how we got to  #
#          our conclusions.                                                    #
################################################################################

devtools::load_all()


# Define the names of all the summary files, which contains the values of the 
# summary statistics for each simulation and condition. This will be used as the 
# basis for the main analyses.
files <- paste0(
    rep(
        c("movement", "fixed"),
        each = 6
    ),
    "_",
    rep(
        c("R10", "T10", "R6R", "T6R", "R6N", "T6N"),
        times = 2
    )
)

data <- lapply(
    files, 
    \(x) data.table::fread(
        file.path("results", "study 2", paste0("summary_", x, ".csv")), 
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
    data, 
    function(x) {
        # Separate the data from before and after the preprocessing
        before <- dplyr::filter(x, preprocessed == "before")
        after <- dplyr::filter(x, preprocessed == "after")

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
    file.path("results", "study 2", "difference_distribution.Rds")
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
            R10 = movement_R10 & movement_R10_effect < 0 & fixed_R10 & fixed_R10_effect < 0,
            R6R = movement_R6R & movement_R6R_effect < 0 & fixed_R6R & fixed_R6R_effect < 0,
            R6N = movement_R6N & movement_R6N_effect < 0 & fixed_R6N & fixed_R6N_effect < 0,
            T10 = movement_T10 & movement_T10_effect < 0 & fixed_T10 & fixed_T10_effect < 0,
            T6R = movement_T6R & movement_T6R_effect < 0 & fixed_T6R & fixed_T6R_effect < 0,
            T6N = movement_T6N & movement_T6N_effect < 0 & fixed_T6N & fixed_T6N_effect < 0
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
    file.path(".", "results", "study 2", "selected_filters.csv")
)





#-------------------------------------------------------------------------------
# Error covariance
#-------------------------------------------------------------------------------

# Read in all trajectory data
types <- paste(
    rep(c("fixed", "movement"), each = 6),
    rep(c("R10", "R6R", "R6N", "T10", "T6R", "T6N"), times = 2),
    sep = "_"
)
data_list <- lapply(
    types,
    \(x) data.table::fread(
        file = file.path("results", "study 2", paste0("trajectory_", x, ".csv")),
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
    data <- data[data$preprocessed == "after", ]

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
    file.path("results", "study 2", "unsystematic error, overall covariance.Rds")
)





#-------------------------------------------------------------------------------
# Temporal error
#-------------------------------------------------------------------------------

# Read in all trajectory data
types <- paste(
    rep(c("fixed", "movement"), each = 6),
    rep(c("R10", "R6R", "R6N", "T10", "T6R", "T6N"), times = 2),
    sep = "_"
)
data_list <- lapply(
    types,
    \(x) data.table::fread(
        file = file.path("results", "study 2", paste0("trajectory_", x, ".csv")),
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
    data <- data[data$preprocessed == "after", ]

    # Define which pipelines you're looking at
    parameters <- lapply(
        unique(data$preprocessing_function),
        function(name) {
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
        }
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
    file.path("results", "study 2", "unsystematic error, autoregression.Rds")
)





#-------------------------------------------------------------------------------
# Visualization
#-------------------------------------------------------------------------------

# Select four cases, namely movement/fixed vs R10/T6N
cases <- paste(
    rep(c("movement", "fixed"), times = 2),
    rep(c("R10", "T6N"), each = 2),
    sep = "_"
)

# Define titles for each case
titles <- list(
    "movement_R10" = "Movement",
    "movement_T6N" = "",
    "fixed_R10" = "Positional",
    "fixed_T6N" = ""
)
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
                    title = titles[[x]],
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
                    panel.grid = ggplot2::element_line(
                        color = "gray75"
                    ),
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
        nrow = 2,
        ncol = 2,
        labels = c("A", "", "B", ""),
        font.label = list(
            size = 30
        )
    )

    ggplot2::ggsave(
        file.path("figures", "study 2", paste0("result_", i, ".png")),
        plt,
        width = 5500, 
        height = 5000,
        unit = "px"
    )
}
