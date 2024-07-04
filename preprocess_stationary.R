################################################################################
# Purpose: Preprocess the stationary data sets.                                #
#                                                                              #
#          General purpose is to find out how successful we are at reducing    #
#          the error related to the stationary data sets. This concerns the    #
#          the systematic and unsystematic error, making use of the assumed    #
#          "real" positions of the tags during these experiments.              #
################################################################################


################################################################################
# PRELIMINARIES

devtools::load_all()

# Get the data that you want to preprocess.
data_files <- c("stationary_14-10-2023", 
                "stationary_21-10-2023",
                "stationary_22-12-2023")
data_list <- lapply(data_files, 
                    \(x) readRDS(file.path("data", "stationary", paste0(x, ".Rds"))))
names(data_list) <- data_files

# Also read in the data on where the anchors were located for these data.
files <- c("anchor_position_14-10-2023", 
           "anchor_position_21-10-2023",
           "anchor_position_22-12-2023")
anchor_list <- lapply(files, 
                      \(x) readRDS(file.path("data", paste0(x, ".Rds"))))
names(anchor_list) <- data_files

# Make a distinction between 6 and 4 anchor data
for(i in c(4, 6)) {
    idx <- paste0(data_files[3], "_", i)
    data_list[[idx]] <- data_list[[3]] %>% 
        dplyr::filter(anchors == i)
    anchor_list[[idx]] <- anchor_list[[3]]
}

# Adjust the data_files vector
data_files <- names(data_list)

# Add the information of the anchors to the datasets themselves. We need this 
# information if we want to correct the distortion. Only retain minimal and 
# maximal x and y coordinates (as this is what leads to its standardization)
for(i in data_files) {
    data_list[[i]] <- data_list[[i]] %>% 
        dplyr::mutate(anchor_min_x = min(anchor_list[[i]]$x), 
                      anchor_max_x = max(anchor_list[[i]]$x),
                      anchor_min_y = min(anchor_list[[i]]$y), 
                      anchor_max_y = max(anchor_list[[i]]$y))
}



################################################################################

# Create several different preprocessing pipelines to be tested. These pipelines
# consist of:
#   - Binning the data to bins of 500msec in size
#   - Handling the systematic distortion with the default parameters in the 
#     `polynomial_distortion` function
#   - Filtering the data with different kinds of filters or not filtering at all
#
# Filtering is based on the performance of the filters in the simulation study, 
# meaning that we will only use one kind of moving window approach (2, average).
# Given that the dynamical and equilibrium filters don't seem to work, I leave 
# them out for now.
columns <- c("X", "Y", "anchor_min_x", "anchor_max_x", "anchor_min_y", "anchor_max_y")
binning <- \(x) nameless::bin(x, 
                              span = 0.5, 
                              \(x) nameless::average(x, cols = columns), 
                              .by = "id")
moving <- \(x) nameless::moving_window(x, 
                                       span = 2, 
                                       fx = \(x) nameless::average(x, cols = columns), 
                                       .by = "id")
distort <- \(x) x %>% 
    dplyr::mutate(x = nameless::normalize_position(x, min_x = anchor_min_x[1], max_x = anchor_max_x[1]), 
                  y = nameless::normalize_position(y, min_x = anchor_min_y[1], max_x = anchor_max_y[1])) %>% 
    nameless::polynomial_distortion() %>% 
    dplyr::mutate(x = nameless::denormalize_position(x, min_x = anchor_min_x[1], max_x = anchor_max_x[1]), 
                  y = nameless::denormalize_position(y, min_x = anchor_min_y[1], max_x = anchor_max_y[1]))

conditions <- list("bin" = list(binning), 

                   # Binning and moving window
                   "mov_bin" = list(moving, binning), 
                   "bin_mov" = list(binning, moving),

                   # Binning and systematic distortion
                   "bin_distort" = list(binning, distort),
                   "distort_bin" = list(distort, binning),

                   # All three combined
                   "distort_mov_bin" = list(distort, moving, binning), 
                   "distort_bin_mov" = list(distort, binning, moving),

                   "mov_distort_bin" = list(moving, distort, binning), 
                   "mov_bin_distort" = list(moving, binning, distort), 
                   
                   "bin_mov_distort" = list(binning, moving, distort),
                   "bin_distort_mov" = list(binning, distort, moving))

saveRDS(conditions, file.path("results", "stationary", "conditions.Rds"))

# Combine the name of the conditions with the names of the data files
data_files <- data.frame(condition = rep(names(conditions), each = length(data_files)), 
                         filename = rep(data_files, times = length(conditions)))





################################################################################
# ANALYSIS

# Create a function that will create a local plot of the data together with the 
# actual positions of the tags (X, Y) and the anchor positions
local_plot <- function(x, title) {
    (nameless::plot(x, per_iteration = FALSE) +
        ggplot2::geom_point(ggplot2::aes(x = X, y = Y), 
                            color = "salmon", 
                            size = 0.25) +
        ggplot2::annotate("point", 
                          x = c(rep(x$anchor_min_x, 2), rep(x$anchor_max_x, 2)),
                          y = rep(c(x$anchor_min_y, x$anchor_max_y), times = 2),
                          color = "cornflowerblue") +
        ggplot2::labs(title = title) +
        ggplot2::lims(x = c(x$anchor_min_x[1] - 0.25, x$anchor_max_x[1] + 0.25), 
                      y = c(x$anchor_min_y[1] - 0.25, x$anchor_max_y[1] + 0.25))) %>% 
        suppressMessages() %>% 
        return()
}

# Create a function in which we will compute the summary statistics of interest
compute_summary_statistics <- function(x) {
    x %>% 
        # Compute the difference between filtered and expected positions. Used 
        # to measure the extent to which systematic error is present in the data
        dplyr::mutate(difference_x = X - x, 
                      difference_y = Y - y) %>% 

        # Center the measured positions around zero. Used to measure the extent 
        # to which unsystematic error is present in the data
        dplyr::group_by(experiment, id) %>%
        dplyr::mutate(centered_x = x - mean(x), 
                      centered_y = y - mean(y)) %>%  
        dplyr::ungroup() %>% 

        # Compute the statistics of interest per experiment and per tag. This 
        # will allow for a more broad view on where it still might go awry
        dplyr::group_by(experiment, id) %>% 
        dplyr::arrange(time) %>% 
        dplyr::summarize(# Statistics with respect to systematic error
                         mean_diff_x = mean(difference_x, na.rm = TRUE), 
                         mean_diff_y = mean(difference_y, na.rm = TRUE), 
                         sd_diff_x = sd(difference_x, na.rm = TRUE), 
                         sd_diff_y = sd(difference_y, na.rm = TRUE), 
                         q025_diff_x = quantile(difference_x, probs = 0.025, na.rm = TRUE),
                         q025_diff_y = quantile(difference_y, probs = 0.025, na.rm = TRUE),
                         q975_diff_x = quantile(difference_x, probs = 0.975, na.rm = TRUE),
                         q975_diff_y = quantile(difference_y, probs = 0.975, na.rm = TRUE), 
                         
                         # Statistics with regard to unsystematic error
                         var_x = var(centered_x), 
                         var_y = var(centered_y), 
                         cov_xy = cov(centered_x, centered_y),
                         cor_xy = cor(centered_x, centered_y),
                         
                         # Autocorrelation in the residuals
                         auto_x = cor(centered_x[2:length(centered_x)], 
                                      centered_x[2:length(centered_x) - 1]), 
                         auto_y = cor(centered_y[2:length(centered_y)], 
                                      centered_y[2:length(centered_y) - 1])) %>% 
        dplyr::ungroup() %>% 
        return()
}

# For the analysis, we should first create a function that will do the 
# preprocessing and check its efficacy. This will allow us to put tidyverse to 
# its maximal use when preprocessing the data.
pipeline_efficacy <- function(x){

    print(paste0(x$filename, ": ", x$condition))

    # Retrieve the data and the pipeline for the condition
    local_data <- data_list[[x$filename]]
    fx <- conditions[[x$condition]]

    # Visualize the funfiltered data. Will be used later to combine with the 
    # filtered data. 
    plots <- list(local_plot(local_data, "Original"))

    # Compute the summary statistics of the data before they are processed 
    # through the pipeline. This will give us values to compare the results 
    # to, which is an overall better approach. Add an indicator that tells us 
    # that this is the original data
    summary_statistics <- compute_summary_statistics(local_data) %>% 
        dplyr::mutate(type = "original")

    # Execute the pipeline for each data set at hand and compare the values to 
    # the presumed "real" positions contained in X and Y. Furthermore compute 
    # some values for the statistics that were computed for the original data 
    # as well, namely autocorrelation and covariance 
    local_data <- local_data %>% 
        # Necessary to group per experiment to ensure that tags are not doubly 
        # assessed
        dplyr::group_by(experiment) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(data = data %>% 
                          as.data.frame() %>% 
                          nameless::execute_pipeline(fx) %>% 
                          list()) %>% 
        tidyr::unnest(data) %>% 
        dplyr::ungroup()

    # Visualize the filtered data. Together with the unfiltered data, this should 
    # help us identify whether the results we are getting make sense
    plots[[2]] <- local_plot(local_data, "Preprocessed")
    ggplot2::ggsave(file.path("figures", 
                              "stationary", 
                              "preprocessed", 
                              paste0(x$filename, "_", x$condition, ".jpg")), 
                    ggpubr::ggarrange(plotlist = plots, 
                                      ncol = 2), 
                    width = 2000, 
                    height = 1200, 
                    unit = "px")

    # Compute the statistics as explained above
    summary_statistics <- compute_summary_statistics(local_data) %>% 
        dplyr::mutate(type = "preprocessed") %>% 
        rbind(summary_statistics)
    
    # Bind the summary statistics to the information in the data files
    summary_statistics <- x %>% 
        dplyr::mutate(data = list(summary_statistics)) %>% 
        tidyr::unnest(data)

    return(summary_statistics)
}

# Define the number of cores to run this on in parallel
n_cores <- max(c(parallel::detectCores() - 1, 1))

# Execute this function for each combination of the data and the pipeline.
results <- data_files %>%
    split(seq_len(nrow(data_files))) %>%
    as.list() %>% 
    parallel::mclapply(pipeline_efficacy,
                       mc.cores = n_cores) %>%
    dplyr::bind_rows() %>% 
    dplyr::relocate(condition:id, type)

data.table::fwrite(results, 
                   file.path("results", "stationary_preprocessing.csv"))





################################################################################
# VISUALIZATION

# Load the needed file
results <- data.table::fread(file.path("results", "stationary_preprocessing.csv"))

# Create a function that will plot distributional histograms for the different 
# statistics and dimensions
make_plot <- function(x, 
                      col) {

    # Extract the different data sets and the different preprocessing strategies
    dataset <- unique(x$filename)
    preprocess <- unique(x$condition)

    # For clarity, get rid of extreme outliers
    y <- x %>% 
        dplyr::select(as.character(col)) %>% 
        unlist() %>% 
        as.numeric()
    idx <- abs(y) <= 3 * sd(y)
    x <- x[idx,]

    # Get the limits of the plots
    limits <- x %>% 
        dplyr::select(as.character(col)) %>% 
        range(na.rm = TRUE)
    limits <- limits + c(-1, 1) * diff(limits) * 0.1

    # Loop over these unique names and create the plots. Put these plots in a 
    # list to be bound together using ggpubr
    plt <- list()
    for(i in seq_along(preprocess)) {
        for(j in seq_along(dataset)) {
            # Select the data
            plot_data <- x %>% 
                dplyr::filter(filename == dataset[j], 
                              condition == preprocess[i]) %>% 
                dplyr::select(as.character(col), type) %>% 
                setNames(c("X", "M"))

            # Create the plot
            idx <- (i - 1) * length(dataset) + j
            plt[[idx]] <- ggplot2::ggplot(data = plot_data, 
                                          ggplot2::aes(x = X, fill = M)) +
                ggplot2::geom_histogram(color = "black",
                                        alpha = 0.5,
                                        bins = 15, 
                                        position = "identity") +
                ggplot2::scale_fill_manual(values = c("original" = "salmon", 
                                                      "preprocessed" = "cornflowerblue")) +
                ggplot2::labs(title = ifelse(i == 1, dataset[j], " "), 
                              y = ifelse(j == 1, preprocess[i], " "), 
                              x = " ") +
                ggplot2::lims(x = limits) +
                ggplot2::theme_minimal() +
                ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, 
                                                                    angle = 90, 
                                                                    vjust = 0.5, 
                                                                    hjust = 0.5), 
                               plot.title = ggplot2::element_text(size = 15, 
                                                                  hjust = 0.5),
                               axis.text = ggplot2::element_text(size = 10, 
                                                                 angle = 45, 
                                                                 hjust = 1))
        }
    }

    # Bind together and return the resulting plot
    plt <- ggpubr::ggarrange(plotlist = plt,
                             nrow = length(preprocess), 
                             ncol = length(dataset),
                             common.legend = TRUE,
                             legend = "bottom")

    return(plt)
}

# Loop over them and save the plots
col <- colnames(results)[-c(1:5)]
for(i in col) {
    ggplot2::ggsave(file.path("figures", "stationary", paste0(i, ".jpg")), 
                    make_plot(results, i), 
                    width = length(unique(results$filename)) * 750, 
                    height = length(unique(results$condition)) * 900, 
                    unit = "px")
}

# Results:
#   - Fixing the distortion also seems to decrease the amount of absolute error
#     in the data. This is true for all datasets, but is least obvious for the 
#     21-10-2023 data. This might be due to the presence of outliers after 
#     fixing the distortions, which in turn can be explained by the fact that we
#     don't have the exact anchor positions for this dataset.
#       - The degree of the absolute error across tags also decreases, as shown
#         by the sd_diff
#   - Using a moving window increases the amount of autocorrelation in the 
#     residuals
#   - The correlation between the x and y dimensions is mostly left untouched
#     by the transformations
#   - Using a moving window seems to decrease the unsystematic error most. 
#     Binning does not really have an effect.
#
#   - With regard to order, the following seem to be true
#       1) When you apply the distortion filter does not seem to matter for any
#          of the statistics
#       2) The moving window seems to decrease the unsystematic error most when
#          you first apply binning.
#       -> Based on these results, pipeline can look like `distort`, `binning`,
#          `moving`
#
#   TO DO
#       - Find a way to quantify this a bit stronger (e.g., significance levels
#         or Bayesian decomposition with BF)