################################################################################
# Purpose: Preprocess the synthetic data sets.                                 #
#                                                                              #
#          General purpose is to find out under what conditions we can         #
#          successfully get rid of unsystematic measurement error. This        #
#          pertains to both the preprocessing strategy (binning, moving        #
#          window, filtering,...) and the nature of the error (random/         #
#          nonrandom, related/unrelated,...).                                  #
#                                                                              #
#          Based on the results of this analysis, we can then decide on a      #
#          strategy to use on the actual data, hopefully ridding it from       #
#          a sufficient amount of unsystematic measurement error.              #
################################################################################

devtools::load_all()

################################################################################
# PRELIMINARIES

# Get the data that you want to preprocess.
data_files <- c("synthetic_unrelated_10",
                "synthetic_related_10",
                "synthetic_temporal_10",
                "synthetic_unrelated_6random", 
                "synthetic_related_6random",
                "synthetic_temporal_6random",
                "synthetic_unrelated_6nonrandom",
                "synthetic_related_6nonrandom",
                "synthetic_temporal_6nonrandom")
data_list <- lapply(data_files, 
                    \(x) data.table::fread(file.path("data", "synthetic", paste0(x, ".csv")), 
                                           data.table = FALSE))
names(data_list) <- data_files

saveRDS(data_list, file.path("results", "synthetic", "data_list.Rds"))

# Create a dataframe that the others will be compared to. This contains the 
# original, non-error-containing binned positions of the agents. In an ideal 
# scenario, the error-containing data sets should provide values of the positions
# that lie as close to this one as possible.
data_original <- data.table::fread(file.path("data", "synthetic", "synthetic_original.csv"), 
                                   data.table = FALSE) %>% 
    dplyr::group_by(nsim) %>% 
    tidyr::nest() %>% 
    dplyr::mutate(new_data = data %>% 
                      as.data.frame() %>% 
                      nameless::bin(span = 0.5, 
                                    fx = \(x) nameless::middle(x),
                                    .by = "id") %>% 
                      tidyr::nest()) %>% 
    dplyr::select(-data) %>% 
    tidyr::unnest(new_data) %>% 
    tidyr::unnest(data) %>% 
    dplyr::ungroup()

data.table::fwrite(data_original, 
                   file.path("results", "synthetic", "preprocessed_original.csv"))

# Create several different preprocessing pipelines to be tested. Created 
# with the following in mind: 
#   - Presence/Absence of a moving-window filtering approach (5 datapoints span)
#   - Counterbalancing moving-window filtering and binning
#   - Weighted average or normal average in moving-window
#   - Weighted averages by index or by actual time
#   - Taking the middle of the bin or the average as the position
binning <- \(x, fx) nameless::bin(x, span = 0.5, fx = fx, .by = "id")
moving <- \(x, fx) nameless::moving_window(x, span = 2, fx = fx, .by = "id")

bin_av <- \(x) binning(x, fx = \(x) nameless::average(x))
bin_mid <-  \(x) binning(x, fx = \(x) nameless::middle(x))

mov_av <- \(x) moving(x, fx = \(x) nameless::average(x))
mov_idx <- \(x) moving(x, fx = \(x) nameless::weighted_average(x, .by = "index"))
mov_time <- \(x) moving(x, fx = \(x) nameless::weighted_average(x, .by = "relative_time"))

conditions <- list(# Only binning
                   "bin_av" = list(bin_av), 
                   "bin_mid" = list(bin_mid),

                   # Moving window with span 2, average
                   "mov_av-bin_av" = list(mov_av, bin_av),
                   "bin_av-mov_av" = list(bin_av, mov_av),
                   "mov_av-bin_mid" = list(mov_av, bin_mid),
                   "bin_mid-mov_av" = list(bin_mid, mov_av),

                   # Moving window with span 2, weighted average on index
                   "mov_idx-bin_av" = list(mov_idx, bin_av),
                   "bin_av-mov_idx" = list(bin_av, mov_idx),
                   "mov_idx-bin_mid" = list(mov_idx, bin_mid),
                   "bin_mid-mov_idx" = list(bin_mid, mov_idx),

                   # Moving window with span 2, weighted average on time
                   "mov_time-bin_av" = list(mov_time, bin_av),
                   "bin_av-mov_time" = list(bin_av, mov_time),
                   "mov_time-bin_mid" = list(mov_time, bin_mid),
                   "bin_mid-mov_time" = list(bin_mid, mov_time))

# Combine the information of the conditions with the information on the data
# itself, matching conditions to data
data_files <- data.frame(filename = rep(data_files, each = length(conditions)), 
                         condition = rep(names(conditions), times = length(data_files)))

data.table::fwrite(data_files, file.path("results", "synthetic", "data_files.csv"))
saveRDS(conditions, file.path("results", "synthetic", "conditions.Rds"))





################################################################################
# ANALYSIS

# Load the needed variables
data_original <- list("preprocessed" = data.table::fread(file.path("results", "synthetic", "preprocessed_original.csv"), 
                                                         data.table = FALSE), 
                      "original" = data.table::fread(file.path("data", "synthetic", "synthetic_original.csv"), 
                                                     data.table = FALSE))
data_files <- data.table::fread(file.path("results", "synthetic", "data_files.csv"), 
                                data.table = FALSE)
conditions <- readRDS(file.path("results", "synthetic", "conditions.Rds"))
data_list <- readRDS(file.path("results", "synthetic", "data_list.Rds"))

binning <- \(x, fx) nameless::bin(x, span = 0.5, fx = fx, .by = "id")
moving <- \(x, fx) nameless::moving_window(x, span = 2, fx = fx, .by = "id")

bin_av <- \(x) binning(x, fx = \(x) nameless::average(x))
bin_mid <-  \(x) binning(x, fx = \(x) nameless::middle(x))

mov_av <- \(x) moving(x, fx = \(x) nameless::average(x))
mov_idx <- \(x) moving(x, fx = \(x) nameless::weighted_average(x, .by = "index"))
mov_time <- \(x) moving(x, fx = \(x) nameless::weighted_average(x, .by = "relative_time"))

# Create a function in which we will compute the summary statistics of interest
compute_summary_statistics <- function(data, kind) {
    data_original[[kind]] %>% 
        # Add the original dataset into the one it is being compared to and 
        # give `X` and `Y` as labels for the real positions (similar to 
        # stationary preprocessing)
        dplyr::rename(X = x, 
                      Y = y) %>% 
        dplyr::full_join(data, by = c("nsim", "time", "id")) %>% 

        # Compute the difference between filtered and expected positions. Used 
        # to measure the extent to which systematic error is present in the data
        dplyr::mutate(difference_x = X - x, 
                      difference_y = Y - y) %>% 

        # Compute the statistics of interest per simulation and id. This 
        # will allow for a more broad view on where it still might go awry
        dplyr::group_by(nsim, id) %>% 
        dplyr::arrange(time) %>% 
        dplyr::summarize(# Statistics about how close we are to the actual 
                         # positions
                         mean_diff_x = mean(difference_x, na.rm = TRUE), 
                         mean_diff_y = mean(difference_y, na.rm = TRUE), 
                         q025_diff_x = quantile(difference_x, probs = 0.025, na.rm = TRUE),
                         q025_diff_y = quantile(difference_y, probs = 0.025, na.rm = TRUE),
                         q975_diff_x = quantile(difference_x, probs = 0.975, na.rm = TRUE),
                         q975_diff_y = quantile(difference_y, probs = 0.975, na.rm = TRUE), 

                         # Statistics about the size of the measurement error
                         # (compared to the actual positions)
                         sd_diff_x = sd(difference_x, na.rm = TRUE), 
                         sd_diff_y = sd(difference_y, na.rm = TRUE), 
                         
                         # Autocorrelation in the residuals
                         auto_x = cor(difference_x[2:length(difference_x)], 
                                      difference_x[2:length(difference_x) - 1],
                                      use = "pairwise.complete.obs"), 
                         auto_y = cor(difference_y[2:length(difference_y)], 
                                      difference_y[2:length(difference_y) - 1], 
                                      use = "pairwise.complete.obs")) %>% 
        dplyr::ungroup() %>% 
        suppressMessages() %>% 
        return()
}

# For the analysis, we should first create a function that will do the 
# preprocessing and check its efficacy. This will allow us to put tidyverse to 
# its maximal use when preprocessing the data.
pipeline_efficacy <- function(x){

    print(paste0(x$filename, ": ", x$condition))

    # Retrieve the data and the pipeline for the condition
    local_data <- data_list[[x$filename]] %>% 
        dplyr::filter(time <= 10)
    fx <- conditions[[x$condition]]

    # Check whether the data have a reference to the simulation number. If not, 
    # add it to the dataframe
    if(is.null(local_data$nsim)) {
        local_data$nsim <- 1
    }

    # Compute the summary statistics of the data before they are processed 
    # through the pipeline. This will give us values to compare the results 
    # to, which is an overall better approach. Add an indicator that tells us 
    # that this is the original data
    summary_statistics <- compute_summary_statistics(local_data, "original") %>% 
        dplyr::mutate(type = "original")

    # Execute the pipeline for each data set at hand   
    local_data <- local_data %>% 
        dplyr::group_by(nsim) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(data = data %>%  
                          as.data.frame() %>% 
                          nameless::execute_pipeline(fx) %>% 
                          list()) %>% 
        tidyr::unnest(data) %>% 
        dplyr::ungroup()

    # Get some summary statistics from this
    summary_statistics <- compute_summary_statistics(local_data, "preprocessed") %>% 
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
results <- data_files[1,] %>%
    split(seq_len(nrow(data_files[1,]))) %>%
    as.list() %>% 
    # parallel::mclapply(pipeline_efficacy,
    #                    mc.cores = n_cores) %>%
    lapply(pipeline_efficacy) %>% 
    dplyr::bind_rows()

data.table::fwrite(results, 
                   file.path("results", "synthetic_preprocessing.csv"))





################################################################################
# VISUALIZATION

# Load the needed file
results <- data.table::fread(file.path("results", "synthetic_preprocessing.csv"))

# Create a function that will plot distributional histograms for the different 
# statistics and dimensions
make_plot <- function(x, 
                      dim, 
                      statistic) {

    # Create the column name to be used when selecting data
    col <- paste0(statistic, 
                  "_diff_", 
                  dim)

    # Extract the different data sets and the different preprocessing strategies
    dataset <- unique(x$filename)
    preprocess <- unique(x$condition)

    # Get the limits of the plots
    limits <- x %>% 
        dplyr::select(as.character(col)) %>% 
        range(na.rm = TRUE)
    limits <- limits + c(-1, 1) * diff(limits) * 0.05

    # Loop over these unique names and create the plots. Put these plots in a 
    # list to be bound together using ggpubr
    plt <- list()
    for(i in seq_along(preprocess)) {
        for(j in seq_along(dataset)) {
            # Select the data
            plot_data <- x %>% 
                dplyr::filter(filename == dataset[j], 
                              condition == preprocess[i]) %>% 
                dplyr::select(as.character(col)) %>% 
                setNames("X")

            # Create the plot
            idx <- (i - 1) * length(dataset) + j
            plt[[idx]] <- ggplot2::ggplot(data = plot_data, 
                                          ggplot2::aes(x = X)) +
                ggplot2::geom_histogram(fill = "cornflowerblue", 
                                        bins = 10) +
                ggplot2::geom_vline(xintercept = 0,
                                    color = "salmon", 
                                    linewidth = 2) +
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
                             ncol = length(dataset))

    return(plt)
}

# Define the dimensions and the statistics
statistics <- c("mean", "sd", "q025", "q975")
dims <- c("x", "y")

# Loop over them and save the plots
for(i in dims) {
    for(j in statistics) {
        ggplot2::ggsave(file.path("figures", "synthetic", paste0(j, "_", i, ".jpg")), 
                        make_plot(results, i, j), 
                        width = 6 * 500, 
                        height = 13 * 600, 
                        unit = "px")
    }
}
