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

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------

# Get the data that you want to preprocess.
data_files <- c("movement_R10",
                "movement_U10",
                "movement_T10",
                "movement_R6R", 
                "movement_U6R",
                "movement_T6R",
                "movement_R6N",
                "movement_U6N",
                "movement_T6N", 
                "fixed_R10",
                "fixed_U10",
                "fixed_T10",
                "fixed_R6R",
                "fixed_U6R",
                "fixed_T6R",
                "fixed_R6N",
                "fixed_U6N",
                "fixed_T6N")
data_list <- lapply(data_files, 
                    \(x) data.table::fread(file.path("data", "simulation_1", paste0(x, ".csv")), 
                                           data.table = FALSE))
names(data_list) <- data_files

saveRDS(data_list, file.path("results", "simulation_1", "data_list.Rds"))





#-------------------------------------------------------------------------------
# Pipelines
#-------------------------------------------------------------------------------

# Create several different preprocessing pipelines to be tested. Created 
# with the following in mind: 
#   - Differences in moving window approach
#       - Span: 3, 5, or 11 data points
#       - Statistic: average, weighted average, local linear or quadratic regression
#       - Variable for weighting: time or index
#       - Distribution for weighting: standard normal or N(0, 0.25^2) (about span of 1.5m)
#   - Kalman filter
#       - Whether to feedback and feedforward the data 
#
# We cross all of these approaches

# Define all moving windows. To allow for the stable creation of different 
# moving windows in for-loops, we will need to create a wrapper-function that 
# takes in the variable arguments and outputs the function to be used in the 
# pipeline. 
fx <- list("av" = \(x) nameless::average(x), 
           "idx" = \(x) nameless::weighted_average(x, .by = "index"),
           "time" = \(x) nameless::weighted_average(x, .by = "relative_time"),
           "lin" = \(x) nameless::linear(x),
           "quad" = \(x) nameless::parabola(x))

# Define the Kalman filters
kalm_rev <- \(x) nameless::kalman_filter(x, reverse = TRUE, .by = "id")
kalm_norev <- \(x) nameless::kalman_filter(x, reverse = FALSE, .by = "id")

# Bind together all conditions exhaustively. With this binding, we keep account 
# of the following:
#   - Counterbalance order of the pipeline
#   - If multiple moving windows combined, no similar ones used (i.e., average 
#     not paired with weighted average, linear not paired with parabola)
#   - If multiple moving windows combined, span of these moving windows is kept
#     the same
#   - If Kalman filter used, either put at the end or at the beginning of the 
#     pipeline
#
# First, create a list for all moving windows in separation
spans <- c(1, 2, 5)
fx_names <- c("av", "idx", "time", "lin", "quad")
combos <- data.frame(spans = rep(spans, each = length(fx_names)), 
                     fx_names = rep(fx_names, times = length(spans)))

mw_alone <- lapply(seq_len(nrow(combos)), 
                   function(i) {
                       # Both steps in this function prevent R from being lazy 
                       # and only compiling the function when we need it. 
                       # Instead, we will ask R to compile the function now, 
                       # preventing problems when the index `i` takes on values 
                       # that fall outside of the range of `combos`

                       # Initialize the function to be performed on the window
                       gx <- fx[[combos$fx[i]]]

                       # Initialize the moving window function itself
                       factory <- \(x) nameless::moving_window(x, span = combos$spans[i], fx = gx, .by = "id")
                       return(list(factory))
                   })
names(mw_alone) <- sapply(seq_len(nrow(combos)), 
                          \(i) paste(combos$spans[i], combos$fx[i], sep = "_"))

# Create a list with the moving windows combined with each other. Same logic
# applies here with regard to lazy R
spans <- c(1, 2, 5)
fx_1 <- c("av", "idx", "time")
fx_2 <- c("lin", "quad")
combos <- data.frame(spans = rep(spans, each = length(fx_1) * length(fx_2)), 
                     fx_1 = rep(rep(fx_1, each = length(fx_2)), times = length(spans)), 
                     fx_2 = rep(fx_2, times = length(spans) * length(fx_1)))

mw_combined <- append(lapply(seq_len(nrow(combos)), 
                             function(i) {
                                 gx_1 <- fx[[combos$fx_1[i]]]
                                 gx_2 <- fx[[combos$fx_2[i]]]

                                 factory_1 <- \(x) nameless::moving_window(x, span = combos$spans[i], fx = gx_1, .by = "id")
                                 factory_2 <- \(x) nameless::moving_window(x, span = combos$spans[i], fx = gx_2, .by = "id")

                                 return(list(factory_1, factory_2))
                             }), 
                      lapply(seq_len(nrow(combos)), 
                             function(i) {
                                 gx_1 <- fx[[combos$fx_1[i]]]
                                 gx_2 <- fx[[combos$fx_2[i]]]

                                 factory_1 <- \(x) nameless::moving_window(x, span = combos$spans[i], fx = gx_1, .by = "id")
                                 factory_2 <- \(x) nameless::moving_window(x, span = combos$spans[i], fx = gx_2, .by = "id")

                                 return(list(factory_2, factory_1))
                             }))
names(mw_combined) <- c(sapply(seq_len(nrow(combos)), 
                               \(i) paste(combos$spans[i], combos$fx_1[i], combos$fx_2[i], sep = "_")), 
                        sapply(seq_len(nrow(combos)), 
                               \(i) paste(combos$spans[i], combos$fx_2[i], combos$fx_1[i], sep = "_")))

# Create a list with the Kalman filters alone
kalman_alone <- list("kr" = list(kalm_rev), 
                     "kn" = list(kalm_norev))

# Add Kalman filter combinations to the moving windows
mw <- append(mw_alone, mw_combined)
mw_names <- names(mw)

kalman_combined <- list()
for(i in seq_along(mw)) {
    # Reversed Kalman filter first and last
    key <- paste0("kr_", mw_names[i])
    kalman_combined[[key]] <- append(list(kalm_rev), mw[[i]])

    key <- paste0(mw_names[i], "_kr")
    kalman_combined[[key]] <- append(mw[[i]], list(kalm_rev))

    # Original Kalman filter first and last
    key <- paste0("kn_", mw_names[i])
    kalman_combined[[key]] <- append(list(kalm_norev), mw[[i]])
    
    key <- paste0(mw_names[i], "_kn")
    kalman_combined[[key]] <- append(mw[[i]], list(kalm_norev))
}

# Combine all lists into one overarching combinations list
conditions <- append(mw_alone,
                     append(kalman_alone, 
                            append(mw_combined, 
                                   kalman_combined)))

# Combine the information of the conditions with the information on the data
# itself, matching conditions to data
data_files <- data.frame(filename = rep(data_files, each = length(conditions)), 
                         condition = rep(names(conditions), times = length(data_files)),
                         original = rep(c("movement", "fixed"), each = length(conditions) * 9))

data.table::fwrite(data_files, file.path("results", "simulation_1", "data_files.csv"))
saveRDS(conditions, file.path("results", "simulation_1", "conditions.Rds"))





################################################################################
# PREPROCESSING

# Load the original data set (for comparison)
original <- list("movement" = data.table::fread(file.path("data", "simulation_1", "movement.csv"), 
                                                data.table = FALSE), 
                 "fixed" = data.table::fread(file.path("data", "simulation_1", "fixed.csv"), 
                                             data.table = FALSE))

# Create a function in which we will compute the summary statistics of interest.
# Key here is to compare the positions that result from preprocessing with the 
# real positions that are contained in the `original` data set.
#
# This function takes in the arguments `data` -- the preprocessed data set -- 
# and `kind` -- the key that defines the original data. It then computes the 
# difference between preprocessed and actual positions and provides some 
# summary statistics.
compute_summary_statistics <- function(data, kind) {
    original[[kind]] %>% 
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

# Create a function that will do the preprocessing and checks its efficacy. 
# It depends on `compute_summary_statistics` and will format the data in a way 
# that is useful to us. Putting this all in a separate function will allow us to 
# put tidyverse to its maximal use when preprocessing the data.
#
# The single argument `x` contains the information on the combination of data set
# and pipeline that are contained in the variable `data_files`.
pipeline_efficacy <- function(x){

    print(paste0(x$filename, ": ", x$condition))

    # Retrieve the data and the pipeline for the condition
    local_data <- data_list[[x$filename]] 
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
    summary_statistics <- compute_summary_statistics(local_data, x$original) %>% 
        dplyr::mutate(preprocessed = "before")

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
    summary_statistics <- compute_summary_statistics(local_data, x$original) %>% 
        dplyr::mutate(preprocessed = "after") %>% 
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
    # lapply(pipeline_efficacy) %>% 
    dplyr::bind_rows()

data.table::fwrite(results, 
                   file.path("results", "simulation_1" "preprocessing_results.csv"))





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
