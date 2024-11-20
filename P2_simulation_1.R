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
# Parallellization
#-------------------------------------------------------------------------------

n_cores <- 3 #max(c(parallel::detectCores() - 1, 1))






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

# Define all moving windows. In these moving windows, we always use `span = 1`, 
# defining a window of 3 observations over which to execute the functions. This 
# number is based on a small tuning study in which performance for `span = 1` 
# outperformed those of `span = 2` and `span = 5` when looking at RMSE.
#
# As to the statistics used in this moving window, we use the simple average, 
# a weighted average over the index of the observation within the window 
# weighted by the standard normal, and a weighted average over the relative time 
# compared to the midpoint in the window weighted by a normal with standard 
# deviation 1/10. This standard deviation is chosen to be similar to the standard
# normal for index, where the time between observations is typically 0.1.
#
# To allow for the stable creation of different moving windows in for-loops, we 
# will need to create a wrapper-function that takes in the variable arguments 
# and outputs the function to be used in the pipeline. 
mw <- list(\(x) nameless::average(x), 
           \(x) nameless::weighted_average(x, 
                                           .by = "index"),
           \(x) nameless::weighted_average(x, 
                                           .by = "relative_time", 
                                           weights = \(x) dnorm(x, mean = 0, sd = 1/10))) %>% 
    lapply(function(x) {
               factory <- \(y) nameless::moving_window(y, 
                                                       span = 1, 
                                                       fx = x, 
                                                       .by = "id")
               return(factory)
           }) %>% 
    `names<-` (c("av", "idx", "time"))

# Define the Kalman filters. For these filters, we provide them with the 
# assumed error variance corresponding to a SD of 3.1cm, which we observed in 
# the calibration experiments. We differentiate between two filters; one that 
# goes in one direction of time (`reverse = FALSE`) and one in which we first 
# train the filter on the reversed data and only then filter the actual data 
# (`reverse = TRUE`).
kalm <- list("kalm" = \(x) nameless::kalman_filter(x, 
                                                   assumed_variance = 0.031^2,
                                                   reverse = FALSE, 
                                                   .by = "id"), 
             "kalm-rev" = \(x) nameless::kalman_filter(x, 
                                                       assumed_variance = 0.031^2, 
                                                       reverse = TRUE, 
                                                       .by = "id"))

# Define the LOESS and LOWESS. Here, we differentiate between the degree of the 
# fitted polynomial (linear or parabolic, `degree = 1` or `degree = 2` resp.)
# and the number of observations accounted for within the fit-window (either 10
# or 15). These settings were based on an initial tuning round.
reg <- list("loess-1-10" = \(x) nameless::local_regression(x, 
                                                           degree = 1, 
                                                           span_obs = 10, 
                                                           surface = "direct"), 
            "loess-1-15" = \(x) nameless::local_regression(x, 
                                                           degree = 1, 
                                                           span_obs = 15, 
                                                           surface = "direct"),
            "loess-2-10" = \(x) nameless::local_regression(x, 
                                                           degree = 2, 
                                                           span_obs = 10, 
                                                           surface = "direct"),
            "loess-2-15" = \(x) nameless::local_regression(x, 
                                                           degree = 2, 
                                                           span_obs = 15, 
                                                           surface = "direct"))

# Create functions that will create all possible pairs and triplets based on the 
# input strings. This will make it easier for us make the combination of the
# different preprocessing functions.
make_pairs <- function(x, y) {
    return(rbind(expand.grid(x, y), 
                 expand.grid(y, x)))
}

make_triplets <- function(x, y, z) {
    return(rbind(expand.grid(x, y, z), 
                 expand.grid(x, z, y), 
                 expand.grid(y, x, z), 
                 expand.grid(y, z, x), 
                 expand.grid(z, x, y), 
                 expand.grid(z, y, x)))
}

create_labels <- function(x) {
    columns <- colnames(x)

    do.call("paste", c(x[columns], sep = "_")) %>% 
        return()
}

# Create a list that contains all of the different combinations between the 
# preprocessing functions. In the end, we get a big list of lists which can be 
# handed to the `execute_pipeline` function.
all_functions <- append(append(mw, kalm), reg)

# Names of the relevant functions
singles <- names(all_functions)
pairs <- rbind(make_pairs(names(mw), names(kalm)), 
               make_pairs(names(mw), names(reg)), 
               make_pairs(names(reg), names(kalm)))
triplets <- make_triplets(names(mw), names(kalm), names(reg))

# Actually putting them in a list
singles <- lapply(singles, 
                  \(x) list(all_functions[[x]])) %>% 
    `names<-` (singles)
pairs <- lapply(seq_len(nrow(pairs)), 
                \(i) lapply(pairs[i,], 
                            \(x) all_functions[[x]])) %>% 
    `names<-` (create_labels(pairs))
triplets <- lapply(seq_len(nrow(triplets)), 
                   \(i) lapply(triplets[i,], 
                               \(x) all_functions[[x]])) %>% 
    `names<-` (create_labels(triplets))

# Combine all lists into one overarching combinations list
conditions <- append(append(singles, pairs), triplets)

# Combine the information of the conditions with the information on the data
# itself, matching conditions to data
data_files <- data.frame(filename = data_files, 
                         original = rep(c("movement", "fixed"), each = length(data_files) / 2))

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
                      difference_y = Y - y, 
                      distance = sqrt((X - x)^2 + (Y - y)^2)) %>% 

        # Compute the statistics of interest per simulation and id. This 
        # will allow for a more broad view on where it still might go awry
        dplyr::group_by(nsim, id) %>% 
        dplyr::arrange(time) %>% 
        dplyr::summarize(# Statistics about how close we are to the actual 
                         # positions
                         mean_diff_x = mean(difference_x, na.rm = TRUE), 
                         mean_diff_y = mean(difference_y, na.rm = TRUE), 
                         mean_dist = mean(distance, na.rm = TRUE),
                         q025_diff_x = quantile(difference_x, probs = 0.025, na.rm = TRUE),
                         q025_diff_y = quantile(difference_y, probs = 0.025, na.rm = TRUE),
                         q025_diff_y = quantile(distance, probs = 0.025, na.rm = TRUE),
                         q975_diff_x = quantile(difference_x, probs = 0.975, na.rm = TRUE),
                         q975_diff_y = quantile(difference_y, probs = 0.975, na.rm = TRUE), 
                         q975_diff_y = quantile(distance, probs = 0.975, na.rm = TRUE), 

                         # Statistics about the size of the measurement error
                         # (compared to the actual positions)
                         rmse_diff_x = sd(difference_x, na.rm = TRUE), 
                         rmse_diff_y = sd(difference_y, na.rm = TRUE), 
                         rmse_dist = sd(distance, na.rm = TRUE),
                         mae_diff_x = mean(abs(difference_x), na.rm = TRUE), 
                         mae_diff_y = mean(abs(difference_y), na.rm = TRUE), 
                         mae_dist = mean(abs(distance), na.rm = TRUE), 
                         
                         # Autocorrelation in the residuals
                         auto_x = cor(difference_x[2:length(difference_x)], 
                                      difference_x[2:length(difference_x) - 1],
                                      use = "pairwise.complete.obs"), 
                         auto_y = cor(difference_y[2:length(difference_y)], 
                                      difference_y[2:length(difference_y) - 1], 
                                      use = "pairwise.complete.obs"), 
                         auto_dist = cor(distance[2:length(distance)], 
                                         distance[2:length(distance) - 1], 
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
    
    # Retrieve the data and the pipeline for the condition
    local_data <- data_list[[x$filename]]
    fx <- names(conditions)

    # Check whether the data have a reference to the simulation number. If not, 
    # add it to the dataframe
    if(is.null(local_data$nsim)) {
        local_data$nsim <- 1
    }

    # Compute the summary statistics of the data before they are processed 
    # through the pipeline. This will give us values to compare the results 
    # to, which is an overall better approach. Add an indicator that tells us 
    # that this is the original data    
    result <- compute_summary_statistics(local_data, x$original) %>% 
        dplyr::mutate(preprocessed = "before", 
                      filename = x$filename, 
                      condition = NA) 

    data.table::fwrite(result, 
                       file.path("results", 
                                 "simulation_1", 
                                 "tmp", 
                                 paste0("tmp0.csv")))

    n <- length(fx)

    # Nest the different simulations in `local_data` so that we can use it 
    # in the mclapply later.
    local_data <- local_data %>% 
        dplyr::group_by(nsim) %>% 
        tidyr::nest()

    # Parallellize the execution of each of the pipelines, given the nested 
    # structure in `local_data`
    parallel::mclapply(seq_along(fx), 
                       function(i) {
                           # Print something so that we know where the function 
                           # is at
                           cat("\rExecuting pipeline", 
                               i, 
                               "of", 
                               n, 
                               "for dataset", 
                               x$filename)
            
                           # Execute the pipeline 
                           result <- lapply(seq_len(nrow(local_data)), 
                                            function(j) {
                                                local_data$data[[j]] %>% 
                                                    as.data.frame() %>% 
                                                    nameless::execute_pipeline(conditions[[fx[i]]], 
                                                                               report = FALSE) %>% 
                                                    dplyr::mutate(nsim = j) %>% 
                                                    return()
                                            })
                           result <- tryCatch(do.call("rbind", result), error = function(e) browser())
            
                           # Compute the summary statistics from the preprocessed 
                           # data and save these results
                           result <- compute_summary_statistics(result, x$original) %>% 
                               dplyr::mutate(preprocessed = "after", 
                                             filename = x$filename, 
                                             condition = fx[i]) 

                           data.table::fwrite(result, 
                                              file.path("results", 
                                                        "simulation_1", 
                                                        "tmp", 
                                                        paste0("tmp", i, ".csv")))
                           
                           rm(list = c("result"))
                           gc()

                           return(NULL)
                       },
                       mc.cores = n_cores)

    cat("\n")

    # Bind all results together
    summary_statistics <- list()
    for(i in 0:length(fx)) {
        summary_statistics[[i + 1]] <- data.table::fread(file.path("results", 
                                                                   "simulation_1", 
                                                                   "tmp", 
                                                                   paste0("tmp", i, ".csv")))
    }
    
    summary_statistics <- tryCatch(do.call("rbind", summary_statistics) %>% 
        dplyr::relocate(filename, condition, preprocessed, nsim),
        error = function(e) browser())

    # Save these results and delete the dataframes created here
    data.table::fwrite(summary_statistics, 
                       file.path("results", "simulation_1", paste0(x$filename, ".csv")))

    rm(list = c("local_data", "summary_statistics"))
    gc()

    return(NULL)
}

# Execute this function for each combination of the data and the pipeline.
for(i in seq_len(nrow(data_files))) {
    # Give us some feedback on which file is being preprocessed now
    cat("\rPreprocessing", 
        data_files$filename[i], 
        ": data file", 
        i, 
        "of", 
        nrow(data_files),
        "\n")

    # Actually preprocess the file
    pipeline_efficacy(data_files[i,])
}





################################################################################
# VISUALIZATION

# # Load the needed files
# filenames <- list.files(path = file.path(".", "results", "simulation_1"), 
#                         pattern = "\\.csv")
# filenames <- filenames[filenames != "data_files.csv"]

# results <- lapply(filenames, 
#                   \(x) data.table::fread(file.path(".", "results", "simulation_1", x), 
#                                          data.table = FALSE))

# # Create a function that takes in a dataframe and creates the plots of interest
# make_plot <- function(x, 
#                       statistics) {

#     # Split data before preprocessing and after preprocessing
#     before <- dplyr::filter(x, preprocessed == "before")
#     after <- dplyr::filter(x, preprocessed == "after")

#     # Get the data of before
#     before <- before %>%
#         dplyr::select(contains(statistics)) %>%
#         setNames("X") %>%
#         dplyr::mutate(M = 1)

#     # Get all conditions out of there
#     conditions <- unique(after$condition)

#     # Fix the limits on the x-axis (within bounds, of course)
#     all_x <- x[, statistics]

#     if(grepl("sd", statistics, fixed = TRUE)) {
#         idx <- all_x < quantile(all_x, probs = 0.95)
#     } else {
#         idx <- all_x < quantile(all_x, probs = 0.975) & all_x > quantile(all_x, probs = 0.025)
#     }

#     xlim <- range(all_x[idx])

#     # Loop over all conditions and create the plot of interest
#     plt <- list()
#     for(i in conditions) {
#         # Get plot data for the condition and the statistic of interest. Bind 
#         # together for before and after
#         plot_data <- after %>%
#             dplyr::filter(condition == i) %>%
#             dplyr::select(contains(statistics)) %>%
#             setNames("X") %>%
#             dplyr::mutate(M = 2) %>%
#             rbind(before) %>%
#             dplyr::mutate(M = factor(M))

#         # Create a histogram as the plot of choice. Include the condition name 
#         # in the plot and make the legend tell us something
#         plt[[i]] <- ggplot2::ggplot(data = plot_data, 
#                                     ggplot2::aes(x = X, fill = M)) +
#             ggplot2::geom_histogram(alpha = 0.5, 
#                                     bins = 15, 
#                                     color = "black", 
#                                     position = "identity") +
#             ggplot2::labs(title = i, 
#                           legend = "Preprocessed") +
#             ggplot2::lims(x = xlim) +
#             ggplot2::scale_fill_manual(labels = c("1" = "Before", 
#                                                   "2" = "After"), 
#                                        values = c("1" = "salmon", 
#                                                   "2" = "cornflowerblue")) +
#             ggplot2::theme_minimal() 
#     }

#     # Bind together and save under figures
#     plt <- ggpubr::ggarrange(plotlist = plt, 
#                              nrow = 17, 
#                              ncol = 17,
#                              common.legend = TRUE, 
#                              legend = "right")

#     ggplot2::ggsave(plt, 
#                     filename = file.path("figures", 
#                                          "simulation_1", 
#                                          "summary_statistics",
#                                          paste0(x$filename[1], "__", statistics, ".png")), 
#                     width = 17 * 600,
#                     height = 17 * 650, 
#                     unit = "px")

#     return(NULL)
# }

# # Create all figures
# for(i in seq_along(results)) {
#     for(j in c("mean_diff_x", "mean_diff_y", "mean_dist", "sd_diff_x", "sd_diff_y", "sd_dist")) {
#         make_plot(results[[i]], j)
#     }
# }
