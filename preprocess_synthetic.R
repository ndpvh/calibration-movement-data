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
data_files <- c("synthetic_related_10", 
                "synthetic_unrelated_10", 
                "synthetic_related_6random", 
                "synthetic_related_6nonrandom",
                "synthetic_unrelated_6random", 
                "synthetic_unrelated_6nonrandom")
data_list <- lapply(data_files, 
                    \(x) data.table::fread(file.path("data", "synthetic", paste0(x, ".csv")), 
                                           data.table = FALSE))
names(data_list) <- data_files

saveRDS(data_list, file.path("results", "synthetic", "data_list.Rds"))

# Add some metadata to these files. This will allow us to quickly distinguish 
# the different files and compare how effective the preprocessing strategies 
# are on them.
data_files <- data.frame(filename = data_files, 
                         related_error = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE), 
                         sampling_rate = c(10, 10, 6, 6, 6, 6), 
                         random_drop = c(NA, NA, TRUE, FALSE, TRUE, FALSE))

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
                                    fx = \(x) nameless::average(x),
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
#   - Presence/Absence of a moving-window filtering approach
#   - Counterbalancing moving-window filtering and binning
#   - Weighted average or normal average in moving-window
#   - Size of the moving-window (5 or 11 data-points, span = 2 or 5)
#   - Weighted averages by index or by actual time
binning <- \(x) nameless::bin(x, span = 0.5, \(x) nameless::average(x), .by = "id")
moving_2 <- \(x, fx) nameless::moving_window(x, span = 2, fx = fx, .by = "id")
moving_5 <- \(x, fx) nameless::moving_window(x, span = 5, fx = fx, .by = "id")

conditions <- list("bin_only" = list(binning), 
                   # Moving window with span 2, average
                   "mov2_av_1" = list(\(x) moving_2(x, fx = \(x) nameless::average(x)), 
                                      binning), 
                   "mov2_av_2" = list(binning, 
                                      \(x) moving_2(x, fx = \(x) nameless::average(x))),
                   # Moving window with span 5, average
                   "mov5_av_1" = list(\(x) moving_5(x, fx = \(x) nameless::average(x)), 
                                      binning), 
                   "mov5_av_2" = list(binning, 
                                      \(x) moving_5(x, fx = \(x) nameless::average(x))),
                   # Moving window with span 2, weighted average on index
                   "mov2_wav_ind_1" = list(\(x) moving_2(x, fx = \(x) nameless::weighted_average(x, .by = "index")), 
                                           binning), 
                   "mov2_wav_ind_2" = list(binning, 
                                           \(x) moving_2(x, fx = \(x) nameless::weighted_average(x, .by = "index"))),
                   # Moving window with span 5, weighted average on index
                   "mov5_wav_ind_1" = list(\(x) moving_5(x, fx = \(x) nameless::weighted_average(x, .by = "index")), 
                                           binning), 
                   "mov5_wav_ind_2" = list(binning, 
                                           \(x) moving_5(x, fx = \(x) nameless::weighted_average(x, .by = "index"))),
                   # Moving window with span 2, weighted average on time
                   "mov2_wav_time_1" = list(\(x) moving_2(x, fx = \(x) nameless::weighted_average(x, .by = "relative_time")), 
                                            binning), 
                   "mov2_wav_time_2" = list(binning, 
                                            \(x) moving_2(x, fx = \(x) nameless::weighted_average(x, .by = "relative_time"))),
                   # Moving window with span 5, weighted average on time
                   "mov5_wav_time_1" = list(\(x) moving_5(x, fx = \(x) nameless::weighted_average(x, .by = "relative_time")), 
                                            binning), 
                   "mov5_wav_time_2" = list(binning, 
                                            \(x) moving_5(x, fx = \(x) nameless::weighted_average(x, .by = "relative_time"))))

# Add some metadata to each of these conditions to be used later
conditions_metadata <- data.frame(condition = names(conditions),
                                  moving_window = c(FALSE, rep(TRUE, 12)), 
                                  span = c(NA, rep(rep(c(2, 5), each = 2), times = 3)), 
                                  weighted = c(NA, rep(FALSE, 4), rep(TRUE, 8)), 
                                  by = c(rep(NA, 5), rep(c("index", "relative_time"), each = 4)))

# Combine the information of the conditions with the information on the data
# itself, matching conditions to data
data_files <- data_files %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(tidyr::nest(conditions_metadata)) %>% 
    tidyr::unnest(data) %>%
    dplyr::ungroup() %>% 
    dplyr::select(filename, 
                  condition, 
                  related_error:random_drop,
                  moving_window:by)

data.table::fwrite(data_files, file.path("results", "synthetic", "data_files.csv"))
saveRDS(conditions, file.path("results", "synthetic", "conditions.Rds"))





################################################################################
# ANALYSIS

# Load the needed variables
data_original <- data.table::fread(file.path("results", "synthetic", "preprocessed_original.csv"))
data_files <- data.table::fread(file.path("results", "synthetic", "data_files.csv"))
conditions <- readRDS(file.path("results", "synthetic", "conditions.Rds"))
data_list <- readRDS(file.path("results", "synthetic", "data_list.Rds"))

binning <- \(x) nameless::bin(x, span = 0.5, \(x) nameless::average(x), .by = "id")
moving_2 <- \(x, fx) nameless::moving_window(x, span = 2, fx = fx, .by = "id")
moving_5 <- \(x, fx) nameless::moving_window(x, span = 5, fx = fx, .by = "id")

# For the analysis, we should first create a function that will do the 
# preprocessing and check its efficacy. This will allow us to put tidyverse to 
# its maximal use when preprocessing the data.
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

    # Check the efficacy by comparing the values obtained after preprocessing to 
    # the positions of the original data set
    local_data <- local_data %>% 
        dplyr::rename(preprocessed_x = x, 
                      preprocessed_y = y) %>% 
        dplyr::full_join(data_original, by = c("nsim", "time", "id")) %>% 
        dplyr::mutate(difference_x = x - preprocessed_x, 
                      difference_y = y - preprocessed_y)

    # Get some summary statistics from this
    summary_statistics <- local_data %>% 
        dplyr::group_by(nsim) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(mean_diff_x = mean(data[[1]]$difference_x, na.rm = TRUE), 
                      mean_diff_y = mean(data[[1]]$difference_y, na.rm = TRUE), 
                      sd_diff_x = sd(data[[1]]$difference_x, na.rm = TRUE), 
                      sd_diff_y = sd(data[[1]]$difference_y, na.rm = TRUE), 
                      q025_diff_x = quantile(data[[1]]$difference_x, probs = 0.025, na.rm = TRUE),
                      q025_diff_y = quantile(data[[1]]$difference_y, probs = 0.025, na.rm = TRUE),
                      q975_diff_x = quantile(data[[1]]$difference_x, probs = 0.975, na.rm = TRUE),
                      q975_diff_y = quantile(data[[1]]$difference_y, probs = 0.975, na.rm = TRUE)) %>% 
        dplyr::select(-data)

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
    dplyr::bind_rows()

data.table::fwrite(results, 
                   file.path("results", "synthetic_preprocessing.csv"))
