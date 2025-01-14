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
library(locfit)

################################################################################
# PRELIMINARIES

#-------------------------------------------------------------------------------
# Parallellization
#-------------------------------------------------------------------------------

n_cores <- 30 #max(c(parallel::detectCores() - 1, 1))






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
mw <- list(\(x) nameless::average(x,
                                  cols = c("x_original", "y_original")), 
           \(x) nameless::weighted_average(x, 
                                           .by = "index",
                                           cols = c("x_original", "y_original")),
           \(x) nameless::weighted_average(x, 
                                           .by = "relative_time", 
                                           weights = \(x) dnorm(x, mean = 0, sd = 1/10),
                                           cols = c("x_original", "y_original"))) %>% 
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
reg <- list("loess-1" = \(x) nameless::local_regression(x, .by = "id", degree = 1), 
            "loess-2" = \(x) nameless::local_regression(x, .by = "id", degree = 2),
            "loess-3" = \(x) nameless::local_regression(x, .by = "id", degree = 3))

# Create functions that will create all possible pairs and triplets based on the 
# input strings. This will make it easier for us make the combination of the
# different preprocessing functions.
make_pairs <- function(x, y) {
    return(rbind(expand.grid(x, y, stringsAsFactors = FALSE), 
                 expand.grid(y, x, stringsAsFactors = FALSE)))
}

make_triplets <- function(x, y, z) {
    return(rbind(expand.grid(x, y, z, stringsAsFactors = FALSE), 
                 expand.grid(x, z, y, stringsAsFactors = FALSE), 
                 expand.grid(y, x, z, stringsAsFactors = FALSE), 
                 expand.grid(y, z, x, stringsAsFactors = FALSE), 
                 expand.grid(z, x, y, stringsAsFactors = FALSE), 
                 expand.grid(z, y, x, stringsAsFactors = FALSE)))
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
                \(i) lapply(as.character(pairs[i,]), 
                            \(x) all_functions[[x]])) %>% 
    `names<-` (create_labels(pairs))
triplets <- lapply(seq_len(nrow(triplets)), 
                   \(i) lapply(as.character(triplets[i,]), 
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
    nameless::pipeline_efficiency(data_list[[data_files$filename[i]]], 
                                  conditions, 
                                  .by = "nsim", 
                                  summary.by = "id",
                                  path = file.path(".", "results", "simulation_1"),
                                  filename = data_files$filename[i], 
                                  metadata = list("filename" = data_files$filename[i]), 
                                  n_cores = n_cores)
}

# Now that we have all results, we will also create overview files containing 
# all fixed and movement results together. Will make interpretation and analysis
# somewhat easier
filenames <- paste(rep(c("fixed", "movement"), each = 9),
                   rep(c("R10", "U10", "T10", "R6R", "U6R", "T6R", "R6N", "U6N", "T6N"), times = 2),
                   sep = "_") %>% 
    paste(".csv", sep = "")

# Merge datafiles together
for(i in c("fixed", "movement")) {
    # Select only those files that are either fixed or movement
    idx <- grepl(filenames, pattern = i, fixed = TRUE)
    selected_files <- filenames[idx]

    # Loop over trajectory or summary
    for(j in c("trajectory", "summary")) {    
        # Load these files and put them in a list
        files <- lapply(selected_files, 
                        \(x) data.table::fread(file.path(".", 
                                                         "results", 
                                                         "simulation_1", 
                                                         paste0(j, "_", x)), 
                                               data.table = FALSE) %>% 
                            dplyr::mutate(error_type = stringr::str_split_i(x, 
                                                                            pattern = "_", 
                                                                            i = 2) %>% 
                                              stringr::str_split_i(pattern = ".csv", 
                                                                   i = 1)))
    
        # Bind these data together
        files <- do.call("rbind", files) %>% 
            dplyr::mutate(movement_type = i)
    
        data.table::fwrite(files, 
                           file.path(".", "results", "simulation_1", paste0(j, "_", i, ".csv")))
    }
}





################################################################################
# VISUALIZATION

#-------------------------------------------------------------------------------
# Per file
#-------------------------------------------------------------------------------

# Load the needed files
filenames <- paste(rep(c("fixed", "movement"), each = 9),
                   rep(c("R10", "U10", "T10", "R6R", "U6R", "T6R", "R6N", "U6N", "T6N"), times = 2),
                   sep = "_") %>% 
    paste(".csv", sep = "")

trajectory_results <- lapply(filenames, 
                             \(x) data.table::fread(file.path(".", "results", "simulation_1", paste0("trajectory_", x)), 
                                                    data.table = FALSE))
summary_results <- lapply(filenames, 
                          \(x) data.table::fread(file.path(".", "results", "simulation_1", paste0("summary_", x)), 
                                                 data.table = FALSE))

# Create all figures
columns <- c("mean_diff_x", 
             "mean_diff_y", 
             "mean_dist", 
             "rmse_diff_x", 
             "rmse_diff_y", 
             "rmse_dist", 
             "mae_diff_x", 
             "mae_diff_y", 
             "mae_dist")

for(i in seq_along(results)) {
    for(j in columns) {
        # Bar plot
        plt <- nameless::barplot(summary_results[[i]], j)
        ggplot2::ggsave(plt, 
                        filename = file.path("figures", 
                                             "simulation_1", 
                                             "barplot summary statistics",
                                             paste0(stringr::str_split_i(filenames[i], 
                                                                         pattern = ".csv", 
                                                                         i = 1), 
                                                    "__", 
                                                    j, 
                                                    ".png")), 
                        width = 15 * 600,
                        height = 15 * 650, 
                        unit = "px")

        # Histograms
        plt <- nameless::histogram(summary_results[[i]], j)

        ggplot2::ggsave(plt, 
                        filename = file.path("figures", 
                                             "simulation_1", 
                                             "histogram summary statistics",
                                             paste0(stringr::str_split_i(filenames[i], 
                                                                         pattern = ".csv", 
                                                                         i = 1), 
                                                    "__", 
                                                    j, 
                                                    ".png")), 
                        width = 15 * 600,
                        height = 15 * 650, 
                        unit = "px")
    }

    # Get all unique preprocessing pipelines and plot the trajectories for 
    # these
    for(j in unique(trajectory_results[[i]]$condition)[-1]) {
        # Trajectories
        plt <- trajectory_results[[i]] %>% 
            dplyr::filter(nsim == 1) %>% 
            dplyr::filter(condition %in% c("", j)) %>% 
            nameless::trajectory()
    
        ggplot2::ggsave(plt, 
                        filename = file.path("figures", 
                                             "simulation_1", 
                                             "trajectory summary statistics",
                                             paste0(stringr::str_split_i(filenames[i], 
                                                                         pattern = ".csv", 
                                                                         i = 1), 
                                                    "__", 
                                                    j, 
                                                    ".png")), 
                        width = 900 * 3, 
                        height = 1000 * 10,
                        unit = "px")
    }
}






#-------------------------------------------------------------------------------
# All files together
#-------------------------------------------------------------------------------

# Load the needed files
filenames <- c("fixed.csv", "movement.csv")
trajectory_results <- lapply(filenames, 
                             \(x) data.table::fread(file.path(".", "results", "simulation_1", paste0("trajectory_", x)), 
                                                    data.table = FALSE))
summary_results <- lapply(filenames, 
                          \(x) data.table::fread(file.path(".", "results", "simulation_1", paste0("summary_", x)), 
                                                 data.table = FALSE))

# Create all figures
columns <- c("mean_diff_x", 
             "mean_diff_y", 
             "mean_dist", 
             "rmse_diff_x", 
             "rmse_diff_y", 
             "rmse_dist", 
             "mae_diff_x", 
             "mae_diff_y", 
             "mae_dist")

for(i in seq_along(filenames)) {
    for(j in columns) {
        # Bar plot
        plt <- nameless::barplot(summary_results[[i]], j)
        ggplot2::ggsave(plt[["plot"]], 
                        filename = file.path("figures", 
                                             "simulation_1", 
                                             "barplot summary statistics",
                                             paste0(stringr::str_split_i(filenames[i], 
                                                                         pattern = ".csv", 
                                                                         i = 1), 
                                                    "__", 
                                                    j, 
                                                    ".png")), 
                        width = 15 * 600,
                        height = 15 * 650, 
                        unit = "px")

        # Histograms
        plt <- nameles::histogram(summary_results[[i]], j)

        ggplot2::ggsave(plt, 
                        filename = file.path("figures", 
                                             "simulation_1", 
                                             "histogram summary statistics",
                                             paste0(stringr::str_split_i(filenames[i], 
                                                                         pattern = ".csv", 
                                                                         i = 1), 
                                                    "__", 
                                                    j, 
                                                    ".png")), 
                        width = 15 * 600,
                        height = 15 * 650, 
                        unit = "px")
    }

    # Get all unique preprocessing pipelines and plot the trajectories for 
    # these
    for(j in unique(trajectory_results[[i]]$condition)[-1]) {
        # Trajectories
        plt <- trajectory_results[[i]] %>% 
            dplyr::filter(nsim == 1) %>% 
            dplyr::filter(condition %in% c("", j)) %>% 
            nameless::trajectory()
    
        ggplot2::ggsave(plt, 
                        filename = file.path("figures", 
                                             "simulation_1", 
                                             "trajectory summary statistics",
                                             paste0(stringr::str_split_i(filenames[i], 
                                                                         pattern = ".csv", 
                                                                         i = 1), 
                                                    "__", 
                                                    j, 
                                                    ".png")), 
                        width = 900 * 3, 
                        height = 1000 * 10,
                        unit = "px")
    }
}    
