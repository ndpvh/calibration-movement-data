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

    # Save the original datafiles and give it a tag of "before"
    local_data %>% 
        dplyr::mutate(preprocessed = "before", 
                      filename = x$filename, 
                      condition = NA) %>% 
        data.table::fwrite(file.path("results", 
                                     "simulation_1", 
                                     "tmp_trajectory", 
                                     paste0("tmp0.csv")))

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
                                 "tmp_summary", 
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

                           # Save this preprocessed trajectory in a temporary 
                           # file
                           result %>% 
                               dplyr::mutate(preprocessed = "after", 
                                             filename = x$filename, 
                                             condition = fx[i]) %>% 
                               data.table::fwrite(file.path("results", 
                                                            "simulation_1", 
                                                            "tmp_trajectory", 
                                                            paste0("tmp", i, ".csv")))
            
                           # Compute the summary statistics from the preprocessed 
                           # data and save these results
                           result <- compute_summary_statistics(result, x$original) %>% 
                               dplyr::mutate(preprocessed = "after", 
                                             filename = x$filename, 
                                             condition = fx[i]) 

                           data.table::fwrite(result, 
                                              file.path("results", 
                                                        "simulation_1", 
                                                        "tmp_summary", 
                                                        paste0("tmp", i, ".csv")))
                           
                           rm(list = c("result"))
                           gc()

                           return(NULL)
                       },
                       mc.cores = n_cores)

    cat("\n")

    # Bind all results together
    summary_statistics <- list()
    trajectories <- list()
    for(i in 0:length(fx)) {
        summary_statistics[[i + 1]] <- data.table::fread(file.path("results", 
                                                                   "simulation_1", 
                                                                   "tmp_summary", 
                                                                   paste0("tmp", i, ".csv")))
        trajectories[[i + 1]] <- data.table::fread(file.path("results", 
                                                             "simulation_1", 
                                                             "tmp_trajectory", 
                                                             paste0("tmp", i, ".csv")))
    }
    
    summary_statistics <- tryCatch(do.call("rbind", summary_statistics) %>% 
        dplyr::relocate(filename, condition, preprocessed, nsim),
        error = function(e) browser())

    # Save these results and delete the dataframes created here
    data.table::fwrite(summary_statistics, 
                       file.path("results", "simulation_1", paste0("summary_", x$filename, ".csv")))
    data.table::fwrite(trajectories, 
                       file.path("results", "simulation_1", paste0("trajectory_", x$filename, ".csv")))

    rm(list = c("local_data", "summary_statistics", "trajectories"))
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





################################################################################
# VISUALIZATION

#-------------------------------------------------------------------------------
# Different plots
#-------------------------------------------------------------------------------

# Create a function that takes in a dataframe and creates the plots of interest
histogram <- function(x, 
                      statistics) {

    # Split data before preprocessing and after preprocessing
    before <- dplyr::filter(x, preprocessed == "before")
    after <- dplyr::filter(x, preprocessed == "after")

    # Get the data of before
    before <- before %>%
        dplyr::select(contains(statistics)) %>%
        setNames("X") %>%
        dplyr::mutate(M = 1)

    # Get all conditions out of there
    conditions <- unique(after$condition)

    # Fix the limits on the x-axis (within bounds, of course)
    all_x <- x[, statistics]

    if(is.na(sd(all_x))) {
        xlim <- c(0, 1)
    } else {
        if(grepl("rmse", statistics, fixed = TRUE) | 
           grepl("dist", statistics, fixed = TRUE) |
           grepl("mae", statistics, fixed = TRUE)) {

            limit <- max(c(quantile(before$X, probs = 0.95), 
                           quantile(after[, statistics], probs = 0.95)))

            limit <- max(c(mean(before$X) + 3 * sd(before$X), 
                           mean(after[, statistics]) + 3 * sd(after[, statistics])))

            idx <- all_x < limit

        } else {
            limits <- c(min(quantile(before$X, probs = 0.025), 
                            quantile(after[, statistics], probs = 0.025)), 
                        max(quantile(before$X, probs = 0.975), 
                            quantile(after[, statistics], probs = 0.975)))

            limits <- c(min(c(mean(before$X) - 3 * sd(before$X), 
                              mean(after[, statistics]) - 3 * sd(after[, statistics]))), 
                        max(c(mean(before$X) + 3 * sd(before$X), 
                              mean(after[, statistics]) + 3 * sd(after[, statistics]))))

            idx <- all_x < limits[2] & all_x > limits[1]
        }

        xlim <- range(all_x[idx]) + 0.05 * c(-1, 1) * diff(range(all_x[idx]))

        if(grepl("mean", statistics, fixed = TRUE) & !grepl("dist", statistics, fixed = TRUE)) {
            xlim <- c(-max(abs(xlim)), max(abs(xlim)))
        }

        
    }

    

    # Loop over all conditions and create the plot of interest
    plt <- list()
    for(i in conditions) {
        # Get plot data for the condition and the statistic of interest. Bind 
        # together for before and after
        plot_data <- after %>%
            dplyr::filter(condition == i) %>%
            dplyr::select(contains(statistics)) %>%
            setNames("X") %>%
            dplyr::mutate(M = 2) %>%
            rbind(before) %>%
            dplyr::mutate(M = factor(M))

        # Create a histogram as the plot of choice. Include the condition name 
        # in the plot and make the legend tell us something
        plt[[i]] <- ggplot2::ggplot(data = plot_data, 
                                    ggplot2::aes(x = X, fill = M)) +
            ggplot2::geom_histogram(alpha = 0.5, 
                                    bins = 15, 
                                    color = "black", 
                                    position = "identity") +
            ggplot2::labs(title = i, 
                          legend = "Preprocessed") +
            ggplot2::lims(x = xlim) +
            ggplot2::scale_fill_manual(labels = c("1" = "Before", 
                                                  "2" = "After"), 
                                       values = c("1" = "salmon", 
                                                  "2" = "cornflowerblue")) +
            ggplot2::theme_minimal() 
    }

    # Bind together and save under figures
    plt <- ggpubr::ggarrange(plotlist = plt, 
                             nrow = 15, 
                             ncol = 15,
                             common.legend = TRUE, 
                             legend = "right")

    return(plt)
}

# Create a function to create a bar plot for each of the conditions
barplot <- function(x, 
                    statistics) {

    # Split data before preprocessing and after preprocessing
    before <- dplyr::filter(x, preprocessed == "before")
    after <- dplyr::filter(x, preprocessed == "after")

    # Get the data of before
    before <- before %>%
        dplyr::select(contains(statistics), condition) %>%
        dplyr::mutate(condition = "before") %>% 
        setNames(c("X", "M"))

    # Fix the limits on the x-axis (within bounds, of course)
    all_x <- x[, statistics]

    if(is.na(sd(all_x))) {
        xlim <- c(0, 1)
    } else {
        if(grepl("rmse", statistics, fixed = TRUE) | 
           grepl("dist", statistics, fixed = TRUE) |
           grepl("mae", statistics, fixed = TRUE)) {

            limit <- max(c(quantile(before$X, probs = 0.95), 
                           quantile(after[, statistics], probs = 0.95)))

            limit <- max(c(mean(before$X) + 3 * sd(before$X), 
                           mean(after[, statistics]) + 3 * sd(after[, statistics])))

            idx <- all_x < limit

        } else {
            limits <- c(min(quantile(before$X, probs = 0.025), 
                            quantile(after[, statistics], probs = 0.025)), 
                        max(quantile(before$X, probs = 0.975), 
                            quantile(after[, statistics], probs = 0.975)))

            limits <- c(min(c(mean(before$X) - 3 * sd(before$X), 
                              mean(after[, statistics]) - 3 * sd(after[, statistics]))), 
                        max(c(mean(before$X) + 3 * sd(before$X), 
                              mean(after[, statistics]) + 3 * sd(after[, statistics]))))

            idx <- all_x < limits[2] & all_x > limits[1]
        }

        xlim <- range(all_x[idx]) + 0.05 * c(-1, 1) * diff(range(all_x[idx]))

        if(grepl("mean", statistics, fixed = TRUE) & !grepl("dist", statistics, fixed = TRUE)) {
            xlim <- c(-max(abs(xlim)), max(abs(xlim)))
        }        
    }

    # Create some plot data that will be used for the barplot
    conditions <- c("before", unique(after$condition))
    plot_data <- after %>% 
        dplyr::select(contains(statistics), condition) %>% 
        setNames(c("X", "M")) %>% 
        rbind(before) %>% 
        dplyr::group_by(M) %>% 
        dplyr::summarize(means = mean(X), 
                         q975 = quantile(X, probs = 0.975),
                         sd = sd(X)) %>% 
        dplyr::ungroup() %>% 
        dplyr::arrange(factor(M, levels = conditions)) %>% 
        dplyr::rename(X = M) %>% 
        dplyr::mutate(M = ifelse(X == "before", 1, 2))

    # Create a barplot using all of this information. The barplot will show 
    # the mean levels of each condition, hopefully providing us with a clearer
    # picture than the histograms
    plt <- ggplot2::ggplot(data = plot_data) +
        # ggplot2::geom_errorbar(ggplot2::aes(x = factor(X), 
        #                                     ymin = 0, 
        #                                     ymax = means + sign(means) * sd)) +
        ggplot2::geom_bar(ggplot2::aes(x = factor(X), 
                                       y = means, 
                                       fill = factor(M)),
                          stat = "identity",
                          color = "black") +
        ggplot2::coord_flip() +
        ggplot2::labs(title = paste("Average performance:", statistics), 
                      legend = "Preprocessed") +
        ggplot2::scale_fill_manual(labels = c("1" = "Before", 
                                              "2" = "After"), 
                                   values = c("1" = "salmon", 
                                              "2" = "cornflowerblue")) +
        ggplot2::theme_minimal() 

    return(plt)
}

# Create a function that will create the wanted plot
trajectory <- function(x) {
    # Retrieve data and bind it together with the preprocessed data. We add the
    # columns x_original and y_original to the dataframe to make sure we can 
    # delete them if present in the preprocessed data before joining with the 
    # original data. Differentially handled by Kalman filters than moving 
    # windows, as the latter needs explicit inclusion of columns while the 
    # former does this by default.
    local_data <- preprocess(x) %>% 
        dplyr::rename(filtered_x = x, 
                      filtered_y = y) %>% 
        dplyr::mutate(x_original = NA, 
                      y_original = NA) %>% 
        dplyr::select(-x_original, -y_original) %>% 
        dplyr::full_join(data_list[[x$filename]], 
                         by = c("nsim", "id", "time"))

    # Loop over each of the id's for a separate plot
    ids <- unique(local_data$id)

    # Create name-plots that denote whatever it is you're seeing
    name_plot <- function(x) {
        return(ggplot2::ggplot() +
            ggplot2::annotate("text", 
                              x = 0, 
                              y = 0,
                              label = x,
                              size = 10,
                              hjust = 0.5, 
                              vjust = 0.5) +
            ggplot2::theme_void())
    }

    plt <- list()
    plt[[1]] <- name_plot(" ")
    plt[[2]] <- name_plot("Unfiltered")
    plt[[3]] <- name_plot("Filtered")

    f <- length(plt) + 1
    for(i in seq_along(ids)) {
        # Get the original data and make them in plot data (x, y, xend, yend)
        original <- to_segments(local_data %>% 
                                    dplyr::filter(nsim == 1), 
                                .vars = c("x_original", "y_original"), 
                                .id = ids[i])            

        # Get filtered and unfiltered data
        other <- list(to_segments(local_data, 
                                  .vars = c("x", "y"), 
                                  .id = ids[i]), 
                      to_segments(local_data, 
                                  .vars = c("filtered_x", "filtered_y"), 
                                  .id = ids[i]))

        # Compute the standard deviations between the actual movement and the 
        # measured movement in both cases. Will be  
        rmse <- c(compute_rmse(local_data, 
                               .vars = c("x", "y"), 
                               .id = ids[i]), 
                  compute_rmse(local_data, 
                               .vars = c("filtered_x", "filtered_y"), 
                               .id = ids[i]))

        # Compute the limits of the plot. Makes sure both plots have the same 
        # limits
        xlims <- c(other[[1]]$x, 
                   other[[2]]$x, 
                   other[[1]]$xend, 
                   other[[2]]$xend) %>% 
            range() 
        ylims <- c(other[[1]]$y, 
                   other[[2]]$y, 
                   other[[1]]$yend, 
                   other[[2]]$yend) %>% 
            range()

        # Make the limits somewhat broader
        xlims <- xlims + c(-1, 1) * diff(xlims) * 0.25
        ylims <- ylims + c(-1, 1) * diff(ylims) * 0.25

        # Create a name plot for the kind of movement
        # Add a name-plot
        plt[[f]] <- name_plot(ids[i])
        f <- f + 1

        # And make the plots for filtered and unfiltered data
        for(j in seq_along(other)) {
            plt[[f]] <- ggplot2::ggplot() +
                # Measured vs real movements
                ggplot2::geom_segment(data = other[[j]], 
                                      ggplot2::aes(x = x, 
                                                   y = y, 
                                                   xend = xend, 
                                                   yend = yend), 
                                      color = "grey75", 
                                      linewidth = 1, 
                                      alpha = 0.1) +
                ggplot2::geom_segment(data = original,
                                      ggplot2::aes(x = x, 
                                                   y = y, 
                                                   xend = xend, 
                                                   yend = yend), 
                                      color = "black", 
                                      linewidth = 1) +
                # Distance between measured and real movements
                ggplot2::annotate("text", 
                                  x = xlims[1] + 0.95 * diff(xlims), 
                                  y = ylims[1] + 0.95 * diff(ylims), 
                                  label = latex2exp::TeX(paste0("$RMSE = ", 
                                                                rmse[j],
                                                                "$")), 
                                  size = 5,
                                  hjust = 1, 
                                  vjust = 1) +
                # Theme, limits, and labels
                ggplot2::labs(x = "x", 
                              y = "y") +
                ggplot2::lims(x = xlims, 
                              y = ylims) +
                ggplot2::theme_minimal() +
                ggplot2::theme(plot.title = ggplot2::element_text(size = 35, hjust = 0.5), 
                               axis.title = ggplot2::element_text(size = 25))
            f <- f + 1
        }
    }

    return(plt)
}






#-------------------------------------------------------------------------------
# Per file
#-------------------------------------------------------------------------------

# Load the needed files
filenames <- paste(rep(c("fixed", "movement"), each = 9),
                   rep(c("R10", "U10", "T10", "R6R", "U6R", "T6R", "R6N", "U6N", "T6N"), times = 2),
                   sep = "_") %>% 
    paste(".csv", sep = "")

results <- lapply(filenames, 
                  \(x) data.table::fread(file.path(".", "results", "simulation_1", x), 
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
        # 
        plt <- histogram(results[[i]], j)

        ggplot2::ggsave(plt, 
                        filename = file.path("figures", 
                                             "simulation_1", 
                                             "histogram summary statistics",
                                             paste0(results[[i]]$filename[1], "__", j, ".png")), 
                        width = 15 * 600,
                        height = 15 * 650, 
                        unit = "px")
    }
}






#-------------------------------------------------------------------------------
# All files together
#-------------------------------------------------------------------------------

# Load the needed files
filenames <- c("fixed.csv", "movement.csv")
results <- lapply(filenames, 
                  \(x) data.table::fread(file.path(".", "results", "simulation_1", x), 
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
        plt <- barplot(results[[i]], j)
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
        plt <- histogram(results[[i]], j)

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
}
