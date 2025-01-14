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

n_cores <- 1 #max(c(parallel::detectCores() - 1, 1))






#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------

# Get the data that you want to preprocess.
data_files <- paste("data", 
                    c("R10", "U10", "T10",
                      "R6R", "U6R", "T6R",
                      "R6N", "U6N", "T6N"),
                    sep = "_")
data_list <- lapply(data_files, 
                    \(x) data.table::fread(file.path("data", "simulation_2", paste0(x, ".csv")), 
                                           data.table = FALSE))
names(data_list) <- data_files

saveRDS(data_list, file.path("results", "simulation_2", "data_list.Rds"))





#-------------------------------------------------------------------------------
# Pipelines
#-------------------------------------------------------------------------------

# Define the pipelines that worked best overall. Based on the results that are 
# described in P3_simulation_1.R, in which we found that the following pipelines
# worked in the majority of the cases:
#   - Kalman filter
#   - Kalman filter (reversed)
#   - LOESS (2nd degree)
#   - LOESS (2nd degree) + Kalman filter
#   - LOESS (2nd degree) + Kalman filter (reversed)
#   - Kalman filter + LOESS (2nd degree)
#   - Kalman filter (reversed) + LOESS (2nd degree)
#   - Kalman filter + LOESS (3rd degree)
#   - Kalman filter (reversed) + LOESS (3rd degree)
#
# These are the pipelines we will use on simulation 2.
kalm <- \(x) nameless::kalman_filter(x, 
                                     assumed_variance = 0.031^2,
                                     reverse = FALSE, 
                                     .by = "id")
kalm_rev <- \(x) nameless::kalman_filter(x, 
                                         assumed_variance = 0.031^2,
                                         reverse = TRUE, 
                                         .by = "id")
loess_2 <- \(x) nameless::local_regression(x, .by = "id", degree = 2)
loess_3 <- \(x) nameless::local_regression(x, .by = "id", degree = 3)

conditions <- list("kalm" = list(kalm), 
                   "kalm-rev" = list(kalm_rev),
                   "loess-2" = list(loess_2),
                   "loess-2_kalm" = list(loess_2, kalm),
                   "loess-2_kalm-rev" = list(loess_2, kalm_rev),
                   "kalm_loess-2" = list(kalm, loess_2),
                   "kalm-rev_loess-2" = list(kalm_rev, loess_2),
                   "kalm_loess-3" = list(kalm, loess_3),
                   "kalm-rev_loess-3" = list(kalm_rev, loess_3))

saveRDS(conditions, file.path("results", "simulation_2", "conditions.Rds"))





################################################################################
# PREPROCESSING

# Load the original data set (for comparison)
original <- data.table::fread(file.path("data", "simulation_2", "data.csv"), 
                              data.table = FALSE)

# Execute the pipeline for each of the datafiles.
for(i in seq_along(data_files)) {
    # Give us some feedback on which file is being preprocessed now
    cat("\rPreprocessing", 
        data_files[i], 
        ": data file", 
        i, 
        "of", 
        length(data_files),
        "\n")

    # Actually preprocess the file
    nameless::pipeline_efficiency(data_list[[data_files[i]]] %>% 
                                      dplyr::filter(nsim %in% 1:5), 
                                  conditions, 
                                  .by = "nsim",
                                  summary.by = "id",
                                  path = file.path(".", "results", "simulation_2"),
                                  filename = data_files[i], 
                                  metadata = list("filename" = data_files[i]), 
                                  n_cores = n_cores)
}

# Now that we have all results, we will also an overview file containing 
# all results together. 
filenames <- paste(data_files, ".csv", sep = "")

# Loop over trajectory and summary
for(i in c("trajectory", "summary")) {
    files <- lapply(filenames, 
                    \(x) data.table::fread(file.path(".", 
                                                     "results", 
                                                     "simulation_2", 
                                                     paste0(i, "_", x)), 
                                           data.table = FALSE) %>% 
                        dplyr::mutate(error_type = stringr::str_split_i(x, 
                                                                        pattern = "_", 
                                                                        i = 2) %>% 
                                          stringr::str_split_i(pattern = ".csv", 
                                                               i = 1)))
    
        # Bind these data together
        files <- do.call("rbind", files)
        data.table::fwrite(files, 
                           file.path(".", "results", "simulation_2", paste0("data_", i, ".csv")))
}





################################################################################
# VISUALIZATION


