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

n_cores <- 3 #max(c(parallel::detectCores() - 1, 1))






#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------

# Get the data that you want to preprocess.
data_files <- paste(
    "data", 
    c("R10", "T10", "R6R", "T6R", "R6N", "T6N"),
    sep = "_"
)
data_list <- lapply(
    data_files, 
    \(x) data.table::fread(
        file.path("data", "study 3", paste0(x, ".csv")), 
        data.table = FALSE
    )
)
names(data_list) <- data_files

saveRDS(data_list, file.path("results", "study 3", "data_list.Rds"))





#-------------------------------------------------------------------------------
# Pipelines
#-------------------------------------------------------------------------------

# Define the pipelines that worked best overall. Based on the results that come 
# from study 2, which were saved in a separate file. Only use those (combinations
# of) filters that performed well on all three summary statistics, which are 
# 19 in total
filters <- readRDS(file.path("results", "study 2", "filters.Rds"))

selection <- data.table::fread(
    file.path("results", "study 2", "selected_filters.csv"),
    data.table = FALSE
)
selection$total <- selection$bias_dist & selection$rmse_dist

# Filter out those combinations of functions that did not perform well, mapping
# the names in the data.frame to the names of the list
selected <- selection$preprocessing_function[selection$total]
filters <- filters[selected]

# Adjust the spans of the LOESS analysis
spans <- seq(
    1 / 1800, 
    1799 / 1800, 
    5 / 1800
)
reg <- list(
    "loess-1" = \(x) nameless::local_regression(x, .by = "id", degree = 1, spans = spans), 
    "loess-2" = \(x) nameless::local_regression(x, .by = "id", degree = 2, spans = spans),
    "loess-3" = \(x) nameless::local_regression(x, .by = "id", degree = 3, spans = spans)
)

for(i in names(filters)) {
    if(stringr::str_detect(i, "loess")) {
        split <- stringr::str_split(i, "_")[[1]]
        idx <- which(split %in% names(reg))
        idy <- which(names(reg) %in% split)

        filters[[i]][[idx]] <- reg[[idy]]
    }
}

saveRDS(
    filters, 
    file.path("results", "study 3", "filters.Rds")
)





################################################################################
# PREPROCESSING

# Load the original data set. Used to compare the filtered data to
original <- data.table::fread(
    file.path("data", "study 3", "data.csv"), 
    data.table = FALSE
)

# Execute this function for each combination of the data and the pipeline.
for(i in seq_along(data_files)) {
    # Give us some feedback on which file is being preprocessed now
    cat(
        "\rPreprocessing", 
        data_files[i], 
        ": data file", 
        i, 
        "of", 
        length(data_files),
        "\n"
    )

    # Actually preprocess the file
    nameless::pipeline_efficiency(
        data_list[[data_files[i]]], 
        filters, 
        .by = "nsim", 
        summary.by = "id",
        path = file.path("results", "study 3"),
        filename = data_files[i], 
        metadata = list("filename" = data_files[i]), 
        n_cores = n_cores
    )
}

# COMMENTED OUT DUE TO MEMORY INTENSITY
#
# # Now that we have all results, we will also create overview files containing 
# # all results together. Will make interpretation and analysis somewhat easier
# filenames <- paste(data_files, ".csv", sep = "")

# # Merge datafiles together. Loop over trajectory or summary
# for(i in c("trajectory", "summary")) {    
#     # Load these files and put them in a list
#     files <- lapply(
#         filenames, 
#         \(x) data.table::fread(
#             file.path("results", "study 3", paste0(i, "_", x)), 
#             data.table = FALSE
#         ) %>% 
#             dplyr::mutate(
#                 error_type = stringr::str_split_i(x, pattern = "_", i = 2) %>% 
#                     stringr::str_split_i(pattern = ".csv", i = 1)
#             )
#     )

#     # Bind these data together and save in a conjoint file
#     files <- do.call("rbind", files) 

#     data.table::fwrite(
#         files, 
#         file.path("results", "study 3", paste0(i, ".csv"))
#     )
# }
