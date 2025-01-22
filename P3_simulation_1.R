################################################################################
# Purpose: Perform the formal analyses for simulation 1 and describe how we    # 
#          got to our inference.                                               #
################################################################################

devtools::load_all()


# Define the names of all the summary files, which contains the values of the 
# summary statistics for each simulation and condition. This will be used as the 
# basis for the main analyses.
movement <- c("fixed", "movement")
errors <- c("_R10", "_R6N", "_R6R",
            "_U10", "_U6N", "_U6R",
            "_T10", "_T6N", "_T6R",
            "")
files <- paste(rep(movement, each = length(errors)), 
               rep(errors, times = length(movement)), 
               sep = "")

data <- lapply(files, 
               \(x) data.table::fread(file.path("results", 
                                                "simulation_1", 
                                                paste0("summary_", x, ".csv")), 
                                      data.table = FALSE)) %>% 
    `names<-` (files)

# Define the columns in which we are interested
columns <- c("bias_dist", "rmse_dist", "mae_dist")

# Loop over the different datafiles and compare the values of the summary 
# statistics for each of the conditions to the values of these same statistics
# before any preprocessing was done. 
distribution <- lapply(data, 
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
                      result <- lapply(columns, 
                                       function(y) {
                                          # Create reference
                                          reference <- before[, y]

                                          # Use mutate on the nested data
                                          after %>% 
                                              dplyr::rowwise() %>% 
                                              dplyr::mutate(data_2 = data %>%
                                                  as.data.frame() %>% 
                                                  dplyr::select(tidyselect::matches(y)) %>% 
                                                  unlist() %>% 
                                                  as.numeric() %>% 
                                                  nameless::compare_distribution(reference) %>% 
                                                  list()) %>% 
                                              dplyr::ungroup() %>% 
                                              dplyr::select(-data) %>% 
                                              tidyr::unnest(data_2) %>% 
                                              return()                                          
                                       })

                      # Bind together into one big dataframe with an additional 
                      # column specifying the statistic that it used
                      statistics <- sapply(seq_along(result), 
                                           \(i) rep(columns[i], nrow(result[[i]]))) %>% 
                          as.character() 
                      result <- do.call("rbind", result)
                      result$statistic <- statistics

                      return(result)
                  }) %>% 
    `names<-` (names(data))

saveRDS(distribution, 
        file.path("results", "simulation_1", "difference_distribution.Rds"))





#-------------------------------------------------------------------------------
# Selecting pipelines
#-------------------------------------------------------------------------------

# Let's check which conditions were significant for the "fixed" and "movement"
# datasets
significant <- data.frame(preprocessing_function = distribution[[1]]$preprocessing_function,
                          statistics = distribution[[1]]$statistic)
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
#   - Negative value of the statistic
#   - Significance of the effect 
#   - Success in at least half of the datasets. This is an arbitrary criterion, 
#     but ensures some pipelines get picked even when all pipelines performed 
#     bad in some datasets (e.g., time-related error has no significant pipelines
#     that reduced the mean distance)
statistics <- c("mae_dist", "rmse_dist")
results <- data.frame(preprocessing_function = significant %>% 
                          dplyr::filter(statistics == "bias_dist") %>% 
                          dplyr::select(preprocessing_function) %>% 
                          unlist() %>% 
                          as.vector())
for(i in statistics) {
    tmp <- significant %>% 
        dplyr::filter(statistics == i) %>% 
        dplyr::mutate(R10 = movement_R10 & movement_R10_effect < 0 & fixed_R10 & fixed_R10_effect < 0,
                      R6R = movement_R6R & movement_R6R_effect < 0 & fixed_R6R & fixed_R6R_effect < 0,
                      R6N = movement_R6N & movement_R6N_effect < 0 & fixed_R6N & fixed_R6N_effect < 0,
                      U10 = movement_U10 & movement_U10_effect < 0 & fixed_U10 & fixed_U10_effect < 0,
                      U6R = movement_U6R & movement_U6R_effect < 0 & fixed_U6R & fixed_U6R_effect < 0,
                      U6N = movement_U6N & movement_U6N_effect < 0 & fixed_U6N & fixed_U6N_effect < 0,
                      T10 = movement_T10 & movement_T10_effect < 0 & fixed_T10 & fixed_T10_effect < 0,
                      T6R = movement_T6R & movement_T6R_effect < 0 & fixed_T6R & fixed_T6R_effect < 0,
                      T6N = movement_T6N & movement_T6N_effect < 0 & fixed_T6N & fixed_T6N_effect < 0) %>%
        dplyr::rowwise() %>% 
        dplyr::mutate(selected = sum(dplyr::across(R10:T6N)) / 9,
                      selected = selected >= 2/3) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(preprocessing_function, selected)

    tmp <- setNames(tmp, c("preprocessing_function", i))

    results <- dplyr::full_join(results, tmp, by = "preprocessing_function")    
}

View(results)
data.table::fwrite(results, 
                   file.path(".", "results", "simulation_1", "selected_pipelines.csv"))
## In general, many pipelines did well according to our standards. The selected
## pipelines usually involve a moving window. When looking at both MAE and RMSE, 
## another 52 functions survive to the next round.
