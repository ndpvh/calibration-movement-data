################################################################################
# Purpose: Perform the formal analyses for simulation 1 and describe how we    # 
#          got to our inference.                                               #
################################################################################

devtools::load_all()


# Define the names of all the summary files, which contains the values of the 
# summary statistics for each simulation and condition. This will be used as the 
# basis for the main analyses.
movement <- c("fixed", "movement")
errors <- c("R10", "R6N", "R6R",
            "U10", "U6N", "U6R",
            "T10", "T6N", "T6R")
files <- paste(rep(movement, each = length(errors)), 
               rep(errors, times = length(movement)), 
               sep = "_")

data <- lapply(files, 
               \(x) data.table::fread(file.path("results", 
                                                "simulation_1", 
                                                paste0("summary_", x, ".csv")), 
                                      data.table = FALSE)) %>% 
    `names<-` (files)

# Also make a general "movement" and "fixed" file that will combine all of the 
# specific files. These will be used to determine which of the preprocessing 
# pipelines can be used in the next step.
data[["movement"]] <- do.call("rbind",
                              data[grepl("movement", files, fixed = TRUE)])
data[["fixed"]] <- do.call("rbind",
                            data[grepl("fixed", files, fixed = TRUE)])

# Define the columns in which we are interested
columns <- c("mean_dist", "rmse_dist", "mae_dist")

# Loop over the different datafiles and compare the values of the summary 
# statistics for each of the conditions to the values of these same statistics
# before any preprocessing was done. 
results <- lapply(data, 
                  function(x) {
                      # Separate the data from before and after the preprocessing
                      before <- dplyr::filter(x, preprocessed == "before")
                      after <- dplyr::filter(x, preprocessed == "after")

                      # Nest the data based on the unique conditions in the 
                      # preprocessed data
                      after <- after %>%
                          dplyr::group_by(condition) %>% 
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

saveRDS(results, 
        file.path("results", "simulation_1", "difference_distribution.Rds"))





#-------------------------------------------------------------------------------
# Selecting pipelines
#-------------------------------------------------------------------------------

# Let's check which conditions were significant for the "fixed" and "movement"
# datasets
significant <- data.frame(condition = results[[1]]$condition,
                          statistics = results[[1]]$statistic)
for(i in names(data)) {
    significant[, i] <- results[[i]]$significant
    significant[, paste0(i, "_effect")] <- results[[i]]$median_diff
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
statistics <- c("mean_dist", "rmse_dist", "mae_dist")
results <- data.frame(condition = significant$condition[significant$statistics == "mean_dist"])
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
                      selected = selected >= 0.5) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(condition, selected)

    tmp <- setNames(tmp, c("condition", i))

    results <- dplyr::full_join(results, tmp, by = "condition")    
}

View(results)
## In general: RMSE did not deliver any pipelines. Mean distance and MAE 
## delivered the same pipelines, namely: LOESS 2 and 3, Kalman filter, and the 
## combination of these two. These will be used in the next step.

## Additional results from inspection: Percentage significant and percentage 
## negative per statistic
##  - Mean distances:
##     - Fixed:
##
##                  significance        negative
##          - R10:  158/158             158/158
##          - R6R:  157/158             158/158
##          - R6N:  158/158             158/158
##          - U10:  158/158             158/158
##          - U6R:  157/158             158/158
##          - U6N:  158/158             158/158
##          - T10:  158/158             158/158
##          - T6R:  157/158             158/158
##          - T6N:  158/158             158/158
##            
##      - Movement
##
##                  significance        negative
##          - R10:  69/158              158/158
##          - R6R:  58/158              59/158
##          - R6N:  15/158              64/158
##          - U10:  69/158              158/158
##          - U6R:  57/158              63/158
##          - U6N:  11/158              64/158
##          - T10:  0/158               86/158
##          - T6R:  134/158             8/158
##          - T6N:  130/158             4/158

##  - RMSE
##     - Fixed
##
##                  significance        negative
##          - R10:  158/158             158/158
##          - R6R:  158/158             158/158
##          - R6N:  154/158             158/158
##          - U10:  158/158             158/158
##          - U6R:  158/158             158/158
##          - U6N:  155/158             158/158
##          - T10:  117/158             158/158
##          - T6R:  91/158              158/158
##          - T6N:  132/158             158/158
##            
##      - Movement:
##
##                  significance        negative
##          - R10:  14/158              158/158
##          - R6R:  67/158              59/158
##          - R6N:  73/158              42/158
##          - U10:  2/158               158/158
##          - U6R:  65/158              59/158
##          - U6N:  73/158              38/158
##          - T10:  18/158              12/158
##          - T6R:  127/158             4/158
##          - T6N:  134/158             7/158
##
##  - MAE:
##      - Fixed:
##
##                  significance        negative
##          - R10:  158/158             158/158
##          - R6R:  158/158             158/158
##          - R6N:  157/158             158/158
##          - U10:  158/158             158/158
##          - U6R:  158/158             158/158
##          - U6N:  157/158             158/158
##          - T10:  158/158             158/158
##          - T6R:  158/158             158/158
##          - T6N:  157/158             158/158
##            
##      - Movement:
##
##                  significance        negative
##          - R10:  69/158              158/158
##          - R6R:  15/158              64/158
##          - R6N:  58/158              59/158
##          - U10:  69/158              158/158
##          - U6R:  11/158              64/158
##          - U6N:  57/158              63/158
##          - T10:  0/158               86/158
##          - T6R:  130/158             4/158
##          - T6N:  134/158             8/158