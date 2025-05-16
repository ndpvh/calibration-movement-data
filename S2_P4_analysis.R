################################################################################
# Purpose: Perform the formal analyses for study 2 and describe how we got to  #
#          our conclusions.                                                    #
################################################################################

devtools::load_all()


# Define the names of all the summary files, which contains the values of the 
# summary statistics for each simulation and condition. This will be used as the 
# basis for the main analyses.
files <- paste0(
    rep(
        c("movement", "fixed"),
        each = 6
    ),
    "_",
    rep(
        c("R10", "T10", "R6R", "T6R", "R6N", "T6N"),
        times = 2
    )
)

data <- lapply(
    files, 
    \(x) data.table::fread(
        file.path("results", "study 2", paste0("summary_", x, ".csv")), 
        data.table = FALSE
    )
) %>% 
    `names<-` (files)

# Define the columns in which we are interested
columns <- c("bias_dist", "rmse_dist", "mae_dist")

# Loop over the different datafiles and compare the values of the summary 
# statistics for each of the conditions to the values of these same statistics
# before any preprocessing was done. 
distribution <- lapply(
    data, 
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
        result <- lapply(
            columns, 
            function(y) {
               # Create reference
               reference <- before[, y]

               # Use mutate on the nested data
               after %>% 
                   dplyr::rowwise() %>% 
                   dplyr::mutate(
                      data_2 = data %>%
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
            }
        )

        # Bind together into one big dataframe with an additional 
        # column specifying the statistic that it used
        statistics <- sapply(
            seq_along(result), 
            \(i) rep(columns[i], nrow(result[[i]]))
        ) %>% 
            as.character() 
        result <- do.call("rbind", result)
        result$statistic <- statistics

        return(result)
    }
) %>% 
    `names<-` (names(data))

saveRDS(
    distribution, 
    file.path("results", "study 2", "difference_distribution.Rds")
)





#-------------------------------------------------------------------------------
# Selecting pipelines
#-------------------------------------------------------------------------------

# Let's check which conditions were significant for the "fixed" and "movement"
# datasets
significant <- data.frame(
    preprocessing_function = distribution[[1]]$preprocessing_function,
    statistics = distribution[[1]]$statistic
)
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
#   - Negative value of the statistic (i.e., measurement became better)
#   - Significance of the effect 
#   - Success in at least 2/3 of the datasets. This is an arbitrary criterion, 
#     but ensures some pipelines get picked even when all pipelines performed 
#     bad in some datasets (e.g., time-related error has no significant pipelines
#     that reduced the mean distance)
statistics <- c("bias_dist", "mae_dist", "rmse_dist")
results <- data.frame(
    preprocessing_function = significant %>% 
        dplyr::filter(statistics == "bias_dist") %>% 
        dplyr::select(preprocessing_function) %>% 
        unlist() %>% 
        as.vector()
)

for(i in statistics) {
    tmp <- significant %>% 
        dplyr::filter(statistics == i) %>% 
        dplyr::mutate(
            R10 = movement_R10 & movement_R10_effect < 0 & fixed_R10 & fixed_R10_effect < 0,
            R6R = movement_R6R & movement_R6R_effect < 0 & fixed_R6R & fixed_R6R_effect < 0,
            R6N = movement_R6N & movement_R6N_effect < 0 & fixed_R6N & fixed_R6N_effect < 0,
            T10 = movement_T10 & movement_T10_effect < 0 & fixed_T10 & fixed_T10_effect < 0,
            T6R = movement_T6R & movement_T6R_effect < 0 & fixed_T6R & fixed_T6R_effect < 0,
            T6N = movement_T6N & movement_T6N_effect < 0 & fixed_T6N & fixed_T6N_effect < 0
        ) %>%
        dplyr::rowwise() %>% 
        dplyr::mutate(
            selected = sum(dplyr::across(R10:T6N)) / 6,
            selected = selected >= 2/3
        ) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(preprocessing_function, selected)

    tmp <- setNames(tmp, c("preprocessing_function", i))

    results <- dplyr::full_join(
        results, 
        tmp, 
        by = "preprocessing_function"
    )    
}

View(results)
data.table::fwrite(
    results, 
    file.path(".", "results", "study 2", "selected_pipelines.csv")
)
