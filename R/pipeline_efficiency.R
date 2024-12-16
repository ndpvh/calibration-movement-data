#' Compare preprocessed with actual positions
#' 
#' This function preprocesses the provided data and then checks how close it gets
#' to the actual data. Importantly, this process assumes that you have the actual 
#' positions under columns "x_actual" and "y_actual" while the measured 
#' positions are under "x" and "y". Preprocessed measures will also end up under
#' "x" and "y", but only after preprocessing.
#' 
#' Depends on the function \code{\link[nameless]{summary_statistics}} which 
#' uses measured/filtered and actual positions to compute mean deviation, 
#' absolute mean deviation, and root mean squared error. 
#' 
#' @param data Dataframe to be preprocessed
#' @param fx List of preprocessing pipelines to be executed. Each pipeline 
#' consists of another list with functiosn to use when preprocessing the data.
#' If names are added to this list, the names will be included in a separate 
#' column \code{function} in the resulting dataframes. Allows one to know which
#' data were created by which function. Otherwise includes number of the pipeline.
#' @param .by Character vector containing columns to group the preprocessing by. 
#' This only concerns columns at the highest level and is not passed on to the 
#' lower-level preprocessing functions. Defaults to \code{NULL}.
#' @param summary List of functions that will be used to summarize the results of
#' preprocessing. Defaults to computation of the RMSE, MAE, and mean difference.
#' @param .vars Character vector or matrix containing the columns for which to 
#' compute the summary statistics. If matrix, each row should contain the columns
#' to compare. Defaults to x, y, and distance.
#' @param path Path under which to save the datafiles that result from this 
#' function. Defaults to a "results" folder at the current directory.
#' @param filename Name to use as an identifier for the results of this analysis.
#' Defaults to an empty string.
#' @param metadata Named list containing metadata that should be added to the 
#' preprocessed datafiles. Consists of column name and value to assign to it.
#' Defaults to an empty list.
#' @param n_cores Integer denoting how many cores to use for the analysis. 
#' Defaults to \code{1}.
#' 
#' 
#' @return Returns \code{NULL} and instead saves several datafiles under the path 
#' specified by the user
#' 
#' @export
pipeline_efficiency <- function(data, 
                                fx, 
                                .by = NULL,
                                summary = list(bias, rmse, mae),
                                .vars = c("X", "Y", "Dist", "x", "y", "dist"),
                                path = file.path(".", "results"), 
                                filename = "",
                                metadata = list(),
                                n_cores = 1) {

    ############################################################################
    # Step 1: Preliminary information

    # Check whether the data have a reference to .by argument. If not, then we 
    # create one in the data, allowing us to just continue with the 
    # preprocessing itself.
    if(!is.null(.by) & is.null(data[, .by])) {
        data[, .by] <- 1
    }

    # Check whether the folders "tmp_trajectory" and "tmp_summary" exist. If not, 
    # create them.
    if(!dir.exists(file.path(path, "tmp_trajectory"))) {
        dir.create(file.path(path, "tmp_trajectory"))
    }

    if(!dir.exists(file.path(path, "tmp_summary"))) {
        dir.create(file.path(path, "tmp_summary"))
    }

    # Save the original datafiles and give it a tag of "before". Can be used for
    # visualization purposes later on.
    data %>% 
        dplyr::mutate(preprocessed = "before",
                      preprocessing_function = NA) %>% 
        add_metadata(metadata = metadata) %>% 
        data.table::fwrite(file.path(path, "tmp_trajectory", "tmp0.csv"))

    # Compute the summary statistics of the data before they are processed 
    # through the pipeline. This will give us values to compare the results 
    # to, which is an overall better approach. Add an indicator that tells us 
    # that this is the original data    
    result <- data %>%
        dplyr::mutate(X = x_original, 
                      Y = y_original, 
                      dist = sqrt(x^2 + y^2),
                      Dist = sqrt(X^2 + Y^2)) %>% 
        summary_statistics(fx = summary, 
                           .vars = .vars, 
                           .by = .by) %>% 
        dplyr::mutate(preprocessed = "before",
                      preprocessing_function = NA) %>% 
        add_metadata(metadata = metadata)

    data.table::fwrite(result, 
                       file.path(path, "tmp_summary", "tmp0.csv"))



    ############################################################################
    # Step 2: Preprocessing

    # Nest the different simulations in `local_data` so that we can use it 
    # in the mclapply later.
    data <- data %>% 
        dplyr::group_by(tidyselect::all_of(.by)) %>% 
        tidyr::nest()

    # Get the names of the functions
    function_names <- if(!is.null(names(fx))) names(fx) else seq_len(length(fx))

    # Create function that will do the preprocessing on the local level, 
    # depending only on an index of the pipeline to use in `fx`
    process <- function(i) {
        # Print something so that we know where the function is at
        cat("\rExecuting pipeline", i, "of", length(fx))

        # Execute the pipeline for each of the nested data structures in `data``
        result <- lapply(seq_len(nrow(local_data)), 
                         function(j) {
                             data$data[[j]] %>% 
                                 as.data.frame() %>% 
                                 nameless::execute_pipeline(fx[[i]], 
                                                            report = FALSE) %>% 
                                 return()
                         })

        result <- do.call("rbind", result)

        # Save this preprocessed trajectory in a temporary file
        result %>% 
            dplyr::mutate(preprocessed = "after", 
                          preprocessing_function = funtion_names[i]) %>% 
            add_metadata(metadata) %>% 
            data.table::fwrite(file.path(path,
                                         "tmp_trajectory", 
                                         paste0("tmp", i, ".csv")))

        # Compute the summary statistics from the preprocessed 
        # data and save these results in a temporary file
        browser()
        result <- result %>% 
            dplyr::mutate(X = x_original, 
                          Y = y_original, 
                          dist = sqrt(x^2 + y^2),
                          Dist = sqrt(X^2 + Y^2)) %>% 
            summary_statistics(fx = summary,
                               .vars = .vars,
                               .by = .by) %>% 
            dplyr::mutate(preprocessed = "after", 
                          preprocessing_function = funtion_names[i]) %>% 
            add_metadata(metadata)

        data.table::fwrite(result, 
                           file.path(path, "tmp_summary", paste0("tmp", i, ".csv")))

        # Remove some of these variables and do garbage collection. Helps in 
        # memory maintenance        
        rm(list = c("result"))
        gc()

        return(NULL)
    }

    # Parallellize the execution of each of the pipelines, given the nested 
    # structure in `local_data`
    parallel::mclapply(seq_along(fx), 
                       process,
                       mc.cores = n_cores)
    cat("\n")



    ############################################################################
    # Step 3: Handling and saving results

    # Bind all results together
    summary_statistics <- list()
    trajectories <- list()
    for(i in 0:length(fx)) {
        summary_statistics[[i + 1]] <- data.table::fread(file.path(path, 
                                                                   "tmp_summary", 
                                                                   paste0("tmp", i, ".csv")))
        trajectories[[i + 1]] <- data.table::fread(file.path(path, 
                                                             "tmp_trajectory", 
                                                             paste0("tmp", i, ".csv")))
    }
    
    summary_statistics <- do.call("rbind", summary_statistics) 
    trajectories <- do.call("rbind", trajectories)

    # Save these results and delete the dataframes created here
    data.table::fwrite(summary_statistics, 
                       file.path(path, paste0("summary_", filename, ".csv")))
    data.table::fwrite(trajectories, 
                       file.path(path, paste0("trajectory_", filename, ".csv")))

    rm(list = c("local_data", "summary_statistics", "trajectories"))
    gc()

    # Also delete all files in the temporary paths. We don't need them anymore.
    files <- list.files(file.path(path, "tmp_summary"), 
                        include.dirs = F, 
                        full.names = T, 
                        recursive = T)
    file.remove(files)

    files <- list.files(file.path(path, "tmp_trajectory"), 
                        include.dirs = F, 
                        full.names = T, 
                        recursive = T)
    file.remove(files)

    # Nothing to return here
    return(NULL)
}
