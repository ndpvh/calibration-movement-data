#' Load the experimental data
#' 
#' This is a basic function that will load in the experimental data that we have
#' in a format that we can analyze. More specifically, it will combine the data
#' from files 'experiments.csv' -- which contains information on the experiments 
#' that were run -- and 'datapoints.csv' -- which contains the positions that 
#' were measured during each experiment.
#' 
#' @export
load_data <- function(){
    # Load the data that contains information on the experiments that were run. 
    # The argument `fill = TRUE` is added to ensure that even if the datafile is
    # corrupt, we get a result.
    exp <- data.table::fread(file.path("data", "raw_data", "experiments.csv"),
                             data.table = FALSE,
                             fill = TRUE) %>%
        # Select only a few columns and rename `id` so that you can later join
        # this dataframe with the one that contains the measured positions
        dplyr::select(id, name, date) %>% 
        dplyr::rename(experiment_id = id)

    # Load the data that contains the measured positions and bind them together
    # with the experiment information
    data.table::fread(file.path("data", "raw_data", "datapoints.csv"),
                      data.table = FALSE,
                      fill = TRUE) %>%
        # Column names from the NUC
        setNames(c("id", "created", "modified", "tag_id", "experiment_id",
                   "timestamp", "x", "y", "person_id")) %>% 
        # Join the position data with the experimental data
        plyr::join(exp, by = "experiment_id") %>% 
        # Delete and rename columns: Only retain timestamps, tag id's, x and y 
        # position, and the name of the experiment
        dplyr::select(timestamp, tag_id, x:y, name) %>% 
        dplyr::rename(tag = tag_id,
                      experiment = name) %>% 
        # Delete data that do not belong to any experiment
        dplyr::filter(!(is.na(experiment))) %>% 
        return()
}

# For future reference: The experiments that have been done are the following:
#   - test experiment:          walking around with no clear goal
#   - my_experiment:
#   - stationary 2, 3, 5, 6:    laying down the tags on the ground without moving 
#                               them any further
#   - movement L X V:           moving with X Volts across the Lth line of the 
#                               Y-axis (movement itself is along x-axis)  
#   - movement Lperp X V:       same principle, but X and Y are inversed
#   - movement diag X V:        same principle, but this time along the diagonals
#   - headbands:                stationary data with humans wearing headbands
#   - tablet batch experiment:  test whether sampling rate drops when using 
#                               many tablets at once