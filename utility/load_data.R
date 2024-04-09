#' Load our experimental data
#' 
#' This is a basic function that will load in the experimental data that we have
#' in a format that we can analyze. More specifically, it will combine the data
#' from files 'experiments.csv' -- which contains information on the experiments 
#' that were run -- and 'datapoints.csv' -- which contains the positions that 
#' were measured during each experiment.
#' 
#' @export
load_data <- function(){
    # Load in both the data about the experiments that were run on both days.
    # The argument `fill = TRUE` is added to ensure that even if the datafile is
    # corrupt, we get a result.
    exp <- fread(file.path("data", "raw_data", "experiments.csv"),
                 data.table = FALSE,
                 fill = TRUE) %>%
        # Select only a few columns and rename `id` so that you can later join
        # this dataframe with the one that contains the measured positions
        select(id, name, date) %>% 
        rename(experiment_id = id)

    # Load and preprocess the data that contains the measured positions
    fread(file.path("data", "raw_data", "datapoints.csv"),
          data.table = FALSE,
          fill = TRUE) %>%
        # Column names from the NUC
        setNames(c("id", "created", "modified", "tag_id", "experiment_id",
                   "timestamp", "x", "y", "person_id")) %>% 
        # Join the position data with the experimental data
        plyr::join(exp, by = "experiment_id") %>% 
        # Delete and rename columns: Only retain timestamps, tag id's, x and y 
        # position, and the name of the experiment
        select(timestamp, tag_id, x:y, name) %>% 
        rename(tag = tag_id,
               experiment = name) %>% 
        # Delete data in which no experiment was running
        filter(!(is.na(experiment))) %>% 
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