################################################################################
# Purpose: Preprocess the calibration data that was gathered on 14/10/2023     #
#          and 21/10/2023. These data are used to estimate the amount of       #
#          systematic and unsystematic noise in the movement data. We          #
#          additionally use these data to find out what the actual sampling    #
#          rate of the data is (expected 5Hz)                                  #
################################################################################

source(file.path("utility", "utility.R"))

# Load in both the data about the experiments that were run on both days
exp <- fread(file.path("data", "experiments.csv"),
             data.table = FALSE) %>%
    # Select only a few columns and rename `id` so that you can later join
    # this dataframe with the one that contains the measured positions
    select(id, name, date) %>% 
    rename(experiment_id = id)

# Load and preprocess the data that contains the measured positions
data <- fread(file.path("data", "datapoints.csv"),
              data.table = FALSE) %>%
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
    filter(!(is.na(experiment))) 

# Experiments that have been done are the following:
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





#-------------------------------------------------------------------------------
# Stationary data
#-------------------------------------------------------------------------------

# Sum all of the stationary data to one datafile. This way, we get measurements
# at all positions of the net that we installed, which we can then compare to 
# the positions these measurements should actually have.
#
# The process below first creates a boolean that indicates whether the
# experiments are the ones that are defined in `stationary_experiments`. Then,
# this boolean is used to filter the data to only contain these data.
stationary_experiments <- paste("stationary", c(2, 3, 5, 6), sep = " ")
stationary <- data %>% 
    group_by(experiment) %>% 
    tidyr::nest() %>% 
    mutate(stationary_boolean = experiment %in% stationary_experiments) %>%
    tidyr::unnest(data) %>% 
    ungroup() %>% 
    filter(stationary_boolean)

# Visualize these data
visualize_positions(stationary)

# Lets make a sum of them all and plot these as well


stationary %>% 
    filter(stationary) %>% 
    ggplot(aes(x = x, y = y)) + 
        geom_point()

########################
# Preparing the data

# Create a function that will create a dotted rectangle of a given length and 
# width. Importantly, both should be natural numbers, as there is one meter 
# apart between each dot
create_rectangle <- function(w, h){
    w <- w + 1 ; h <- h + 1
    cbind(X = rep((1:w) - 1, each = h),
          Y = rep((1:h) - 1, times = w)) %>% 
        as_tibble() %>% 
        mutate(tag = row_number()) %>% 
        return()
}

# Create a function that will classify the measurements into several discrete
# groups, based on their continuous value.
# Assumptions are that each of the tags are enough apart so that the measurement 
# error is not severe enough to not measure their locations accurately
classify_dimension <- function(x, numcat){
    (x - min(x)) %>% 
        `/` (max(x) - min(x)) %>%
        `*` (numcat) %>% 
        round() %>% 
        return()
}

# Delete the original tags and replace them with the ones created by the grid
# Furthermore assign the "real" grid coordinates to each of the tags as to be 
# able to estimate the bias and variation of the measurement
stationary <- stationary %>% 
    # Only get data from the stationary measurements and delete tags
    filter(stationary) %>% 
    select(-tag) %>% 
    # Correct the measurements to begin at the origin for the best of your
    # abilities
    mutate(x = x,
           y = y) %>%
    # Assign each of the measurements to discrete categories as defined in 
    # classify_dimension
    mutate(X = classify_dimension(x, 10),
           Y = classify_dimension(y, 8))

# Quick test of whether this works: Apparently it does! (Change X to Y for the 
# other dimension)
ggplot(data = stationary, aes(x = x, y = y, color = factor(X))) +
    geom_point()

# Join the stationary dataframe with a rectangle that is created with the same 
# specifications
stationary <- stationary %>% 
    plyr::join(create_rectangle(10, 8), 
               by = c("X", "Y"))

# Finally, correct the "real positions" by adding the minimum of x and y to it.
# While this is not going to be a perfect fit, it does give us a close approximation
# to the actual positions of X and Y, which were not measured in our first 
# calibration experiment.
stationary <- stationary %>% 
    mutate(X = X + min(x),
           Y = Y + min(y))

# Check whether it worked with some plots. It actually seems to work, and the
# relationship is quite okay
ggplot(data = stationary, aes(x = x, y = X)) +
    geom_point() +
    geom_smooth()
ggplot(data = stationary, aes(x = y, y = Y)) +
    geom_point() +
    geom_smooth()

# Also visualize the relationship between x and y to make sure there are no
# interactions happening: Otherwise 2D estimation
ggplot(data = stationary, aes(x = x, y = y)) +
    geom_point() + 
    geom_smooth()

# Save these data as being preprocessed
saveRDS(stationary,
        file.path("Calibration experiments", 
                  "Data", 
                  "preprocessed_stationary.Rds"))
