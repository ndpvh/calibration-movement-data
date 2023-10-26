library(tidyverse)
library(data.table)

# Load in both the data about the experiments
exp <- fread(file.path("Calibration experiments", "Data", "experiments.csv"),
             data.table = FALSE) %>% 
    # Select only a few columns and rename `id` so that you can later join
    # this dataframe with the other one
    select(id, name) %>% 
    rename(experiment_id = id)

# Load and preprocess the data about the positions
pos <- fread(file.path("Calibration experiments", "Data", "datapoints.csv"),
             data.table = FALSE) %>%
    # Column names from NUC
    setNames(c("id", "created", "modified", "tag_id", "experiment_id",
               "timestamp", "x", "y", "person_id")) %>% 
    # Join the position data with the experimental data
    plyr::join(exp, by = "experiment_id") %>% 
    # Delete and rename columns
    select(timestamp, tag_id, x:y, name) %>% 
    rename(tag = tag_id,
           experiment = name) %>% 
    # Delete data in which no experiment was running
    filter(!(is.na(experiment))) 

# Experiments that have been done are the following:
#   - test experiment: walking around
#   - my_experiment:
#   - stationary 2, 3, 5, 6: laying down the tags on the net and letting them
#                                    sit there
#   - movement L X V: moving with X Volts across the Lth line of the Y-axis
#                     (movement itself is along X-axis)  
#   - movement Lperp X V: same principle, but X and Y are inversed
#   - movement diag X V: same principle, but this time along the diagonals

###################################################
# Stationary data
###################################################

# Just visualizing the data a bit
visualize_positions <- function(name){
    plot_data <- pos %>% 
        filter(experiment == name) %>% 
        select(x, y)

    plot(plot_data$x, plot_data$y)
}

visualize_positions("stationary 2")
visualize_positions("stationary 3")
visualize_positions("stationary 5")
visualize_positions("stationary 6")

# Lets make a sum of them all and plot these as well
to_sum <- c("stationary 2", 
            "stationary 3",
            "stationary 5",
            "stationary 6")
stationary <- pos %>% 
    group_by(experiment) %>% 
    tidyr::nest() %>% 
    mutate(stationary = experiment %in% to_sum) %>% 
    tidyr::unnest(data) %>% 
    ungroup()

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
