################################################################################
# Purpose: Preprocess the calibration data that was gathered on 14/10/2023     #
#          and 21/10/2023. These data are used to estimate the amount of       #
#          systematic and unsystematic noise in the movement data. We          #
#          additionally use these data to find out what the actual sampling    #
#          rate of the data is (expected 5Hz)                                  #
################################################################################

source(file.path("utility", "utility.R"))

# Load the data of interest
data <- load_data()

#-------------------------------------------------------------------------------
# Stationary data: 14/10/2023
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
    filter(experiment %in% stationary_experiments)

# Visualize these data
visualize_positions(stationary)

# Create a function that will assign a measured row number that corresponds to 
# the row numbers given in the `rectangle` function. It simply works by standardizing 
# the measured positions, multiplying it by the number of rows in that specific
# direction (x or y), and then rounding of the value.
assign_row <- function(x, n_rows){
    (x - min(x)) %>% 
        `/` (max(x) - min(x)) %>%
        `*` (n_rows) %>% 
        round() %>% 
        return()
}

# Use the `assign_row` function to determine in which row each measured position 
# `x` and `y` are located. Then remove the pre-assigned tag numbers and replace 
# them with tag numbers that are created by the `rectangle` function. To make 
# the measured positions `x` and `y` similar to the supposed real positions 
# `X` and `Y`, we add the minimal values of `x` and `y` to `X` and `Y` resp.
#
# The reason why tag numbers should be replaced is because there are less tags
# (60) than measured positions (90). We later want to group by tags, and using 
# their original labels would then lead us into trouble.
stationary <- stationary %>% 
    mutate(X = assign_row(x, 10),
           Y = assign_row(y, 8)) %>% 
    select(-tag) %>% 
    plyr::join(rectangle(c(10, 8)), 
               by = c("X", "Y")) %>% 
    mutate(X = X + min(x), 
           Y = Y + min(y))

# Given the visualization below, an extra correction is needed. More specifically,
# we find that the measured size of the y-direction is smaller than the
# actual size, falling somewhat short of 8 meters. To allow us to correct for
# this, we should try to center the measured space into the actual space, which
# is what we do in the following bit of code.
#
# Importantly, no such correction seems to be needed in the x-direction.
stationary <- stationary %>% 
    mutate(Y = Y - abs(max(Y) - max(y))/2)

# Visualize whether our approximation of the real positions is adequate
visualize_positions(stationary) +
    geom_hline(yintercept = unique(stationary$Y),
               color = "red") +
    geom_vline(xintercept = unique(stationary$X),
               color = "red")

# Save these data as being preprocessed
saveRDS(stationary,
        file.path("data", 
                  "calibration", 
                  "stationary", 
                  "preprocessed_stationary_14-10-2023.Rds"))





#-------------------------------------------------------------------------------
# Stationary data: 21/10/2023
#-------------------------------------------------------------------------------

# Get the stationary data for the second calibration experiment. Here, there is 
# less data than in the previous one, and therefore less to preprocess.
# Importantly, this only uses the data for the tags that were put on the floor: 
# Stationary data when tags where in headbands is filtered out.
#
# The process is the same as before.
stationary <- data %>% 
    filter(experiment == "STATIONARY 2")

# Visualize these data
visualize_positions(stationary)

# Use the same ideas as before to preprocess these data and assign them positions
# that we can assume are close to the actual positions.
stationary <- stationary %>% 
    mutate(X = assign_row(x, 10),
           Y = assign_row(y, 8)) %>% 
    select(-tag) %>% 
    plyr::join(rectangle(c(10, 8)), 
               by = c("X", "Y")) %>% 
    mutate(X = X + min(x), 
           Y = Y + min(y))

# Visualize whether our approximation of the real positions is adequate.
#
# It seems that this time, we don't really need a correction. Unfortunately, 
# however, there does seem to be some systematic push to the right in the 
# x-direction due to outliers on row 3, column 1. This cannot be corrected in 
# a standardized way, and we will leave it at that.
visualize_positions(stationary) +
    geom_hline(yintercept = unique(stationary$Y),
               color = "red") +
    geom_vline(xintercept = unique(stationary$X),
               color = "red")

# Save these data as being preprocessed
saveRDS(stationary,
        file.path("data", 
                  "calibration", 
                  "stationary", 
                  "preprocessed_stationary_21-10-2023.Rds"))





#-------------------------------------------------------------------------------
# Stationary data: 22/12/2023
#-------------------------------------------------------------------------------

# Get the stationary data for the third calibration experiment. This time, we 
# did the calibration in the JK building. In this calibration round, we also 
# measured with either 4 or 6 anchors.
stationary <- data %>% 
    filter(str_sub(experiment, 1, 12) == "stationarity") %>%
    filter(str_sub(experiment, -10, -1) == "22-12-2023")

# Visualize these data
visualize_positions(stationary)

# Create a function that will assign a measured row number that corresponds to 
# the row numbers given in the `rectangle` function. It simply works by standardizing 
# the measured positions, multiplying it by the number of rows in that specific
# direction (x or y), and then rounding of the value.
assign_row <- function(x, n_rows){
    (x - min(x)) %>% 
        `/` (max(x) - min(x)) %>%
        `*` (n_rows) %>% 
        round() %>% 
        return()
}

# Use the same ideas as before to preprocess these data and assign them positions
# that we can assume are close to the actual positions. Here, we need to take 
# some extra precautions, as there are some datapoints that skew the assignment
# of the grid. Therefore done in two separate steps.
stationary <- stationary %>% 
    # First mutate the x-direction: Anomalies are located in this direction. 
    #
    # First, we filter based on the positions of the anchors.
    # Then we filter based on the row that is assigned to the data. 
    # Importantly, x should be corrected so that it starts at 0, allowing us to 
    # adequately filter out the anomalies.
    filter(between(x, 2.348, 14.35)) %>% 
    mutate(X = assign_row(x, 10)) %>% 
    filter(between(x - min(x), min(X), max(X))) %>% 
    # Now, reasign the rows for the x-direction and also handle the y-direction.
    # The reasignment is done to correct for any distortions that the now-
    # deleted anomalies might caused in the assingment of the rows.
    mutate(X = assign_row(x, 10), 
           Y = assign_row(y, 7)) %>% 
    # Delete the tags and use the function `rectangle` to produce new "tags" 
    # that correspond to individual positions within the grid. 
    select(-tag) %>% 
    plyr::join(rectangle(c(10, 7)), 
               by = c("X", "Y")) %>% 
    # Correct `X` and `Y` to start at the positions of `x` and `y`
    mutate(X = X + min(x), 
           Y = Y + min(y))

# Visualize whether our approximation of the real positions is adequate.
#
# It seems that this time, we don't really need a correction. Unfortunately, 
# however, there does seem to be some systematic push to the right in the 
# x-direction due to outliers on row 3, column 1. This cannot be corrected in 
# a standardized way, and we will leave it at that.
visualize_positions(stationary) +
    geom_hline(yintercept = unique(stationary$Y),
               color = "red") +
    geom_vline(xintercept = unique(stationary$X),
               color = "red")

# Save these data as being preprocessed
saveRDS(stationary,
        file.path("data",  
                  "calibration", 
                  "stationary", 
                  "preprocessed_stationary_22-12-2023.Rds"))

# As an additional thing, we will also save two kinds of these data, depending 
# on whether 4 anchors were used for the measurement or not. Then we can use 
# the code in `calibration_stationary` as it currently is.
six_anchors <- paste0("stationarity ", 1:4, " - 22-12-2023")

stationary %>% 
    filter(experiment %in% six_anchors) %>% 
    saveRDS(file.path("data",  
                      "calibration", 
                      "stationary", 
                      "preprocessed_stationary_6_anchors_22-12-2023.Rds"))

stationary %>% 
    filter(!(experiment %in% six_anchors)) %>% 
    saveRDS(file.path("data",  
                      "calibration", 
                      "stationary", 
                      "preprocessed_stationary_4_anchors_22-12-2023.Rds"))





#-------------------------------------------------------------------------------
# Moving data
#-------------------------------------------------------------------------------

# Get the data in which movement was introduced and measured.
#
# Make use of the fact that all data of these calibrations starts with the word 
# 'movement', on which you can select.
moving <- data %>% 
    mutate(partial_experiment = stringr::str_sub(experiment, end = 8)) %>% 
    filter(partial_experiment == "movement") %>% 
    select(-partial_experiment)

# Check whether there are any "real" experiments and let these replace their 
# "fake" counterparts. This is important to let our method of extracting the 
# correct information work.
real <- moving %>% 
    mutate(real_experiment = stringr::str_sub(experiment, start = -4)) %>% 
    filter(real_experiment == "real") %>% 
    select(-real_experiment) %>% 
    mutate(experiment_to_replace = stringr::str_sub(experiment, end = -6))

moving <- moving %>% 
    # Make sure you can adequately match the data in `real` to the data in this 
    # dataframe before filtering it out
    mutate(partial_experiment = stringr::str_sub(experiment, end = 14)) %>% 
    filter(partial_experiment %notin% unique(real$experiment_to_replace)) %>% 
    select(-partial_experiment)

moving <- real %>% 
    select(-experiment) %>% 
    rename(experiment = experiment_to_replace) %>% 
    rbind(moving)

# Extract the necessary information from the label we gave each of the movement 
# experiments by splitting the strings. Then convert this information to a 
# dataframe and bind it to the `moving` data. Afterwards, shuffle the columns
# so that the dataframe has a similar format to the `stationary` data.
moving <- stringr::str_split(moving$experiment, pattern = " ") %>% 
    as.data.frame() %>% 
    t() %>% 
    as_tibble(rownames = NULL) %>% 
    setNames(c("movement", "row_measured", "volts", "V")) %>% 
    select(row_measured, volts) %>% 
    cbind(moving) %>% 
    select(timestamp:experiment, row_measured, volts) %>% 
    mutate(volts = as.numeric(volts))

# Determine the expermental speed at which the tag is moving. More specifically, 
# for each of the volts and rows that were measured, append the measured time 
# it took for the tag to move from one side of the grid to the other.
moving <- fread(file.path("data", "moving_times_14-10-2023.txt"),
                data.table = FALSE) %>% 
    plyr::join(moving) %>% 
    select(timestamp:experiment, row_measured:seconds)

# Delete the data that fell outside of the grid. For this, we can use the
# `stationary` dataset that we measured on the same day, it being the one of
# October 14th.
# 
# This is required to have an adequate idea of speed, as it was measured from
# one side of the grid to the other. 
stationary <- readRDS(file.path("data",  
                                "calibration", 
                                "stationary", 
                                "preprocessed_stationary_14-10-2023.Rds"))

# Determine the limits in the y and x direction
ylims <- stationary %>% 
    summarize(min_y = min(Y), 
              max_y = max(Y))
xlims <- stationary %>% 
    summarize(min_x = min(X), 
              max_x = max(X))

# Delete the data that falls outside of these limits
moving <- moving %>% 
    filter(x > xlims$min_x,
           x < xlims$max_x, 
           y > ylims$min_y, 
           y < ylims$max_y)

# TO DO: Filter the data based on the computations done on the stationary data.

# Visualize these data together with the approximated bounds in which they fall.
#
# Does not seem to completely fit, so find out what is going wrong here.
visualize_positions(moving) +
    geom_hline(yintercept = unique(unlist(ylims)),
               color = "red") +
    geom_vline(xintercept = unique(unlist(xlims)),
               color = "red")

# Save these data as being preprocessed
saveRDS(moving,
        file.path("data", 
                  "calibration", 
                  "movement",
                  "preprocessed_movement_14-10-2023.Rds"))





#-------------------------------------------------------------------------------
# Headband data
#-------------------------------------------------------------------------------

# Get the data in which we measured stationary tags with headbands.
#
# Make use of the fact that all data of these calibrations starts with the word 
# 'headbands', on which you can select.
headband <- data %>% 
    mutate(partial_experiment = stringr::str_sub(experiment, end = 9)) %>% 
    filter(partial_experiment == "headbands") %>% 
    select(-partial_experiment)

# TO DO: Filter the data based on the computations done on the stationary data.

# Save these data as being preprocessed
saveRDS(headband,
        file.path("data",
                  "calibration",
                  "headbands",
                  "preprocessed_headbands_21-10-2023.Rds"))



