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
        file.path("data", "preprocessed_stationary_14-10-2023.Rds"))





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
        file.path("data", "preprocessed_stationary_21-10-2023.Rds"))
