################################################################################
# Purpose: Prepare the calibration data that was gathered on 14/10/2023        #
#          and 21/10/2023 for analysis. These data are used to estimate        #
#          the amount of systematic and unsystematic noise in the movement     #
#          data. We additionally use these data to find out what the actual    #
#          sampling rate of the data is (expected 5Hz)                         #
################################################################################

devtools::load_all()
data <- load_data()

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

#-------------------------------------------------------------------------------
# Stationary data: 14/10/2023
#-------------------------------------------------------------------------------

# EXPERIMENTAL INFORMATION: 
#
# During this experiment, several tags were laid down on the intersections of 
# several 1m x 1m grids. There were more intersections than tags, so we had to 
# measure the tag's positions in four different batches, each covering a 
# different part of the grid

# Combine all of the stationary data to one datafile. This way, we get measurements
# at all positions of the net that we installed, which we can then compare to 
# the positions these measurements should actually have.
#
# The process below first creates a boolean that indicates whether the
# experiments are the ones that are defined in `stationary_experiments`. Then,
# this boolean is used to filter the data to only contain these data.
stationary_experiments <- paste("stationary", c(2, 3, 5, 6), sep = " ")
stationary <- data %>% 
    dplyr::filter(experiment %in% stationary_experiments) %>% 
    dplyr::mutate(time = withr::with_options(list(digits = 16),
                                             as.numeric(timestamp)), 
                  time = time - min(time)) %>% 
    dplyr::rename(id = tag)

# Use the `assign_row` function to determine in which row each measured position 
# `x` and `y` are located. Then remove the pre-assigned tag numbers and replace 
# them with new id's. Finally, To make we add the mean minimal values of `x` and 
# `y` to `X` and `Y` resp. to make the measured positions `x` and `y` similar to 
# the supposed real positions.
#
# The reason why tag numbers should be replaced is because there are less tags
# (60) than measured positions (99). We later want to group by tags, and using 
# their original labels would then lead us into trouble.
stationary <- stationary %>% 
    # Assign X, Y, and id's
    dplyr::mutate(X = assign_row(x, 10),
                  Y = assign_row(y, 8),
                  id = paste0(X, Y) %>% 
                      factor() %>% 
                      as.numeric()) %>% 
    # Group by the id numbers, compute the mean X and Y values per id, and add 
    # the measured center of the space to X and Y. This will bring their values 
    # closer to the real positions of the tags
    dplyr::mutate(X = X - mean(X) + mean(x), 
                  Y = Y - mean(Y) + mean(y)) %>% 
    # Delete data that consists of only a few observation (within the same 
    # measuring session)
    dplyr::group_by(experiment, id) %>% 
    tidyr::nest() %>% 
    dplyr::mutate(rows = nrow(as.data.frame(data))) %>% 
    dplyr::filter(rows > 10) %>% 
    dplyr::select(-rows) %>% 
    tidyr::unnest(data) %>% 
    dplyr::ungroup()

# Visualize the measured position and the real positions
plot(stationary, 
     as_points = TRUE, 
     per_iteration = FALSE) +
    ggplot2::geom_point(ggplot2::aes(x = X, y = Y), 
                        color = "red") +
    ggplot2::lims(x = range(stationary$X) + c(-1, 1) * 0.5,
                  y = range(stationary$Y) + c(-1, 1) * 0.5)

# Save these data as being preprocessed
saveRDS(stationary,
        file.path("data", 
                  "stationary", 
                  "stationary_14-10-2023.Rds"))





#-------------------------------------------------------------------------------
# Stationary data: 21/10/2023
#-------------------------------------------------------------------------------

# EXPERIMENTAL INFORMATION: 
#
# Same procedure as the experiment at the 14th.

# Get the stationary data for the second calibration experiment. Here, there is 
# less data than in the previous one, and therefore less to preprocess.
# Importantly, this only uses the data for the tags that were put on the floor: 
# Stationary data when tags where in headbands is filtered out.
#
# The process is the same as before.
stationary <- data %>% 
    dplyr::filter(experiment == "STATIONARY 2") %>% 
    dplyr::mutate(time = withr::with_options(list(digits = 16),
                                             as.numeric(timestamp)), 
                  time = time - min(time)) %>% 
    dplyr::rename(id = tag)

# Use the `assign_row` function to determine in which row each measured position 
# `x` and `y` are located. Then remove the pre-assigned tag numbers and replace 
# them with new id's. Finally, To make we add the mean minimal values of `x` and 
# `y` to `X` and `Y` resp. to make the measured positions `x` and `y` similar to 
# the supposed real positions.
stationary <- stationary %>% 
    # Assign X, Y, and id's
    dplyr::mutate(X = assign_row(x, 10),
                  Y = assign_row(y, 8),
                  id = paste0(X, Y) %>% 
                      factor() %>% 
                      as.numeric()) %>% 
    # Group by the id numbers, compute the mean X and Y values per id, and add 
    # the measured center of the space to X and Y. This will bring their values 
    # closer to the real positions of the tags
    dplyr::mutate(X = X - mean(X) + mean(x), 
                  Y = Y - mean(Y) + mean(y)) %>% 
    # Delete data that consists of only a few observation (within the same 
    # measuring session)
    dplyr::group_by(experiment, id) %>% 
    tidyr::nest() %>% 
    dplyr::mutate(rows = nrow(as.data.frame(data))) %>% 
    dplyr::filter(rows > 10) %>% 
    dplyr::select(-rows) %>% 
    tidyr::unnest(data) %>% 
    dplyr::ungroup()

# Visualize the measured position and the real positions
plot(stationary, 
     as_points = TRUE, 
     per_iteration = FALSE) +
    ggplot2::geom_point(ggplot2::aes(x = X, y = Y), 
                        color = "red") +
    ggplot2::lims(x = range(stationary$X) + c(-1, 1) * 0.5,
                  y = range(stationary$Y) + c(-1, 1) * 0.5)

# Save these data as being preprocessed
saveRDS(stationary,
        file.path("data", 
                  "stationary", 
                  "stationary_21-10-2023.Rds"))





#-------------------------------------------------------------------------------
# Stationary data: 22/12/2023
#-------------------------------------------------------------------------------

# EXPERIMENTAL INFORMATION: 
#
# This one took place in the JK building. In this calibration round, we varied
# the number of anchors with which we measured the position (4 or 6 anchors).
# We distinguish between having all data (4 and 6 together) and having only the 
# data for one of the conditions.

# Get the stationary data for the third calibration experiment. 
stationary <- data %>% 
    dplyr::filter(stringr::str_sub(experiment, 1, 12) == "stationarity") %>%
    dplyr::filter(stringr::str_sub(experiment, -10, -1) == "22-12-2023") %>% 
    dplyr::mutate(time = withr::with_options(list(digits = 16),
                                             as.numeric(timestamp)), 
                  time = time - min(time)) %>% 
    dplyr::rename(id = tag)

# There are some data outside of the grid in this data set, which was probably
# due to some tags remaining active at the control booth.
#
# Delete these data first (otherwise all goes wrong)
stationary <- stationary %>% 
    # Delete data outside of a given bound
    dplyr::filter(x < 13.5) 

# Assign X, Y, and id's
stationary <- stationary %>% 
    dplyr::mutate(X = assign_row(x, 10),
                  Y = assign_row(y, 7),
                  id = paste0(X, Y) %>% 
                      factor() %>% 
                      as.numeric()) %>% 
    # Group by the id numbers, compute the mean X and Y values per id, and add 
    # the measured center of the space to X and Y. This will bring their values 
    # closer to the real positions of the tags
    dplyr::mutate(X = X - mean(X) + mean(x), 
                  Y = Y - mean(Y) + mean(y)) %>% 
    # Delete data that consists of only a few observation (within the same 
    # measuring session)
    dplyr::group_by(experiment, id) %>% 
    tidyr::nest() %>% 
    dplyr::mutate(rows = nrow(as.data.frame(data))) %>% 
    dplyr::filter(rows > 10) %>% 
    dplyr::select(-rows) %>% 
    tidyr::unnest(data) %>% 
    dplyr::ungroup()

# Add in information about how many anchors were used to measure these positions.
six_anchors <- paste0("stationarity ", 1:4, " - 22-12-2023")
stationary <- stationary %>% 
    dplyr::mutate(anchors = ifelse(experiment %in% six_anchors, 6, 4)) 

# Visualize the measured position and the real positions
plot(stationary, 
     as_points = TRUE, 
     per_iteration = FALSE) +
    ggplot2::geom_point(ggplot2::aes(x = X, y = Y), 
                        color = "red") +
    ggplot2::lims(x = range(stationary$X) + c(-1, 1) * 0.5,
                  y = range(stationary$Y) + c(-1, 1) * 0.5)

# Save these data as being preprocessed
saveRDS(stationary,
        file.path("data",
                  "stationary", 
                  "stationary_22-12-2023.Rds"))

# NOTE:
#
# For some reason, some tags didn't produce a lot of data when only 4 anchors 
# were used and have, consequentially, been deleted. Might be useful to do the 
# experiment again at some point