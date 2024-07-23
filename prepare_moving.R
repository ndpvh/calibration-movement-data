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



