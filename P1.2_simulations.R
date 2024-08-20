################################################################################
# Purpose: This is Part 1.2 for simulating data. Here, we add measurement      #
#          error to the data sets generated in Part 1.1 for Simulation 1 and   # 
#          2. Both data sets will thus contain the same amount of noise.       #
#                                                                              #
#          RQ: Which method performs best in filtering the data and ridding    #
#              it of unsystematic error.                                       #
################################################################################

devtools::load_all()

################################################################################
# PRELIMINARIES

# Create a dataframe containing the necessary information to load the data sets
# when asked for
filenames <- data.frame(folder = c("simulation_1", "simulation_1", "simulation_2"), 
                        name = c("movement", "fixed", "data"))





################################################################################
# MEASUREMENT ERROR

# Introduce measurement error into the data. This is meant to mimic the kind of 
# measurement error that we find in the data. We distinguish between three types:
# (a) random and independent across dimensions, (b) random and dependent across
# dimensions, (c) time-dependent. Of this, the latter is closest to the actual 
# observations in the calibration data.

#-------------------------------------------------------------------------------
# Independence across dimensions
#-------------------------------------------------------------------------------

# Transform the dataset to include random error. This error will be equal to 
# a few centimeters, based on earlier estimates of the measurement error, 
# specifically coming from the calibration test of 22-12-2023, where the upper 
# bound of the 99%CI was about 8cm. Using this information, we can derive the 
# standard deviation as follows:
#    - Find qnorm(0.995, 0, 1), which is equal to about 2.58. This value 
#      communicates for what value of x the probability of falling below it is 
#      equal to 99.5%.
#    - Realizing that for our distribution, this value should be equal to 0.08, 
#      we can derive the needed standard deviation as being equal to 0.08 / 2.58, 
#      which is approximately 0.031
#    - To check our result, we can get the probability of falling below 8 for 
#      the specified distribution by using pnorm(0.08, 0, 0.031). This value 
#      should be equal to 99.5%, and indeed it is approximately so
set.seed(425) # Call it Karma - Silverstein
for(i in seq_len(nrow(filenames))) {
    # Add normally distributed measurement error to the data
    data <- data.table::fread(file.path("data", 
                                        filenames$folder[i], 
                                        paste0(filenames$name[i], ".csv")), 
                              data.table = FALSE) %>% 
    dplyr::rename(x_original = x, 
                  y_original = y) %>% 
    dplyr::mutate(x = rnorm(length(x_original), x_original, 0.031), 
                  y = rnorm(length(y_original), y_original, 0.031))

    # Save the data
    data.table::fwrite(data, 
                       file.path("data", 
                                 filenames$folder[i], 
                                 paste0(filenames$name[i], "_U10.csv")))
}





#-------------------------------------------------------------------------------
# Dependence across dimensions
#-------------------------------------------------------------------------------

# Transform the dataset to include random error. We use a similar approach as 
# earlier, but now also include correlations. To do this, we first create a 
# matrix that contains standard deviations on its diagonal, and another matrix 
# that contains the correlations on its off-diagonal while having 1's on its 
# diagonal. We then use these matrices to create the covariance matrix, from
# which we then generate the residuals to be added to the x and y coordinates
SD <- diag(rep(0.031, 2))
COR <- c(1, 0.25, 0.25, 1) %>% 
    matrix(ncol = 2)
COV <- SD %*% COR %*% SD

set.seed(7244) # Falling on Deaf Ears - Hail the Sun
for(i in seq_len(nrow(filenames))) {
    data <- data.table::fread(file.path("data", 
                                        filenames$folder[i], 
                                        paste0(filenames$name[i], ".csv")), 
                              data.table = FALSE)

    # Create the residuals 
    residuals <- MASS::mvrnorm(nrow(data), c(0, 0), COV)

    # Add normally distributed measurement error to the data
    data <- data %>% 
        dplyr::rename(x_original = x, 
                      y_original = y) %>% 
        dplyr::mutate(x = x_original + residuals[,1], 
                      y = y_original + residuals[,2])

    # Save the data
    data.table::fwrite(data, 
                       file.path("data", 
                                 filenames$folder[i], 
                                 paste0(filenames$name[i], "_R10.csv")))
}





#-------------------------------------------------------------------------------
# Time-dependence
#-------------------------------------------------------------------------------

# Here, we will need to create a vector autoregressive model that will account 
# both for contemporaneous and lagged measurement error. The parameters that 
# are used here are taken from estimations we did on the stationary calibration 
# data. The measurement error is added for each experiment and id separately.
#
# Start by creating the function that will create the measurement error to be 
# added to the observations.
add_residuals <- function(x) {
    # Create the parameters to be used
    B <- c(7.77e-1, -1.37e-3, 4.98e-3, 7.85e-1) %>% 
        matrix(nrow = 2, ncol = 2)
    S <- c(9.66e-5, 7.25e-7, 7.25e-7, 9.32e-5) %>% 
        matrix(nrow = 2, ncol = 2)

    # Create a matrix that will contain the residuals and add the initial 
    # condition
    residuals <- matrix(0, nrow = nrow(x), ncol = 2)
    residuals[1,] <- MASS::mvrnorm(1, c(0, 0), S)

    # Loop over and create the other residuals
    for(i in seq_len(nrow(x) - 1)) {
        residuals[i + 1,] <- B %*% residuals[i,] + MASS::mvrnorm(1, c(0, 0), S)
    }

    # Add the residuals to the dataframe and return
    x %>% 
        dplyr::arrange(time) %>% 
        dplyr::mutate(x = x_original + residuals[,1], 
                      y = y_original + residuals[,2]) %>% 
        return()
}

# Actually add the measurement error to the data. As mentioned above, done per 
# simulation and per tag separately
set.seed(55) # Sober Exit(s) - Static Dress
for(i in seq_len(nrow(filenames))) {
    # Add normally distributed measurement error to the data
    data <- data.table::fread(file.path("data", 
                                        filenames$folder[i], 
                                        paste0(filenames$name[i], ".csv")), 
                              data.table = FALSE) %>% 
        dplyr::rename(x_original = x, 
                      y_original = y) %>% 
        dplyr::group_by(nsim, id) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(data = data %>% 
                          as.data.frame() %>% 
                          add_residuals() %>% 
                          list()) %>% 
        tidyr::unnest(data) %>% 
        dplyr::arrange(nsim, time, id)

    # Save the data
    data.table::fwrite(data, 
                       file.path("data", 
                                 filenames$folder[i], 
                                 paste0(filenames$name[i], "_T10.csv")))
}





################################################################################
# SAMPLING FREQUENCY

# Our data often does not have an actual sampling frequency of 10Hz, but it is 
# ususally lower. Account for this kind of missing data, and the unequal time 
# intervals that this creates in the data (and which may make our job 
# difficult). We will retain a sampling frequency of 6Hz.

# Create new filenames
filenames <- data.frame(folder = c("simulation_1", "simulation_1", "simulation_1",
                                   "simulation_1", "simulation_1", "simulation_1",
                                   "simulation_2", "simulation_2", "simulation_2"), 
                        name = c("movement_R10", "movement_T10", "movement_U10", 
                                 "fixed_R10", "fixed_T10", "fixed_U10", 
                                 "data_R10", "data_T10", "data_U10"))

#-------------------------------------------------------------------------------
# Random missingness
#-------------------------------------------------------------------------------

# Determine which iterations you will delete. Number of deletions will impose a 
# mean sampling rate of 6Hz across different people
sample_idx <- \(x) sample(seq_along(x), 
                          round(0.6 * length(x)))

set.seed(313) # 3's & 7's - Queens of the Stone Age
for(i in seq_len(nrow(filenames))) {
    data <- data.table::fread(file.path("data", 
                                        filenames$folder[i], 
                                        paste0(filenames$name[i], ".csv")), 
                              data.table = FALSE) %>% 
        dplyr::group_by(nsim) %>% 
        dplyr::mutate(idx = dplyr::row_number() %in% sample_idx(y)) %>% 
        dplyr::filter(idx) %>% 
        dplyr::select(-idx)

    data.table::fwrite(data, 
                       file.path("data",
                                 filenames$folder[i], 
                                 paste0(stringr::str_sub(filenames$name[i], 0, -3), 
                                        "6R.csv")))
}





#-------------------------------------------------------------------------------
# Nonrandom missingness
#-------------------------------------------------------------------------------

# The approach we take here is an easy one, where we draw random time points at 
# which no position is measured. Once chosen, we then either keep it at this one
# time point, or we make the period at which no measurements are obtained longer 
# (5 observations, 10 observations, or 15 observations, it being 500msec, 1sec,
# or 1.5sec long). We try to approximate each time as being as long as the other.
#
# In practice, we delete the indices in two waves: One in which we take care of 
# the longer problems, then one in which we delete only a single data point 
# randomly up until as there are as many missing observations for the nonrandom 
# and the random missing data (created earlier)
impute_missing <- function(x) {
    N <- nrow(x)

    # Create relative indices per person for the blocked missing data
    idx <- cbind(sample(1:(length(unique(x$time)) - 15),           # Sample time points to delete (-15 so that you don't get overshoots in deletion)
                        round(0.3 * N / 10), replace = TRUE),      # As many as to have about 75% of missing data being in longer blocks
                 sample(unique(x$id),                              # Sample participants whose date to delete in unison
                        round(0.3 * N / 10), replace = TRUE)) %>%  # As many as there are blocks to delete         
        as.data.frame() %>% 
        setNames(c("from", "participant")) %>% 
        dplyr::mutate(from = as.numeric(from)) %>% 
        dplyr::mutate(to = from + rep(c(5, 10, 15) - 1, each = round(length(from) / 3)))

    # Explicate all indices to be deleted
    idx <- idx %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(indices = seq(from, to) %>% 
                          as.vector() %>% 
                          data.frame() %>% 
                          setNames("indices") %>% 
                          tidyr::nest()) %>%
        tidyr::unnest(indices) %>% 
        tidyr::unnest(data) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(-from, -to)

    # Delete these indices already in the way that was previously used. 
    # Importantly, the deletion is relative to the participant, so we have to 
    # account for the participant in this deletion.
    x <- x %>% 
        dplyr::group_by(id) %>% 
        tidyr::nest() %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(data = data %>%
            as.data.frame() %>%
            dplyr::mutate(index = dplyr::row_number() %in% idx$indices[idx$participant == id]) %>%
            dplyr::filter(!index) %>% 
            dplyr::select(-index) %>%
            list()) %>% 
        tidyr::unnest(data) %>%
        dplyr::ungroup()

    # Now sample the remaining time points to be deleted from the remaining data 
    # points
    sample_idx <- \(x) sample(seq_along(x), 
                              round(0.6 * N))
    x <- x %>% 
        dplyr::mutate(index = dplyr::row_number() %in% sample_idx(y)) %>% 
        dplyr::filter(index) %>% 
        dplyr::select(-index)

    return(x)
}

# Delete these time points from the data for each individual separately
set.seed(9410) # Newsstand Rock (exposition) - Rx Bandits
for(i in seq_len(nrow(filenames))) {
    data <- data.table::fread(file.path("data", 
                                        filenames$folder[i], 
                                        paste0(filenames$name[i], ".csv")), 
                              data.table = FALSE) %>% 
        dplyr::group_by(nsim) %>% 
        tidyr::nest() %>%
        dplyr::mutate(data = data %>%  
            as.data.frame() %>%
            impute_missing() %>%
            list()) %>%
        tidyr::unnest(data)

    data.table::fwrite(data, 
                       file.path("data",
                                 filenames$folder[i], 
                                 paste0(stringr::str_sub(filenames$name[i], 0, -3), 
                                        "6N.csv")))
}
