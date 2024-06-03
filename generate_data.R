################################################################################
# Purpose: Generate some synthetic data that can be used to initially test our # 
#          preprocessing strategies. For this, we will use the predped package #
#          to simulate pedestrian movement. Then, we add different kinds of    #
#          measurement error to these data.                                    #
#                                                                              #
#          The primary reason for generating these data is to find out which   #
#          methods are best to use for filtering the data and ridding it of    #
#          unsystematic error.                                                 #
################################################################################

devtools::load_all()

################################################################################
# SIMULATION

# Create a setting in which the agents will walk around. This will just be 
# a very simple environment. The room will be 10m x 5m, and the coordinate system
# will start in the lower-left corner at (0, 0). There will be 10 objects to 
# interact with
setting <- predped::background(shape = predped::rectangle(center = c(5, 2.5),
                                                          size = c(10, 5)),
                               objects = list(predped::rectangle(center = c(2.1, 1.1),
                                                                 size = c(0.2, 0.2),
                                                                 interactable = TRUE),
                                              predped::rectangle(center = c(2.1, 2.5),
                                                                 size = c(0.2, 0.2),
                                                                 interactable = TRUE),
                                              predped::rectangle(center = c(2.1, 3.9),
                                                                 size = c(0.2, 0.2),
                                                                 interactable = TRUE),
         
                                              predped::rectangle(center = c(7.9, 1.1),
                                                                 size = c(0.2, 0.2),
                                                                 interactable = TRUE),
                                              predped::rectangle(center = c(7.9, 2.5),
                                                                 size = c(0.2, 0.2),
                                                                 interactable = TRUE),
                                              predped::rectangle(center = c(7.9, 3.9),
                                                                 size = c(0.2, 0.2),
                                                                 interactable = TRUE),
         
                                              predped::rectangle(center = c(4.2, 1.1),
                                                                 size = c(0.2, 0.2),
                                                                 interactable = TRUE),
                                              predped::rectangle(center = c(5.8, 1.1),
                                                                 size = c(0.2, 0.2),
                                                                 interactable = TRUE),
         
                                              predped::rectangle(center = c(4.2, 3.9),
                                                                 size = c(0.2, 0.2),
                                                                 interactable = TRUE),
                                              predped::rectangle(center = c(5.8, 3.9),
                                                                 size = c(0.2, 0.2),
                                                                 interactable = TRUE)),
                                entrance = c(0, 2.5))

# Create a model. Here, I will just keep one archetype of people in there to 
# simplify everything a bit.
model <- predped::predped(id = "synthetic",
                          setting = setting, 
                          archetypes = "BaselineEuropean")

# Simulate synthetic data. Some things to note here: 
#   - Sampling rate is set to 10Hz, from which we will then randomly downsample
#     (done to create unequal spacing in time). The argument `time_step` is 
#     thus set to 1/10
#   - Number of iterations is set to have about 15 minutes of simulation, i.e.
#     to 15 * 60 * 10 = 9000. The is the amount of data that we will have per
#     participant 
#   - Number of participants is set to 10, with no way to get more participants
#   - Number of goals is set to an impossible amount to reach
#
# We run these simulations 100 times, allowing us to estimate the variability 
# around the efficacy of a given preprocessing pipeline. One simulation should 
# take about 22min 4sec, so do give yourself some time to run this (about 37 
# hours total, not in parallel).

# Create the different seeds that will be used to create the data.
set.seed(74327) # Knowing That You've Arrived - Tides of Man
seeds <- sample(1:10000, 100)

# Create a function that will take in the seeds and output the data
simulate_data <- function(x) {
    print(x)

    set.seed(seeds[x])

    # Create the trace using the predped package
    trace <- predped::simulate(model, 
                               initial_number_agents = 10, 
                               max_agents = 10, 
                               iterations = 9000, 
                               time_step = 1 / 10,
                               goal_number = 1000, 
                               goal_duration = 10, 
                               order_goal_stack = FALSE)
    
    # Load the data and transform to a dataframe. Save this dataframe as being the 
    # "real data". These data will have the following columns: 
    #      - x, y: Coordinates at which the pedestrian was standing
    #      - time: A variable denoting the time that has passed in seconds
    #      - id: The name of the pedestrian
    data <- lapply(seq_along(trace), 
                   \(i) sapply(trace[[i]]$agents, 
                               \(y) c(predped::id(y), 
                                      predped::position(y))) %>% 
                       t() %>% 
                       as.data.frame() %>% 
                       setNames(c("id", "x", "y")) %>% 
                       cbind(time = (i - 1)/10) %>% 
                       dplyr::select(time, id:y))
    data <- do.call("rbind", data)
    
    # Make sure everything is in the correct format
    data <- data %>% 
        dplyr::mutate(x = as.numeric(x), 
                      y = as.numeric(y), 
                      time = as.numeric(time))

    return(data)
}

# Carry out the simulation in parallel
n_cores <- parallel::detectCores()
data <- parallel::mclapply(seq_along(seeds), 
                           \(x) simulate_data(x),
                           mc.cores = getOption("mc.cores", n_cores - 1))

data <- do.call("rbind", data)

# Save these data
data.table::fwrite(data, file.path("data", "synthetic", "synthetic_original.csv"))





################################################################################
# MEASUREMENT ERROR

####################################
# Create a dataset with random measurement error (no relation between the 
# different dimensions)

# Load the original dataset
data <- data.table::fread(file.path("data", "synthetic", "synthetic_original.csv"), 
                          data.table = FALSE)

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
set.seed(546) # Aaahh!!! Real Spiders - Wind Walkers
data <- data %>% 
    dplyr::rename(x_original = x, 
                  y_original = y) %>% 
    dplyr::mutate(x = rnorm(length(x_original), x_original, 0.031), 
                  y = rnorm(length(y_original), y_original, 0.031))

# Save these data
data.table::fwrite(data, file.path("data", "synthetic", "synthetic_unrelated_10.csv"))



####################################
# Create a dataset with random measurement error (relationship between error 
# in both dimensions, correlation of 0.25, which is much more than observed)

# Load the original dataset
data <- data.table::fread(file.path("data", "synthetic", "synthetic_original.csv"), 
                          data.table = FALSE)

# Transform the dataset to include random error. We use a similar approach as 
# earlier, but no also include correlations. To do this, we first create a 
# matrix that contains standard deviations on its diagonal, and another matrix 
# that contains the correlations on its off-diagonal while having 1's on its 
# diagonal. We then use these matrices to create the covariance matrix, from
# which we then generate the residuals to be added to the x and y coordinates
SD <- diag(rep(0.031, 2))
COR <- c(1, 0.25, 0.25, 1) %>% 
    matrix(ncol = 2)
COV <- SD %*% COR %*% SD

set.seed(7244) # Falling on Deaf Ears - Hail the Sun
residuals <- MASS::mvrnorm(nrow(data), c(0, 0), COV)

data <- data %>% 
    dplyr::rename(x_original = x, 
                  y_original = y) %>% 
    dplyr::mutate(x = x_original + residuals[,1], 
                  y = y_original + residuals[,2])

# Save these data
data.table::fwrite(data, file.path("data", "synthetic", "synthetic_related_10.csv"))



####################################
# Create unequal time bins by random missingness. Keep a sampling rate of about
# 6Hz (observed in data)

# Load the two measurement error datasets
data <- list(data.table::fread(file.path("data", "synthetic", "synthetic_unrelated_10.csv"),
                               data.table = FALSE),
             data.table::fread(file.path("data", "synthetic", "synthetic_related_10.csv"),
                               data.table = FALSE))

# Determine which iterations you will delete. Number of deletions will impose a 
# mean sampling rate of 6Hz across different people
sample_idx <- \(x) sample(seq_along(x), 
                          round(0.6 * length(x)))

set.seed(313) # 3's & 7's - Queens of the Stone Age
data <- lapply(data, 
               \(x) x %>% 
                   dplyr::group_by(nsim) %>% 
                   dplyr::mutate(idx = dplyr::row_number() %in% sample_idx(y)) %>% 
                   dplyr::filter(idx) %>% 
                   dplyr::select(-idx))

# And save the results
data.table::fwrite(data[[1]], 
                   file.path("data", "synthetic", "synthetic_unrelated_6random.csv"))
data.table::fwrite(data[[2]], 
                   file.path("data", "synthetic", "synthetic_related_6random.csv"))



####################################
# Create unequal time bins by nonrandom missingness. Keep a sampling rate of 
# about 6Hz (observed in the data)

# Load the two measurement error datasets
data <- list(data.table::fread(file.path("data", "synthetic", "synthetic_unrelated_10.csv"),
                               data.table = FALSE),
             data.table::fread(file.path("data", "synthetic", "synthetic_related_10.csv"),
                               data.table = FALSE))

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
N <- nrow(data[[1]])
impute_missing <- function(x) {
    # Create relative indices per person for the blocked missing data
    idx <- cbind(sample(1:length(unique(x$time)),                  # Sample time points to delete
                        round(0.3 * N / 10), replace = TRUE),      # As many as to have about 75% of missing data being in longer blocks
                 sample(unique(data[[1]]$id),                      # Sample participants whose date to delete in unison
                        round(0.3 * N / 10), replace = TRUE)) %>%  # As many as there are blocks to delete         
        as.data.frame() %>% 
        setNames(c("from", "participant")) %>% 
        dplyr::mutate(from = as.numeric(from)) %>% 
        dplyr::mutate(to = from + rep(c(5, 10, 15) - 1, each = round(length(from) / 3)))

    # Explicate all indices to be deleted
    idx <- idx %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(indices = multi_seq(from, to) %>% 
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
        dplyr::mutate(index = dplyr::row_number() %in% idx$indices[idx$participant == id]) %>%
        dplyr::filter(!index) %>% 
        dplyr::select(-index) %>% 
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
data <- lapply(data, 
               \(x) x %>% 
                   dplyr::group_by(nsim) %>% 
                   impute_missing())

# And save the results
data.table::fwrite(data[[1]], 
                   file.path("data", "synthetic", "synthetic_unrelated_6nonrandom.csv"))
data.table::fwrite(data[[2]], 
                   file.path("data", "synthetic", "synthetic_related_6nonrandom.csv"))





################################################################################
# VISUALIZATION

# Names of all the files to make gifs of off
names <- c("synthetic_original", 
           "synthetic_unrelated_10",
           "synthetic_related_10",
           "synthetic_unrelated_6random", 
           "synthetic_related_6random",
           "synthetic_unrelated_6nonrandom",
           "synthetic_related_6nonrandom")

# Loop over all of these files
for(i in names) {
    data <- data.table::fread(file.path("data", 
                                        "synthetic",
                                        paste0(i, ".csv")), 
                              data.table = FALSE)

    plt <- plot(data)
    gifski::save_gif(lapply(plt, \(x) print(x)), 
                     file.path("data", 
                             "synthetic", 
                             paste0(i, ".gif")),
                     width = 2000, 
                     height = 1000,
                     delay = 1/10)
}
