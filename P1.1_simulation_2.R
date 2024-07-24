################################################################################
# Purpose: Generate some synthetic data that can be used to initially test our # 
#          preprocessing strategies. For this, we will use the predped package #
#          to simulate pedestrian movement.                                    #
#                                                                              #
#          RQ: Which method performs best in filtering the data and ridding    #
#              it of unsystematic error.                                       #
#                                                                              #
#          This is Part 1.1 of Simulation 2, meaning that we only generate the #
#          "real" positions of the pedestrians at each time point. Measurement #
#          error is added in tandem with Simulation 1 in Part 1.2.             #
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
data$nsim <- rep(seq_along(seeds), each = 9001 * 10)

# Save these data
data.table::fwrite(data, file.path("data", "simulation_2", "data.csv"))