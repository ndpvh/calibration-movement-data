################################################################################
# Purpose: Generate some synthetic data that can be used to initially test our # 
#          preprocessing strategies. For this, we will use the predped package #
#          to simulate pedestrian movement. Then, we add different kinds of    #
#          measurement error to these data.                                    #
################################################################################

library(predped)
library(tidyverse)

################################################################################
# SIMULATION

# Create a setting in which the agents will walk around. This will just be 
# a very simple environment. The room will be 10m x 5m, and the coordinate system
# will start in the lower-left corner at (0, 0). There will be 10 objects to 
# interact with
setting <- background(shape = predped::rectangle(center = c(5, 2.5),
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
set.seed(74327) # Tides of Man - Knowing That You've Arrived
trace <- predped::simulate(model, 
                           initial_number_agents = 10, 
                           max_agents = 10, 
                           iterations = 9000, 
                           time_step = 1 / 10,
                           goal_number = 1000, 
                           goal_duration = 10, 
                           order_goal_stack = FALSE)

# Make into plots and a gif for inspection
#   - Get size of rectangle of setting
#   - Create plots
#   - Create gif
points <- shape(setting)@points
poly_size <- c(max(points[,1] - min(points[,1])),
               max(points[,2] - min(points[,2])))

plt <- plot(trace, trace = TRUE, linewidth = 1.5)
plt <- lapply(plt, 
              \(x) x + 
                   ggplot2::theme(legend.position = "none",
                                  plot.title = ggplot2::element_text(size = 5 * max(poly_size),
                                                                     hjust = 0.5),
                                  axis.text = ggplot2::element_text(size = 2.5 * max(poly_size))))

gifski::save_gif(lapply(plt, \(x) print(x)),                                     
                 file.path("data", "synthetic", "synthetic_trace.gif"),
                 delay = 1/50, 
                 width = poly_size[1] * 200, 
                 height = poly_size[2] * 200)

# Save the results
saveRDS(trace, file.path("data", "synthetic", "synthetic_trace.Rds"))





################################################################################
# MEASUREMENT ERROR

# Load the data and transform to a dataframe. Save this dataframe as being the 
# "real data"
trace <- readRDS(file.path("data", "synthetic", "synthetic_trace.Rds"))

data <- lapply(seq_along(trace), 
               \(i) sapply(trace[[i]]$agents, 
                           \(y) c(predped::id(y), 
                                  predped::position(y))) %>% 
                   t() %>% 
                   as.data.frame() %>% 
                   setNames(c("id", "x", "y")) %>% 
                   cbind(time = i) %>% 
                   select(time, id:y))
data <- do.call("rbind", data)

saveRDS(data. file.path("data", "synthetic", "synthetic_original.Rds"))
