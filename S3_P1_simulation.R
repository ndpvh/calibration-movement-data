################################################################################
# Purpose: Generate some more complex synthetic data that can be used to test  # 
#          some strategies of filtering the unsystematic error from the data.  #
#          Here, we use the predped package to simulate pedestrian movement    #
#          in a realistic scenario.                                            #
#                                                                              #
#          Structure of the file is as follows:                                #
#              L18-150:    Simulating positions                                #
#              L156-417:   Adding noise                                        #
#                                                                              #
#          Note that when adding noise, we add both time-independent noise     #
#          and a time-dependent noise. We use the parameters estimated in the  #
#          S1_P1_analysis script for this purpose.                             #
################################################################################

devtools::load_all()

################################################################################
# SIMULATION

# Create a setting in which the agents will walk around. This will just be 
# a very simple environment. The room will be 10m x 5m, and the coordinate system
# will start in the lower-left corner at (0, 0). There will be 10 objects to 
# interact with
setting <- predped::background(
    # Rectangular environment to move into
    shape = predped::rectangle(
        center = c(10, 5),
        size = c(20, 10)
    ),

    # Objects within the environment
    objects = list(
        predped::rectangle(
            id = "table",
            center = c(10, 5),
            size = c(0.5, 0.5)
        ),
        predped::rectangle(
            id = "tablet 1",
            center = c(2.1, 1.1),
            size = c(0.1, 0.2),
            orientation = pi/4
        ),
        predped::rectangle(
            id = "tablet 2",
            center = c(2.1, 5),
            size = c(0.1, 0.2)
        ),
        predped::rectangle(
            id = "tablet 3",
            center = c(2.1, 8.9),
            size = c(0.2, 0.1),
            orientation = pi/4
        ),
        predped::rectangle(
            id = "tablet 4",
            center = c(7.4, 8.9),
            size = c(0.2, 0.1)
        ),
        predped::rectangle(
            id = "tablet 5",
            center = c(12.6, 8.9),
            size = c(0.2, 0.1)
        ),
        predped::rectangle(
            id = "tablet 6",
            center = c(17.9, 8.9),
            size = c(0.1, 0.2),
            orientation = pi/4
        ),
        predped::rectangle(,
            id = "tablet 7",
            center = c(17.9, 5),
            size = c(0.1, 0.2)
        ),
        predped::rectangle(
            id = "tablet 8",
            center = c(17.9, 1.1),
            size = c(0.2, 0.1),
            orientation = pi/4
        ),
        predped::rectangle(
            id = "tablet 9",
            center = c(7.4, 1.1),
            size = c(0.2, 0.1)
        ),
        predped::rectangle(
            id = "tablet 10",
            center = c(12.6, 1.1),
            size = c(0.2, 0.1)
        )
    ),

    # Define the entrance to the environment
    entrance = c(0, 2.5)
)

# Plot the environment
predped::plot(setting, dark_mode = TRUE)

# Knit the setting together with the parameters to use within the simulation
model <- predped::predped(
    setting = setting, 
    archetypes = "BaselineEuropean"
)

# Create a list of goals that the agents will be assigned at each time. 
goals <- list(
    predped::goal(
        id = "table 1", 
        position = c(10, 5) + c(0, 0.55)/2,
        counter = 10
    ),
    predped::goal(
        id = "table 2", 
        position = c(10, 5) + c(0, -0.55)/2,
        counter = 10
    ),
    predped::goal(
        id = "table 3", 
        position = c(10, 5) + c(0.55, 0)/2,
        counter = 10
    ),
    predped::goal(
        id = "table 4", 
        position = c(10, 5) + c(-0.55, 0)/2,
        counter = 10
    ),
    predped::goal(
        id = "tablet 1", 
        position = c(2.1, 1.1) + c(0.15, 0.15)/2,
        counter = 10
    ),
    predped::goal(
        id = "tablet 2", 
        position = c(2.1, 5) + c(0.15, 0)/2,
        counter = 10
    ),
    predped::goal(
        id = "tablet 3", 
        position = c(2.1, 8.9) + c(0.15, -0.15)/2,
        counter = 10
    ),
    predped::goal(
        id = "tablet 4", 
        position = c(7.4, 8.9) + c(0, -0.15)/2,
        counter = 10
    ),
    predped::goal(
        id = "tablet 5", 
        position = c(12.6, 8.9) + c(0, -0.15)/2,
        counter = 10
    ),
    predped::goal(
        id = "tablet 6", 
        position = c(17.9, 8.9) + c(-0.15, -0.15)/2,
        counter = 10
    ),
    predped::goal(
        id = "tablet 7", 
        position = c(17.9, 5) + c(-0.15, 0)/2,
        counter = 10
    ),
    predped::goal(
        id = "tablet 8", 
        position = c(17.9, 1.1) + c(-0.15, 0.15)/2,
        counter = 10
    ),
    predped::goal(
        id = "tablet 9", 
        position = c(7.4, 1.1) + c(0, 0.15)/2,
        counter = 10
    ),
    predped::goal(
        id = "tablet 10", 
        position = c(12.6, 1.1) + c(0, 0.15)/2,
        counter = 10
    )
)

# Visualize the background together with the goals within it.
goals_xy <- t(sapply(goals, predped::position))
predped::plot(setting, dark_mode = TRUE) +
    ggplot2::annotate(
        "point", 
        x = goals_xy[, 1], 
        y = goals_xy[, 2],
        color = 'cornflowerblue'
    )


# Create a function that will handle the goals. At each time, we need to see 
# which tablets are free and which ones aren't, updating a list containing free
# an non-free goals. 
#
# As per the requirements of predped, this function takes in a state and returns 
# another state.
update_goals <- function(state) {
    # Check at which iteration we currently are. If at the first iteration, we 
    # need to change the goals of the agents and create some variables to 
    # account for.
    if(predped::iteration(state) == 0) {
        # Create some variables that will help us with goal handling
        goal_id <- sapply(goals, predped::id)
        free_goals <- goal_id
        goal_agent <- c()

        # Loop over all agents and assign them their goals
        agent_list <- predped::agents(state)
        for(i in seq_along(agent_list)) {
            # Remove all previously assigned goals
            predped::goals(agent_list[[i]]) <- list()

            # Change the current goal and make the agent plan their route
            selected <- sample(free_goals, 1)
            predped::current_goal(agent_list[[i]]) <- goals[goal_id == selected][[1]]
            predped::status(agent_list[[i]]) <- "plan"

            # Update the free_goals list 
            free_goals <- free_goals[free_goals != selected]

            # Update the goal_agent list
            goal_agent[predped::id(agent_list[[i]])] <- selected
        }

    # If at any other iteration, we need to check whether any of the agents 
    # currently has an exit goal. If so, we know that we need to assign them 
    # a new goal.
    } else {
        # Extract the variables of interest
        vars <- predped::variables(state)

        goal_id <- vars[["goal_id"]]
        free_goals <- vars[["free_goals"]]
        goal_agent <- vars[["goal_agent"]]

        # Find out whether any agents have an exit goal at this moment
        agent_list <- predped::agents(state)

        exiting <- sapply(agent_list, \(x) predped::id(predped::current_goal(x)))
        exiting <- exiting == "goal exit"

        # If no-one is exiting, then we don't need to do anything. If someone is 
        # exiting, however, then we need to provide them with a new goal drawn
        # from the goal list
        if(any(exiting)) {
            # Loop over all agents
            for(i in seq_along(agent_list)) {
                if(!exiting[i]) {
                    next
                }

                # Sample one of the free goals
                selected <- sample(free_goals, 1)
                predped::current_goal(agent_list[[i]]) <- goals[goal_id == selected][[1]]
                predped::status(agent_list[[i]]) <- "plan"

                # Update the free_goals list: Delete the new goal and add the old
                # one
                free_goals <- c(
                    free_goals[free_goals != selected],
                    goal_agent[predped::id(agent_list[[i]])]
                )

                # Update the goal_agent list
                goal_agent[predped::id(agent_list[[i]])] <- selected
            }
        }
    }

    # Update the agent list
    predped::agents(state) <- agent_list

    # Save the variables in the slot of the state
    predped::variables(state)[["goal_id"]] <- goal_id 
    predped::variables(state)[["free_goals"]] <- free_goals 
    predped::variables(state)[["goal_agent"]] <- goal_agent

    return(state)
}

# Simulate synthetic data. Some things to note here: 
#   - Number of iterations is set to have about 15 minutes of simulation, i.e.
#     to 15 * 60 * 2 = 1800 The is the amount of data that we will have per
#     participant 
#   - Number of participants is set to 10, with no way to get more participants
#   - Number of goals is set to an impossible amount to reach
#
# We run these simulations 100 times, allowing us to estimate the variability 
# around the efficacy of a given preprocessing pipeline. Note that these 
# simulations may take a while
set.seed(74327) # Knowing That You've Arrived - Tides of Man
seeds <- sample(1:10000, 100)

data <- lapply(
    seq_along(seeds), 
    function(i) {
        print(i)
        set.seed(seeds[i])
        
        # Create a trace using predped
        trace <- predped::simulate(
            model,
            initial_number_agents = 10,
            max_agents = 10,
            iterations = 1800,
            fx = update_goals
        )
        saveRDS(
            trace, 
            file.path("data", "study 3", paste0("trace_", i, ".gif"))
        )

        # Save the GIF, allowing you to inspect the actual data
        plt <- predped::plot(
            trace,
            dark_mode = TRUE
        )
        gifski::save_gif(
            lapply(plt, print),
            file.path("data", "study 3", paste0("data_", i, ".gif")),
            delay = 1/15
        )

        # Transform to a dataset
        data <- predped::time_series(trace) %>% 
            dplyr::select(time, id, x, y) %>% 
            dplyr::mutate(
                nsim = i,
                x = as.numeric(x),
                y = as.numeric(y),
                time = as.numeric(time)
            )

        return(data)
    }
)
data <- do.call("rbind", data)

# Save these data
data.table::fwrite(
    data, 
    file.path("data", "study 3", "data.csv")
)





################################################################################
# ADDING NOISE

# An important thing to note here: We used a delay 500msec, while the temporal
# error has been estimated on smaller time-scales. We know that the parameters
# of the VAR depend on the time interval, but one could argue that at smaller 
# time-scales, we can expect error to be more temporally correlated than at 
# larger time-scales, meaning this test is more stringent.
#
# Furthermore note that we use the same type of error as for the first 
# simulation. Therefore less well documented in this version of the code.



# Overall error ################################################################

# Read in the data
data <- data.table::fread(
    file.path("data", "study 3", "data.csv"),
    data.table = FALSE
)

# Transform the dataset to include random error. 
covariances <- do.call(
    "rbind",
    readRDS(file.path("results", "study 1", "unsystematic error, overall covariance.Rds"))
)
max_var <- max(covariances$ub[covariances$covariance != "cov_xy"])
max_cov <- max(covariances$ub[covariances$covariance == "cov_xy"])
S <- matrix(
    c(max_var, max_cov, max_cov, max_var),
    nrow = 2,
    ncol = 2
)

# Create the residuals and add them to the data
set.seed(6622537) # Sesame Street is no Place for Romance - A Lot Like Birds
residuals <- MASS::mvrnorm(
    nrow(data), 
    c(0, 0), 
    S
)

data <- data %>% 
    dplyr::rename(
        x_actual = x, 
        y_actual = y
    ) %>% 
    dplyr::mutate(
        x = x_actual + residuals[,1], 
        y = y_actual + residuals[,2]
    )

# Save the data
data.table::fwrite(
    data, 
    file.path("data", "study 3", "data_R10.csv")
)



# Temporal error ###############################################################

# Read in the data
data <- data.table::fread(
    file.path("data", "study 3", "data.csv"),
    data.table = FALSE
)

# Here, we will need to create a vector autoregressive model that will account 
# both for contemporaneous and lagged measurement error. 
params <- readRDS(file.path("results", "study 1", "unsystematic error, autoregression parameters.Rds"))
B <- matrix(
    params[3:6, 3],
    nrow = 2, 
    ncol = 2
)
S <- matrix(
    c(
        max(params[c(7, 10), 3]), 
        max(params[8:9, 3]), 
        max(params[8:9, 3]), 
        max(params[c(7, 10), 3])
    ),
    nrow = 2, 
    ncol = 2
)

# Start by creating the function that will create the measurement error to be 
# added to the observations.
add_residuals <- function(x) {
    # Create a matrix that will contain the residuals and add the initial 
    # condition
    residuals <- MASS::mvrnorm(
        nrow(x), 
        c(0, 0), 
        S
    )

    y <- matrix(0, nrow = nrow(x), ncol = 2)
    y[1,] <- residuals[1,]

    # Loop over and create the other residuals
    for(i in seq_len(nrow(x) - 1)) {
        y[i + 1,] <- y[i,] %*% B + residuals[i + 1,]
    }

    # Add the residuals to the dataframe and return
    x %>% 
        dplyr::arrange(time) %>% 
        dplyr::mutate(
            x = x_actual + residuals[, 1], 
            y = y_actual + residuals[, 2]
        ) %>% 
        return()
}

# Add the error to the data
set.seed(64710) # Sharon Tate, Despite Everything - The Sound of Animals Fighting
data <- data %>% 
    dplyr::rename(
        x_actual = x, 
        y_actual = y
    ) %>% 
    dplyr::group_by(nsim, id) %>% 
    tidyr::nest() %>% 
    dplyr::mutate(
        data = data %>% 
            as.data.frame() %>% 
            add_residuals() %>% 
            list()
    ) %>% 
    tidyr::unnest(data) %>% 
    dplyr::arrange(nsim, time, id)

# Save the data
data.table::fwrite(
    data, 
    file.path("data", "study 3", "data_T10.csv")
)



# Missing data #################################################################

# Create missing data that gets rid of 30% of the data, just as we did for the 
# data in Study 2.
#
# Create the two functions needed for this
random_missing <- function(x) {
    idx <- sample(
        seq_len(nrow(x)), 
        round(0.7 * nrow(x))
    )

    x[idx, c("x", "y")] <- NA 
    return(x)
}

nonrandom_missing <- function(x) {
    N <- nrow(x)

    # Create relative indices per person for the blocked missing data
    idx <- data.frame(
        from = sample(
            1:(length(unique(x$time)) - 15), 
            round(0.3 * N / 10), 
            replace = TRUE
        ), 
        participants = sample(
            unique(x$id), 
            round(0.3 * N / 10), 
            replace = TRUE
        )
    ) %>%  
        dplyr::mutate(
            from = as.numeric(from),
            to = from + rep(c(5, 10, 15) - 1, each = round(length(from) / 3))
        ) 

    # Explicate all indices to be deleted
    idx <- idx %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(
            indices = seq(from, to) %>% 
                as.vector() %>% 
                data.frame() %>% 
                setNames("indices") %>% 
                tidyr::nest()
        ) %>%
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
        dplyr::mutate(
            data = data %>%
                as.data.frame() %>%
                dplyr::mutate(
                    number = dplyr::row_number(),
                    index = number %in% idx$indices[idx$participants == id]
                ) %>%
                dplyr::filter(!index) %>% 
                dplyr::select(-index, -number) %>%
                list()
        ) %>% 
        tidyr::unnest(data) %>% 
        dplyr::ungroup()

    # Now sample the remaining time points to be deleted from the remaining data 
    # points
    x <- x %>% 
        random_missing()

    return(x)
}

# Define the files for which to impute the missing data, loop over them and 
# impose
filenames <- paste0("data", c("_R", "_T"))

set.seed(443) # Pale Blue Dot - The Receiving End of Sirens
for(i in seq_along(filenames)) {
    data <- data.table::fread(
        file.path("data", "study 3", paste0(filenames[i], "10.csv")),
        data.table = FALSE
    )
    
    # Random missing
    tmp <- data %>% 
        dplyr::group_by(nsim) %>% 
        random_missing() %>% 
        dplyr::filter(!is.na(x))

    data.table::fwrite(
        tmp, 
        file.path("data", "study 3", paste0(filenames[i], "6R.csv"))
    )

    # Nonrandom missing
    tmp <- data %>% 
        dplyr::group_by(nsim) %>% 
        tidyr::nest() %>%
        dplyr::mutate(data = data %>%  
            as.data.frame() %>%
            nonrandom_missing() %>%
            list()) %>%
        tidyr::unnest(data) %>% 
        dplyr::filter(!is.na(x))

    data.table::fwrite(
        tmp, 
        file.path("data", "study 3", paste0(filenames[i], "6N.csv"))
    )
}