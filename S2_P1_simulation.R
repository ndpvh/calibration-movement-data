################################################################################
# Purpose: Generate some simplistic synthetic data that can be used to test    # 
#          some strategies of filtering the unsystematic error from the data.  #
#          In this first step, we generate data according to the shape of      #
#          some geometric functions, namely a points -- representing           #
#          stationary data -- and either a circle, a square, or a spiral --    #
#          representing for movement data.                                     #
#                                                                              #
#          Structure of the file is as follows:                                #
#              L21-177:    Creating the functions necessary for simulation     #
#              L183-295:   Simulating positions                                #
#              L301-580:   Adding noise                                        #
#              L587-717:   Visualization                                       #
#                                                                              #
#          Note that when adding noise, we add both time-independent noise     #
#          and a time-dependent noise. We use the parameters estimated in the  #
#          S1_P1_analysis script for this purpose.                             #
################################################################################

devtools::load_all()

################################################################################
# FUNCTIONS

# Create a bunch of utility functions that will help generate the data, mostly 
# for the movement data. Each function will take in a given number of data points
# to generate (N) and a partition number that defines how many partitions 
# should be taken of the circumference of the geometric shape (p)

# For circular movement, smooth and counterclockwise
circle <- function(N, p) {
    # Define the parameters of the circle
    radius <- 1
    center <- c(0, 0)

    # Get the movement angle of the pedestrian after each time step and compute 
    # the vector of all angles the pedestrian will be at 
    angle <- 2 * pi / p
    angle <- seq(1, N, 1) * angle 

    # Ceate the (x, y) coordinates using these angles
    locations <- data.frame(
        x = center[1] + radius * cos(angle), 
        y = center[2] + radius * sin(angle)
    )

    return(locations)
}

# For rectangular movement, abrupt and counterclockwise
rectangle <- function(N, p) {
    # Define the parameters of the rectangle
    width <- 1 
    height <- 1
    center <- c(0, 0)

    # Compute the circumference of the rectangle and define the distance that 
    # the pedestrians will walk on along this circumference. Furthermore define
    # the distance the pedestrian has travelled at each time point.
    circumference <- 2 * width + 2 * height
    d <- circumference / p
    d <- seq(1, N, 1) * d

    # Define the sides on which the pedestrian walks, starting from one of the 
    # corners of the rectangle
    sides <- data.frame(
        x = center[1] + c(-1, 1, 1, -1) * width / 2,
        y = center[1] + c(-1, -1, 1, 1) * height / 2,
        side = 1:4,
        adjustment_x = c(1, 0, -1, 0),
        adjustment_y = c(0, 1, 0, -1)
    )

    # Now for the actual movement, for which we need a dirty trick: 
    #   - Define the distances that are bigger than the circumference and 
    #     divide them up in how many times the pedestrian has already made a 
    #     complete loop around the rectangle. Correct these distances for the 
    #     size of the circumference (i.e., the distances for each loop start 
    #     at 0)
    #   - Assign the distances to a given side and correct these so that they 
    #     start at 0 for each side
    #   - Add the distances to the respective dimension that they should be 
    #     assigned to
    locations <- data.frame(
        d = d, 
        loop = ceiling(d / circumference)
    ) %>% 
        # Correct for circumference
        dplyr::mutate(d = d - (loop - 1) * circumference) %>% 
        # Assign sides and correct for sides. Assumption: You start in the lower 
        # left corner
        dplyr::mutate(
            side = 1 + (d >= 2 * width + height) + (d >= width + height) + (d >= width), 
            w = ifelse(side == 1, 0, ifelse(side == 4, 2, 1)), 
            h = ifelse(side >= 3, 1, 0), 
            d = d - width * w - height * h
        ) %>% 
        # Use the information on the sides to define when to add and to delete 
        # a distance in a given dimension
        plyr::join(sides, by = "side") %>% 
        dplyr::mutate(
            d_x = d * adjustment_x, 
            d_y = d * adjustment_y
        ) %>% 
        dplyr::mutate(
            x = ifelse(adjustment_x != 0, x + d_x, x), 
            y = ifelse(adjustment_y != 0, y + d_y, y)
        ) %>% 
        # Select only the (x, y) coordinates
        dplyr::select(x, y)
                      
    return(locations)
}

# For spiral movement, smooth, outward, and counterclockwise
spiral <- function(N, p) {
    # Define the parameters to be used for the spiral
    velocity <- 1
    distance <- sin(2 * pi / p) / sin(pi / 2 - pi / p)
    radius <- distance

    # Find the first angle that satisfies the distance d to the second point 
    # (starting at the center)
    alpha <- distance / velocity

    # Create a locations matrix and loop over all data points (starting at the 
    # origin makes sure we can start at index 2)
    locations <- matrix(0, nrow = N, ncol = 2) 
    for(i in 2:N) {
        # Approximate the points with a circle with radius r
        locations[i,] <- radius * c(cos(alpha), sin(alpha))

        # Adjust the angle and radius of the approximating circle using the 
        # distance to cross
        alpha <- alpha + distance / radius 
        radius <- velocity * alpha
    }

    # Additional check: Due to approximations, it might be that the distance 
    # between two points is not equal to the distance to cross. Try to flag this
    # and generate additional data points when this is found
    locations <- locations %>% 
        # Compute the distances between each of the points
        as.data.frame() %>% 
        setNames(c("x", "y")) %>% 
        dplyr::mutate(distance = c(NA, sqrt((x[2:N] - x[2:N - 1])^2 + (y[2:N] - y[2:N - 1])^2))) %>% 
        # Check whether any of these is an outlier
        dplyr::mutate(
            outlier = distance > mean(distance, na.rm = TRUE) + 4 * sd(distance, na.rm = TRUE), 
            outlier = ifelse(is.na(outlier), FALSE, outlier)
        )

    # Delete the data points that fall before the outlier and add some additional
    # points at the end of the spiral (no outliers should be observed there)
    if(any(locations$outlier)) {
        # Delete
        idx <- c(1:N)[locations$outlier]
        locations <- locations[-c(1:idx), c("x", "y")]

        # Replace
        for(i in 1:idx) {
            # Approximate the points with a circle with radius r
            locations <- rbind(
                locations, 
                radius * c(cos(alpha), sin(alpha))
            )

            # Adjust the parameters
            alpha <- alpha + distance / radius 
            radius <- velocity * alpha
        }
    } else {
        locations <- locations %>% 
            dplyr::select(x, y)
    }

    return(locations)
}





################################################################################
# SIMULATION

# Now that the functions are ready, actually simulate some data with them. 
# Some things to note: 
#   - For each data set, about 100 data points simulated
#   - For the movement data, partitions `p` are taken to be {20, 40, 80}
#   - For the fixed data, just coordinates placed on 9 locations within a 2 x 2 
#     square and repeated for the duration of the experiment
#   - The time variable is created such that the data was sampled at a 10Hz rate
#     (i.e., data were collected for about 10sec)
#   - We repeat these same data 100 times, allowing us to aggregate across 
#     different instantiations of simulated measurement error (allowing greater
#     generality)

# Movement #####################################################################

N <- 100
fx <- list(
    "circle" = \(x) circle(N, x), 
    "rectangle" = \(x) rectangle(N, x), 
    "spiral" = \(x) spiral(N, x)
)
fx_names <- names(fx)

set.seed(33211) # The "Fun" in Dysfunction - Hail the Sun
data <- lapply(
    fx_names,
    function(name) {
        tmp <- lapply(
            c(20, 40, 80),
            function(p) {
                # Generate the data
                tmp <- fx[[name]](p) %>% 
                    dplyr::mutate(
                        time = (dplyr::row_number() - 1) / 10,
                        id = paste0(name, "_", p)
                    )

                # Rescale the data so that the distance between each point is 
                # about 0.14, which agrees with the average speed of 1.4 m/s of 
                # a human. This is easily achieved by scaling the (x, y) space 
                # with a constant factor 0.14 / current distance between points.
                tmp <- tmp %>% 
                    dplyr::mutate(
                        distance = c(
                            NA, 
                            sqrt((x[2:N] - x[2:N - 1])^2 + (y[2:N] - y[2:N - 1])^2)
                        ),
                        x = 0.14 * x / mean(distance, na.rm = TRUE), 
                        y = 0.14 * y / mean(distance, na.rm = TRUE)
                    ) %>% 
                    dplyr::select(x, y, time, id)

                return(tmp)
            }
        )

        return(do.call("rbind", tmp))
    }
)
data <- do.call("rbind", data)

# Replicate these same data 100 times. This makes sure that the randomness
# of the added error has no part in the results 
data <- do.call(
    "rbind", 
    replicate(
        100, 
        data, 
        simplify = FALSE
    )
) %>% 
    dplyr::mutate(nsim = rep(1:100, each = nrow(data)))

# Save the results
data.table::fwrite(
    data, 
    file.path("data", "study 2", "movement.csv")
)



# Stationary ###################################################################

# Create stationary data. This is a bit easier, as its just a repeated point
data <- data.frame(
    x = rep(-1:1, each = N * 3), 
    y = rep(
        rep(-1:1, each = N), 
        times = 3
    )
) %>% 
    dplyr::mutate(id = rep(paste0("fixed_", 1:9), each = N)) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(time = (dplyr::row_number() - 1) / 7) %>% 
    dplyr::ungroup()

# Same repeat will be used here
data <- do.call(
    "rbind", 
    replicate(
        100, 
        data, 
        simplify = FALSE
    )
) %>% 
    dplyr::mutate(nsim = rep(1:100, each = nrow(data)))

data.table::fwrite(
    data, 
    file.path("data", "study 2", "fixed.csv")
)





################################################################################
# ADDING NOISE

# Overall error ################################################################

# Transform the dataset to include random error. This error will be equal to 
# a few centimeters, based on earlier estimates of the measurement error, 
# specifically coming from the calibration test of 22-12-2023, where the upper 
# bound of the 99%CI was about 6cm. 
#
# We use the maximal variance for both x- and y-coordinates, which is equal to 
# 0.004468. We also take the covariance into account, again selecting the 
# maximal one. This one is equal to 0.0280662.
filenames <- c("movement", "fixed")

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

set.seed(425) # Call it Karma - Silverstein
for(i in seq_along(filenames)) {
    # Read in the data
    data <- data.table::fread(
        file.path("data", "study 2", paste0(filenames[i], ".csv")),
        data.table = FALSE
    )

    # Create the residuals and add them to the data
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
        file.path("data", "study 2", paste0(filenames[i], "_R10.csv"))
    )
}



# Temporal error ###############################################################

# Here, we will need to create a vector autoregressive model that will account 
# both for contemporaneous and lagged measurement error. The parameters that 
# are used here are taken from estimations we did on the stationary calibration 
# data. The measurement error is added for each experiment and id separately.
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
set.seed(55) # Sober Exit(s) - Static Dress
for(i in seq_along(filenames)) {
    # Add normally distributed measurement error to the data
    data <- data.table::fread(
        file.path("data", "study 2", paste0(filenames[i], ".csv")),
        data.table = FALSE
    ) %>% 
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
        file.path("data", "study 2", paste0(filenames[i], "_T10.csv"))
    )
}



# Missing data #################################################################

# Our data often does not have an actual sampling frequency of 10Hz, but it is 
# ususally lower. Account for this kind of missing data, and the unequal time 
# intervals that this creates in the data (and which may make our job 
# difficult). We will retain a sampling frequency of 7Hz.
#
# Two types: Random missingness and nonrandom missingness. Create two functions
# that will impute these missings.
random_missing <- function(x, 
                           N = round(0.3 * nrow(x))) {
    idx <- sample(
        seq_len(nrow(x)), 
        N
    )

    x[idx, c("x", "y")] <- NA 
    return(x)
}

# We draw random time points at which no position is measured. Once chosen, we 
# then either keep it at this one time point, or we make the period at which no 
# measurements are obtained longer (5 observations, 10 observations, or 15 
# observations, it being 500msec, 1sec, or 1.5sec long). We try to approximate 
# each time as being as long as the other.
#
# In practice, we delete the indices in two waves: One in which we take care of 
# the longer problems, then one in which we delete only a single data point 
# randomly up until as there are as many missing observations for the nonrandom 
# and the random missing data (created earlier)
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
        random_missing(N = round(0.3 * N) - (N - nrow(x)))

    return(x)
}

# Define the files for which to impute the missing data, loop over them and 
# impose
filenames <- paste0(
    rep(c("fixed", "movement"), each = 2),
    rep(c("_R", "_T"), times = 2)
)

set.seed(9410) # Newsstand Rock (exposition) - Rx Bandits
for(i in seq_along(filenames)) {
    data <- data.table::fread(
        file.path("data", "study 2", paste0(filenames[i], "10.csv")),
        data.table = FALSE
    )
    
    # Random missing
    tmp <- data %>% 
        dplyr::group_by(nsim) %>% 
        tidyr::nest() %>%
        dplyr::mutate(data = data %>%
            as.data.frame() %>%
            random_missing() %>%
            list()) %>%
        tidyr::unnest(data) %>%
        dplyr::filter(!is.na(x))

    data.table::fwrite(
        tmp, 
        file.path("data", "study 2", paste0(filenames[i], "6R.csv"))
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
        file.path("data", "study 2", paste0(filenames[i], "6N.csv"))
    )
}





################################################################################
# VISUALIZATION

# Visualize the types of data that we are looking at here. Only make the 
# distinction between fixed and movement, and between the different types of 
# error.
fixed <- rbind(
    data.table::fread(
        file.path("data", "study 2", "fixed_R10.csv"),
        data.table = FALSE
    ) %>% 
        dplyr::mutate(type = "R"),
    data.table::fread(
        file.path("data", "study 2", "fixed_T10.csv"),
        data.table = FALSE
    ) %>% 
        dplyr::mutate(type = "T")
)
movement <- rbind(
    data.table::fread(
        file.path("data", "study 2", "movement_R10.csv"),
        data.table = FALSE
    ) %>% 
        dplyr::mutate(type = "R"),
    data.table::fread(
        file.path("data", "study 2", "movement_T10.csv"),
        data.table = FALSE
    ) %>% 
        dplyr::mutate(type = "T")
)

# Select only a few of the ids to showcase
data <- rbind(fixed, movement) %>% 
    dplyr::filter(id %in% c("fixed_1", "circle_40", "rectangle_40", "spiral_40")) %>% 
    dplyr::filter(nsim == 1) 

# Make the plots themselves
plots <- lapply(
    c("actual", "R", "T"),
    function(x) {
        # Select the data of interest, based on the type of the error. For the 
        # "actual" data, we need to do something special
        if(x == "actual") {
            selected_data <- data[data$type == "R", ] %>% 
                dplyr::mutate(
                    x = x_actual, 
                    y = y_actual
                ) %>% 
                dplyr::select(x, y, id)
        } else {
            selected_data <- data[data$type == x, ] %>% 
                dplyr::select(x, y, id)
        }

        # Loop over all of the id's and plot their positions
        plt <- lapply(
            unique(selected_data$id),
            function(y) {
                plot_data <- selected_data[selected_data$id == y, ]

                title <- list(
                    "fixed_1" = "Point",
                    "circle_40" = "Circle",
                    "rectangle_40" = "Rectangle",
                    "spiral_40" = "Spiral"
                )
                limits <- data[data$id == y, ] %>% 
                    dplyr::select(x, y) %>% 
                    unlist() %>% 
                    as.numeric() %>% 
                    range()

                return(
                    ggplot2::ggplot(plot_data, 
                                    ggplot2::aes(x = x, 
                                                 y = y)) +
                        ggplot2::geom_point(size = 3, 
                                            color = "cornflowerblue",
                                            shape = 19) +
                        ggplot2::labs(x = "x",
                                      y = "y",
                                      title = ifelse(x == "actual", 
                                                     title[[y]], 
                                                     " ")) +
                        ggplot2::lims(x = limits, 
                                      y = limits) +
                        ggplot2::theme_minimal() +
                        ggplot2::theme(plot.title = ggplot2::element_text(size = 32,
                                                                          hjust = 0.5),
                                       axis.title = ggplot2::element_text(size = 25),
                                       axis.text = ggplot2::element_text(size = 15), 
                                       legend.title = ggplot2::element_text(size = 25),
                                       legend.text = ggplot2::element_text(size = 15),
                                       panel.background = ggplot2::element_rect(fill = NA, 
                                                                                linewidth = 1.5))
                )
            }
        )

        return(
            ggpubr::ggarrange(
                plotlist = plt,
                nrow = 1
            )
        )
    }
)

name <- ggpubr::ggarrange(
    plotlist = lapply(
        c("Actual positions", "Random error", "Temporal error"),
        \(x) nameless::name_plot(x, size = 10)
    ),
    ncol = 1
)

ggplot2::ggsave(
    file.path("figures", "study 2", "data.png"),
    ggpubr::ggarrange(
        name,
        ggpubr::ggarrange(
            plotlist = plots,
            nrow = 3
        ),
        ncol = 2,
        widths = c(0.15, 0.75)
    ),
    width = 5600,
    height = 3800,
    unit = "px"
)



