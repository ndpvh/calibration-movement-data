################################################################################
# Purpose: Generate some simplistic synthetic data that can be used to test    # 
#          a lot of preprocessing strategies. For this, we will generate data  #
#          according to some geometric functions, namely a circle, a square,   # 
#          and a spiral (for movement data) and points (for fixed data).       #
#                                                                              #
#          RQ: Which method performs best in filtering the data and ridding    #
#              it of unsystematic error.                                       #
#                                                                              #
#          This is Part 1.1 of Simulation 1, meaning that we only generate the #
#          "real" positions of the pedestrians at each time point. Measurement #
#          error is added in tandem with Simulation 2 in Part 1.2.             #
################################################################################                                     #
################################################################################

devtools::load_all()

################################################################################
# PRELIMINARIES

# Create a bunch of utility functions that will help generate the data, mostly 
# for the movement data. Each function will take in a given number of data points
# to generate (N) and a partition number that defines how many partitions 
# should be taken of the circumference of the geometric shape (p)
circle <- function(N, p) {
    # Define the parameters of the circle
    r <- 1                  # Radius 
    c <- c(0, 0)            # Center

    # Get the movement angle of the pedestrian after each time step and compute 
    # the vector of all angles the pedestrian will be at 
    angle <- 2 * pi / p     # Angle of movement for each time step
    angle <- seq(1, N, 1) * angle 

    # Ceate the (x, y) coordinates using these angles
    locations <- cbind(x = c[1] + r * cos(angle), 
                       y = c[2] + r * sin(angle)) %>% 
        as.data.frame() %>% 
        setNames(c("x", "y"))

    return(locations)
}

rectangle <- function(N, p) {
    # Define the parameters of the rectangle
    w <- 1                  # Width
    h <- 1                  # Height
    c <- c(0, 0)            # Center

    # Compute the circumference of the rectangle and define the distance that 
    # the pedestrians will walk on along this circumference. Furthermore define
    # the distance the pedestrian has travelled at each time point.
    circumference <- 2 * w + 2 * h 
    d <- circumference / p
    d <- seq(1, N, 1) * d

    # Define the sides on which the pedestrian walks, starting from one of the 
    # corners of the rectangle
    sides <- rbind(c(c[1] - w / 2, c[2] - h / 2, 1, 1, 0), 
                   c(c[1] + w / 2, c[2] - h / 2, 2, 0, 1), 
                   c(c[1] + w / 2, c[2] + h / 2, 3, -1, 0),
                   c(c[1] - w / 2, c[2] + h / 2, 4, 0, -1)) %>% 
        as.data.frame() %>% 
        setNames(c("x", "y", "side", "adjustment_x", "adjustment_y"))

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
    locations <- cbind(d = d, 
                       loop = ceiling(d / circumference)) %>% 
        as.data.frame() %>% 
        # Correct for circumference
        dplyr::mutate(d = d - (loop - 1) * circumference) %>% 
        # Assign sides and correct for sides. Assumption: You start in the lower 
        # left corner
        dplyr::mutate(side = 1 + (d >= 2 * w + h) + (d >= w + h) + (d >= w), 
                      width = ifelse(side == 1, 0, ifelse(side == 4, 2, 1)), 
                      height = ifelse(side >= 3, 1, 0), 
                      d = d - width * w - height * h) %>% 
        # Use the information on the sides to define when to add and to delete 
        # a distance in a given dimension
        plyr::join(sides, by = "side") %>% 
        dplyr::mutate(d_x = d * adjustment_x, 
                      d_y = d * adjustment_y) %>% 
        dplyr::mutate(x = ifelse(adjustment_x != 0, x + d_x, x), 
                      y = ifelse(adjustment_y != 0, y + d_y, y)) %>% 
        # Select only the (x, y) coordinates
        dplyr::select(x, y)
                      
    return(locations)
}

spiral <- function(N, p) {
    # Define the parameters to be used for the spiral
    v <- 1                                        # Velocity
    d <- sin(2 * pi / p) / sin(pi / 2 - pi / p)   # Distance to next point
    r <- d                                        # Radius of the approximating circle

    # Find the first angle that satisfies the distance d to the second point 
    # (starting at the center)
    alpha <- d / v

    # Create a locations matrix and loop over all data points (starting at the 
    # origin makes sure we can start at index 2)
    locations <- matrix(0, nrow = N, ncol = 2) 
    for(i in 2:N) {
        # Approximate the points with a circle with radius r
        locations[i,] <- c(r * cos(alpha), r * sin(alpha))

        # Adjust the angle and radius of the approximating circle using the 
        # distance to cross
        alpha <- alpha + d / r 
        r <- v * alpha
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
        dplyr::mutate(outlier = distance > mean(distance, na.rm = TRUE) + 4 * sd(distance, na.rm = TRUE), 
                      outlier = ifelse(is.na(outlier), FALSE, outlier))

    # Delete the data points that fall before the outlier and add some additional
    # points at the end of the spiral (no outliers should be observed there)
    if(any(locations$outlier)) {
        # Delete
        idx <- c(1:N)[locations$outlier]
        locations <- locations[-c(1:idx), c("x", "y")]

        # Replace
        for(i in 1:idx) {
            # Approximate the points with a circle with radius r
            locations <- rbind(locations, 
                               c(r * cos(alpha), r * sin(alpha)))

            # Adjust the parameters
            alpha <- alpha + d / r 
            r <- v * alpha
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
#     (i.e., data were collected for 10sec)
#   - We repeat these same data 100 times, allowing us to aggregate across 
#     different instantiations of simulated measurement error (allowing greater
#     generality)
#
# Movement data first
N <- 100
fx <- list("circle" = \(x) circle(N, x), 
           "rectangle" = \(x) rectangle(N, x), 
           "spiral" = \(x) spiral(N, x))
fx_names <- names(fx)

data <- NULL
set.seed(33211) # The "Fun" in Dysfunction - Hail the Sun
for(i in seq_along(fx)) {
    for(j in c(20, 40, 80)) {
        # Generate the data
        tmp <- fx[[i]](j) %>% 
            dplyr::mutate(time = (dplyr::row_number() - 1) / 10, 
                          id = paste0(fx_names[i], "_", j))

        # Rescale the data so that the distance between each point is about 
        # 0.14: Agrees with the average speed of 1.4 m/s of a human, furthermore 
        # allowing us to straightforwardly plug in empirical measures of the 
        # measurement error.
        #
        # This is easily achieved by scaling the (x, y) space with a constant 
        # factor 0.14 / current distance between points
        tmp <- tmp %>% 
            dplyr::mutate(distance = c(NA, sqrt((x[2:N] - x[2:N - 1])^2 + (y[2:N] - y[2:N - 1])^2))) %>% 
            dplyr::mutate(x = 0.14 * x / mean(distance, na.rm = TRUE), 
                          y = 0.14 * y / mean(distance, na.rm = TRUE)) %>% 
            dplyr::select(x, y, time, id)

        # Bind it together with the data, if it already exists
        if(is.null(data)) {
            data <- tmp
        } else {
            data <- rbind(data, tmp)
        }
    }
}

# We will use these same data 100 times. This is to make sure that the randomness
# of the added error has no part in the results 
row.names(data) <- NULL
data <- do.call("rbind", replicate(100, data, simplify = FALSE)) %>% 
    dplyr::mutate(nsim = rep(1:100, each = nrow(data)))

data.table::fwrite(data, file.path("data", "simulation_1", "movement.csv"))

# Now for the fixed data
data <- cbind(rep(-1:1, each = N * 3), 
              rep(rep(-1:1, each = N), times = 3)) %>% 
    as.data.frame() %>% 
    setNames(c("x", "y")) %>% 
    # Add id
    dplyr::mutate(id = rep(paste0("fixed_", 1:9), each = N)) %>% 
    # Add timing
    dplyr::group_by(id) %>% 
    dplyr::mutate(time = (dplyr::row_number() - 1) / 10) %>% 
    dplyr::ungroup()

# Same repeat will be used here
row.names(data) <- NULL
data <- do.call("rbind", replicate(100, data, simplify = FALSE)) %>% 
    dplyr::mutate(nsim = rep(1:100, each = nrow(data)))

data.table::fwrite(data, file.path("data", "simulation_1", "fixed.csv"))
