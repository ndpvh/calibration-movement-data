################################################################################
# Purpose: Generate some simplistic synthetic data that can be used to test    # 
#          some strategies of filtering the unsystematic error from the data.  #
#          In this first step, we generate data according to the shape of      #
#          some geometric functions, namely a points -- representing           #
#          stationary data -- and either a circle, a square, or a spiral --    #
#          representing for movement data.                                     #
#                                                                              #
#          Structure of the file is as follows:                                #
#              L21-177: Creating the functions necessary for simulation        #
#              Lx-x: Simulating positions                                      #
#              Lx-x: Adding noise                                              #
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
#   - The time variable is created such that the data was sampled at a 7Hz rate
#     (i.e., data were collected for about 14sec)
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
                        time = (dplyr::row_number() - 1) / 7,
                        id = paste0(name, "_", p)
                    )

                # Rescale the data so that the distance between each point is 
                # about 0.2, which agrees with the average speed of 1.4 m/s of 
                # a human. This is easily achieved by scaling the (x, y) space 
                # with a constant factor 0.14 / current distance between points.
                tmp <- tmp %>% 
                    dplyr::mutate(
                        distance = c(
                            NA, 
                            sqrt((x[2:N] - x[2:N - 1])^2 + (y[2:N] - y[2:N - 1])^2)
                        ),
                        x = 0.2 * x / mean(distance, na.rm = TRUE), 
                        y = 0.2 * y / mean(distance, na.rm = TRUE)
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
