#' Plot moving positions
#' 
#' Plot the pedestrian's movements based on the positions in a dataframe. This 
#' dataframe should consist of several columns, namely (a) a time variable that
#' orders the positions (`time`), (b) a grouping variable that denotes which 
#' positions belong to which pedestrian (`id`), and (c) the x and y positions 
#' themselves (`x` and `y` resp.).
#' 
#' Disclaimer: For this function to work properly, the time-variable should be 
#' an increasing number: It should not be in datatime format.
#' 
#' @param x Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' @param as_points Logical denoting whether to plot the positions as a single 
#' point (`TRUE`) or as a circle with a given orientation (see also `predped`, 
#' `FALSE`). Defaults to `FALSE`.
#' @param per_iteration Logical denoting whether to make multiple plots of 
#' positions at each time point (`TRUE`) or to just plot all measured locations 
#' at once (`FALSE`). If `FALSE`, `as_points` is ignored (will automatically be
#' `TRUE`). Defaults to `TRUE`
#' @param ... Additional arguments to pass to ggplot2.
#' 
#' @return List of plots, one for each time points in the dataframe.
#' 
#' @export
plot <- function(x, 
                 as_points = FALSE, 
                 per_iteration = TRUE,
                 limits = NULL,
                 ...) {

    # Get limits of the plot
    if(is.null(limits)) {
        limits <- cbind(x = range(x$x) + 0.25 * c(-1, 1),
                        y = range(x$y) + 0.25 * c(-1, 1)) %>% 
            as.data.frame()
    }
    
    # Transform the time-variable to integers denoting the order of the 
    # positions
    x <- x %>% 
        dplyr::group_by(id) %>% 
        dplyr::mutate(time = dplyr::min_rank(time)) %>% 
        dplyr::ungroup()

    # Add the orientation of the pedestrians at each iteration. This is done by 
    # taking the angle between a pedestrian's previous position and the current
    # one
    x <- x %>% 
        # Create references to the previous position for each pedestrian
        dplyr::group_by(id) %>% 
        dplyr::mutate(x_prev = c(NA, x[2:length(x) - 1]),
                      y_prev = c(NA, y[2:length(y) - 1])) %>% 
        dplyr::ungroup() %>% 
        # Find the orientation needed to go from the previous position to the 
        # current one
        dplyr::mutate(orientation = atan2(y - y_prev, x - x_prev)) %>% 
        # Delete the references to the prevous positions
        dplyr::select(-x_prev, -y_prev)

    # Make a nested data structure and order the time-variable again (just to 
    # be sure)
    x <- x %>% 
        dplyr::group_by(time) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(time = sort(time))

    # Now that we have this, we should try to visualize the positions at which 
    # the pedestrians were standing at each time point. Depending on which type 
    # of plot the user wanted to have, this will be a bit different.
    if(as_points | !per_iteration) {
        fx <- plot_point
    } else {
        fx <- plot_pedestrian
    } 

    # Check whether to plot per iteration, or whether not to. If plotted per 
    # iteration, we will have to loop over these iterations. If not, then we 
    # just use all locations at once.
    if(per_iteration) {
        cat("\n")

        # Now apply the plots for each iteration
        plt <- list()
        for(i in seq_len(nrow(x))) {            
            cat(paste0("\rMaking plot for time point ", i))

            # Select the data for a given time point
            plot_data <- x %>% 
                dplyr::filter(time == i) %>% 
                tidyr::unnest(data) %>% 
                dplyr::ungroup()

            # Make the plot
            plt[[i]] <- fx(plot_data, ...) +
                ggplot2::labs(title = paste0("time point ", i)) +
                ggplot2::lims(x = limits$x, 
                              y = limits$y)
        }
        cat("\n")

    } else {
        x <- x %>% 
            tidyr::unnest(data) %>% 
            dplyr::ungroup()

        plt <- fx(x, ...) +
            ggplot2::lims(x = limits$x, 
                          y = limits$y)
    }
    
    return(plt)
}

# Plot that will use the position data to plot points at those location in the 
# dataframe.
plot_point <- function(x, ...) {
    View(x)
    plt <- ggplot2::ggplot(data = x, 
                           ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point(...) +
        ggplot2::xlab("x") + 
        ggplot2::ylab("y") +
        ggplot2::theme_minimal() +
        ggplot2::coord_equal()

    return(plt)
}

# Plot that will use the positions and orientations in the data to plot 
# pedestrian-like entities that are walking around. 
plot_pedestrian <- function(x, ...) {
    # For each person denoted by `id`, create a list of points that can be used 
    # to segmentize circles with. Happens in several steps:
    #   - First create all angles at which you would like to have a point
    #   - Then group by the `id`s
    #   - Use the measured positions and the fixed radius 0.25 to create all 
    #     points that you will combine
    #   - Make segments instead of just coordinates by putting some in x, y, and 
    #     xend, yend (done within `id`)
    #   - Use the orientation of the movement to make an additional segment that 
    #     will denote the direction in which the pedestrian is walking. This is 
    #     based on the same principles
    angles <- seq(0, 2 * pi, length.out = 25)

    # Circle points around the center
    y <- x %>% 
        dplyr::group_by(id) %>% 
        dplyr::mutate(data = cbind(x_circ = x + cos(angles) * 0.25, 
                                   y_circ = y + sin(angles) * 0.25) %>% 
                          data.frame() %>% 
                          tidyr::nest()) %>%
        dplyr::ungroup() %>%  
        tidyr::unnest(data) %>% 
        tidyr::unnest(data) 

    # Segmentation of these points
    y <- y %>% 
        dplyr::group_by(id) %>% 
        dplyr::mutate(x = c(x_circ[2:length(x_circ) - 1], NA),
                      xend = c(x_circ[2:length(x_circ)], NA),
                      y = c(y_circ[2:length(y_circ) - 1], NA), 
                      yend = c(y_circ[2:length(y_circ)], NA)) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(-x_circ, -y_circ) %>% 
        dplyr::filter(!is.na(x))

    # Add orientation segments
    y <- x %>% 
        dplyr::mutate(xend = x + cos(orientation) * 0.25, 
                      yend = y + sin(orientation) * 0.25) %>% 
        rbind(y)

    # Make the plot itself
    plt <- ggplot2::ggplot(data = y, 
                           ggplot2::aes(x = x, 
                                        y = y, 
                                        xend = xend, 
                                        yend = yend)) +
        ggplot2::geom_segment(...) +
        ggplot2::xlab("x") + 
        ggplot2::ylab("y") +
        ggplot2::theme_minimal() +
        ggplot2::coord_equal()

    return(plt)
}
