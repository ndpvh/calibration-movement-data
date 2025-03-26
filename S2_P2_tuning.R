################################################################################
# Purpose: Tune and test out some of the building blocks that we will use to   #
#          filter the simulated datasets.                                      #
#                                                                              #
#          Structure of the file is as follows:                                #
#              L21-177:    Reading in the data                                 #
#              L183-295:   Simulating positions                                #
#              L301-580:   Adding noise                                        #
#              L587-717:   Visualization                                       #
################################################################################

devtools::load_all()

################################################################################
# DATA

# Read in all of the datasets to be used at this tuning stage. Put them all in 
# a list that you can loop over when needed.
filenames <- paste0(
    rep(c("movement", "fixed"), each = 6),
    rep(c("_R10", "_T10", "_R6R", "_R6N", "_T6R", "_T6N"), times = 2)
)

data_list <- lapply(
    filenames, 
    \(name) data.table::fread(
        file.path("data", "study 2", paste0(name, ".csv")),
        data.table = FALSE
    ) %>% 
        dplyr::filter(nsim %in% 1:20) %>% 
        dplyr::mutate(
            kind = stringr::str_split(name, "_")[[1]][1],
            type = stringr::str_split(name, "_")[[1]][2]
        )
)
data <- do.call("rbind", data_list)





################################################################################
# FILTER THE DATA

# Filter functions #############################################################

# Define all of the filters to be used. Done in two steps. First, we define the
# Kalman filter and LOESS/LOWESS. Only then do we define the different moving 
# windows.
filters <- list(
    # Use the variance as found in Study 1
    "kalm" = list(
        \(x) nameless::kalman_filter(
            x, 
            reverse = FALSE, 
            .by = "id", 
            assumed_variance = 0.004468
        )
    ),

    # Define LOESSs of different degrees
    "loess_1" = list(
        \(x) nameless::local_regression(
            x,
            .by = "id",
            degree = 1
        )
    ),
    "loess_2" = list(
        \(x) nameless::local_regression(
            x,
            .by = "id",
            degree = 2
        )
    ),
    "loess_3" = list(
        \(x) nameless::local_regression(
            x,
            .by = "id",
            degree = 3
        )
    ),
    "loess_4" = list(
        \(x) nameless::local_regression(
            x,
            .by = "id",
            degree = 4
        )
    )
)

# Define the different statistics to compute/use in the moving windows for 
# filtering. Mostly consists of (weighted) averages and polynomial approximations.
#
# We need to specify that we want to include the additional columns containing 
# the real positions of the data.
fx <- list(
    # Normal average
    "av" = \(x) nameless::average(
        x, 
        cols = c("x_actual", "y_actual")
    ), 

    # Weighted averages based on the index of the observations within the 
    # window and based on the relative time compared to the middle observation
    # within the window
    "idx" = \(x) nameless::weighted_average(
        x, 
        .by = "index", 
        cols = c("x_actual", "y_actual")
    ),
    "time" = \(x) nameless::weighted_average(
        x, 
        .by = "relative_time", 
        weights = \(x) dnorm(x, mean = 0, sd = 1/10), 
        cols = c("x_actual", "y_actual")
    ),

    # Polynomial approximations, namely linear and quadratic. Somewhat related
    # to the LOESS, but more dumbed down.
    "lin" = \(x) nameless::linear(
        x, 
        cols = c("x_actual", "y_actual")
    ),
    "quad" = \(x) nameless::parabola(
        x, 
        cols = c("x_actual", "y_actual")
    )
)

# Make the combination of the spans and function names for the moving windows.
spans <- c(1, 2, 5)
fx_names <- names(fx)
combos <- data.frame(
    spans = rep(spans, each = length(fx_names)),  
    fx_names = rep(fx_names, times = length(spans))
)

# Create the moving windows themselves
moving_windows <- lapply(
    seq_len(nrow(combos)), 
    function(i) {
        # Initialize the function to be performed on the window
        gx <- fx[[combos$fx[i]]]

        # Initialize the moving window function itself
        factory <- \(x) nameless::moving_window(
            x, 
            span = combos$spans[i], 
            fx = gx, 
            .by = "id"
        )

        return(list(factory))
    }
)
names(moving_windows) <- sapply(
    seq_len(nrow(combos)), 
    \(i) paste(combos$fx[i], combos$spans[i], sep = "_")
)

# Add the moving windows to the filters-list
filters <- append(
    filters,
    moving_windows
)



# General functions ############################################################

# Create a function that will execute the pipeline and output the adjusted 
# movement data.
preprocess <- function(fx) {
    # Check whether the data have a reference to the simulation number. If not, 
    # add it to the dataframe
    if(is.null(data$nsim)) {
        data$nsim <- 1
    }

    # Preprocess the data through a group_by of simulation, error type, and 
    # movement type.
    result <- data %>% 
        dplyr::group_by(kind, type, nsim) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(
            data = data %>% 
                as.data.frame() %>% 
                execute_pipeline(fx = fx, report = FALSE) %>% 
                list()
        ) %>% 
        tidyr::unnest(data)

    return(result)


    # # Execute the pipeline. Changed from doing this in a nested data structure 
    # # to ease debugging when necessary. Also is less susceptible to errors 
    # # triggered in glue (which I couldn't resolve).
    # result <- list() ; f <- 1
    # for(i in unique(local_data$nsim)) {
    #     result[[f]] <- local_data %>% 
    #         dplyr::filter(nsim == i) %>% 
    #         execute_pipeline(pipelines[[fx]], report = FALSE) %>% 
    #         dplyr::mutate(nsim = i)

    #     f <- f + 1
    # }

    # return(do.call("rbind", result))
}

# Filter the data using the different filters that were defined.
results <- lapply(
    filters,
    \(x) preprocess(x)
)
names(results) <- names(filters)

saveRDS(
    results,
    file.path("results", "study 2", "tuning.Rds")
)





################################################################################
# VISUALIZATION

# Trajectories #################################################################

# Create a function that will transform a dataframe to plot_data containing 
# information on the segments that were walked between locations.
segmentize <- function(x, 
                       .vars, 
                       .id) {
    x  %>% 
        # Filter out the id of interest and select the specific range of variables
        # you are interested in
        dplyr::filter(id == .id) %>% 
        dplyr::rename(
            X = .data[[.vars[1]]], 
            Y = .data[[.vars[2]]]
        ) %>% 
        dplyr::select(nsim, time, X, Y) %>% 
        dplyr::arrange(nsim, time) %>% 

        # Group by each simulation, nest, and then transform the points we have
        # to segments by defining the x, y, xend, and yend as used in ggplot2
        dplyr::group_by(nsim) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(
            data = data %>% 
                as.data.frame() %>% 
                dplyr::mutate(
                    x = c(X[2:length(X) - 1], NA), 
                    y = c(Y[2:length(X) - 1], NA), 
                    xend = c(X[2:length(X)], NA), 
                    yend = c(Y[2:length(X)], NA),
                    time = c(diff(time), NA)
                ) %>% 
                dplyr::filter(abs(time) < 0.15) %>% 
                dplyr::select(-X, -Y, -time) %>% 
                dplyr::filter(!is.na(x)) %>% 
                list()
        ) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(data) %>% 
        tidyr::unnest(data) %>% 
        return()
}

# Create a function that will create the wanted plot
trajectory <- function(data) {
    # Loop over each of the id's for a separate plot
    ids <- unique(data$id)

    # Create nameplots already
    plots <- lapply(
        ids, 
        function(x) {
            # Distinguish original data from filtered data when creating plot
            # data. Only the first simulation is necessary for the actual 
            # positions
            original <- data[data$nsim == 1, ] %>% 
                segmentize(
                    .vars = c("x_actual", "y_actual"),
                    .id = x
                ) 

            other <- list(
                segmentize(
                    data, 
                    .vars = c("x", "y"),
                    .id = x
                ),
                segmentize(
                    data, 
                    .vars = c("x_filtered", "y_filtered"),
                    .id = x
                )
            )

            # Compute the MAD for each of the datasets 
            mad <- c(
                data[data$id == x, ] %>% 
                    dplyr::mutate(
                        x = x - x_actual,
                        y = y - y_actual,
                        dist = sqrt(x^2 + y^2)
                    ) %>% 
                    dplyr::select(dist) %>% 
                    unlist() %>% 
                    mean(), 
                data[data$id == x, ] %>% 
                    dplyr::mutate(
                        x = x_filtered - x_actual,
                        y = y_filtered - y_actual,
                        dist = sqrt(x^2 + y^2)
                    ) %>% 
                    dplyr::select(dist) %>% 
                    unlist() %>% 
                    mean()
            )

            # Compute the limits of the plot. Makes sure both plots have the same 
            # limits
            xlims <- c(other[[1]]$x, 
                       other[[2]]$x, 
                       other[[1]]$xend, 
                       other[[2]]$xend) %>% 
                range() 
            ylims <- c(other[[1]]$y, 
                       other[[2]]$y, 
                       other[[1]]$yend, 
                       other[[2]]$yend) %>% 
                range()

            # Make the limits somewhat broader
            xlims <- xlims + c(-1, 1) * diff(xlims) * 0.25
            ylims <- ylims + c(-1, 1) * diff(ylims) * 0.25

            # Create a name plot for the kind of movement add a name-plot
            plt <- lapply(
                other, 
                \(y) ggplot2::ggplot() +
                    ggplot2::geom_segment(data = y, 
                                          ggplot2::aes(x = x,
                                                       y = y,
                                                       xend = xend,
                                                       yend = yend),
                                          color = "grey75",
                                          linewidth = 1,
                                          alpha = 0.1) +
                    ggplot2::geom_segment(data = original, 
                                          ggplot2::aes(x = x,
                                                       y = y,
                                                       xend = xend,
                                                       yend = yend),
                                          color = "black",
                                          linewidth = 1,
                                          alpha = 0.1) +
                    # Distance between measured and real movements
                    ggplot2::annotate("text", 
                                      x = xlims[1] + 0.95 * diff(xlims), 
                                      y = ylims[1] + 0.95 * diff(ylims), 
                                      label = latex2exp::TeX(paste0("$MAD = ", 
                                                                    mad[j],
                                                                    "$")), 
                                      size = 5,
                                      hjust = 1, 
                                      vjust = 1) +
                    # Theme, limits, and labels
                    ggplot2::labs(x = "x", 
                                  y = "y") +
                    ggplot2::lims(x = xlims, 
                                  y = ylims) +
                    ggplot2::theme_minimal() +
                    ggplot2::theme(plot.title = ggplot2::element_text(size = 35, hjust = 0.5), 
                                   axis.title = ggplot2::element_text(size = 25))
            )

            # Add a nameplot to this list
            plt <- list(
                nameless::name_plot(x, size = 20),
                plt
            )

            return(
                ggpubr::ggarrange(
                    plotlist = plt,
                    nrow = 1
                )
            )
        }
    )

    plots <- append(
        ggpubr::ggarrange(
            nameless::name_plot(" ", size = 20),
            nameless::name_plot("Unfiltered", size = 20),
            nameless::name_plot("Filtered", size = 20),
            nrow = 1
        ),
        plots    
    )

    return(
        ggpubr::ggarrange(
            plotlist = plots,
            ncol = 1
        )
    )
}

# Create a dataset that contains all filtered and unfiltered data
data <- lapply(
    results,
    \(x) x %>% 
        # Make sure that you can distinguish real positions, measured positions,
        # and filtered positions through selecting (x, y) from the filtered
        # data and binding it to the unfiltered data
        dplyr::rename(
            x_filtered = x, 
            y_filtered = y
        ) %>% 
        dplyr::mutate(
            x_actual = NA, 
            y_actual = NA
        ) %>% 
        dplyr::select(-x_actual, -y_actual) %>% 
        dplyr::full_join(
            data, 
            by = c("kind", "type", "nsim", "id", "time")
        )
)
names(data) <- names(results)

# Get the unique values for the kinds and the errors. Then you loop over them 
# and create the plot for that unique combination
for(i in names(data)) {
    data_i <- data[[i]]
    for(j in unique(data_i$kind)) {
        for(k in unique(data_i$type)) {
            # Select the relevant data
            idx <- data_i$kind == j & data_i$type == k
            selected_data <- data_i[idx, ]

            # Once you have these data, you can create a trajectory plot and save
            # it
            ggplot2::ggsave(
                file.path("figures", "study 2", "tuning", paste0(i, "_", j, "_", k, ".png")),
                trajectory(selected_data),
                width = 10000,
                height = 3000,
                unit = "px"
            )
        }
    }
}





# Looking at the MAD ###########################################################

# Compute the MADs for each of the results
summarized <- lapply(
    seq_along(results),
    \(i) results[[i]] %>% 
        dplyr::group_by(kind, type) %>% 
        dplyr::mutate(
            x = x - x_actual, 
            y = y - y_actual, 
            dist = sqrt(x^2 + y^2)
        ) %>% 
        dplyr::summarize(
            mad = mean(dist),
            sd_mad = sd(dist)
        ) %>% 
        dplyr::mutate(
            filter = names(results)[i]
        ) %>% 
        dplyr::ungroup()
)
summarized <- do.call("rbind", summarized)

# Check the mean MAD for each filtering technique. These results are used to 
# determine which filters to use in the next step.
summarized %>% 
    dplyr::group_by(filter) %>% 
    dplyr::summarize(
        mad = mean(mad),
        sd_mad = mean(sd_mad)
    )

# Choices made:
#   - Span 5 seems best for the summary statistics
#   - The "linear" and "quadratic" approximations have no advantage versus 
#     over the LOESS, so those are left out
#   - LOESS of the 4th degree seems to have decreased performance compared to 
#     other LOESS