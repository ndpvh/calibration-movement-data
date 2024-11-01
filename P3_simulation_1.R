################################################################################
# Purpose: Do formal analyses of the results retrieved in P2.                  #
#                                                                              #
#          First, we visualize the results of preprocessing the data in a      #
#          nice figure displaying actual movement and comparing it to measured #
#          movement before and after preprocessing.                            #
#                                                                              #
#          Then, we move on to analyzing the results of the simulation,        #
#          showcasing which of the preprocessing pipelines was the best        #
#          overall (across fixed and movement datasets).                       #
################################################################################

devtools::load_all()

################################################################################
# PRELIMINARIES

#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------

# Define the measurement error and the kind of data. This way, we can read in 
# all version for movement and fixed data and bind them together in one big 
# datafile. This will greatly enhance the comphrensibility and information in 
# our plots
technique <- c("R10", "U10", "T10", "R6R", "U6R", "T6R", "R6N", "U6N", "T6N")
kind <- c("movement", "fixed")

# Read in the datafiles (used for plotting the figures)
data_list <- lapply(kind, 
                    function(x) {
                        # Create filenames
                        filename <- paste(rep(x, each = length(technique)), 
                                          technique, 
                                          sep = "_")

                        # Read in all datafiles
                        data <- lapply(filename, 
                                       \(i) data.table::fread(file.path("data", 
                                                                        "simulation_1", 
                                                                        paste0(i, ".csv"))))

                        # Don't use all data, but only a selection
                        data <- lapply(data, 
                                       \(x) dplyr::filter(x, nsim %in% 1:20))

                        # Add an indicator to each of the datafiles and change 
                        # nsim so that it adds up across all datasets. 
                        # Otherwise might give issues
                        start <- seq(0, (length(filename) - 1) * 100, by = 100)
                        for(i in seq_along(start)) {
                            data[[i]]$nsim <- data[[i]]$nsim + start[i]
                        }

                        # Bind together 
                        return(do.call("rbind", data))
                    })
names(data_list) <- kind





################################################################################
# VISUALIZATION

#-------------------------------------------------------------------------------
# Pipelines
#-------------------------------------------------------------------------------

# Define all moving windows. To allow for the stable creation of different 
# moving windows in for-loops, we will need to create a wrapper-function that 
# takes in the variable arguments and outputs the function to be used in the 
# pipeline. 
fx <- list("av" = \(x) nameless::average(x), 
           "idx" = \(x) nameless::weighted_average(x, .by = "index"),
           "time" = \(x) nameless::weighted_average(x, .by = "relative_time", weights = \(x) dnorm(x, mean = 0, sd = 1/10)),
           "lin" = \(x) nameless::linear(x),
           "quad" = \(x) nameless::parabola(x))

# Make the combination of the spans and function names for the moving windows.
spans <- c(1, 2, 5)
fx_names <- names(fx)
combos <- data.frame(spans = rep(spans, each = length(fx_names)), 
                     fx_names = rep(fx_names, times = length(spans)))

# Create the moving windows themselves
pipelines <- lapply(seq_len(nrow(combos)), 
                    function(i) {
                       # Initialize the function to be performed on the window
                       gx <- fx[[combos$fx[i]]]

                       # Initialize the moving window function itself
                       factory <- \(x) nameless::moving_window(x, span = combos$spans[i], fx = gx, .by = "id")
                       return(list(factory))
                   })
names(pipelines) <- sapply(seq_len(nrow(combos)), 
                           \(i) paste(combos$spans[i], combos$fx[i], sep = "_"))

# Define the Kalman filters and put them in the pipelines
pipelines[["kalm_rev_cv"]] <- list(\(x) nameless::kalman_filter(x, reverse = TRUE, .by = "id", assumed_variance = 0.02^2))
pipelines[["kalm_norev_cv"]] <- list(\(x) nameless::kalman_filter(x, reverse = FALSE, .by = "id", assumed_variance = 0.02^2))

# Also add a loess of varying degrees
pipelines[["loess_1_75"]] <- list(\(x) nameless::local_regression(x, degree = 1, span = 0.75))
pipelines[["loess_2_75"]] <- list(\(x) nameless::local_regression(x, degree = 2, span = 0.75))
pipelines[["loess_1_50"]] <- list(\(x) nameless::local_regression(x, degree = 1, span = 0.50))
pipelines[["loess_2_50"]] <- list(\(x) nameless::local_regression(x, degree = 2, span = 0.50))
pipelines[["loess_1_25"]] <- list(\(x) nameless::local_regression(x, degree = 1, span = 0.25))
pipelines[["loess_2_25"]] <- list(\(x) nameless::local_regression(x, degree = 2, span = 0.25))

# Define the link between pipelines and data
data_files <- data.frame(filename = rep(names(data_list), each = length(pipelines)), 
                         condition = rep(names(pipelines), times = length(data_list)))





#-------------------------------------------------------------------------------
# Functions for plotting and executing
#-------------------------------------------------------------------------------

# Create a function that will execute the pipeline and output the adjusted 
# movement data.
preprocess <- function(x) {
    # Retrieve the data and the pipeline for the condition
    local_data <- data_list[[x$filename]]
    fx <- x$condition

    # Check whether the data have a reference to the simulation number. If not, 
    # add it to the dataframe
    if(is.null(local_data$nsim)) {
        local_data$nsim <- 1
    }

    # Execute the pipeline 
    result <- local_data %>% 
        dplyr::group_by(nsim) %>%
        tidyr::nest() %>% 
        dplyr::mutate(data = data %>% 
                          as.data.frame() %>% 
                          nameless::execute_pipeline(pipelines[[fx]], 
                                                     report = FALSE) %>% 
                          list()) %>% 
        tidyr::unnest(data) %>% 
        dplyr::ungroup()

    return(result)
}

# Create a function that will transform a dataframe to plot_data containing 
# information on the segments that were walked between locations.
to_segments <- function(x, 
                        .vars, 
                        .id) {
    x  %>% 
        dplyr::filter(id == .id) %>% 
        dplyr::rename(X = tidyselect::all_of(.vars[1]), 
                      Y = tidyselect::all_of(.vars[2])) %>% 
        dplyr::select(nsim, time, X, Y) %>% 
        dplyr::arrange(nsim, time) %>% 
        dplyr::group_by(nsim) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(data = data %>% 
                          as.data.frame() %>% 
                          dplyr::mutate(x = c(X[2:length(X) - 1], NA), 
                                        y = c(Y[2:length(X) - 1], NA), 
                                        xend = c(X[2:length(X)], NA), 
                                        yend = c(Y[2:length(X)], NA),
                                        time = c(diff(time), NA)) %>% 
                          dplyr::filter(abs(time) < 0.15) %>% 
                          dplyr::select(-X, -Y, -time) %>% 
                          dplyr::filter(!is.na(x)) %>% 
                          list()) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(data) %>% 
        tidyr::unnest(data) %>% 
        return()
}

# Create a function that will compute the RMSE for a given set of data
compute_rmse <- function(x, 
                         .vars, 
                         .id) {
    x %>% 
        dplyr::filter(id == .id) %>% 
        dplyr::rename(X1 = tidyselect::all_of(.vars[1]), 
                      Y1 = tidyselect::all_of(.vars[2]), 
                      X2 = x_original, 
                      Y2 = y_original) %>% 
        dplyr::mutate(distance = sqrt((X1 - X2)^2 + (Y1 - Y2)^2)) %>% 
        dplyr::select(distance) %>% 
        unlist() %>% 
        as.numeric() %>% 
        sd() %>% 
        round(digits = 4) %>% 
        format(nsmall = 4) %>% 
        return()
}

# Create a function that will create the wanted plot
trajectory <- function(x) {
    # Retrieve data and bind it together with the preprocessed data. We add the
    # columns x_original and y_original to the dataframe to make sure we can 
    # delete them if present in the preprocessed data before joining with the 
    # original data. Differentially handled by Kalman filters than moving 
    # windows, as the latter needs explicit inclusion of columns while the 
    # former does this by default.
    local_data <- preprocess(x) %>% 
        dplyr::rename(filtered_x = x, 
                      filtered_y = y) %>% 
        dplyr::mutate(x_original = NA, 
                      y_original = NA) %>% 
        dplyr::select(-x_original, -y_original) %>% 
        dplyr::full_join(data_list[[x$filename]], 
                         by = c("nsim", "id", "time"))

    # Loop over each of the id's for a separate plot
    ids <- unique(local_data$id)

    # Create name-plots that denote whatever it is you're seeing
    name_plot <- function(x) {
        return(ggplot2::ggplot() +
            ggplot2::annotate("text", 
                              x = 0, 
                              y = 0,
                              label = x,
                              size = 10,
                              hjust = 0.5, 
                              vjust = 0.5) +
            ggplot2::theme_void())
    }

    plt <- list()
    plt[[1]] <- name_plot(" ")
    plt[[2]] <- name_plot("Unfiltered")
    plt[[3]] <- name_plot("Filtered")

    f <- length(plt) + 1
    for(i in seq_along(ids)) {
        # Get the original data and make them in plot data (x, y, xend, yend)
        original <- to_segments(local_data %>% 
                                    dplyr::filter(nsim == 1), 
                                .vars = c("x_original", "y_original"), 
                                .id = ids[i])            

        # Get filtered and unfiltered data
        other <- list(to_segments(local_data, 
                                  .vars = c("x", "y"), 
                                  .id = ids[i]), 
                      to_segments(local_data, 
                                  .vars = c("filtered_x", "filtered_y"), 
                                  .id = ids[i]))

        # Compute the standard deviations between the actual movement and the 
        # measured movement in both cases. Will be  
        rmse <- c(compute_rmse(local_data, 
                               .vars = c("x", "y"), 
                               .id = ids[i]), 
                  compute_rmse(local_data, 
                               .vars = c("filtered_x", "filtered_y"), 
                               .id = ids[i]))

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

        # Create a name plot for the kind of movement
        # Add a name-plot
        plt[[f]] <- name_plot(ids[i])
        f <- f + 1

        # And make the plots for filtered and unfiltered data
        for(j in seq_along(other)) {
            plt[[f]] <- ggplot2::ggplot() +
                # Measured vs real movements
                ggplot2::geom_segment(data = other[[j]], 
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
                                      linewidth = 1) +
                # Distance between measured and real movements
                ggplot2::annotate("text", 
                                  x = xlims[1] + 0.95 * diff(xlims), 
                                  y = ylims[1] + 0.95 * diff(ylims), 
                                  label = latex2exp::TeX(paste0("$RMSE = ", 
                                                                rmse[j],
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
            f <- f + 1
        }
    }

    return(plt)

}

# Now that we have this, we can create plots for each of the different pipelines
for(i in seq_len(nrow(data_files))) {
    cat(paste0("\rCreating plot ", i, " of ", nrow(data_files)))

    # Can comment out if you want to look at a specific analysis
    # if(!grepl("kalm_rev_cv", data_files$condition[i], fixed = TRUE)) {
    #     next
    # }

    # Create plots
    plt <- trajectory(data_files[i, ]) 

    # Bind them together and save them somewhere
    plt <- ggpubr::ggarrange(plotlist = plt, 
                             nrow = 10, 
                             ncol = 3)
    ggplot2::ggsave(plt, 
                    filename = file.path("figures", 
                                         "simulation_1", 
                                         "preprocessed", 
                                         paste0(data_files$filename[i], "__", data_files$condition[i], ".png")),
                    width = 900 * 3, 
                    height = 1000 * 10,
                    unit = "px")

    if(i == nrow(data_files)) {
        print(" ")
    }
}






















