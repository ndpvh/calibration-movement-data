#' Create a pane with a chosen title in it
#' 
#' @param x Title to be plotted in an empty plot.
#' @param ... Additional arguments for \code{annotate}.
#' 
#' @return Ggplot object
#' 
#' @export 
name_plot <- function(x, ...) {
    plt <- ggplot2::ggplot() +
        ggplot2::annotate(
            "text",
            label = x, 
            x = 0, 
            y = 0,
            ...
        ) +
        ggplot2::theme_void()

    return(plt)
}


# # Create a function that takes in a dataframe and creates the plots of interest
# histogram <- function(x, 
#                       statistics) {

#     # Split data before preprocessing and after preprocessing
#     before <- dplyr::filter(x, preprocessed == "before")
#     after <- dplyr::filter(x, preprocessed == "after")

#     # Get the data of before
#     before <- before %>%
#         dplyr::select(contains(statistics)) %>%
#         setNames("X") %>%
#         dplyr::mutate(M = 1)

#     # Get all conditions out of there
#     conditions <- unique(after$preprocessing_function)

#     # Fix the limits on the x-axis (within bounds, of course)
#     all_x <- x[, statistics]

#     if(is.na(sd(all_x))) {
#         xlim <- c(0, 1)
#     } else {
#         if(grepl("rmse", statistics, fixed = TRUE) | 
#            grepl("dist", statistics, fixed = TRUE) |
#            grepl("mae", statistics, fixed = TRUE)) {

#             limit <- max(c(quantile(before$X, probs = 0.95), 
#                            quantile(after[, statistics], probs = 0.95)))

#             limit <- max(c(mean(before$X) + 3 * sd(before$X), 
#                            mean(after[, statistics]) + 3 * sd(after[, statistics])))

#             idx <- all_x < limit

#         } else {
#             limits <- c(min(quantile(before$X, probs = 0.025), 
#                             quantile(after[, statistics], probs = 0.025)), 
#                         max(quantile(before$X, probs = 0.975), 
#                             quantile(after[, statistics], probs = 0.975)))

#             limits <- c(min(c(mean(before$X) - 3 * sd(before$X), 
#                               mean(after[, statistics]) - 3 * sd(after[, statistics]))), 
#                         max(c(mean(before$X) + 3 * sd(before$X), 
#                               mean(after[, statistics]) + 3 * sd(after[, statistics]))))

#             idx <- all_x < limits[2] & all_x > limits[1]
#         }

#         xlim <- range(all_x[idx]) + 0.05 * c(-1, 1) * diff(range(all_x[idx]))

#         if(grepl("mean", statistics, fixed = TRUE) & !grepl("dist", statistics, fixed = TRUE)) {
#             xlim <- c(-max(abs(xlim)), max(abs(xlim)))
#         }

        
#     }

    

#     # Loop over all conditions and create the plot of interest
#     plt <- list()
#     for(i in conditions) {
#         # Get plot data for the condition and the statistic of interest. Bind 
#         # together for before and after
#         plot_data <- after %>%
#             dplyr::filter(preprocessing_function == i) %>%
#             dplyr::select(contains(statistics)) %>%
#             setNames("X") %>%
#             dplyr::mutate(M = 2) %>%
#             rbind(before) %>%
#             dplyr::mutate(M = factor(M))

#         # Create a histogram as the plot of choice. Include the condition name 
#         # in the plot and make the legend tell us something
#         plt[[i]] <- ggplot2::ggplot(data = plot_data, 
#                                     ggplot2::aes(x = X, fill = M)) +
#             ggplot2::geom_histogram(alpha = 0.5, 
#                                     bins = 15, 
#                                     color = "black", 
#                                     position = "identity") +
#             ggplot2::labs(title = i, 
#                           legend = "Preprocessed") +
#             ggplot2::lims(x = xlim) +
#             ggplot2::scale_fill_manual(labels = c("1" = "Before", 
#                                                   "2" = "After"), 
#                                        values = c("1" = "salmon", 
#                                                   "2" = "cornflowerblue")) +
#             ggplot2::theme_minimal() 
#     }

#     # Bind together and save under figures
#     plt <- ggpubr::ggarrange(plotlist = plt, 
#                              nrow = 15, 
#                              ncol = 15,
#                              common.legend = TRUE, 
#                              legend = "right")

#     return(plt)
# }

# # Create a function to create a bar plot for each of the conditions
# barplot <- function(x, 
#                     statistics) {

#     # Split data before preprocessing and after preprocessing
#     before <- dplyr::filter(x, preprocessed == "before")
#     after <- dplyr::filter(x, preprocessed == "after")

#     # Get the data of before
#     before <- before %>%
#         dplyr::select(contains(statistics), preprocessing_function) %>%
#         dplyr::mutate(preprocessing_function = "before") %>% 
#         setNames(c("X", "M"))

#     # Fix the limits on the x-axis (within bounds, of course)
#     all_x <- x[, statistics]

#     if(is.na(sd(all_x))) {
#         xlim <- c(0, 1)
#     } else {
#         if(grepl("rmse", statistics, fixed = TRUE) | 
#            grepl("dist", statistics, fixed = TRUE) |
#            grepl("mae", statistics, fixed = TRUE)) {

#             limit <- max(c(quantile(before$X, probs = 0.95), 
#                            quantile(after[, statistics], probs = 0.95)))

#             limit <- max(c(mean(before$X) + 3 * sd(before$X), 
#                            mean(after[, statistics]) + 3 * sd(after[, statistics])))

#             idx <- all_x < limit

#         } else {
#             limits <- c(min(quantile(before$X, probs = 0.025), 
#                             quantile(after[, statistics], probs = 0.025)), 
#                         max(quantile(before$X, probs = 0.975), 
#                             quantile(after[, statistics], probs = 0.975)))

#             limits <- c(min(c(mean(before$X) - 3 * sd(before$X), 
#                               mean(after[, statistics]) - 3 * sd(after[, statistics]))), 
#                         max(c(mean(before$X) + 3 * sd(before$X), 
#                               mean(after[, statistics]) + 3 * sd(after[, statistics]))))

#             idx <- all_x < limits[2] & all_x > limits[1]
#         }

#         xlim <- range(all_x[idx]) + 0.05 * c(-1, 1) * diff(range(all_x[idx]))

#         if(grepl("mean", statistics, fixed = TRUE) & !grepl("dist", statistics, fixed = TRUE)) {
#             xlim <- c(-max(abs(xlim)), max(abs(xlim)))
#         }        
#     }

#     # Create some plot data that will be used for the barplot
#     conditions <- c("before", unique(after$preprocessing_function))
#     plot_data <- after %>% 
#         dplyr::select(contains(statistics), preprocessing_function) %>% 
#         setNames(c("X", "M")) %>% 
#         rbind(before) %>% 
#         dplyr::group_by(M) %>% 
#         dplyr::summarize(means = mean(X), 
#                          q025 = quantile(X, probs = 0.025),
#                          q975 = quantile(X, probs = 0.975),
#                          sd = sd(X)) %>% 
#         dplyr::ungroup() %>% 
#         dplyr::arrange(factor(M, levels = conditions)) %>% 
#         dplyr::rename(X = M) %>% 
#         dplyr::mutate(M = ifelse(X == "before", 1, 2))

#     # Create a barplot using all of this information. The barplot will show 
#     # the mean levels of each condition, hopefully providing us with a clearer
#     # picture than the histograms
#     plt <- ggplot2::ggplot(data = plot_data) +
#         ggplot2::geom_errorbar(ggplot2::aes(x = reorder(factor(X), -means), 
#                                             ymin = q025, 
#                                             ymax = q975)) +
#         ggplot2::geom_bar(ggplot2::aes(x = reorder(factor(X), -means),
#                                        y = means, 
#                                        fill = factor(M)),
#                           stat = "identity",
#                           color = "black") +
#         ggplot2::coord_flip() +
#         ggplot2::labs(title = paste("Average performance:", statistics), 
#                       legend = "Preprocessed") +
#         ggplot2::scale_fill_manual(labels = c("1" = "Before", 
#                                               "2" = "After"), 
#                                    values = c("1" = "salmon", 
#                                               "2" = "cornflowerblue")) +
#         ggplot2::theme_minimal() 

#     return(list("plot" = plt, "data" = plot_data))
# }

# # Create a function that will create the wanted plot
# trajectory <- function(x) {
#     # Create a function that will transform a dataframe to plot_data containing 
#     # information on the segments that were walked between locations.
#     to_segments <- function(x, 
#                             .vars, 
#                             .id) {
#         x  %>% 
#             dplyr::filter(id == .id) %>% 
#             dplyr::rename(X = tidyselect::all_of(.vars[1]), 
#                           Y = tidyselect::all_of(.vars[2])) %>% 
#             dplyr::select(nsim, time, X, Y) %>% 
#             dplyr::arrange(nsim, time) %>% 
#             dplyr::group_by(nsim) %>% 
#             tidyr::nest() %>% 
#             dplyr::mutate(data = data %>% 
#                               as.data.frame() %>% 
#                               dplyr::mutate(x = c(X[2:length(X) - 1], NA), 
#                                             y = c(Y[2:length(X) - 1], NA), 
#                                             xend = c(X[2:length(X)], NA), 
#                                             yend = c(Y[2:length(X)], NA),
#                                             time = c(diff(time), NA)) %>% 
#                               dplyr::filter(abs(time) < 0.15) %>% 
#                               dplyr::select(-X, -Y, -time) %>% 
#                               dplyr::filter(!is.na(x)) %>% 
#                               list()) %>% 
#             dplyr::ungroup() %>% 
#             dplyr::select(data) %>% 
#             tidyr::unnest(data) %>% 
#             return()
#     }

#     # Loop over each of the id's for a separate plot
#     ids <- unique(x$id)

#     # Create name-plots that denote whatever it is you're seeing
#     name_plot <- function(x) {
#         return(ggplot2::ggplot() +
#             ggplot2::annotate("text", 
#                               x = 0, 
#                               y = 0,
#                               label = x,
#                               size = 10,
#                               hjust = 0.5, 
#                               vjust = 0.5) +
#             ggplot2::theme_void())
#     }

#     plt <- list()
#     plt[[1]] <- name_plot(" ")
#     plt[[2]] <- name_plot("Unfiltered")
#     plt[[3]] <- name_plot("Filtered")

#     f <- length(plt) + 1
#     for(i in seq_along(ids)) {
#         # Get the original data and make them in plot data (x, y, xend, yend)
#         original <- to_segments(dplyr::filter(x, preprocessed == "before"), 
#                                 .vars = c("x_original", "y_original"), 
#                                 .id = ids[i])

#         # Get filtered and unfiltered data
#         other <- list(to_segments(dplyr::filter(x, preprocessed == "before"), 
#                                   .vars = c("x", "y"), 
#                                   .id = ids[i]), 
#                       to_segments(dplyr::filter(x, preprocessed == "after"), 
#                                   .vars = c("x", "y"), 
#                                   .id = ids[i]))

#         # Compute the standard deviations between the actual movement and the 
#         # measured movement in both cases. Will be  
#         RMSE <- c(x %>% 
#                       dplyr::filter(preprocessed == "before") %>% 
#                       dplyr::mutate(dist = (x - x_original)^2 + (y - y_original)^2,
#                                     dist = sqrt(dist)) %>% 
#                       dplyr::summarize(rmse = rmse(dist)) %>% 
#                       dplyr::select(rmse) %>% 
#                       unlist() %>% 
#                       as.numeric(), 
#                   x %>% 
#                       dplyr::filter(preprocessed == "after") %>% 
#                       dplyr::mutate(dist = (x - x_original)^2 + (y - y_original)^2,
#                                     dist = sqrt(dist)) %>% 
#                       dplyr::summarize(rmse = rmse(dist)) %>% 
#                       dplyr::select(rmse) %>% 
#                       unlist() %>% 
#                       as.numeric())

#         # Compute the limits of the plot. Makes sure both plots have the same 
#         # limits
#         xlims <- c(other[[1]]$x, 
#                    other[[2]]$x, 
#                    other[[1]]$xend, 
#                    other[[2]]$xend) %>% 
#             range() 
#         ylims <- c(other[[1]]$y, 
#                    other[[2]]$y, 
#                    other[[1]]$yend, 
#                    other[[2]]$yend) %>% 
#             range()

#         # Make the limits somewhat broader
#         xlims <- xlims + c(-1, 1) * diff(xlims) * 0.25
#         ylims <- ylims + c(-1, 1) * diff(ylims) * 0.25

#         # Create a name plot for the kind of movement
#         # Add a name-plot
#         plt[[f]] <- name_plot(ids[i])
#         f <- f + 1

#         # And make the plots for filtered and unfiltered data
#         for(j in seq_along(other)) {
#             plt[[f]] <- ggplot2::ggplot() +
#                 # Measured vs real movements
#                 ggplot2::geom_segment(data = other[[j]], 
#                                       ggplot2::aes(x = x, 
#                                                    y = y, 
#                                                    xend = xend, 
#                                                    yend = yend), 
#                                       color = "grey75", 
#                                       linewidth = 1, 
#                                       alpha = 0.1) +
#                 ggplot2::geom_segment(data = original,
#                                       ggplot2::aes(x = x, 
#                                                    y = y, 
#                                                    xend = xend, 
#                                                    yend = yend), 
#                                       color = "black", 
#                                       linewidth = 1) +
#                 # Distance between measured and real movements
#                 ggplot2::annotate("text", 
#                                   x = xlims[1] + 0.95 * diff(xlims), 
#                                   y = ylims[1] + 0.95 * diff(ylims), 
#                                   label = latex2exp::TeX(paste0("$RMSE = ", 
#                                                                 RMSE[j],
#                                                                 "$")), 
#                                   size = 5,
#                                   hjust = 1, 
#                                   vjust = 1) +
#                 # Theme, limits, and labels
#                 ggplot2::labs(x = "x", 
#                               y = "y") +
#                 ggplot2::lims(x = xlims, 
#                               y = ylims) +
#                 ggplot2::theme_minimal() +
#                 ggplot2::theme(plot.title = ggplot2::element_text(size = 35, hjust = 0.5), 
#                                axis.title = ggplot2::element_text(size = 25))
#             f <- f + 1
#         }
#     }

#     return(ggpubr::ggarrange(plotlist = plt, 
#                              nrow = 10, 
#                              ncol = 3))
# }