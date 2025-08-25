################################################################################
# Purpose: Perform an application of the filtering pipeline and check the      #
#          nonfiltered vs filtered speeds of the sled, for which we have real  #
#          speeds available.                                                   #
################################################################################

devtools::load_all()

################################################################################
# DATAFILES

# Load in the datafiles in which all results are saved. 
data <- data.table::fread(
    file.path("data", "raw_data", "datapoints.csv"),
    data.table = FALSE,
    fill = TRUE
)

# Load the experiment names and join together with the bigger datafile. 
data <- data.table::fread(
    file.path("data", "raw_data", "experiments.csv"),
    data.table = FALSE
) %>% 
    dplyr::rename(experiment_id = id) %>% 
    dplyr::select(experiment_id, name) %>% 
    dplyr::full_join(data, by = "experiment_id") %>% 
    dplyr::select(-experiment_id) %>% 
    dplyr::rename(experiment = name)

# Select the data that we need and differentiate between the actual walking data
# and the measurements of the points participants had to walk towards
data <- data %>% 
    dplyr::filter(
        stringr::str_detect(experiment, "REC E") & 
        stringr::str_detect(experiment, "walking")
    ) 

boundaries <- data %>% 
    dplyr::filter(stringr::str_detect(experiment, "centerpoint")) %>% 
    dplyr::group_by(tag_id) %>% 
    dplyr::summarize(
        x = mean(x, na.rm = TRUE), 
        y = mean(y, na.rm = TRUE)
    )

data <- data %>% 
    dplyr::filter(stringr::str_detect(experiment, "centerpoint", negate = TRUE))

# Delete all data that fall outside of the boundaries defined by the centerpoint 
# data. Additionally, define all data that defines whether a person has already 
# turned or not. 
data <- data %>% 
    dplyr::filter(
        x <= max(boundaries$x) & x >= min(boundaries$x), 
        y <= max(boundaries$y) & y >= min(boundaries$y)
    ) %>% 
    dplyr::mutate(
        turned = y <= boundaries$y[boundaries$tag_id == 15]
    ) %>% 
    dplyr::filter(
        !(tag_id %in% c(12, 17))
    )

# Add a time-variable to the data (expressed in sec) and remove unnecessary 
# information
data <- data %>% 
    dplyr::mutate(
        time = withr::with_options(
            list(digits = 16),
            as.numeric(timestamp)
        ), 
        time = time - min(time)
    ) %>% 
    dplyr::select(tag_id, time, x, y, turned)

# Add information on the condition, both based on whether a new round was started
# and where the final observation took place. Also filter out those people of 
# whom we don't have both turning and non-turning data, needing at least two 
# datapoints in order to compute speeds
center <- boundaries[boundaries$tag_id == 15, ]
data <- data %>% 
    dplyr::group_by(tag_id) %>% 
    tidyr::nest() %>% 
    dplyr::mutate(
        data = data[[1]] %>% 
            # Define the different rounds of walking that a participant did
            dplyr::arrange(time) %>% 
            dplyr::mutate(
                diff_time = c(0, diff(time)),
                round = ifelse(
                    diff_time > 30,
                    1, 
                    0
                ),
                round = cumsum(round) + 1
            ) %>% 
            dplyr::select(-diff_time) %>% 

            # Per round, we also need to do some additional stuff
            dplyr::group_by(round) %>% 
            tidyr::nest() %>% 
            dplyr::mutate(
                data = data[[1]] %>% 
                    # Delete data that do not have enough information to estimate
                    # speeds before and after turning
                    dplyr::mutate(
                        n_turned = sum(turned),
                        n_not_turned = sum(!turned)
                    ) %>% 
                    dplyr::mutate(complete = n_turned >= 2 & n_not_turned >= 2) %>% 
                    dplyr::select(-n_turned, -n_not_turned) %>% 

                    # Define the angles at which participants walked according 
                    # to their frame of reference, where walking straight ahead
                    # is equal to 0 degrees. Angles are transformed from radians
                    # to degrees.
                    dplyr::mutate(
                        angle = atan2(
                            y[length(y)] - center$y, 
                            x[length(x)] - center$x
                        ), 
                        angle = 180 * angle / pi + 90, 
                    ) %>% 
                    list()
            ) %>% 
            tidyr::unnest(data) %>% 
            dplyr::ungroup() %>% 
            dplyr::filter(complete) %>% 
            dplyr::select(-complete) %>% 
            list()
    ) 

# Link the walked angles to the experimental conditions. To make sure each 
# condition happens only once, we update the exp_angles list multiple times.
# We create a function for this and then apply it within a mutate statement
exp_angles <- c(72.5, 50, 32.5, 20, 10, 0, -10, -20, -32.5, -50, -72.5)
apply_condition <- function(x) {
    # Create a data.frame containing all the angles per round
    rounds <- unique(x$round)
    angles <- lapply(
        rounds, 
        \(x) matrix(
            c(x, exp_angles), 
            nrow = 1
        )
    )
    angles <- do.call("rbind", angles) %>% 
        as.data.frame() %>% 
        setNames(c("round", paste0("a_", exp_angles)))

    # Subtract the observed angles from this per round
    angles <- x %>% 
        dplyr::group_by(round, angle) %>% 
        tidyr::nest() %>% 
        dplyr::select(-data) %>% 
        dplyr::ungroup() %>% 
        dplyr::rename(obs = angle) %>% 
        dplyr::full_join(angles, by = "round") %>% 
        dplyr::mutate(
            a_72.5 = abs(a_72.5 - obs), 
            a_50 = abs(a_50 - obs),
            a_32.5 = abs(a_32.5 - obs),
            a_20 = abs(a_20 - obs),
            a_10 = abs(a_10 - obs),
            a_0 = abs(a_0 - obs),
            `a_-10` = abs(`a_-10` - obs),
            `a_-20` = abs(`a_-20` - obs),
            `a_-32.5` = abs(`a_-32.5` - obs),
            `a_-50` = abs(`a_-50` - obs),
            `a_-72.5` = abs(`a_-72.5` - obs),
        )

    # Retrieve whichever index was observed to be the lowest for each row
    angles$idx <- 0
    for(i in seq_len(nrow(angles))) {
        angles$idx[i] <- which.min(angles[i, 3:13]) %>% 
            as.numeric()
    }

    # Create a function to change indices when there are doubles
    correct_idx <- function(y) {
        for(i in 1:nrow(y)) {
            if(y$diff_idx[i] == 0) {
                y$idx[i - 1] <- ifelse(
                    (y$idx[i - 1] + 1) %in% unique(y$idx), 
                    y$idx[i - 1] - 1, 
                    y$idx[i - 1] + 1
                )
            }
        }

        return(y)
    } 

    # Add information on the manipulated speed and angles, finally 
    # differentiating between different conditions
    angles <- angles %>% 
        dplyr::arrange(round) %>% 
        dplyr::mutate(
            diff_idx = c(1, diff(idx))
        ) %>% 

        # Manipulated speed
        dplyr::mutate(
            exp_speed = ifelse(
                diff_idx < 0, 
                1, 
                0
            ),
            exp_speed = cumsum(exp_speed) + 1
        ) %>% 

        # Manipulated angle: First we need to do some preprocessing, ensuring 
        # that each index only occurs once per round for the speeds
        dplyr::group_by(exp_speed) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(
            data = data[[1]] %>% 
                correct_idx() %>% 
                list()
        ) %>% 
        tidyr::unnest(data) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(
            exp_angle = exp_angles[idx]
        ) 

    # Clean up and add this information to the actual data.frame
    angles %>% 
        dplyr::select(exp_speed, exp_angle, round) %>% 
        dplyr::full_join(x, by = "round") %>% 
        return()
}

data <- data %>% 
    dplyr::filter(tag_id != 40) %>% 
    dplyr::mutate(
        data = data[[1]] %>% 
            apply_condition() %>% 
            list()
    ) %>% 
    tidyr::unnest(data) %>% 
    dplyr::ungroup() 

# Add information on the anchor locations to these data: Needed to filter out 
# the systematic error
anchors <- readRDS(file.path("data", "anchor_position_16-11-2024.Rds"))
data <- data %>%
    dplyr::mutate(
        anchor_xmin = min(anchors[, 2]),
        anchor_xmax = max(anchors[, 2]),
        anchor_ymin = min(anchors[, 3]),
        anchor_ymax = max(anchors[, 3]),
    )

# Finally, there was one additional round where one of the participants redid 
# the experiment whilst having a high pace. For this analysis, we leave these
# data out of the equation
data <- data %>% 
    dplyr::filter(exp_speed != 4)





################################################################################
# FILTERS

# Generate a data_list containing two copies of the same data; One having the 
# original dataset, the other having the filtered dataset
data_list <- list(
    "unfiltered" = data, 
    "filtered" = data 
)



############################
# Systematic distortion: Use the parameters estimated on this day to restore 
# the observed distortion

# Create a function that will correct the data based on a polynomial estimated 
# on these data
correct <- function(x) {
    # Retrieve the parameters of the polynomial as estimated in Study 1
    params <- readRDS(file.path("results", "study 1", paste0("polynomial_16-11-2024.Rds")))
    B <- params$B 

    # Create the X matrix for these data
    X <- matrix(
        1, 
        nrow = nrow(x),
        ncol = 1
    )

    for(i in 1:9) {
        X <- cbind(X, x$x^i, x$y^i)
    }

    # Compute the result Y
    Y <- X %*% B

    # Remove expected instances that fall outside of the bounds
    Y <- ifelse(
        (Y > 1) | (Y < -1), 
        NA, 
        Y
    )
    idx <- is.na(Y[, 1]) | is.na(Y[, 2])
    Y[idx, ] <- NA

    # Add to data.frame
    x$x <- Y[, 1]
    x$y <- Y[, 2]

    return(x)
}

# Correct the data and transform the positions back to their original locations.
data_list[["filtered"]] <- data_list[["filtered"]] %>% 
    # Transform to -1, 1 range based on anchor positions
    dplyr::mutate(
        x = 2 * (x - anchor_xmin) / (anchor_xmax - anchor_xmin) - 1,
        y = 2 * (y - anchor_ymin) / (anchor_ymax - anchor_ymin) - 1
    ) %>% 
    # Use the polynomial to correct the distortion
    correct() %>% 
    # Transform back to original scale
    dplyr::filter(!is.na(x)) %>% 
    dplyr::mutate(
        x = (anchor_xmax - anchor_xmin) * (x + 1) / 2 + anchor_xmin,
        y = (anchor_ymax - anchor_ymin) * (y + 1) / 2 + anchor_ymin
    ) 



############################
# Unsystematic error: Use the Kalman filter to filter out unsystematic error
# using the observed error variance of this day

# Use the Kalman filter to filter the unsystematic error of the data
filter <- list(
    \(x) nameless::kalman_filter(
        x, 
        assumed_variance = 0.00303601,
        reverse = FALSE,
        .by = c("round", "tag_id")
    )
)

data_list[["filtered"]] <- nameless::execute_pipeline(
    data_list[["filtered"]],
    filter, 
    report = FALSE 
)


############################
# Compare filtered and unfiltered results visually and save the result
base::plot(data$x, data$y)
base::plot(data_list[["filtered"]]$x, data_list[["filtered"]]$y)

saveRDS(
    data_list, 
    file.path("results", "study 4", "data_list.Rds")
)





################################################################################
# ANALYSIS

data_list <- readRDS(file.path("results", "study 4", "data_list.Rds"))

# For each of the datasets, compute the observed distances that have been 
# travelled between consecutive points in time and the speed required for this, 
# and this for each line and voltage separately
results <- lapply(
    data_list, 
    \(x) x %>% 
        dplyr::group_by(exp_speed, exp_angle, tag_id) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(
            data = data[[1]] %>% 
                dplyr::arrange(time) %>% 
                dplyr::mutate(
                    diff_x = c(NA, x[2:length(x)] - x[2:length(x) - 1]), 
                    diff_y = c(NA, y[2:length(y)] - y[2:length(y) - 1]),
                    distance = diff_x^2 + diff_y^2,
                    distance = sqrt(distance),
                    diff_time = c(NA, time[2:length(time)] - time[2:length(time) - 1]),
                    speed = distance / diff_time
                ) %>% 
                dplyr::select(-diff_x, -diff_y, -diff_time) %>% 
                list()
        ) %>% 
        tidyr::unnest(data) %>% 
        dplyr::ungroup() %>% 
        dplyr::filter(!is.na(speed)) %>% 
        dplyr::select(exp_speed, exp_angle, tag_id, time, x, y, speed, distance, turned) %>% 
        return()
) %>% 
    `names<-` (names(data_list))

# Manipulation check: What is the range of speeds for each speed manipulation
model <- aov(
    data = results[[1]],
    speed ~ factor(exp_speed)
)
summary(model)
lsr::etaSquared(model)

results[[1]] %>% 
    dplyr::mutate(exp_angle = abs(exp_angle)) %>% 
    dplyr::group_by(exp_speed, exp_angle) %>% 
    dplyr::summarize(
        mean_speed = mean(speed), 
        sd_speed = sd(speed),
        q025_speed = quantile(speed, 0.025),
        q975_speed = quantile(speed, 0.975)
    )


model <- aov(
    data = results[[2]],
    speed ~ factor(exp_speed)
)
summary(model)
lsr::etaSquared(model)

results[[2]] %>% 
    dplyr::mutate(exp_angle = abs(exp_angle)) %>% 
    dplyr::group_by(exp_speed, exp_angle) %>% 
    dplyr::summarize(
        mean_speed = mean(speed), 
        sd_speed = sd(speed),
        q025_speed = quantile(speed, 0.025),
        q975_speed = quantile(speed, 0.975)
    )

# Having these speeds, we will now look at amplitudes or changes in speed. 
# We will do this in a typical way: We define the baseline speed in the period 
# before turning as the mean of all observations in that period and define the 
# amplitude as the difference between this baseline and the minimum speed 
# observed within a 2sec window after the turning point
amplitudes <- lapply(
    results, 
    \(x) x %>% 
        # First, center the time variable at 0, allowing us to use relative 
        # times per round of walking
        dplyr::group_by(tag_id, exp_speed, exp_angle) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(
            data = data[[1]] %>% 
                dplyr::mutate(
                    diff_turned_1 = c(diff(turned), 0), 
                    diff_turned_2 = c(0, diff(turned)),
                    time_turned = mean(
                        c(
                            time[diff_turned_1 == 1], 
                            time[diff_turned_2 == 1]
                        )
                    ),
                    time = time - time_turned
                ) %>% 
                dplyr::select(-diff_turned_1, -diff_turned_2, -time_turned) %>% 
                list()
        ) %>% 
        tidyr::unnest(data) %>% 
        dplyr::ungroup() %>% 
        
        # Compute the amplitudes of the speed. Importantly, we use all available
        # data for the baseline, but only 2sec for the minimum speed, for which
        # we filter the data first
        dplyr::filter(time <= 1) %>% 
        dplyr::group_by(tag_id, exp_speed, exp_angle, turned) %>% 
        dplyr::summarize(
            mean_speed = mean(speed), 
            min_speed = min(speed)
        ) %>% 
        dplyr::ungroup() %>% 
        dplyr::group_by(tag_id, exp_speed, exp_angle) %>% 
        tidyr::pivot_wider(
            values_from = c(mean_speed, min_speed), 
            names_from = turned
        ) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(
            amplitude_speed = min_speed_TRUE - mean_speed_FALSE
        ) %>% 
        dplyr::select(tag_id, exp_speed, exp_angle, amplitude_speed)
) 

# We expect no differences between turning left or right, so we will take the 
# absolute value of exp_angle
amplitudes <- lapply(
    amplitudes, 
    \(x) x %>% 
        dplyr::mutate(exp_angle = abs(exp_angle))
) %>% 
    `names<-` (names(results))

# Let's look at some things: Are there differences between speed and angle on 
# the observed amplitude?
#
# Both seem to agree on a main effect of speed, but not of the other variables
model <- aov(
    data = amplitudes[[1]],
    amplitude_speed ~ exp_speed * exp_angle
)
summary(model)
lsr::etaSquared(model)

amplitudes[[1]] %>% 
    dplyr::group_by(exp_speed, exp_angle) %>% 
    dplyr::summarize(
        mean = mean(amplitude_speed), 
        sd = sd(amplitude_speed),
        q025 = quantile(amplitude_speed, 0.025),
        q975 = quantile(amplitude_speed, 0.975)
    )

model <- aov(
    data = amplitudes[[2]],
    amplitude_speed ~ exp_speed * exp_angle
)
summary(model)
lsr::etaSquared(model)

amplitudes[[2]] %>% 
    dplyr::group_by(exp_speed, exp_angle) %>% 
    dplyr::summarize(
        mean = mean(amplitude_speed), 
        sd = sd(amplitude_speed),
        q025 = quantile(amplitude_speed, 0.025),
        q975 = quantile(amplitude_speed, 0.975)
    )





################################################################################
# VISUALIZATION

# Plot of the speeds for each participant at each time
for(i in seq_along(results)) {
    plt <- lapply(
        unique(results[[i]]$exp_speed),
        function(x) {
            tmp <- results[[i]]
            tmp <- tmp[tmp$exp_speed == x, ]

            plt <- lapply(
                unique(tmp$exp_angle),
                function(y) {
                    tmp <- tmp[tmp$exp_angle == y, ]

                    # Center the time variable at 0 around the turning point 
                    # per person in this dataset. Will make everything easier
                    # to interpret
                    tmp <- tmp %>% 
                        dplyr::group_by(tag_id, exp_speed, exp_angle) %>% 
                        tidyr::nest() %>% 
                        dplyr::mutate(
                            data = data[[1]] %>% 
                                dplyr::mutate(
                                    diff_turned_1 = c(diff(turned), 0), 
                                    diff_turned_2 = c(0, diff(turned)),
                                    time_turned = mean(
                                        c(
                                            time[diff_turned_1 == 1], 
                                            time[diff_turned_2 == 1]
                                        )
                                    ),
                                    time = time - time_turned
                                ) %>% 
                                dplyr::select(-diff_turned_1, -diff_turned_2, -time_turned) %>% 
                                list()
                        ) %>% 
                        tidyr::unnest(data) %>% 
                        dplyr::ungroup()

                    # Creation of actual plot
                    plt <- ggplot2::ggplot(
                        data = tmp, 
                        ggplot2::aes(
                            x = time, 
                            y = speed,
                            color = factor(tag_id)
                        )
                    ) + 
                        ggplot2::geom_line() +
                        ggplot2::geom_vline(
                            xintercept = 0, 
                            color = "red"
                        ) +
                        ggplot2::labs(
                            title = paste0(
                                "speed = ", 
                                tmp$exp_speed[1], 
                                ", angle = ", 
                                y
                            )
                        ) +
                        ggplot2::theme(
                            legend.position = "none"
                        )

                    return(plt)
                }
            )

            plt <- ggpubr::ggarrange(
                plotlist = plt, 
                nrow = 1
            )
        }
    )

    plt <- ggpubr::ggarrange(
        plotlist = plt, 
        ncol = 1
    )

    ggplot2::ggsave(
        file.path("figures", "study 4", paste0("turning - velocity ", names(results)[i], ".png")),
        plt,
        width = 15000, 
        height = 3500, 
        unit = "px",
        limitsize = FALSE
    )
} 

# Create histograms of the amplitudes per condition
plotlist <- list()
limits <- do.call("rbind", amplitudes) %>% 
    dplyr::summarize(
        min = min(amplitude_speed), 
        max = max(amplitude_speed)
    ) %>% 
    as.numeric()

for(i in names(amplitudes)) {
    plot_data <- amplitudes[[i]]

    plt <- lapply(
        c(0, 10, 20, 32.5, 50, 72.5),
        \(x) ggplot2::ggplot(
            data = plot_data[plot_data$exp_angle == x, ],
            ggplot2::aes(
                x = amplitude_speed, 
                fill = factor(exp_speed),
                color = factor(exp_speed)
            )
        ) +
            ggplot2::geom_density(
                alpha = 0.2
            ) +
            ggplot2::geom_vline(
                xintercept = 0, 
                color = "black",
                linetype = "dotted",
                linewidth = 2
            ) +
            ggplot2::labs(
                x = "Amplitude", 
                y = "Density",
                title = x,
                fill = "Speed:  "
            ) + 
            ggplot2::lims(
                x = limits
            ) +
            ggplot2::scale_fill_manual(
                labels = c(
                    "1" = "Low", 
                    "2" = "Medium",
                    "3" = "High"
                ),
                values = c(
                    "1" = "cornflowerblue", 
                    "2" = "salmon",
                    "3" = "goldenrod"
                )
            ) +
            ggplot2::scale_color_manual(
                values = c(
                    "1" = "cornflowerblue", 
                    "2" = "salmon",
                    "3" = "goldenrod"
                )
            ) +
            ggplot2::guides(color = "none") +
            ggplot2::theme(
                panel.background = ggplot2::element_rect(
                    fill = "white"
                ),
                panel.border = ggplot2::element_rect(
                    fill = NA,
                    color = "black",
                    linewidth = 1.5
                ),
                panel.grid.major = ggplot2::element_line(
                    color = "gray75"
                ),
                panel.grid.minor = ggplot2::element_blank(),
                plot.title = ggplot2::element_text(
                    hjust = 0.5, 
                    size = 30
                ),
                axis.text.y = ggplot2::element_text(size = 8),
                axis.title = ggplot2::element_text(size = 25),
                legend.title = ggplot2::element_text(size = 25),
                legend.text = ggplot2::element_text(size = 20)
            )
    )

    plotlist <- append(plotlist, plt)

    plt <- ggpubr::ggarrange(
        plotlist = plt, 
        nrow = 1,
        common.legend = TRUE, 
        legend = "bottom"
    )

    ggplot2::ggsave(
        file.path("figures", "study 4", paste0("amplitude ", i, ".png")),
        plt, 
        width = 10000, 
        height = 1900, 
        unit = "px"
    )
}

plt <- ggpubr::ggarrange(
    plotlist = plotlist, 
    nrow = 2, 
    ncol = 6,
    labels = c("A", rep(" ", 5), "B", rep(" ", 5)),
    font.label = list(size = 30),
    common.legend = TRUE, 
    legend = "bottom"
)

ggplot2::ggsave(
    file.path("figures", "study 4", "amplitude.png"),
    plt, 
    width = 10000, 
    height = 3800, 
    unit = "px"
)
