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

# Load the experiment names and join together with the bigger datafile. This 
# will make it easier for us to select on specific dates on which experiments 
# were ran, as well as on different conditions for each of these experiments.
data <- data.table::fread(
    file.path("data", "raw_data", "experiments.csv"),
    data.table = FALSE
) %>% 
    dplyr::rename(experiment_id = id) %>% 
    dplyr::select(experiment_id, name) %>% 
    dplyr::full_join(data, by = "experiment_id") %>% 
    dplyr::select(-experiment_id) %>% 
    dplyr::rename(experiment = name)

# We will use the data of 14-10-2023 and select only the moving data from them. 
# Then try to get information on the line and voltage that were used for the 
# movement. Finally find out which experiments were retaken and delete their 
# wrong counterpart.
data <- data %>% 
    dplyr::filter(
        stringr::str_detect(timestamp, "2023-10-14"), 
        stringr::str_detect(experiment, "movement")
    ) %>% 
    dplyr::group_by(experiment) %>% 
    tidyr::nest() %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(
        line = stringr::str_split(experiment, " ")[[1]][2],
        voltage = stringr::str_split(experiment, " ")[[1]][3]
    ) %>% 
    dplyr::filter(experiment != "movement 9 3 V real") %>% 
    tidyr::unnest(data) %>% 
    dplyr::ungroup()

# Add information on the anchor locations to these data
anchors <- readRDS(file.path("data", "anchor_position_14-10-2023.Rds"))
data <- data %>%
    dplyr::mutate(
        anchor_xmin = min(anchors[, 2]),
        anchor_xmax = max(anchors[, 2]),
        anchor_ymin = min(anchors[, 3]),
        anchor_ymax = max(anchors[, 3]),
    )

# Add a time-variable to the data_list. Is a transformation of the timestamp 
# provided in the actual datapoints, but now as a numeric in sec
data <- data %>% 
    dplyr::group_by(experiment) %>% 
    dplyr::mutate(
        time = withr::with_options(
            list(digits = 16),
            as.numeric(timestamp)
        ), 
        time = time - min(time)
    ) %>% 
    dplyr::ungroup()

# Finally, add the observed time it took for each sled to reach the end in each 
# condition. Combine this with the information on the distance that the sled 
# travelled and compute observed speed scores from this (distance is measured in
# meter and speed in meters per second). For simplicity, transform the time 
# measure from milliseconds to seconds
data <- data.table::fread(
    file.path("data", "raw_data", "moving_times_14-10-2023.txt"),
    data.table = FALSE
) %>% 
    dplyr::rename(
        line = row_measured,
        voltage = volts
    ) %>% 
    dplyr::mutate(voltage = as.character(voltage)) %>% 
    dplyr::full_join(
        data, 
        by = c("line", "voltage")
    ) %>% 
    dplyr::mutate(
        perp = stringr::str_detect(line, "perp"),
        diag = stringr::str_detect(line, "diag"),
        total_distance = ifelse(
            perp, 
            8,
            ifelse(
                diag, 
                sqrt(10^2 + 8^2),
                10
            )
        ),
        speed_actual = total_distance / seconds
    ) %>% 
    dplyr::select(-perp, -diag)

# Remove data that fall outside of the anchor positions, as we're sure this 
# shouldn't occur. Additionally, remove data that fall outside of the bounds of 
# the measured lines
lines <- data %>% 
    dplyr::group_by(line) %>% 
    dplyr::summarize(
        mean_x = mean(x), 
        mean_y = mean(y), 
        median_x = median(x), 
        median_y = median(y),
        q025_x = quantile(x, 0.025),
        q025_y = quantile(y, 0.025),
        q975_x = quantile(x, 0.975),
        q975_y = quantile(y, 0.975)
    )

data <- data %>% 
    dplyr::filter(
        (x <= anchor_xmax) & (x >= anchor_xmin), 
        (y <= anchor_ymax) & (y >= anchor_ymin)
    ) %>% 
    dplyr::filter(
        (x <= lines$q975_x[lines$line == "11perp"]) & (x >= lines$q025_x[lines$line == "1perp"]), 
        (y <= lines$q975_y[lines$line == "9"]) & (y >= lines$q025_y[lines$line == "1"])
    )

# Unfortunately, the experimental names do not contain all needed information, 
# as there are repetitions of each movement within the data. To snuff these out, 
# we realize that the sled always started at the same position, meaning that 
# we should see a rather discrete jump in either x- or y-coordinates, which 
# would indicate when we have restarted. We just have to find these moments
data <- data %>% 
    dplyr::group_by(experiment) %>% 
    tidyr::nest() %>% 
    dplyr::mutate(
        data = data[[1]] %>% 
            dplyr::arrange(time) %>% 
            dplyr::mutate( 
                diff_x = c(NA, diff(x)),
                diff_y = c(NA, diff(y)),
                sign_x = sign(diff_x),
                sign_y = sign(diff_y)
            ) %>% 
            list()
    ) %>% 
    tidyr::unnest(data) %>% 
    dplyr::ungroup() %>% 
    # Only retain the observations that move in the correct direction: For x, 
    # this is a negative sign, for y, this is a , and for diag, this is a ...
    dplyr::mutate(
        okay_x = ifelse(
            stringr::str_detect(line, "perp"), 
            TRUE,
            ifelse(
                stringr::str_detect(line, "diag"), 
                sign_x == 1, 
                sign_x == -1
            )
        ), 
        okay_y = ifelse(
            stringr::str_detect(line, "perp") | stringr::str_detect(line, "diag") , 
            sign_y == -1,
            TRUE
        )
    ) %>% 
    dplyr::filter(okay_x & okay_y) %>% 
    dplyr::select(-diff_x, -diff_y, -sign_x, -sign_y, -okay_x, -okay_y)

# Delete data for which there is no movement. A small simulation shows that 
# under the assumption that the two dimensions have the same error S and are 
# independent of each other, and when we assume that there is no movement (means
# remain the same over time), then the average distance we get out of it is 
# about lognormally distributed. 
#
# We use this strategy for our observed variance 0.0021025 and take the cutoff 
# point to be the 95% quantile of the resulting distribution (100,000 samples), 
# finding the cutoff point to be 0.1584984. Looking at the median value for the 
# distance at the lowest speed (voltage = 3), we find that this distance is 
# 0.2011218, giving us a sanity check for not removing data that may be 
# potentially interesting.
# data <- data %>% 
#     dplyr::group_by(line, voltage) %>% 
#     tidyr::nest() %>% 
#     dplyr::mutate(
#         data = data[[1]] %>% 
#             dplyr::arrange(time) %>% 
#             dplyr::mutate(
#                 measured_distance = c(
#                     NA, 
#                     (x[2:length(x)] - x[2:length(x) - 1])^2 + (y[2:length(y)] - y[2:length(y) - 1])^2
#                 ), 
#                 measured_distance = sqrt(measured_distance)
#             ) %>% 
#             list()        
#     ) %>% 
#     tidyr::unnest(data) %>% 
#     dplyr::ungroup() %>% 
#     dplyr::filter(measured_distance > 0.1584984) %>% 
#     dplyr::select(-measured_distance)





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
    params <- readRDS(file.path("results", "study 1", paste0("polynomial_14-10-2023.Rds")))
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
        assumed_variance = 0.0021025,
        reverse = FALSE,
        .by = "experiment"
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
        dplyr::group_by(line, voltage) %>% 
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
        dplyr::select(line, voltage, time, x, y, speed_actual, speed, distance) %>% 
        return()
) %>% 
    `names<-` (names(data_list))

# Check some summary statistics: What are the mean speed and the variance around
# this: Ideally, the observed mean lies closer to the real speed, and the 
# observed variance is smaller in the filtered data.
#
# Let's bootstrap these to have a good overview
bootstrap <- function(x, 
                      iterations = 10000, 
                      vectorized_iterations = 100) {
    
    # Get the sample size of the data. Needed to ensure that each of the samples
    # has an equal size to the actual data
    N <- nrow(x)

    # Determine how many times you will have to run the `vectorized_iterations` 
    # to attain the `iterations`
    whole_number <- floor(iterations / vectorized_iterations) 
    iters <- c(
        rep(vectorized_iterations, each = whole_number), 
        iterations %% vectorized_iterations
    )

    # Remove iterations that are equal to 0 (only the case if 
    # vectorized_iterations) is a diviser of iterations
    iters <- iters[iters != 0]

    # Do a mix of vectorized and unvectorized bootstrapping to spare your system's 
    # memory.
    results <- list() ; f <- 1
    for(i in seq_along(iters)) {
        # Sample a number of indices for x that is equal to the sample size times 
        # the number of samples one wants to draw
        idx <- sample(
            1:N, 
            N * iters[i], 
            replace = TRUE
        )

        # Extend the dataframe to account for these values and bind them with an 
        # identity number that conveys the sample they are in
        results[[i]] <- x[idx,] %>% 
            dplyr::mutate(sample_id = rep(f:(f + iters[i] - 1), each = N)) %>% 
            # Compute the covariances based on the corrected x- and y-positions. 
            # Importantly, this is done for each separate bootstrapped sample.
            dplyr::group_by(sample_id) %>% 
            dplyr::summarize(
                speed_mean = mean(speed),
                speed_var = var(speed),
                total_distance = sum(distance),
                speed_actual = speed_actual[1]
            ) %>% 
            dplyr::ungroup() %>% 
            # Delete all other information: Only keep variances and covariance
            dplyr::select(
                sample_id, 
                speed_mean, 
                speed_var, 
                total_distance, 
                speed_actual
            )

        f <- f + iters[i]
    } 
    return(do.call("rbind", results))
}

summarized <- lapply(
    results, 
    \(x) x %>% 
        dplyr::group_by(line, voltage) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(
            data = data[[1]] %>% 
                bootstrap() %>% 
                list()
        ) %>% 
        tidyr::unnest(data) %>% 
        dplyr::summarize(
            m_mean_speed = mean(speed_mean),
            s_mean_speed = sd(speed_mean),
            q025_mean_speed = quantile(speed_mean, 0.025),
            q975_mean_speed = quantile(speed_mean, 0.975),

            m_var_speed = mean(speed_var),
            s_var_speed = sd(speed_var),
            q025_var_speed = quantile(speed_var, 0.025),
            q975_var_speed = quantile(speed_var, 0.975),

            total_distance = mean(total_distance),
            speed_actual = speed_actual[1]
        ) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(
            speed_contained = speed_actual >= q025_mean_speed & speed_actual <= q975_mean_speed
        )
) %>% 
    `names<-` (names(data_list))

saveRDS(
    summarized, 
    file.path("results", "study 4", "moving.Rds")
)

# Let's check the coverage percentage for both types of data
do.call(
    "rbind",
    lapply(
        names(summarized), 
        \(x) summarized[[x]] %>% 
            dplyr::group_by(line) %>% 
            tidyr::nest() %>% 
            dplyr::mutate(
                data = data[[1]] %>% 
                    dplyr::summarize(
                        coverage = sum(speed_contained) / length(speed_contained),

                        m_mean_speed = mean(m_mean_speed),
                        s_mean_speed = mean(s_mean_speed),
                        q025_mean_speed = mean(q025_mean_speed),
                        q975_mean_speed = mean(q975_mean_speed),

                        m_var_speed = mean(m_var_speed),
                        s_var_speed = mean(s_var_speed),
                        q025_var_speed = mean(q025_var_speed),
                        q975_var_speed = mean(q975_var_speed)

                    )
            ) %>% 
            tidyr::unnest(data) %>% 
            dplyr::ungroup() %>% 
            dplyr::mutate(type = x) %>% 
            dplyr::relocate(type)
    )
) %>% 
    View()

do.call(
    "rbind",
    lapply(
        names(summarized), 
        \(x) summarized[[x]] %>% 
            dplyr::group_by(voltage) %>% 
            tidyr::nest() %>% 
            dplyr::mutate(
                data = data[[1]] %>% 
                    dplyr::summarize(
                        coverage = sum(speed_contained) / length(speed_contained),

                        m_mean_speed = mean(m_mean_speed),
                        s_mean_speed = mean(s_mean_speed),
                        q025_mean_speed = mean(q025_mean_speed),
                        q975_mean_speed = mean(q975_mean_speed),

                        m_var_speed = mean(m_var_speed),
                        s_var_speed = mean(s_var_speed),
                        q025_var_speed = mean(q025_var_speed),
                        q975_var_speed = mean(q975_var_speed)
                    )
            ) %>% 
            tidyr::unnest(data) %>% 
            dplyr::ungroup() %>% 
            dplyr::mutate(type = x) %>% 
            dplyr::relocate(type)
    )
) %>% 
    View()

# Let's check the precision with which speed was measured
sapply(summarized, \(x) mean(x$m_var_speed))
sapply(summarized, \(x) mean(x$q975_var_speed - x$q025_var_speed))

# Bootstrap the ability to coverage for the two datasets
coverage <- lapply(
    names(summarized), 
    function(x) {
        data <- summarized[[x]]

        # Get the number of datapoints in the dataset. This is the amount
        # we'll resample
        N <- nrow(data)

        # Sample a number of indices for x that is equal to the sample size times 
        # the number of samples one wants to draw
        idx <- sample(
            1:N, 
            N * 10000, 
            replace = TRUE
        )

        results <- data[idx,] %>% 
            dplyr::mutate(
                sample_id = rep(
                    1:10000, 
                    each = N
                )
            ) %>% 
            # Compute the covariances based on the corrected x- and y-positions. 
            # Importantly, this is done for each separate bootstrapped sample.
            dplyr::group_by(sample_id) %>% 
            dplyr::summarize(
                coverage = sum(speed_contained) / length(speed_contained),
                type = x
            ) %>% 
            dplyr::ungroup() 

        return(results)
    }
)
coverage <- do.call("rbind", coverage) 
coverage <- data.frame(
    sample_id = coverage$sample_id[coverage$type == "filtered"], 
    filtered = coverage$coverage[coverage$type == "filtered"],
    unfiltered = coverage$coverage[coverage$type == "unfiltered"]
) %>% 
    dplyr::mutate(diff = filtered - unfiltered) %>% 
    dplyr::summarize(
        mean_filtered = mean(filtered), 
        sd_filtered = sd(filtered),
        q025_filtered = quantile(filtered, 0.025),
        q975_filtered = quantile(filtered, 0.975),

        mean_unfiltered = mean(unfiltered), 
        sd_unfiltered = sd(unfiltered),
        q025_unfiltered = quantile(unfiltered, 0.025),
        q975_unfiltered = quantile(unfiltered, 0.975),

        mean_diff = mean(diff), 
        sd_diff = sd(diff),
        q025_diff = quantile(diff, 0.025),
        q975_diff = quantile(diff, 0.975)
    ) %>% 
    View()





################################################################################
# VISUALIZATION

plt <- lapply(
    names(data_list),
    function(x) {
        data <- data_list[[x]]

        plt <- ggplot2::ggplot(
            data = data, 
            ggplot2::aes(
                x = x, 
                y = y
            )
        ) +
            ggplot2::geom_point(
                size = 1, 
                color = "black",
                fill = NA,
                shape = 21
            ) +
            ggplot2::labs(
                title = ifelse(x == "filtered", "Filtered", "Unfiltered"),
                x = "x", 
                y = "y"
            ) +
            ggplot2::lims(
                x = c(3.5, 14.25),
                y = c(1, 9.75)
            ) +
            ggplot2::coord_equal() +
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
                legend.position = "none"
            )
        
        return(plt)
    }
)

plt <- ggpubr::ggarrange(
    plotlist = plt, 
    nrow = 1
)

ggplot2::ggsave(
    file.path("figures", "study 4", "moving.png"),
    plt,
    width = 2500,
    height = 1250, 
    unit = "px"
)    
