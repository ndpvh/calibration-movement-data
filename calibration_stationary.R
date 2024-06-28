################################################################################
# Purpose: Analyze the preprocessed stationary data gathered on 14/10/2023     #
#          and 21/10/2023. The main goal of this analysis is to inform us on   #
#          ways to deal with systematic and unsystematic noise in our data.    #
#          Consequently, the results from this analysis were used to create    #
#          the `filter_measurements` function in the "utility" folder.         #
################################################################################

################################################################################
# PRELIMINARIES

devtools::load_all()

# Get the data that you want to preprocess.
data_files <- c("stationary_14-10-2023", 
                "stationary_21-10-2023",
                "stationary_22-12-2023")
data_list <- lapply(data_files, 
                    \(x) readRDS(file.path("data", "stationary", paste0(x, ".Rds"))))
names(data_list) <- data_files

# Make a distinction between 6 and 4 anchor data
for(i in c(4, 6)) {
    idx <- paste0(data_files[3], "_", i)
    data_list[[idx]] <- data_list[[3]] %>% 
        dplyr::filter(anchors == i)
}

# Adjust the data_files vector
data_files <- names(data_list)





################################################################################
# UNSYSTEMATIC ERROR

# PURPOSE: Estimate the error covariance of the stationary data. For this, we 
#          use a nonparametric bootstrap of the standard deviation around the 
#          central measured position (mean, mode, or median).
#
#          We have two variations: One for the overall variation across tags, 
#          and one that is tag-specific. For the second one, we will check whether
#          the error relates to the position of the tag in space.
#
#          Importantly, the analysis makes the assumption of independence of 
#          error across time. This assumption is tested beforehand.





#------------------------------------------------------------------------------#
# Assumption of time-independence

# Create a function that will compute the autocorrelation for a given variable
autocorr <- function(x) {
    # If too few data points, skip
    if(nrow(x) < 3) {
        return(NA)
    }

    # Arrange data according to time
    x <- x %>% 
        dplyr::arrange(time) %>% 
        dplyr::select(x) %>% 
        unlist() %>% 
        as.numeric()

    return(cor(x[2:length(x)], x[2:length(x) - 1]))
}

# Test the time-independence of the measured positions per experiment per tag
corrs <- lapply(data_list, 
                \(x) x %>% 
                    dplyr::arrange(experiment, id, time) %>% 
                    dplyr::group_by(experiment, id) %>% 
                    tidyr::nest() %>% 
                    dplyr::mutate(auto_x = purrr::map(data, 
                                                      \(X) X %>% 
                                                          dplyr::select(time, x) %>% 
                                                          setNames(c("time", "x")) %>% 
                                                          autocorr()), 
                                  auto_y = purrr::map(data, 
                                                      \(X) X %>% 
                                                          dplyr::select(time, y) %>% 
                                                          setNames(c("time", "x")) %>% 
                                                          autocorr())) %>% 
                    dplyr::select(-data))
names(correlations) <- data_files

# Results: 
#   There is quite a high autocorrelation, meaning that errors are not 
#   independent over time. 
#
# Consequence: 
#   In our estimation of the error covariances, we need to account for this 
#   time dependence
#
# How: 
#   In our approach, we will use the analytic least-squares method to estimate 
#   the mean, autoregressive component, and the covariance matrix of a simple 
#   VAR(1) on the data, defined as:
#
#       \bm{y}_t = \bm{\alpha} + B \bm{y}_{t - 1} + \bm{\epsilon}_t
#       \bm{\epsilon}_t \sim N(\bm{0}, \Sigma).
#
#   We will then take the mean value for the transition and 
#   covariance matrix as an approximation for the actual time-dependence and 
#   random error. This info can then be used to estimate the actual positions 
#   of the tags at time t as: 
#
#       \bm{\mu}_t = \bm{y}_t - B \bm{y}_{t - 1} - \bm{\epsilon}_t
#       \bm{\epsilon}_t \sim N(\bm{0}, \Sigma),
#
#   where the values of \bm{\mu}_t are estimated using the min-log-likelihood.
#   Note that we use the data to estimate an equal amount of parameters. This 
#   may not be optimal, but that's left to practice to figure out. Furthermore 
#   note that while we use the same data sets to estimate the parameters of the 
#   VAR(1) and the values of \bm{\mu}_t, in the actual data we won't do this, 
#   but rather use the values that we get from these stationary data to estimate
#   \bm{\mu}_t

# Let's create a function to estimate the parameters of the VAR(1)
autoregression <- function(x) {
    # Arrange the variables according to time
    x <- x %>% 
        dplyr::arrange(time)

    # Prepare the variables
    Y <- x %>% 
        dplyr::mutate(x = ifelse(time == max(time), NA, x), 
                      y = ifelse(time == max(time), NA, y)) %>% 
        dplyr::filter(!is.na(x)) %>% 
        dplyr::select(x, y) %>% 
        as.matrix() %>% 
        t()

    X <- x %>% 
        dplyr::mutate(x = ifelse(time == min(time), NA, x), 
                      y = ifelse(time == min(time), NA, y),
                      intercept = 1) %>% 
        dplyr::filter(!is.na(x)) %>% 
        dplyr::select(intercept, x, y) %>% 
        as.matrix() %>% 
        t()

    # Do least-squares
    B <-  Y %*% t(X) %*% solve(X %*% t(X))

    # Compute the residuals of the model and use them to estimate the covariance
    # matrix.
    e <- Y - B %*% X
    S <- cov(t(e))

    # Extract the other parameters and return in a vector
    return(c(as.vector(B[,2:3]), as.vector(S)))
}

# Apply the function to each of the data sets
results <- lapply(data_list, 
                  \(x) x %>% 
                      dplyr::group_by(experiment, id) %>% 
                      dplyr::summarize(experiment = experiment[1], 
                                       id = id[1],
                                       data = cbind(x, y, time) %>% 
                                           as.data.frame() %>% 
                                           setNames(c("x", "y", "time")) %>% 
                                           autoregression() %>% 
                                           t() %>% 
                                           as.data.frame()) %>% 
                      tidyr::unnest(data) %>% 
                      dplyr::ungroup() %>% 
                      setNames(c("experiment", "id", 
                                 "auto_x", "cross_yx", "cross_xy", "auto_y", 
                                 "var_x", "cov_yx", "cov_xy", "var_y")))

# Get the overall quantiles 0.025, 0.50, and 0.975 for each of the parameters
results <- do.call("rbind", results) %>% 
    dplyr::select(auto_x:var_y) %>% 
    as.matrix() %>% 
    matrixStats::colQuantiles(probs = c(0.025, 0.5, 0.975))

saveRDS(results, 
        file.path("results", "stationary", "var_params.Rds"))

# As a small (quick) test, check whether the dynamic filtering approach would 
# work on these stationary data (recovery done in tests, but gives some
# weird values in the simulated data: This thus serves as a sanity check)
n_cores <- parallel::detectCores()
# results <- parallel::mclapply(data_list, 
results <- lapply(list(data_list[[1]]),
                              \(x) dynamic_filter(x, 
                                                  maxeval = 1e3, 
                                                  print_level = 1) %>% 
                                                #   itermax = 1e3,
                                                #   NP = 2000,
                                                #   trace = TRUE) %>% 
                                  suppressWarnings())#,
                            #   mc.cores = n_cores - 1)

# Currently seems to be infeasible due to time and computational constraints, 
# unfortunately. Way to go seems to be through simple reduction of measurement
# error. Maybe you could try something else than DEoptim?
#
# -> At this moment, filtering done on the measurement level, but maybe I should 
#    only filter the residuals. Might give different results. Flagged this 
#    because values of y_hat are weird in the optimizer





#------------------------------------------------------------------------------#
# Overall measurement error

# Create a function that will create bootstrapped data in an efficient way. 
# Assumption here is that x is a dataframe that contains, among other, the 
# x and y coordinates that we want to bootstrap.
#
# Importantly, it only bootstraps rows: The relationship between x and y remains
# untouched by this function. This will allow us to not only estimate the 
# variance components for both dimensions, but also the covariance between them.
bootstrapped_covariance <- function(x, 
                                    iterations,
                                    vectorized_iterations = 100) {

    # Get the sample size of the data. Needed to ensure that each of the samples
    # has an equal size to the actual data
    N <- nrow(x)

    # Determine how many times you will have to run the `vectorized_iterations` 
    # to attain the `iterations`
    whole_number <- floor(iterations / vectorized_iterations) 
    iters <- c(rep(vectorized_iterations, each = whole_number), 
               iterations %% vectorized_iterations)

    # Remove iterations that are equal to 0 (only the case if 
    # vectorized_iterations) is a diviser of iterations
    iters <- iters[iters != 0]

    # Do a mix of vectorized and unvectorized bootstrapping to spare your system's 
    # memory.
    results <- list() ; f <- 1
    for(i in seq_along(iters)) {
        # Sample a number of indices for x that is equal to the sample size times 
        # the number of samples one wants to draw
        idx <- sample(1:N, 
                      N * iters[i], 
                      replace = TRUE)

        # Extend the dataframe to account for these values and bind them with an 
        # identity number that conveys the sample they are in
        results[[i]] <- x[idx,] %>% 
            dplyr::mutate(sample_id = rep(f:(f + iters[i] - 1), each = N)) %>% 
            # Compute the covariances based on the corrected x- and y-positions. 
            # Importantly, this is done for each separate bootstrapped sample.
            dplyr::group_by(sample_id) %>% 
            dplyr::mutate(var_x = var(x), 
                          var_y = var(y), 
                          cov_xy = cov(x, y)) %>% 
            dplyr::ungroup() %>% 
            # Delete all other information: Only keep variances and covariance
            dplyr::group_by(sample_id, var_x, var_y, cov_xy) %>% 
            tidyr::nest() %>% 
            dplyr::select(-data) %>% 
            dplyr::ungroup()

        f <- f + iters[i]
    } 
    return(do.call("rbind", results))
}

# Loop over the different dates and do all your estimation 
set.seed(39) # The Messenger - Thrice
results <- list()
for(i in seq_along(data_list)){
    print(data_files[i])

    # Load data from specific date and change it somewhat
    data <- data_list[[i]] %>% 
        # Compute several center statistics, namely mean, median, and mode per
        # tag
        dplyr::group_by(id) %>%         
        dplyr::mutate(mu_x = mean(x),
                      mu_y = mean(y),
                      median_x = median(x), 
                      median_y = median(y), 
                      mode_x = mean(modeest::mlv(x, method = "mfv")), 
                      mode_y = mean(modeest::mlv(y, method = "mfv"))) %>% 
        dplyr::ungroup() %>% 
        # Create a "corrected" version of the x- and y-positions using the mean.
        # This one will be used in the estimation of the error covariances.
        #
        # Can be replaced with median or mode too, if desired.
        dplyr::mutate(x = x - mu_x, 
                      y = y - mu_y)

    # Bootstrap the data using this function and immediately compute the necessary
    # summary statistics: 2 variances and 1 covariance.
    covariances <- bootstrapped_covariance(data, 10000, vectorized_iterations = 250) 

    # Create summary statistics for each of the covariances, and more specifically 
    # given quantiles of the bootstrapped distribution. This should give us an 
    # idea of how badly off we are.
    #
    # Here, we use a 99% CI, just to be sure
    results[[data_files[i]]] <- covariances %>% 
        # Summarize the different variables into CI and mean
        dplyr::summarize(lb_var_x = quantile(var_x, 0.005), 
                         lb_var_y = quantile(var_y, 0.005),
                         lb_cov_xy = quantile(cov_xy, 0.005), 
                         m_var_x = mean(var_x), 
                         m_var_y = mean(var_y), 
                         m_cov_xy = mean(cov_xy), 
                         ub_var_x = quantile(var_x, 0.995), 
                         ub_var_y = quantile(var_y, 0.995),
                         ub_cov_xy = quantile(cov_xy, 0.995)) %>% 
        # Restructure the dataframe to be more useful: Put each of the 
        # covariances as a row, and the lower bounds, mean, and upper bounds 
        # as the columns 
        unlist() %>% 
        matrix(nrow = 3, ncol = 3) %>% 
        as.data.frame() %>% 
        setNames(c("lb", "mean", "ub")) %>% 
        cbind(covariance = c("var_x", "var_y", "cov_xy")) 

    # Release the memory that is held up by the bootstrapped data and the data 
    # itself.
    rm(data, covariances)
}

# Save the results
saveRDS(results, 
        file.path("results", "stationary", "unsystematic_error.Rds"))

# Interpretation of the results: 
#   - Covariances between x and y are as good as 0, so not relationship in the 
#     error between both dimensions
#   - Variances are about: 
#         - 14-10-2023: x: 99%CI = [0.00020, 0.00022], mean = 0.00021
#                       y: 99%CI = [0.00018, 0.00020], mean = 0.00019
#         - 21-10-2023: x: 99%CI = [0.00046, 0.00050], mean = 0.00048
#                       y: 99%CI = [0.00031, 0.00032], mean = 0.00031
#         - 22-12-2023: x: 99%CI = [0.00734, 0.00762], mean = 0.00748
#                       y: 99%CI = [0.00804, 0.00824], mean = 0.00814
#   - Variances are lowest for the first calibration period (14-10-2023), twice 
#     as high for the second calibration period (21-10-2023), and highest in the 
#     third calibration period (22-12-2023). This has some consequences:
#         - There may be more error on the sides than in the center, as the 
#           second calibration period only measured the sides
#         - Either the precision with which you measure depends on the room, or 
#           something went wrong in the third calibration session
#   - In the worst case scenario -- which is calibration period 3, upper bound 
#     of the 99%CI -- there is about 8.73cm of standard error in the x-direction, 
#     and 9.08cm of standard error in the y-direction 





#------------------------------------------------------------------------------#
# Measurement error per tag

# Loop over the different dates and do all your estimation, but this time 
# dispatching on the id
set.seed(7244) # Falling on Deaf Ears - Hail the Sun
results <- list()
for(i in seq_along(data_list)){
    print(data_files[i])

    # Load data from specific date and change it somewhat
    data <- data_list[[i]] %>% 
        # Compute several center statistics, namely mean, median, and mode per
        # tag
        dplyr::group_by(id) %>%         
        dplyr::mutate(mu_x = mean(x),
                      mu_y = mean(y),
                      median_x = median(x), 
                      median_y = median(y), 
                      mode_x = mean(modeest::mlv(x, method = "mfv")), 
                      mode_y = mean(modeest::mlv(y, method = "mfv"))) %>% 
        dplyr::ungroup() %>% 
        # Create a "corrected" version of the x- and y-positions using the mean.
        # This one will be used in the estimation of the error covariances.
        #
        # Can be replaced with median or mode too, if desired.
        dplyr::mutate(x = x - mu_x, 
                      y = y - mu_y)

    # Bootstrap the data using this function and immediately compute the necessary
    # summary statistics: 2 variances and 1 covariance.
    covariances <- data %>% 
        dplyr::group_by(experiment, id) %>% 
        tidyr::nest() %>% 
        dplyr::mutate(covariance = bootstrapped_covariance(as.data.frame(data), 
                                                           10000, 
                                                           vectorized_iterations = 10000) %>% 
                          list()) %>% 
        dplyr::select(-data) %>% 
        tidyr::unnest(covariance) %>% 
        dplyr::ungroup() 

    # Create summary statistics for each of the covariances, and more specifically 
    # given quantiles of the bootstrapped distribution. This should give us an 
    # idea of how badly off we are.
    #
    # Here, we use a 99% CI, just to be sure
    results[[data_files[i]]] <- covariances %>% 
        # Summarize the different variables into CI and mean
        dplyr::group_by(experiment, id) %>% 
        dplyr::summarize(lb_var_x = quantile(var_x, 0.005), 
                         lb_var_y = quantile(var_y, 0.005),
                         lb_cov_xy = quantile(cov_xy, 0.005), 
                         m_var_x = mean(var_x), 
                         m_var_y = mean(var_y), 
                         m_cov_xy = mean(cov_xy), 
                         ub_var_x = quantile(var_x, 0.995), 
                         ub_var_y = quantile(var_y, 0.995),
                         ub_cov_xy = quantile(cov_xy, 0.995)) %>% 
        dplyr::ungroup()

    # Release the memory that is held up by the bootstrapped data and the data 
    # itself.
    rm(data, covariances)
}

# Add the X and Y locations to the dataframes. Needed to be able to infer when 
# standard errors are the worst
results <- lapply(seq_along(data_list), 
                  \(x) data_list[[x]] %>% 
                      dplyr::group_by(experiment, id) %>% 
                      dplyr::summarize(X = mean(X), 
                                       Y = mean(Y), 
                                       x = mean(x), 
                                       y = mean(y)) %>% 
                      dplyr::inner_join(results[[x]]))
names(results) <- data_files

# Save the results
saveRDS(results, 
        file.path("results", "stationary", "unsystematic_error_per_tag.Rds"))

# Visualize the location of the coordinates together with their 99%CI based on 
# the mean error variance
error_plot <- function(x, 
                       title = "",
                       linewidth = 1,
                       ...) {
    # Create polygons that define a circle with center equal to the measured 
    # position and a radius equal to the distance from the mean to the bounds 
    # of a 99%CI (computed with the mean error standard deviation)
    compute_normal <- function(vx, vy, cv, mx, my) {
        co <- cbind(rep(seq(mx - 5 * sqrt(vx), mx + 5 * sqrt(vx), length.out = 250), 
                        each = 250),
                    rep(seq(my - 5 * sqrt(vy), my + 5 * sqrt(vy), length.out = 250), 
                        times = 250))
        
        S <- matrix(c(vx, cv, cv, vy), 
                    nrow = 2, 
                    ncol = 2)

        return(cbind(co, 
                     mvtnorm::dmvnorm(co,
                                      mean = c(mx, my),
                                      sigma = S)))
    }

    plot_data <- x %>% 
        dplyr::group_by(experiment, id) %>% 
        dplyr::mutate(data = compute_normal(m_var_x, 
                                            m_var_y, 
                                            m_cov_xy, 
                                            x, 
                                            y) %>% 
                          as.data.frame() %>% 
                          setNames(c("x", "y", "z")) %>% 
                          dplyr::mutate(z = z / sum(z)) %>% 
                          list()) %>% 
        dplyr::select(experiment, id, data) %>% 
        tidyr::unnest(data) %>% 
        dplyr::mutate(experiment_id = paste0(experiment, "_", id)) %>% 
        dplyr::ungroup()

    # Create the plot itself. For some reason, it doesn't want to make the plot
    # by using geom_contour (most probably because I don't have continuous 
    # values for z across the whole range of x and y). Therefore all contours 
    # added in an annotate
    expid <- unique(plot_data$experiment_id)

    plt <- ggplot2::ggplot()
    for(i in expid) {
        tmp <- dplyr::filter(plot_data, experiment_id == i)
        plt <- plt + 
            ggplot2::stat_contour(data = tmp, 
                                  ggplot2::aes(x = x, y = y, z = z),
                                  color = "black",
                                  linewidth = linewidth,
                                  alpha = 1)
    }
    plt <- plt +
        ggplot2::coord_equal() +
        ggplot2::labs(x = "x", 
                      y = "y", 
                      title = title) +
        ggplot2::theme_minimal() +
        ggplot2::theme(...)

    return(plt)
}

plots <- lapply(data_files, 
                \(x) error_plot(results[[x]], 
                                title = x, 
                                linewidth = 0.01,
                                axis.title = ggplot2::element_text(size = 20),
                                plot.title = ggplot2::element_text(size = 30,
                                                                   hjust = 0.5),
                                panel.border = ggplot2::element_rect(color = "black", 
                                                                     fill = NA,
                                                                     linewidth = 2)))

ggplot2::ggsave(file.path("figures", "stationary", "unsystematic_error_per_tag.jpg"), 
                ggpubr::ggarrange(plotlist = plots,
                                  ncol = 1),
                width = 5000, 
                height = 20000, 
                unit = "px",
                limitsize = FALSE)

# Let's do an additionaly analysis: Examine the relationship between the distance
# of a tag to the center and the error that we observe
results <- lapply(results, 
                  \(x) x %>% 
                      dplyr::mutate(distance_idealized = sqrt((X - mean(X))^2 + (Y - mean(Y))^2), 
                                    distance_idealized_x = abs(X - mean(X)),
                                    distance_idealized_y = abs(Y - mean(Y)),
                                    distance_measured = sqrt((x - mean(x))^2 + (y - mean(y))^2), 
                                    distance_measured_x = abs(x - mean(x)),
                                    distance_measured_y = abs(y - mean(y))) %>% 
                      dplyr::rowwise() %>% 
                      dplyr::mutate(lb_var = mean(c(lb_var_x, lb_var_y)), 
                                    m_var = mean(c(m_var_x, m_var_y)),
                                    ub_var = mean(c(ub_var_x, ub_var_y))) %>% 
                      dplyr::ungroup())

# Save the results
saveRDS(results, 
        file.path("results", "stationary", "unsystematic_error_per_tag.Rds"))

# Create plots of the relationship between each
dists <- c("distance_idealized", "distance_idealized_x", "distance_idealized_y", 
           "distance_measured", "distance_measured_x", "distance_measured_y")
spread <- c("m_var", "m_var_x", "m_var_y", 
            "m_var", "m_var_x", "m_var_y")

plots <- list() ; f <- 1
for(i in data_files) {
    for(j in seq_along(dists)) {
        plot_data <- results[[i]] %>% 
            dplyr::select(dists[j], spread[j]) %>% 
            setNames(c("x", "y"))
        plots[[f]] <- ggplot2::ggplot(data = plot_data, 
                                      ggplot2::aes(x = x, 
                                                   y = y)) +
            ggplot2::geom_point(size = 4, 
                                shape = 21, 
                                color = "cornflowerblue") +
            ggplot2::labs(title = i,
                          x = dists[j], 
                          y = spread[j])
        f <- f + 1
    }
}

ggplot2::ggsave(file.path("figures", "stationary", "unsystematic_error_distance.jpg"), 
                ggpubr::ggarrange(plotlist = plots, 
                                  ncol = length(dists), 
                                  nrow = length(data_files)), 
                width = 6000, 
                height = 5500, 
                unit = "px")

# Interpretation of the results: 
#   - There is some variation in the unsystematic error per tag
#   - There is evidence that this variation is related to the position of the 
#     tag in the grid: The closer to the edges, the greater the measurement 
#     error (in all datasets), and in the middle there seems to be a problem for
#     the 22-12-2023 data
#       - This only seems to be the case slightly, and primarily in the 
#         y-direction

#-------------------------------------------------------------------------------
# Sampling rate
#-------------------------------------------------------------------------------

# Here, we estimate the attained sampling rate for our stationary data.

# Create a function that will take in the duration, order them according to 
# size, and then compute the mean difference between each of the durations
sampling_rate <- function(x){
    x %>% 
        as.numeric() %>% 
        sort() %>% 
        diff() %>% 
        mean() %>% 
        return()
}

# Again loop over each of the stationary datasets
for(i in seq_along(stationary)){
    # Load the stationary data for a given date and convert the timestamps to 
    # milliseconds
    data <- load_stationary(stationary[i]) %>% 
        mutate(duration = convert_to_millisecond(timestamp))

    # Compute the mean sampling rate for each tag separately and make the 
    # dataframe somewhat easier to interpret
    result <- data %>% 
        # Get the sampling rate in msec
        group_by(tag) %>% 
        mutate(msec = sampling_rate(duration)) %>% 
        ungroup() %>% 
        # Convert the sampling rate to Hz: 1 / s -> 1000 / msec
        mutate(Hz = 1000 / msec) %>% 
        # Only retain those variables that you want to interpret
        group_by(tag, msec, Hz) %>% 
        tidyr::nest() %>% 
        select(-data)

    # Save this result
    save_result(result, 
                "sampling_rate", 
                stationary[i])

    # Make another histogram to visualize the result across tags 
    plt <- ggplot(result, 
                  aes(x = Hz)) +
        geom_histogram(color = "black", 
                       fill = "cornflowerblue") +
        labs(title = "Sampling rate per tag", 
             x = "Hz") +
        geom_vline(xintercept = 5, 
                   color = "red")

    ggsave(file.path("figures", 
                     "calibration", 
                     "stationary",
                     paste0("sampling_rate_", stationary[i], ".png")),
           plot = plt, 
           units = "px", 
           width = 1000, 
           height = 1100)
}

# Interpretation of the result:
#   - On the 14th of October, sampling rates varied substantially between 2 and 
#     4Hz instead of the expected 5Hz. Sampling rates should thus be increased 
#     for our future experiments.
#   - On the 21st of October, we found a higher than expected sampling rate. This 
#     might suggest that our initial attempts of changing the sampling rate from 
#     5Hz to 8Hz on that day might have been successful, leading to an observed 
#     sampling rate of 6-7Hz 
#   - Some tags seem to perform worse and only send out responses every once in 
#     a while.
#
# Additional comments after calibration on 22-12-2023
#   - On this day, sampling with the 6 anchors remained relatively similar to 
#     the sampling rate on 21-10-2023 (although with a slight loss of frequency).
#   - Sampling with 4 anchors had an effect on the sampling frequency. Need to 
#     find out how to increase it again.





#-------------------------------------------------------------------------------
# Systematic distortions
#-------------------------------------------------------------------------------

# Here, we estimate a function to filter out systematic distortions in the data. 
# We do this by estimating a polynomial on the stationary data and checking
# the differences between the functions estimated on both days. We have two 
# major assumption in this analysis, which is: 
#   - The actual measured position is equal to the average of all measured 
#     positions for a given tag
#   - The systematic error in the x- and y-dimensions is independent


# Create a function that will output a formula for a polynomial of the degree^th
# degree with `dependent` as dependent variable and `independent` as the 
# independent variable.
polynomial_formula <- function(independent, dependent, degree){
    # Start with putting the dependent variable against the independent variable
    polyn <- paste(dependent, "~", independent)

    # Add each term of the polynomial to this string. Should be of the format 
    # I(x^degree)
    for(i in 2:degree){
        polyn <- paste(polyn, "+", "I(",
                       paste0(independent, "^", i), ")")
    }

    # Return this string as a formula to be used in `lm`
    polyn %>% 
        formula() %>% 
        return()
}

# Create a formula that will do the actual estimation, based on data, a dependent 
# variable, an independent variable, and a degree
fit_polynomial <- function(data, independent, dependent, degree){
    polynomial_formula(independent, dependent, 10) %>% 
        lm(data = data) %>% 
        summary() %>% 
        return()
}

# Create a function that will extract all of the coefficients
extract_coefficients <- function(x){
    x$coefficients[,1] %>% 
        as.numeric() %>% 
        return()
}

# As one can predict already: Loop over all stationary data again.
for(i in seq_along(stationary)){
    # Load the data for a given date and standardize all data: Measured positions
    # and supposed real positions. Standardization is performed to make the
    # function scalable to other data in which the absolute positions are not 
    # the same as these of the calibration phase.
    #
    # Importantly, the kind of standardization that is performed is not the 
    # standard way of doing this. Rather, we want to transform all measured positions 
    # that fall within the measurement space (i.e., the space bounded by the 
    # anchors) to fall within -1 and 1. This will allow us to more readily 
    # translate the calibration on one day to the measurements on another, as 
    # you are not dependent on the actual measurements for the standardization
    # (in contrast to the scaled equivalent of the Z-score). Importantly, the 
    # `minmax_standardize` function does this transformation, as defined in the 
    # utility functions.
    #
    # Unfortunately, we don't have the anchors' positions, so we will assume that 
    # they lie at the sides of the "real positions" grid. In reality, this will 
    # probably be off by a few tens of centimeters, but we will have to do 
    # another, more precisely done calibration to combat the issues that this 
    # brings.
    data <- load_stationary(stationary[i]) %>% 
        mutate(standardized_x = minmax_standardize(x, 
                                                   min_x = min(X), 
                                                   max_x = max(X)), 
               standardized_y = minmax_standardize(y, 
                                                   min_x = min(Y), 
                                                   max_x = max(Y)), 
               standardized_X = minmax_standardize(X), 
               standardized_Y = minmax_standardize(Y))

    # Estimate the 10th degree polynomial on the standardized data
    pars <- list("x" = fit_polynomial(data, 
                                      "standardized_x", 
                                      "standardized_X", 
                                      10),
                 "y" = fit_polynomial(data, 
                                      "standardized_y", 
                                      "standardized_Y", 
                                      10))

    # Save these results
    save_result(pars, 
                "polynomial", 
                stationary[i])
                 
    # Create a dataframe to be used in `correct_distortion` which will 
    # contain all of the parameters for the separate dimensions
    result <- cbind("x" = extract_coefficients(pars[["x"]]), 
          "y" = extract_coefficients(pars[["y"]])) %>% 
        as.data.frame()

    # Save these parameters as well 
    save_result(result, 
                "parameters_polynomial", 
                stationary[i])

    # Check how well this polynomial is able to correct for the distortion 
    # in the data.
    coordinates <- data %>% 
        ungroup() %>% 
        select(standardized_x, standardized_y) %>% 
        correct_distortion(result)

    # Create a plot of the locations for the distorted and the corrected data.
    #
    # First, create a function for this
    location_plot <- function(data, x, y, title){
        # Create plot data with the string column names `x` and `y`
        plot_data <- data[, c(x, y)] %>% 
            as.data.frame() %>% 
            setNames(c("X", "Y"))

        # Create the plot itself
        plt <- ggplot(plot_data, 
                      aes(x = X, y = Y)) +
            geom_point(size = 2, 
                       color = "black") +
            labs(x = "Standardized X", 
                 y = "Standardized Y", 
                 title = title)

        return(plt)
    }
    
    # Create and bind the plots of the distorted and corrected data
    plt <- ggarrange(location_plot(data, 
                                   "standardized_x", 
                                   "standardized_y", 
                                   "Original"), 
                     location_plot(coordinates, 
                                   "standardized_x", 
                                   "standardized_y", 
                                   "Corrected"),
                     nrow = 1)

    # Save this plot 
    ggsave(file.path("figures", 
                     "calibration", 
                     "stationary",
                     paste0("correct_distortion_", stationary[i], ".png")), 
           plot = plt, 
           units = "px", 
           width = 2000, 
           height = 1200)
}

# Interpretation of the results: 
#   - For both datasets, the distortions seems to disappear somewhat
#   - However, they do still show some distortion of the distances between each 
#     of the measured tags
#   - Question: is this a problem in our measurements, or a problem in where we 
#     put the tags (i.e., is our assumption of the supposed coordinates X and Y 
#     correct, or did we make a human error here)
#   - Parameters between the two datasets are similar, though do not match 
#     exactly. Given that their are more measured positions in the data from 
#     the 14th of October, we will assume these estimates have a higher power, 
#     and will use those for the correction of our distortion.
#
# Additional comments after calibration on 22-12-2023
#   - The 4-anchor data is corrected better than the 6-anchor data, but both 
#     are worse than the 14-10-2023 counterpart.