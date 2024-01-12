#' Preprocess the measurements
#' 
#' This function will take in the raw measurements and then convert them to a 
#' useable, preprocessed datafile. The preprocessing consists of the following 
#' steps. First, timestamps are converted to a vector of durations, starting 
#' with 0 as the first measurement occasion. Then, time is binned in steps of 
#' 500msec and the median of the measurements within such a bin is taken as the 
#' measured position for that bin. Finally, distortions in the measured positions
#' are corrected for through a polynomial of the 10th degree. 
#' 
#' Importantly, each of these steps is done for each tag and experiment separately. 
#' It is furthermore useful to note that the data that come into this function are 
#' the ones that are read in by the `load_data` function. In that function, some 
#' rudementary preprocessing is also done, but this mostly pertains to binding 
#' dataframes together.
#' 
#' This function depends on several utility functions, namely 
#' `convert_to_millisecond` and `correct_distortion`. It furthermore depends 
#' on the dataframe format that is outputted by `load_data`.
#' 
#' @param x Dataframe that contains timestamps under column `timestamp` and the 
#' measures x- and y-positions under `x` and `y` respectively. Additionally, the 
#' tag number should be contained under `tag` and the experiment name under
#' `experiment`. As mentioned in the description, the dataframe that is loaded 
#' through `load_data` conforms to all these requirements.
#' @param anchor_positions Dataframe that contains the x- and y-positions of the 
#' anchors under the columns `x` and `y`. These are used to correct the
#' measurements in the dataframe for the systematic distortions. If this is not 
#' provided, the measurements are corrected for through the use of the 
#' measured positions, which is less ideal.
#' 
#' @export 
preprocess_measurements <- function(x, anchor_positions = NULL){
    # Create a duration column that starts at 0 for each of the different tags 
    # and experiments
    x <- x %>% 
        mutate(duration = convert_to_millisecond(timestamp)) %>% 
        group_by(experiment, tag) %>% 
        mutate(duration = duration - min(duration)) %>% 
        ungroup()
        
    # Bin the data in bins of 500msec and take the median measurements of x and 
    # y as the measured position in that bin. Median is chosen for increased 
    # robustness against outliers, which happen sometimes (see examples in the 
    # created figures).
    x <- x %>% 
        # Create an integer that denotes a time bin. This integer is created by 
        # dividing the duration in msec by 500 msec and then rounding up to the 
        # nearest integer. To make sure the duration for 0 is also contained 
        # within the correct time bin (1), we add a trivial number to the 
        # durations when they are 0, making it so that they round towards the 
        # integer 1.
        mutate(duration = ifelse(duration == 0, duration + 0.0001, duration)) %>%
        mutate(time_bin = ceiling(duration / 500)) %>% 
        # Group by time bins, tags, and experiments and take the median position
        # within this time bin. Also get a measure for the number of measurements 
        # on which the medians are based. This should give us an idea of the
        # sampling rate per tag, as well as the robustness of our result.
        group_by(experiment, tag, time_bin) %>% 
        mutate(median_x = median(x), 
               median_y = median(y), 
               measurements_per_bin = length(x)) %>% 
        ungroup() %>% 
        # Keep only those variables that you want to keep in the data
        group_by(experiment, tag, time_bin, 
                 median_x, median_y, measurements_per_bin) %>% 
        tidyr::nest() %>% 
        select(-data) %>% 
        ungroup() %>% 
        # Rename the position columns to `x` and `y` again for ease 
        rename(x = median_x, 
               y = median_y)
    
    # Correct the distortion in the measurements using parameters that were 
    # estimated on the October 14th calibration data. This is done in several 
    # steps:
    #   - Load in these parameters
    #   - Determine the positions of the anchors relative to which the measured
    #     positions are standardized. If not provided, use the minimal and 
    #     maximal values of the data itself as an alternative
    #   - Do the correction for the distortion
    #   - Transform the data back to their previous range
    #
    # The calibration data of October 14th is used as it contains more 
    # measurements than the data of October 21st, increasing the former's power.
    # Might change with more calibration experiments happening.
    parameters <- readRDS(file.path("results", 
                                    "calibration",
                                    "stationary",
                                    "parameters_polynomial_14-10-2023.Rds"))

    if(is.null(anchor_positions)){
        tmp <- x
    } else {
        tmp <- anchor_positions
    }
    xlim <- c(min(tmp$x), max(tmp$x))
    ylim <- c(min(tmp$y), max(tmp$y))

    corrected <- x %>% 
        # Select `x` and `y` and standardize them (parameters are based on 
        # standardized positions)
        select(x, y) %>% 
        mutate(x = minmax_standardize(x, min_x = xlim[1], max_x = xlim[2]), 
               y = minmax_standardize(y, min_x = ylim[1], max_x = ylim[2])) %>% 
        # Correct the distortions using the estimated parameters
        correct_distortion(parameters) %>% 
        # Transform to a dataframe for ease of manipulation 
        as.data.frame() %>% 
        setNames(c("x", "y"))

    x <- x %>% 
        # Change the values of x and y to the corrected values, and get rid of 
        # the standardization (get them back to their absolute units)
        mutate(x = minmax_destandardize(corrected$x, min_x = xlim[1], max_x = xlim[2]), 
               y = minmax_destandardize(corrected$y, min_x = ylim[1], max_x = ylim[2]))

    # Impute NAs in time bins that are empty. This is done to make it very
    # explicit whenever a given tag does not have a measurement in a time bin, 
    # provoking explicit thought about how to handle such missing data.
    x <- x %>% 
        group_by(experiment) %>% 
        tidyr::nest() %>% 
        mutate(data = map(data, 
                          function(x) impute_missings(x))) %>%
        tidyr::unnest(data)

    return(x)
}

# Create a function that will impute missing values based on whether a time bin 
# is present or not in the data. The assumption here is that each tag should 
# have an equal number of measurements within a given experiment. Based on this
# assumption, we can then couple missing time bins for a given tag to an NA.
impute_missings <- function(x){
    # Convert x to a dataframe if it isn't already one
    if(!(is.data.frame(x))){
        x <- x %>% 
            as.data.frame()
    }

    # Get the maximal time bin, communicating how many measurements there 
    # should be for each tag
    max_bin <- x %>%
        select(time_bin) %>% 
        unlist() %>% 
        as.numeric() %>% 
        max()

    # Given this maximal time bin, we can generate a numeric containing all 
    # time bins that should be present for each of the tags.
    bins <- 1:max_bin 

    # Create a dataframe that contains only NAs for `x` and `y`, 0's for 
    # `measurements_per_bin`, and the tag numbers and bin numbers. This 
    # dataframe will be used to match agains the actual dataframe, replacing all 
    # nonvalues with the actual values if they can be found in the original
    # dataframe.
    tags <- unique(x$tag)
    null_frame <- cbind(rep(bins, each = length(tags)), 
                        rep(tags, times = length(bins))) %>% 
        as.data.frame() %>% 
        setNames(c("time_bin", "tag"))

    # Join the `null_frame` together with the actual dataframe `x`. If a given 
    # combination of time bin and tag is not found in `x`, the result will give 
    # a missing value at those combinations. Afterwards, return this dataframe
    null_frame %>% 
        plyr::join(x, by = c("tag", "time_bin")) %>% 
        # Replace the NAs in `measurements_per_bin` with 0, so that it is clear 
        # that you don't have any measurements in this period
        mutate(measurements_per_bin = ifelse(is.na(measurements_per_bin), 
                                             0, 
                                             measurements_per_bin)) %>% 
        return()
}
