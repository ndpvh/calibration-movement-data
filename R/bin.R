#' Bin together data
#' 
#' Bin data together within a given window. While related, it is different from
#' `moving_window` in that the data quantity is reduced. Specifically, while 
#' `moving_window` keeps the same amount of data, `bin` will reduce the amount 
#' of data by taking all data within a bin together and summarizing it to one 
#' value.
#' 
#' A lot of this function is the same as for `moving window`: It only differs in 
#' its discrete jump to another bin instead of going from data point to data 
#' point. 
#' 
#' @param data Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' @param span Numeric denoting the size of the bins. Will pertain to the `time` 
#' variable, and should be specified in seconds (e.g., 0.5 means 500msec).
#' @param fx Function to execute on the data that falls within the bin. Should 
#' take in a similar data structure to `data` and output a single-row dataframe
#' with columns `x` and `y`.
#' @param .by String denoting whether the moving window should be taken with 
#' respect to a given grouping variable. Defaults to `NULL`
#' 
#' @return Binned dataframe with a similar structure as `data`
#' 
#' @export 
bin <- function(data, 
                span, 
                fx,
                .by = NULL) {

    # Check whether there is a grouping variable to account for. If not, then 
    # we just use a dummy
    if(!is.null(.by)) {
        group <- unique(data[,.by])
    } else {
        group <- 1
    }

    # Go over each of the data points and smooth the data using the moving 
    # window. We dispatch/loop over all the different possibilities in the 
    # grouping variable, for each of which we will create an individual 
    # dataframe to be smoothed.
    #
    # Within this loop, we do the following:
    #   - Create a group-specific dataframe
    #   - Determine which time points to contain within each bin based on the 
    #     specified span
    #   - Apply the specified function to the data
    #   - Put these data in a list
    all_data <- list()
    for(i in seq_along(group)) {
        # Get the data specific to the group (if .by is not NULL)
        if(!is.null(.by)) {
            individual_data <- data[data[, .by] == group[i],]
        } else {
            individual_data <- data
        }

        # Add an integer index to the data. This will allow researchers some 
        # flexibility with regard to how to assign weights to the different
        # observations. Furthermore specify the bin number to which each of the 
        # data points will be assigned.
        individual_data <- individual_data %>% 
            dplyr::mutate(index = dplyr::row_number(),
                          bin_number = ceiling(time / span)) %>% 
            # If time starts at 0, we want this data point to be added to the 
            # first bin
            dplyr::mutate(bin_number = ifelse(bin_number == 0, 1, bin_number))

        # Add nested tibbles that contain the data for the specified span and 
        # apply the specified function to the data
        all_data[[i]] <- individual_data %>% 
            # Group by bin number
            dplyr::group_by(bin_number) %>% 
            tidyr::nest() %>% 
            # Apply the function to the nested data
            dplyr::rowwise() %>% 
            dplyr::mutate(data = data %>% 
                              as.data.frame() %>% 
                              dplyr::mutate(relative_time = time - (bin_number - 0.5) * span,
                                            index = index - mean(index)) %>% 
                              fx()) %>% 
            tidyr::unnest(data) %>% 
            dplyr::ungroup() %>% 
            # Add the variables time and id to the mix again. As time, take the 
            # middle of the bin
            dplyr::rowwise() %>% 
            dplyr::mutate(time = (bin_number - 0.5) * span, 
                          id = group[i]) %>% 
            dplyr::ungroup()
    }

    # Bind all data together and sort based on time
    data <- do.call("rbind", all_data) %>% 
        dplyr::relocate(time:id, x:y) %>% 
        dplyr::arrange(time) %>% 
        as.data.frame()    

    return(data)
}