#' Converts a substring to a numeric
#' 
#' Is used in the `convert_to_millisecond` function by lack of current better 
#' alternative. This function is definitely not meant for more extensive use.
#' 
#' @param x A string of which you want to extract a number
#' @param first An integer denoting the first character in the string that 
#' contains the number
#' @param last An integer denoting the last character in the string that 
#' contains the number
#' 
#' @export
subnumeric <- function(x, first, last){
    substring(x, first, last) %>% 
        as.numeric() %>% 
        return()
}

#' Convert timestamps into milliseconds
#' 
#' Takes in the timestamps of our measured data and converts it to an integer 
#' number that denotes the milliseconds that correspond to this timestamp.
#' 
#' @param x Vector of times. Will be converted to a POSIXlt timestamp.
#' 
#' @export
convert_to_millisecond <- function(x){
    # Only get hours, minutes, seconds, and milliseconds
    times <- x %>% 
        as.POSIXlt(format = "%H:%M:%OS") %>% 
        format("%H:%M:%OS3")

    # Unfortunately, need a dirty hack to be able to work with this further.
    # Assume this format and define hours, minutes, seconds, and milliseconds.
    times %>% 
        (function(x){
            return(3600000 * subnumeric(x, 1, 2) +
                   60000 * subnumeric(x, 4, 5) +
                   1000 * subnumeric(x, 7, 8) +
                   subnumeric(x, 10, 11))
        }) %>% 
        return()
}