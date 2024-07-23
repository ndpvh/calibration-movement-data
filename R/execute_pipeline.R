#' Perform a preprocessing pipeline
#' 
#' This function takes in movement data and outputs the resulting data when a 
#' given preprocessing pipeline is used on it. The user will have to specify 
#' which functions they want to use on the data and in what order through the 
#' `fx` argument. The function will then output the resulting data after 
#' performing all functions in `fx`, or the result of preprocessing the data 
#' in the way specified.
#' 
#' The dataframe and the functions applied to the dataframe are free for the 
#' user to use. However, within the scope of this project, I will always make 
#' sure that functions both require a dataframe with columns `time`, `id`, `x`, 
#' and `y` as input and output a dataframe with the same specifications. This way,
#' I will be free to specify whichever preprocessing steps in whichever order.
#' 
#' @param data Dataframe that contains the columns `time`, `id`, `x`, and `y`.
#' @param fx List of functions to execute on the data. Make sure that these 
#' functions take in only a single argument, namely the data. If you want to 
#' specify optional arguments, you will have to do so beforehand.
#' 
#' @return Dataframe of a similar structure as `data` on which each of the 
#' functions in `fx` have been executed
#' 
#' @export 
execute_pipeline <- function(data, 
                             fx) {
    for(i in seq_along(fx)) {
        cat(paste0("\rExecuting function ", i))
        
        data <- fx[[i]](data)
    }
    # cat("\n")
    return(data)
}