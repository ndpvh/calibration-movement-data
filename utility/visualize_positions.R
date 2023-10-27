#' Visualize measured positions of a given experiment
#' 
#' This is a very basic plotting function in which the different measured 
#' positions of a given experiment are plotted. This is useful when one wants 
#' to have an initial idea of what the data look like, but cannot provide more
#' information that this initial idea. For this, more specialized functions are
#' needed.
#' 
#' @param data A dataframe with measured positions and an experimental name.
#' @param exp A string or vector of strings defining the name of the experiment 
#' in which one is interested. If this is `NULL`, it will visualize all of the 
#' measured positions across all of the experiments.
#' @param exp_col A string defining the name of the column in which to find the
#' names of the experiments. Defaults to "experiment".
#' @param x_col A string defining the name of the column in which to find the 
#' measured positions in the x-axis. Defaults to "x".
#' @param y_col A string defining the name of the column in which to find the 
#' measured positions in the y-axis. Defaults to "y".
#' @param ... Arguments to provide to `geom_point`.
#' 
#' @export
visualize_positions <- function(data, 
                                exp = NULL,
                                exp_col = "experiment", 
                                x_col = "x", 
                                y_col = "y",
                                ...){
    # If `exp` is undefined, then it should contain all of the experiments
    suppressWarnings(if(is.null(exp)){
        exp <- data[,exp_col] %>% 
            unique() %>% 
            unlist()
    })
    print(exp)
    
    # Use ggplot to plot the data of interest  
    data %>% 
        # Rename the columns to the names that are expected
        setnames(old = c(exp_col, x_col, y_col), 
                 new = c("experiment", "x", "y")) %>% 
        # Only retrieve the data of a given experiment
        filter(experiment %in% exp) %>% 
        # And finally plot the data
        ggplot(aes(x = x, y = y)) +
            geom_point(...) +
            theme_bw()
}
