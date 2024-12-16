# Imports from other packages
#' @importFrom magrittr %>%
#' @import locfit

# A vectorized sequence function
multi_seq <- Vectorize(seq.default, 
                       vectorize.args = c("from", "to", "by", "length.out"))

# A negated %in% function
`%notin%` <- Negate(`%in%`)

#' Utility function to add metadata to dataframe
#' 
#' @param data Data.frame to add the metadata to.
#' @param metadata Named list containing column name and value to assign.
#' 
#' @return Data.frame containing the metadata.
#' 
#' @export 
add_metadata <- function(data,
                         metadata) {

    columns <- names(metadata)
    for(i in columns) {
        data[, i] <- metadata[[i]]
    }

    return(data)
}