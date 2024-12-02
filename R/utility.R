# Imports from other packages
#' @importFrom magrittr %>%
#' @import locfit

# A vectorized sequence function
multi_seq <- Vectorize(seq.default, 
                       vectorize.args = c("from", "to", "by", "length.out"))

# A negated %in% function
`%notin%` <- Negate(`%in%`)
