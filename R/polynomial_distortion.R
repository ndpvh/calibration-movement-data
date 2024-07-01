#' Correct systematic error through polynomial
#' 
#' Use a polynomial of the n^th degree with the provided parameters to correct 
#' the systematric error in the measured locations `x` and `y`.
#' 
#' @param data Dataframe containing the measured locations of x and y.
#' @param params Parameters as outputted by the least-squares procedure (k x 2, 
#' where k is the number of parameters)
#' @param n Degree of the polynomial. Defaults to `1`
#' @param simple Whether to only use main effects of `x` and `y` (`TRUE`) or to 
#' also include interaction effects of varying degrees (`FALSE`). Defaults to 
#' `TRUE`.
#' 
#' @return Dataframe with corrected values of x and y
#' 
#' @export
polynomial_distortion <- function(data,
                                  params = NULL, 
                                  n = 1,
                                  simple = TRUE){

    # Check whether parameters are provided. If not, then change to the default,
    # which is the mean of all estimated parameters across all stationary data
    if(is.null(params)) {
        params <- readRDS(file.path("results", "stationary", "polynomial_params.Rds"))
        params <- Reduce("+", params) / length(params)

        # Adjust the other variables so that X is correctly defined
        n = 5
        simple = FALSE
    }

    # Create the independent variables of the polynomial. Here, we dispatch on 
    # the kind of polynomial you want to create, as there are two options to 
    # move forward.
    #
    # When `simple == TRUE`, we will create a simple polynomial in which the 
    # dimensions `x` and `y` are always treated together, so that: 
    #
    #   [X, Y] = \sum_{j = 0} B_j [x, y]^j
    #
    # Interactions between the dimensions are then captured by the off-diagonal 
    # elements in the 2 x 2 matrix B_j, which is unique for each degree j.
    #
    # Another approach generalizes this polynomial and includes explicit 
    # interaction effects for each degree, so that: 
    #
    #   [X, Y] = \sum_{j = 0, i = 0, i + j \leq n} [\beta_1, \beta_2]_j x^i y^j
    #
    # where you use the parameter vector [\beta_1, \beta_2] to scale the scalar
    # value of the polynomial. This is closer to an actual polynomial approach, 
    # but comes with a larger number of parameters to estimate
    if(simple) {
        # Just select the data `x` and `y` and take them both to the desired power
        X <- data %>% 
            dplyr::select(x, y)
        X <- lapply(seq_len(n), 
                    \(i) X^i)
        X <- append(list(as.data.frame(rep(1, nrow(data)))), 
                    X)
    } else {
        # Here, things are somewhat more complicated. First, we make the 
        # combinations for taking something to a given power and delete those 
        # combinations that would lead to a degree greater than n. Then, we 
        # loop over these combinations and take each of the independent variables
        # to that power.
        degrees <- cbind(x = rep(0:n, each = n + 1), 
                         y = rep(0:n, times = n + 1)) %>% 
            as.data.frame() %>% 
            dplyr::mutate(degree = x + y) %>% 
            dplyr::filter(degree <= n) %>% 
            dplyr::select(-degree)

        X <- lapply(seq_len(nrow(degrees)), 
                    \(i) data$x^degrees$x[i] * data$y^degrees$y[i])
    }    

    X <- do.call("cbind", X) %>% 
        as.matrix()

    # Solve the equation with the provided parameters
    data[,c("x", "y")] <- X %*% params
    return(data)
}