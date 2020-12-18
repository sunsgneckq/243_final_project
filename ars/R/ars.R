###########################################################
## STAT 243 Final Project
## Adaptive Rejection Sampling Algorithm
##
## Date: 12/15/2020
## Authors: Keqin Cao, Mark Campmier, Colleen Sun
###########################################################


# Load Required Packages
library("assertthat")
# Source the helper functions
#source("R/helper.R")
#load_all("/Users/sunsgne/Desktop/243_final_project/ars/R")


#' @title Adaptive Rejection Sampling
#' @description Generates samples from a log-concave distribution via adaptive
#' rejection sampling.
#' @param f vectorized function that's log-concave
#' @param N positive numeric; number of samples to generate
#' @param x0 numeric vector; bounds of the sampling domain
#' @param bounds numeric vector of length 2; boundaries of the underlying sampling
#' distribution f(x)
#' @param ... additional arguments to be passed to \code{f}
#' @return A vector of N samples generated from the f(x) distribution.
#' @details 
#' If \code{x0} is not specified, it assumes the default value of \code{c(-1.0, 1.0)}. 
#' If \code{bounds} is not specified, it assumes the default value of \code{c(-Inf, Inf)}. 
#' 
#' The argument \code{f} must be a vectorized, log-concave function.
#' 
#' \code{x0} can be a numeric vector of any length, as long as the values are 
#' valid for the given \code{f}.
#' The \code{bounds} must be defined for \code{f}. For example, setting 
#' \code{bounds = c(-Inf, Inf)} will cause an error for exponential distributions 
#' because they are not defined on (-Inf, 0).
#' 
#' If the value of \code{f(x)} is too small for the range of x given by \code{bounds}, 
#' the program will throw an error.
#' 
#' @references 
#' Gilks and Wild.Adaptive Rejection Sampling for Gibbs Sampling(1992)
#' 
#' Paciorek.STAT243 Unit 9 - Simulation(2020)
#' 
#' @examples 
#' ars(dnorm, 100)
#' ## generating 100 samples from Exp(0.5) distribution
#' ars(dexp, 100, x0 = 0.5, bounds = c(0, Inf), rate = 0.5)
#' @export ars
ars <-
  function (f, N, x0 = c(-1.0, 1.0), bounds = c(-Inf, Inf), ...) {
    ## Check if there is any missing input arguments
    assert_that(!missing(f), !missing(N), msg = "Missing input arguments")
    ## Check if the input is a function
    is_function = function(x) {
      assert_that(is.function(x))
    }
    ## Check if the input f is a function
    assert_that(is_function(f), msg = "f must be a function")
    ## Check if the N sample size is integer
    assert_that(is.numeric(N), msg = "N must be a numeric input")
    ## Check if the N sample size is larger than 0
    assert_that(N > 0, msg = "N must be larger than 0")
    ## Check if the x0 has all numeric input
    assert_that(is.numeric(x0), msg = "x0 must be numeric input")
    ## Check that the bound has the length of 2
    assert_that(length(bounds) == 2, msg = "Length of bounds must be 2")
    ## Check if bounds is numeric input
    assert_that(is.numeric(bounds), msg = "bounds must be numeric values")
    ## Raise warning message and replace if they are the same
    if (bounds[1] == bounds[2]) {
      warning("Upper bound and lower bound cannot be the same, switch to default values")
      bounds <- c(-Inf, Inf)
    }
    ## Check if lower bound is smaller than upper bound
    assert_that(bounds[1] < bounds[2], msg = "Lower bound must be smaller than upper bound")
    ## Check if x0 is within bounds
    assert_that(all((x0 > bounds[1]) &
                      (x0 < bounds[2])), msg = "x0 must be inside the bounds")

    ## Round N to the next smallest integer
    N <- ceiling(N)

    ## Define log function
    h <- function(x) {
      return (log(f(x, ...)))
    }


    ## Sort x0 in order
    x0 <- sort(x0)
    ## Define and init some variables
    x_k <- c()
    h_k <- c()
    dh_k <- c()
    x <- c()
    ## Take the log function of x0 as defined
    h0 <- h(x0)
    ## Check if length of h0 is 0
    assert_that(length(h0) != 0, msg = "log function of x0 cannot be NaN/infinite")
    ## Define a finite dx value
    dx <- 1E-8
    ## Cauculate dh0 finite difference
    dh0 <- (h(x0 + dx) - h0) / dx

    ## Exclude invalid x0
    x0 <- x0[finite_check(h0, dh0)]
    ## Exclude invalid h0
    h0 <- h0[finite_check(h0, dh0)]
    ## Exclude invalid dh0
    dh0 <- dh0[finite_check(h0, dh0)]


    #### Main Function ####
    ## Set the starter i to be 0
    i <- 0
    ## Set the maximum iteration before stops to be 15000
    max_iter <- 15000
    ## Define the boolean
    check_boundary <- FALSE
    ## Loop until we have N samples
    while (length(x) < N) {
      i <- i + 1
      max_iter_check(i, max_iter)
      chunk_size_vectorized <- min(c(N, i ** 2))

    #### Init ####

    ## Length of x0 cannot be 0
      if (0 < length(x0)) {
        ## Update
        ## Append x_k to x0 vector
        x_k <- append(x_k, x0)
        ## Append h_k to h0 vector
        h_k <- append(h_k, h0)
        ## Append dh_k to dh0 vector
        dh_k <- append(dh_k, dh0)
        x0 <- c()
        h0 <- c()
        dh0 <- c()
        ## Sort x_k's in ascending order
        sorted_string_index <- sort(x_k, index.return = TRUE)
        ## Filter sorted x_k value
        x_k <- sorted_string_index$x
        ## Filter sorted h_k value
        h_k <- h_k[sorted_string_index$ix]
        ## Filter sorted dh_k value
        dh_k <- dh_k[sorted_string_index$ix]
        ## Compute the length of dh_k
        length_dh_k <- length(dh_k)
        ## Define the boundary
        eps <- 1E-7
        if (length_dh_k > 1) {
          while (!all((abs(dh_k[1:length_dh_k - 1] - dh_k[2:length_dh_k]) > eps) &
                      ((x_k[2:length_dh_k] - x_k[1:length_dh_k - 1]) > dx))) {
            ## Check for duplication and remove based on condition
            x_k <- x_k[duplication_check(dh_k, x_k, length_dh_k, eps, dx)]
            h_k <- h_k[duplication_check(dh_k, x_k, length_dh_k, eps, dx)]
            dh_k <- dh_k[duplication_check(dh_k, x_k, length_dh_k, eps, dx)]
            length_dh_k <- length(dh_k)
            if (length_dh_k == 1) {
              break
            }
          }
          ## Log concavity check
          log_concavity_check(eps, length_dh_k, dh_k)
        }
        ## Use helper function to compute z_k, u_k(x), l_k(x), s_k(x)
        z_k <- intercept_z_j(x_k, h_k, dh_k, bounds)
        u_k <- slope_intercept_z_j(x_k, h_k, dh_k)
        l_k <- slope_intercept_l_j(x_k, h_k)
        s_k <- beta_u_x(u_k$m, u_k$b, z_k)
      }

      ##### SAMPLING STEP

      ## draw x from exp(u(x))
      sampling_x <-
      sampling_x(chunk_size_vectorized, s_k$beta, u_k$m, s_k$w, z_k)
      x_s <- sampling_x$x
      J <- sampling_x$J

      ## random uniform for rejection sampling
      w <- runif(chunk_size_vectorized)

      ## Check how close it is to the distribution
      if (!all(boundary_check(x_s, bounds))) {
        ## Check if the samples are inside bounds
        boundary_warning(check_boundary)
        ## Filter out values that are outside bounds
        J <- J[boundary_check(x_s, bounds)]
        x_s <- x_s[boundary_check(x_s, bounds)]
        w <- w[boundary_check(x_s, bounds)]
      }

      ## only use x-values where l_k(x) > -Inf
      boundary <- (J - ifelse(x_s < x_k[J], 1, 0) >= 1) &
        (J - ifelse(x_s < x_k[J], 1, 0) < length(x_k))

      # New Index for J
      new_index <- J - ifelse(x_s < x_k[J], 1, 0)
      ## Perform first rejection test on u(x)/length_dh_k(x)
      y <- exp(x_s[boundary] * (l_k$m[new_index[boundary]] -
                                  u_k$m[J[boundary]]) +
                 l_k$b[new_index[boundary]] -
                 u_k$b[J[boundary]])
      first_rejection_accept <- (w[boundary] <= y)

      ## Append accepted values to x
      x <- append(x, x_s[boundary][first_rejection_accept])
      ## Append rejected values to x0
      x0 <- append(x_s[!boundary], x_s[boundary][!first_rejection_accept])

      if (length(x0) > 0) {
        ## Evaluate h(x0)
        h0 <- h(x0)
        ## Finite difference approximation to h'(x)
        dh0 <- (h(x0 + dx) - h0) / dx

        ## Check for NaNs and infinities
        x0 <- x0[finite_check(h0, dh0)]
        h0 <- h0[finite_check(h0, dh0)]
        dh0 <- dh0[finite_check(h0, dh0)]
        ## Append first rejection accept value with finite check ready
        w <- append(w[!boundary], w[boundary]
                    [!first_rejection_accept])[finite_check(h0, dh0)]
        J <- append(J[!boundary], J[boundary]
                    [!first_rejection_accept])[finite_check(h0, dh0)]
        ## Append accepted values of the second rejection test on the u(x)/h(x)
        x <- append(x, x0[(w <= exp(h0 - x0 * u_k$m[J] - u_k$b[J]))])
      }
    }
    # Return all N samples in x
    res <- x[1:N]
    return(res)
  }
