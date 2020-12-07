###########################################################
## STAT 243 Final Project
## Adaptive Rejection Sampling Algorithm
##
## Date: 12/15/2020
## Authors: Keqin Cao
###########################################################


# Load Required Packages
library("assertthat")
source("helper.R")

ars <- function (f, N, x0 = c(-1.0 , 1.0), bounds = c(-Inf, Inf), ...) {
  ## Check if there is any missing input arguments
  assertthat::assert_that(!missing(f), !missing(N), msg = "Missing input arguments")
  ## Check if the input is a function
  is_function = function(x) {
    assert_that(is.function(x))
  }
  assert_that(is_function(f), msg = "f must be a function")
  ## Check that the bound has the length of 2
  assert_that(length(bounds) == 2)
  bound_length_2 <- function(x) {
    length(x) == 2
  }
  ## Failure if the length of the bound is not equal to 2
  on_failure(bound_length_2) == function(call, env) {
    return("The length of the bounds (lower + upper) must be 2")
  }
  ## Check that the lower bound is not equal to upper bound
  assert_that(bounds[1] != bounds[2])
  ## Raise warning message and replace if they are the same
  if (bounds[1] == bounds[2]) {
    warning("Upper bound and lower bound cannot be the same, switch to default values")
    bounds <- c(-Inf, Inf)
  }
  ## Check if lower bound is smaller than upper bound
  assert_that(bounds[1] < bounds[2], msg = "Lower bound must be smaller than upper bound")
  ## Check if lower bound is numeric input
  assert_that(is.numeric(bounds[1]), msg = "lower bound must be numeric values")
  ## Check if upper bound is numeric input
  assert_that(is.numeric(bounds[2]), msg = "upper bound must be numeric values")
  ## Check if x0 is within bounds
  assert_that( all((x0 > bounds[1]) & (x0 < bounds[2])), msg = "x0 must be inside the bounds")



  ## Define log function
  h <- function(x){
    return(log(f(x, ...elt())))
  }

  ## Define variables
  x0<- sort(x0)
  h0 <- h(x0)
  assert_that(length(h0) != 0, msg = "log function of x0 cannot be NaN/infinite")
  # h'(x) finite difference approx
  x_finite <- 1E-8
  ## Finite difference approximation to h'(x)
  dh0 <- (h(x0 + x_finite) - h0)/x_finite
  # Extract for valid x0
  x0 <- x0[is.finite(h0)&is.finite(dh0)]
  # Extract for valid h0
  h0 <- h0[is.finite(h0)&is.finite(dh0)]
  # Extract for valid dh0
  dh0 <- dh0[is.finite(h0)&is.finite(dh0)]

  ## helper variables
  dx <- 1E-8  # mesh size for finite difference approximation to h'(x)
  eps <- 1E-7  # tolerance for non-log-concavity test
  max_iters <- 10000  # prevent infinite loop
  current_iter <- 0
  bds_warn <- FALSE  # Flag for boundary warning in main while loop

  ## initialize vectors
  x0 <- sort(x0)  # Ensures x0 is in ascending order
  x_j <- c()  # evaluated x values
  h_j <- c()  # evaluated h(x) values
  dh_j <- c()  # evaluated h'(x) values
  x <- c()  # accepted sample vector

  ## Evaluate h(x0)
  h0 <- h(x0)
  ## Finite difference approximation to h'(x)
  dh0 <- (h(x0 + dx) - h0)/dx

  ## Check for NaNs and infinities
  isnum <- is.finite(h0)&is.finite(dh0)
  x0 <- x0[isnum]
  h0 <- h0[isnum]
  dh0 <- dh0[isnum]

}
