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

ars <-
  function (f,N,x0 = c(-1.0 , 1.0),bounds = c(-Inf, Inf),...) {
    ## Check if there is any missing input arguments
    assert_that(!missing(f),!missing(N), msg = "Missing input arguments")
    ## Check if the input is a function
    is_function = function(x) {
      assert_that(is.function(x))
    }
    assert_that(is_function(f), msg = "f must be a function")
    assert_that(is.numeric(N), msg = "N must be a numeric input")
    ## Check that the bound has the length of 2
    assert_that(length(bounds) == 2, msg = "length of bounds must be 2")
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
    assert_that(all((x0 > bounds[1]) &
                      (x0 < bounds[2])), msg = "x0 must be inside the bounds")



    ## Define log function
    h <- function(x) {
      return (log(f(x, ...)))
    }


    ## Define variables
    x0 <- sort(x0)
    x_k <- c()
    h_k <- c()
    dh_k <- c()
    x <- c()

    h0 <- h(x0)
    assert_that(length(h0) != 0, msg = "log function of x0 cannot be NaN/infinite")
    ## Finite difference approximation to h'(x)
    dx <- 1E-8
    dh0 <- (h(x0 + dx) - h0) / dx

    ## Exclude invalid x0
    x0 <- x0[finite_check(h0, dh0)]
    ## Exclude invalid h0
    h0 <- h0[finite_check(h0, dh0)]
    ## Exclude invalid dh0
    dh0 <- dh0[finite_check(h0, dh0)]

    #### Main Function ####
    i <- 0
    max_iter <- 15000
    check_boundary <- FALSE
    ## Loop until we have N samples
    while (length(x) < N) {
      i <- i + 1
      max_iter_check(i, max_iter)
      chunk_size_vectorized <- min(c(N, i ** 2))
    #### INITIALIZATION ####
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
        x_k <- sorted_string_index$x
        h_k <- h_k[sorted_string_index$ix]
        dh_k <- dh_k[sorted_string_index$ix]

        length_dh_k <- length(dh_k)
        eps <- 1E-7
        if (length_dh_k > 1) {
          while (!all((abs(dh_k[1:length_dh_k - 1] - dh_k[2:length_dh_k]) > eps) &
                      ((x_k[2:length_dh_k] - x_k[1:length_dh_k - 1]) > dx))) {
            ## Check for duplication and remove
            x_k <- x_k[duplication_check(dh_k, x_k, length_dh_k, eps, dx)]
            h_k <- h_k[duplication_check(dh_k, x_k, length_dh_k, eps, dx)]
            dh_k <- dh_k[duplication_check(dh_k, x_k, length_dh_k, eps, dx)]

            length_dh_k <- length(dh_k)
            if (length_dh_k == 1) {
              break
            }
          }
          #log concavity check
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

      ## Warn if samples were outside of bounds
      ## This happens if bounds given don't reflect
      ## the actual bounds of the distribution
      if (!all(boundary_check(x_s, bounds))) {
        ## Flag so warning only happens once
        boundary_warning(check_boundary)
        ## Filter out values that are outside bounds
        J <- J[boundary_check(x_s, bounds)]
        x_s <- x_s[boundary_check(x_s, bounds)]
        w <- w[boundary_check(x_s, bounds)]
      }

      ## only use x-values where l_k(x) > -Inf
      boundary <- (J - ifelse(x_s < x_k[J], 1, 0) >= 1) &
        (J - ifelse(x_s < x_k[J], 1, 0) < length(x_k))

      ## Perform first rejection test on u(x)/length_dh_k(x)
      # New Index for J
      new_index <- J - ifelse(x_s < x_k[J], 1, 0)
      y <- exp(x_s[boundary] * (l_k$m[new_index[boundary]] -
                                  u_k$m[J[boundary]]) +
                 l_k$b[new_index[boundary]] -
                 u_k$b[J[boundary]])
      first_rejection_accept <- (w[boundary] <= y)

      ## Append accepted values to x
      x <- append(x, x_s[boundary][first_rejection_accept])
      ## Append rejected values to x0
      x0 <- append(x_s[!boundary],
                   x_s[boundary][!first_rejection_accept])

      if (length(x0) > 0) {
        ## Evaluate h(x0)
        h0 <- h(x0)
        ## Finite difference approximation to h'(x)
        dh0 <- (h(x0 + dx) - h0) / dx

        ## Check for NaNs and infinities
        x0 <- x0[finite_check(h0, dh0)]
        h0 <- h0[finite_check(h0, dh0)]
        dh0 <- dh0[finite_check(h0, dh0)]


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
