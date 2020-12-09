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
  on_failure(bound_length_2) = function(call, env) {
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
  h <- function(x) {
    return (log(f(x, ...)))
  }

  ## Define variables
  dx <- 1E-8
  eps <- 1E-7
  x0<- sort(x0)
  x_k <- c()  # evaluated x values
  h_k <- c()  # evaluated h(x) values
  dh_k <- c()  # evaluated h'(x) values
  x <- c()  # accepted sample vector

  h0 <- h(x0)
  assert_that(length(h0) != 0, msg = "log function of x0 cannot be NaN/infinite")


  ## initialize vectors

  ## Finite difference approximation to h'(x)
  dh0 <- (h(x0 + dx) - h0)/dx
  # Extract for valid x0
  x0 <- x0[is.finite(h0)&is.finite(dh0)]
  # Extract for valid h0
  h0 <- h0[is.finite(h0)&is.finite(dh0)]
  # Extract for valid dh0
  dh0 <- dh0[is.finite(h0)&is.finite(dh0)]

  print(dh0)


  i <- 0
  eps<-1E-7
  bound_warning<- FALSE
  while(length(x) < N){
    i <- i + 1
    if(i > 10000) {
      stop("To prevent infinite loop, max number of iteration is set")
    }
      chunk_size <- min(c(N, i**2))

      ### Initizlization and update
      if (0 <length(x0)) {

        x_k <- append(x_k, x0)
        h_k <- append(h_k, h0)
        dh_k <- append(dh_k, dh0)
        x0 <- c()
        h0 <- c()
        dh0 <- c()
        sorted_string <- sort(x_k, index.return = TRUE)
        x_k <- sorted_string$x
        h_k <- h_k[sorted_string$ix]
        dh_k <- dh_k[sorted_string$ix]

        length_of_dh_k <- length(dh_k)

        if (length_of_dh_k > 1) {
          while (!all((abs(dh_k[1:length_of_dh_k-1] - dh_k[2:length_of_dh_k]) > eps) & ((x_k[2:length_of_dh_k] - x_k[1:length_of_dh_k-1]) > dx))) {

            ## Only keep values with dissimilar neighbors
            ## Always keep first index (one is always unique)
            duplication_check <- append(TRUE, ((abs(dh_k[1:length_of_dh_k-1] - dh_k[2:length_of_dh_k]) > eps) & ((x_k[2:length_of_dh_j] - x_k[1:length_of_dh_j-1]) > dx)))
            x_k <- x_k[duplication_check]
            h_k <- h_k[duplication_check]
            dh_k <- dh_k[duplication_check]

            length_of_dh_k <- length(dh_k)
            if (length_of_dh_k == 1) {
              break
            }
          }

          if (length_of_dh_k > 1) {
            ## Ensure log-concavity of function
            if(!all((dh_k[1:length_of_dh_k-1] - dh_k[2:length_of_dh_k]) >= eps)) {
              stop('Input function f not log-concave.')
            }
          }
        }
        z_k <- intercept_z_j(x_k, h_k, dh_k, bounds)
        u_k <- slope_intercept_z_j(x_k, h_k, dh_k)
        m_u <- u_k$m
        b_u <- u_k$b
        l_k <- slope_intercept_l_j(x_k, h_k)


        s_k <- beta_u_x(u_k$m, u_k$b, z_k)
        beta_k <- s_k$beta
        weight_k <- s_k$w
      }

     ### Sampling step
      sampling <- sampling_x(chunk_size, beta_k, u_k$m,weight_k, z_k)
      x_s <- sampling$x
      J<- sampling$J
      i <- runif(chunk_size)
      print(x_s)
     # assert_that(all(x_s > bounds[1]) , msg = "x* not within bounds")

      if (!all((x_s > bounds[1]) & (x_s < bounds[2]))) {
        if (!bound_warning<- FALSE) {
          warning('x* not inside the bounds')
          bound_warning <- TRUE
        }
        check_bound <- ((x_s > bounds[1]) & (x_s < bounds[2]))
        print(check_bound)
        J<- J[check_bound]
        i<- i[check_bound]
        x_s <-x_s[check_bound]
      }

      ## only use x-values where l_j(x) > -Inf
      boundary <- (J - ifelse(x_s < x_k[J], 1, 0) >= 1) &
        (J - ifelse(x_s < x_k[J], 1, 0) < length(x_k))
      # rejection on u(x)/l(x)
      J_l <- J - ifelse(x_s < x_k[J], 1, 0)
      y <- exp(x_s[boundary]*(l_k$m[J_l[boundary]] -u_k$m[J[boundary]]) +  l_k$b[J_l[boundary]] - u_k$b[J[boundary]])
      counts <- (i[boundary] <= y)
      ## Append accepted values to x
      x <- append(x, x_s[boundary][i])
      ## Append rejected values to x0
      x0 <- append(x_s[!boundary], x_s[boundary][!counts])


      if (length(x0) > 0) {

        ## Evaluate h(x0)
        h0 <- h(x0)
        ## Finite difference approximation to h'(x)
        dh0 <- (h(x0 + dx) - h0)/dx

        ## Check for NaNs and infinities
        check_NA <- is.finite(h0)&is.finite(dh0)
        x0 <- x0[check_NA]
        h0 <- h0[check_NA]
        dh0 <- dh0[check_NA]

        w <- append(i[boundary], i[boundary][!counts])[check_NA]
        J <- append(J[!boundary], J[boundary][counts])[check_NA]

        ## Perform second rejection test on u(x)/h(x)
        counts <- (w <= exp(h0 - x0*u_k$k[J] - u_k$b[J]))
        ## Append accepted values to x
        x <- append(x, x0[counts])
      }
    }
    ## Only return N samples (vectorized operations makes x sometimes larger)
    return (x[1:N])
  }

