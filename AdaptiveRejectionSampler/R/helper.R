

#' @description calculate the intercepts of the piecewise function u_j
#' @param x evaluated x
#' @param h h(x)
#' @param dh derivative of h(x)
#' @description Calculates where the tangents to \eqn{u_k(x)} at neighboring x
#' values intersect. \eqn{u_k} is a piecewise linear upper hull formed by the
#' tangents to h(x).
#' @param x numeric vector of x values to evaluate the tangents at
#' @param h concave function h(x) evaluated at x
#' @param dh derivative of h(x) evaluated at x

#' @return Numeric vector of where the tangents to \eqn{u_k} intersect,
#' including bounds values.

intercept_z_j <- function(x, h, dh, bounds) {
  length_h <- length(h)
  return (append(bounds[1],

                 append((h[2:length_h] - h[1:length_h - 1] -
                           x[2:length_h] * dh[2:length_h] +
                           x[1:length_h - 1] * dh[1:length_h - 1]) / (dh[1:length_h -
                                                                           1] -
                                                                        dh[2:length_h]), bounds[2]
                 )))
}

#' @description Calculates the slopes and intercepts of piecewise function
#' \eqn{u_k} at given x values.
#' @param x numeric vector of x values to evaluate at
#' @param h h(x) evaluated at x
#' @param dh derivative of h(x) evaluated at x
#' @return A list of length 2. \code{m} numeric vector of slopes, \code{b}
#' numeric vector of intercepts.

slope_intercept_z_j <- function(x, h, dh) {
  return (list(m = dh, b = h - x * dh))
}


#' @description calculate the intercepts of the piece-wise function l_j
#' @param x evaluated x
#' @param h h(x)
#' @output intercept and slope of piece-wise l_j(x)

slope_intercept_l_j <- function(x, h) {
  length_h <- length(h)
  return (list(
    m =  (h[2:length_h] - h[1:length_h - 1]) /
      (x[2:length_h] - x[1:length_h -1]),
    b = (x[2:length_h] * h[1:length_h - 1] -
           x[1:length_h - 1] * h[2:length_h]) / (x[2:length_h] -
                                                   x[1:length_h - 1])
  ))
}


integration_check <- function(x, y) {
  if (x + y <= 0.0) {
    stop('Area of s(x)=exp(u(x)) <= 0')
  }
}

#' @description calculate the normalized beta and weights for sampling u(x)
#' @param m slope
#' @param b intercept
#' @param bounds bounds for function
#' @output normalized beta
#' @output w: area weights

beta_u_x <- function (m, b, bounds) {
  length_bounds <- length(bounds)
  b_selection_finite <- ifelse(is.finite(exp(b)), exp(b), 0.0)

  ## treat m=0 case
  non_zero_boolean <- (m != 0.0)
  non_zero_boolean_summation <- sum((b_selection_finite[non_zero_boolean]
                                     / m[non_zero_boolean]) *
                                      (exp(m[non_zero_boolean] *
                                             bounds[2:length_bounds][non_zero_boolean]) -
                                         exp(m[non_zero_boolean] *
                                               bounds[1:length_bounds -
                                                        1][non_zero_boolean])))
  summation_z <- sum(b_selection_finite[!non_zero_boolean] *
                       (bounds[2:length_bounds][!non_zero_boolean] -
                          bounds[1:length_bounds - 1][!non_zero_boolean]))

  ## Enforce that the integral of s(x) > 0
  integration_check(non_zero_boolean_summation, summation_z)
  ## Normalize beta
  beta <- b_selection_finite / (non_zero_boolean_summation + summation_z)
  ## Calculate weights for both m!=0 and m=0 cases
  weight <-
    ifelse(non_zero_boolean, (beta / m) * (exp(m * bounds[2:length_bounds]) -
                                             exp(m * bounds[1:length_bounds - 1])),
           beta * (bounds[2:length_bounds] - z[1:length_bounds -1]))
  return (list(beta = beta, w = ifelse((weight > 0) & is.finite(weight), weight, 0.0)))
}

  #' @description Calculates the slopes and intercepts of piecewise function
  #' \eqn{l_j} at given x values. \eqn{l_j} is a piecewise linear lower hull
  #' formed by adjacent chords.
  #' @param x numeric vector of x values to evaluate at
  #' @param h h(x) evaluated at x
  #' @return A list of length 2. \code{m} numeric vector of slopes, \code{b}
  #' numeric vector of intercepts.

  slope_intercept_l_j <- function(x, h) {
    length_h <- length(h)
    return (list(m = (h[2:length_h] - h[1:length_h-1]) /
                   (x[2:length_h] - x[1:length_h-1]),
                 b = (x[2:length_h] * h[1:length_h-1] - x[1:length_h-1] *
                        h[2:length_h]) / (x[2:length_h] - x[1:length_h-1])))
  }


  #' @description Calculates the normalized beta and weights for sampling \eqn{u_k}.
  #' @param m numeric vector of slopes
  #' @param b numeric vector of intercepts
  #' @param bounds numeric vector of abscissae bounds for function
  #' @return A list of length 2. \code{beta} numeric vector of normalized beta,
  #' \code{w} numeric vector of weights.

  beta_u_x <- function(m, b, bounds) {
    # Number of abscissae bounds on the domain
    length_bounds <- length(bounds)
    # Assign infinity and NaN to 0 for easy filter later
    exp_b_selection <- ifelse(is.finite(exp(b)), exp(b), 0.0)
    # Filter out non-zero slope
    none_zero <- (m != 0.0)
    # Check all the selection of m do not contain 0 for the slope
    assert_that(all(none_zero != 0.0))
    # treat m = 0 case
    sum_none_zero <- sum((exp_b_selection[none_zero] / m[none_zero]) *
                           (exp(m[none_zero] * bounds[2:length_bounds][none_zero]) -
                              exp(m[none_zero] * bounds[1:length_bounds-1][none_zero])))
    sum_z <- sum(exp_b_selection[!none_zero] *
                   (bounds[2:length_bounds][none_zero] -
                      bounds[1:length_bounds-1][!none_zero]))

    ## Check integration larger than 0
    #assert_that(sum_z + sum_none_zero <= 0.0, msg = " s(x)=exp(u(x)) <= 0")

    # Normalize beta
    beta_normalization <- exp_b_selection/(sum_none_zero + sum_z)
    # Calculate weights for all m, with separation for zero and none zero cases
    weights <- ifelse(none_zero,
                      (beta_normalization/m)*(exp(m*bounds[2:length_bounds]) -
                                                exp(m*bounds[1:length_bounds-1])),
                      beta_normalization*(bounds[2:length_bounds] -
                                            bounds[1:length_bounds-1]))
    return(list(beta = beta_normalization,
                w = ifelse((weights > 0) & is.finite(weights), weights, 0.0)))
  }


  #' @description Samples from \eqn{s_k(x)} function distribution, select segment
  #' of sampling based on the areas and use inverse CDF for sampling
  #' @param N Number of samples
  #' @param beta normalized beta given by \code{beta_u_x} function
  #' @param m slope given by \code{slope_intercept_z_j} function
  #' @param weights weights given by \code{beta_u_x} function
  #' @param bounds intercepts given by \code{intercept_z_j} function
  #' @return A list of two vectors of length N of sampled x and J.

  sampling_x <- function(N, beta, m, weights, bounds) {
    sampling <-
      sample(length(weights), N, replace = TRUE, prob = weights)
    random_N <- runif(N)  # uniform random samples

    ## Inverse CDF sampling of x
    x_CDF <-
      ifelse(m[sampling] != 0,
             (1 / m[sampling]) * log((m[sampling] * weights[sampling] / beta[sampling]) *
                                       random_N + exp(m[sampling] * bounds[sampling])),
             bounds[sampling] + weights[sampling] * random_N / beta[sampling])

    return (list(x = x_CDF, J = sampling))
  }




  max_iter_check <- function(i, max_iter) {
    if (i > max_iter) {
      stop('You have reached the maximum number of iteration to prevent overflow')
    }
  }
  #' @description Checks whether the function h is log-concave.
  #' @param bound numeric
  #' @param L numeric
  #' @param x numeric vector
  #' @return Gives an error message if the function h is not log-concave.

      log_concavity_check <- function(bound, L, x) {
        if (L > 1) {
          ## Check log-concavity of function
          if(!all((x[1:L-1] - x[2:L]) >= bound)) {
            stop('Input function f not log-concave.')
          }
        }
      }


      #' @description Checks whether values are within given bounds.
      #' @param x numeric vector
      #' @param bounds numeric vector of length 2, lower bound and upper bound
      #' @return \code{TRUE} if \code{x} is within the given bounds, \code{FALSE} otherwise.
      boundary_check<- function(x, bound) {
        ((x > bound[1]) & (x < bound[2]))
      }


      #' @description Whether the boundary condition is being checked.
      #' @param boolean_check boolean
      #' @return Gives a warning message if the \code{boolean_check} is \code{FALSE}.
      boundary_warning <- function(boolean_check) {
        if (!boolean_check) {
          warning('Check bounds! The sampled x is outside bounds')
          boolean_check <- TRUE
        }
      }


      log_concavity_check <- function(bound, L, x) {
        if (L > 1) {
          ## Check log-concavity of function
          if (!all((x[1:L - 1] - x[2:L]) >= bound)) {
            stop('Input function f not log-concave.')
          }
        }
      }


      finite_check <- function(x, y) {
        is.finite(x) & is.finite(y)
      }

      boundary_check <- function(x, bound) {
        ((x > bound[1]) & (x < bound[2]))
      }


        #' @description Checks whether the algorithm has reached a maximum number of
        #' iterations.
        #' @param i current iteration
        #' @param max_iter maximum number of iterations
        #' @return Gives an error message if the algorithm has reached the maximum number
        #' of iterations.
        max_iter_check <- function(i,max_iter){
          if (i > max_iter) {
            stop('You have reached the maximum number of iteration to prevent overflow')
          }
        }


        #' @description Checks for NaNs and infinities.
        #' @param x numeric vector
        #' @param y numeric vector
        #' @return \code{TRUE} if \code{x} and \code{y} are finite, \code{FALSE} otherwise.
        finite_check<- function(x,y) {
          is.finite(x) && is.finite(y)
        }


        #' @description Breaks out of the current loop if the length of \code{x} is 1.
        #' @param x numeric vector
        length_stopper<- function(x){
          length_dh_k <- length(x)
          if (length_dh_k == 1) {
            break
          }
        }


        duplication_check <- function(x, y, L, eps, dx) {
          dup <-
            append(TRUE, ((abs(x[1:L - 1] - x[2:L]) > eps) &
                            ((y[2:L] - y[1:L - 1]) > dx)))
        }


