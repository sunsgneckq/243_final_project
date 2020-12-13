
#' @description calculate the intercepts of the piecewise function u_j
#' @param x evaluated x
#' @param h h(x)
#' @param dh derivative of h(x)
#' @param bounds pre-defined bounds for h(x)
#' @output intercept of piecewise u_j(x)

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

#' @description calculate the slope and intercept of the function u(x)
#' @param x evaluated x
#' @param h h(x)
#' @param dh derivative of h(x)
#' @output slope of piece-wise u_j(x),

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

#' @description sample from u(x) function distribution, select segmment of sampling based on the areas and use inverse CDF for sampling
#' @param N Number of samples
#' @param beta normalized beta from above function
#' @param m slope
#' @param weights weights
#' @param bounds bounds of intercepts
#' @output sampled x and J

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

boundary_warning <- function(boolean_check) {
  if (!boolean_check) {
    warning('Check bounds! The Sampled x is outside bounds ')
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

duplication_check <- function(x, y, L, eps, dx) {
  dup <-
    append(TRUE, ((abs(x[1:L - 1] - x[2:L]) > eps) &
                    ((y[2:L] - y[1:L - 1]) > dx)))
}
