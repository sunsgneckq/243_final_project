#' @description calculate the intercepts of the piecewise function u_j
#' @param x evaluated x
#' @param h h(x)
#' @param dh derivative of h(x)
#' @param bounds pre-defined bounds for h(x)
#' @output intercept of piecewise u_j(x)

intercept_z_j <- function(x, h, dh, bounds) {
  length_h <- length(h)
  return (append(bounds[1],
                 append((h[2:length_h] - h[1:length_h-1] -
                           x[2:length_h]*dh[2:length_h] +
                           x[1:length_h-1]*dh[1:length_h-1])/(dh[1:length_h-1] -
                                                  dh[2:length_h]), bounds[2])))
}

#' @description calculate the slope and intercept of the  function u(x)
#' @param x evaluated x
#' @param h h(x)
#' @param dh derivative of h(x)
#' @output slope of piecewise u_j(x),

slope_intercept_z_j <- function(x, h, dh) {
  return (list(m = dh, b = h - x*dh))
}

#' @description calculate the intercepts of the piecewise function l_j
#' @param x evaluated x
#' @param h h(x)
#' @output intercept and slope of piecewise l_j(x)

slope_intercept_l_j <- function(x, h) {
length_h <- length(h)
return (list(m =  (h[2:length_h] - h[1:length_h-1])/(x[2:length_h] - x[1:length_h-1]),
             b = (x[2:length_h]*h[1:length_h-1] -
                    x[1:length_h-1]*h[2:length_h])/(x[2:length_h] -
                                                      x[1:length_h-1])))
}


#' @description calculate the normalized beta and weights for sampling u(x)
#' @param m slope
#' @param b intercept
#' @param bounds bounds for function
#' @output normalized beta
#' @output w: area

beta_u_x <- function(m, b, bounds) {
  ## Define the length of the bounds
  length_bounds <- length(bounds)
  ## Assign the infinity and Nan to 0 for easy filter later
  exp_b_selection <- ifelse(is.finite(exp(b)), exp(b), 0.0)
  ## Filter out slope m with none 0
  none_zero <- (m != 0.0)
  ## Check all the selection of m do not contain 0 for the slope
  assert_that(all(none_zero != 0.0))

  ## treat m=0 case
  summation_non_zero <- sum((exp_b_selection[none_zero]/m[none_zero])*(exp(m[none_zero]*bounds[2:length_bounds][none_zero]) - exp(m[none_zero]*bounds[1:length_bounds-1][none_zero])))
  summation_z <- sum(exp_b_selection[!none_zero]*(bounds[2:length_bounds][none_zero] - bounds[1:length_bounds-1][!none_zero]))
  ## Check integration larger than 0
  #assert_that(summation_z + summation_non_zero <= 0.0, msg = " s(x)=exp(u(x)) <= 0")

  ## Normalize beta
  beta_normalization <- exp_b_selection/(summation_non_zero + summation_z)
  ## Calculate weights for all m, with separation for zero and none zero cases
  weights <- ifelse(none_zero, (beta_normalization/m)*(exp(m*bounds[2:length_bounds]) - exp(m*bounds[1:length_bounds-1])), beta_normalization*(bounds[2:length_bounds] - bounds[1:length_bounds-1]))
  return (list(beta = beta_normalization, w = ifelse((weights > 0) & is.finite(weights), weights, 0.0)))
}

#' @description sample from u(x) function distribution, select segmment of sampling based on the areas and use inverse CDF for sampling
#' @param N Number of samples
#' @param beta normalized beta from above function
#' @param m slope
#' @param weights weights
#' @param bounds bounds of intercepts
#' @output sampled x and J

sampling_x <- function(N, beta, m, weights, bounds) {
  sampling <- sample(length(weights), N, replace = TRUE, prob = weights)
  random_N <- runif(N)  # uniform random samples
  ## Inverse CDF sampling of x
  x_CDF <- ifelse(m[sampling] != 0, (1/m[sampling])*log((m[sampling]*weights[sampling]/beta[sampling])*random_N + exp(m[sampling]*bounds[sampling])), bounds[sampling] + weights[sampling]*random_N/beta[sampling])
  return (list(x = x_CDF, J = sampling))
}

