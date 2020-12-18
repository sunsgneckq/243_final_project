test_that("intercept_z_j works as expected", {
  dx <- 1E-8
  x <- runif(10,0,2)
  l <- length(x)
  hx <- dnorm(x, log=T)
  dhx <- (dnorm(x+dx, log=T) - dnorm(x, log=T))/dx
  bounds <- c(-Inf, Inf)
  
  res <- c(-Inf, (hx[2:l] - hx[1:l-1] - x[2:l] * dhx[2:l] + x[1:l-1] * dhx[1:l-1]) / 
             (dhx[1:l-1] - dhx[2:l]), Inf)
  
  expect_length(intercept_z_j(x,hx,dhx,bounds), length(x) + 1)
  expect_equal(intercept_z_j(x, hx, dhx, bounds), res)
  expect_true(all(intercept_z_j(x, hx, dhx, bounds) >= 
                    intercept_z_j(x, hx, dhx, bounds)[1]) &
                all(intercept_z_j(x, hx, dhx, bounds) <= 
                      intercept_z_j(x, hx, dhx, bounds)[11])) #intercepts are between bounds
})


test_that("slope_intercept_l_j works as expected", {
  x <- runif(10,0,2)
  l <- length(x)
  hx <- dnorm(x, log=T)
  m <- (h[2:l] - h[1:l-1]) / (x[2:l] - x[1:l-1])
  b <- (x[2:l] * h[1:l-1] - x[1:l-1] * h[2:l]) / (x[2:l] - x[1:l-1])
  expect_equal(slope_intercept_l_j(x,h)[['m']], m)
  expect_equal(slope_intercept_l_j(x,h)[['b']], b)
})


test_that("beta_u_x works as expected", {
  bounds <- c(-4:1)
  l <- length(bounds)
  b <- runif(5, -1, 1)
  m <- runif(5, -2, 0)
  b_finite <- ifelse(is.finite(exp(b)), exp(b), 0.0)
  non_zero <- (m != 0.0)
  non_zero_sum <- sum((b_finite[non_zero] / m[non_zero]) * 
                        (exp(m[non_zero] *bounds[2:l][non_zero]) - 
                           exp(m[non_zero] * bounds[1:l-1][non_zero])))
  
  z_sum <- sum(b_finite[!non_zero]*(bounds[2:l][!non_zero]-bounds[1:l-1][!non_zero]))
  beta <- b_finite / (non_zero_sum + z_sum)
  weight <- ifelse(non_zero, 
                   (beta / m) * (exp(m * bounds[2:l]) - exp(m * bounds[1:l-1])),
                   beta * (bounds[2:l] - bounds[1:l-1]))
  w <- ifelse((weight > 0) & is.finite(weight), weight, 0.0)

  expect_equal(beta_u_x(m, b, bounds)[['beta']], beta)
  expect_equal(beta_u_x(m, b, bounds)[['w']], w)
  expect_length(beta_u_x(m, b, bounds)[['beta']], length(b))
  expect_length(beta_u_x(m, b, bounds)[['w']], length(b))
})


test_that("log_concavity_check works as expected", {
  x <- c(0:10)
  expect_error(log_concavity_check(bound = 2, L = 2, x), "Input function f not log-concave.")
  expect_error(log_concavity_check(bound = 1.1, L = 2, x = runif(10)), "Input function f not log-concave.")
  expect_error(log_concavity_check(bound = 2, L = 1, x), NA)
})


test_that("duplication_check works as expected", {
  eps <- 1E-7
  dx <- 1E-8
  expect_length(duplication_check(1:5, 1:5, 5, eps, dx), 5)
  expect_equal(duplication_check(1:5, 2:6, 5, eps, dx), c(T,T,T,T,T))
  expect_equal(duplication_check(1:5, c(1,2,7,2,2), 5, eps, dx), c(T,T,T,F,F))
})



