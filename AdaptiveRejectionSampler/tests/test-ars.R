require(EnvStats) # for Pareto distribution
require(extraDistr) # for Laplace distribution

test_that("Missing arguments are dealt with", {
  expect_error(ars(), "Missing input arguments")
  expect_error(ars(N = 100), "Missing input arguments")
})


test_that("Input checks work properly", {
  expect_error(ars(100, 100), "f must be a function")
  expect_error(ars(dnorm, "10"), "N must be a numeric input")
  expect_error(ars(dnorm, -10), "N must be larger than 0")
  expect_error(ars(dnorm, 100, x0 = "x0"), "x0 must be numeric input")
  
  expect_error(ars(dnorm, 100, bounds = 2), "Length of bounds must be 2")
  expect_error(ars(dnorm, 100, bounds = c("2", "-1.5")), "bounds must be numeric values")
  expect_error(ars(dnorm, 100, bounds = c(3.5, 1)), "Lower bound must be smaller than upper bound")
  expect_warning(ars(dnorm, 100, bounds=c(1,1)), "Upper bound and lower bound cannot be the same, switch to default values")
  expect_error(ars(dnorm, 100, x0 = c(-1, 5), bounds = c(0, 0.1)), "x0 must be inside the bounds")
})


test_that("Inputting f(x) non-log-concave throws error", {
  expect_error(ars(dt, 100, df = 2), "Input function f not log-concave.")
  expect_error(ars(dcauchy, 100), "Input function f not log-concave.")
  expect_error(ars(dlnorm, 100, 1, c(0, Inf)), "Input function f not log-concave.")
  expect_error(ars(df, 100, 0.5, c(0, Inf), df1 = 5, df2 = 2), "Input function f not log-concave.")
  expect_error(ars(EnvStats::dpareto, 100, 1.5, c(1, Inf), location = 1), "Input function f not log-concave.")
})


test_that("Possible errors for inputting f(x) log-concave", {
  expect_error(ars(dexp, 100, c(0.5))) # dexp not defined on (-Inf,0)
  expect_error(ars(dnorm, 100, x0 = c(-0.5))) # invalid x0
  expect_error(ars(dexp, 100, 0.5, c(0, Inf), rate = 50)) # f(x) too small for the given range of x"
})


test_that("Result as expected for normal distribution", {
  N <- 500
  set.seed(123)
  norm_output <- ars(dnorm, N, x0 = c(10,13), mean = 12, sd = 3)
  true_mean <- 12
  true_var <- 9

  mean_CI_95 <- t.test(norm_output)$conf.int
  var_CI_95 <- (N - 1) * var(norm_output) / qchisq(c(0.975, 0.025), N-1)
  
  expect_length(norm_output, N)
  expect_true(mean_CI_95[1] <= true_mean & mean_CI_95[2] >= true_mean)
  expect_true(var_CI_95[1] <= true_var & var_CI_95[2] >= true_var)
})


test_that("Result as expected for exponential distribution", {
  N <- 500
  set.seed(123)
  exp_output <- ars(dexp, N, x0 = 0.5, bounds = c(0, Inf), rate = 0.5)
  true_mean = 1/0.5
  mean_CI_95 <- t.test(exp_output)$conf.int
  
  expect_length(exp_output, N)
  expect_true(mean_CI_95[1] <= true_mean & mean_CI_95[2] >= true_mean)
})


test_that("Result as expected for uniform distribution", {
  N <- 500
  set.seed(123)
  unif_output <- ars(dunif, N, x0 = 0.5, bounds = c(0,1))
  true_mean <- 1/2
  mean_CI_95 <- t.test(unif_output)$conf.int
  
  expect_length(unif_output, N)
  expect_true(mean_CI_95[1] <= true_mean & mean_CI_95[2] >= true_mean)
})


test_that("Result as expected for gamma distribution", {
  N <- 500
  set.seed(123)
  gamma_output <- ars(dgamma, N, x0 = 1, bounds = c(0, Inf), shape = 2, rate = 2)
  true_mean <- 1
  mean_CI_95 <- t.test(gamma_output)$conf.int
  
  expect_length(gamma_output, N)
  expect_true(mean_CI_95[1] <= true_mean & mean_CI_95[2] >= true_mean)
})


test_that("Result as expected for beta distribution", {
  N <- 500
  set.seed(123)
  beta_output <- ars(dbeta, N, c(0.2, 0.6), c(0,1), shape1 = 2, shape2 = 5)
  true_mean <- 2/(2+5)
  mean_CI_95 <- t.test(beta_output)$conf.int
  
  expect_length(beta_output, N)
  expect_true(mean_CI_95[1] <= true_mean & mean_CI_95[2] >= true_mean)
})


test_that("Result as expected for laplace distribution", {
  N <- 500
  set.seed(123)
  laplace_output <- ars(extraDistr::dlaplace, N, c(-10, 5, 0, 5, 10))
  true_mean <- 0
  mean_CI_95 <- t.test(laplace_output)$conf.int
  
  expect_length(laplace_output, N)
  expect_true(mean_CI_95[1] <= true_mean & mean_CI_95[2] >= true_mean)
})

