
test_that("Missing arguments are dealt with", {
  expect_error(ars(), "Missing input arguments")
  expect_error(ars(N = 100), "Missing input arguments")
})


test_that("Invalid inputs are dealt with", {
  expect_error(ars(100, 100), "f must be a function")
  expect_error(ars(dnorm, dnorm), "N must be a numeric input") 
  expect_error(ars(dnorm, 100, bounds = 2), "The length of the bounds (lower + upper) must be 2")
  expect_warning(ars(dnorm, 100, bounds=c(1,1)), 
                 "Upper bound and lower bound cannot be the same, switch to default values")
  expect_error(ars(dnorm, 100, x0 = c(-1, 5), bounds = c(0, 0.1)), "x0 must be inside the bounds")
})


test_that("Inputting f(x) log-concave", {
  expect_error(ars(dexp, 100, c(0.5), c(0, Inf)), NA)
  expect_error(ars(dexp, 100, c(0.5)), "Invalid bounds") # dexp not defined on (-Inf,0)
  expect_error(ars(dnorm, 100, x0 = c(-0.5)), "Invalid x0") ## Invalid bounds and x0 can be 
                                                            # hard to check for
  expect_error(ars(dexp, 100, 0.5, c(0, Inf), rate = 50), "f(x) too small for the given range of x") 
                                                            # X~exp(50) gets really small as x increases
})


test_that("Inputting f(x) not log-concave", {
  expect_error(ars(dlnorm, 100, 1, c(0,Inf)), "Input function f not log-concave.")
  expect_error(ars(dt, 100, df = 2), "Input function f not log-concave.")
})


test_that("Output as expected for an exponential distribution", {
  N <- 100
  ars_output <- ars(dexp, N, x0 = 0.5, bounds = c(0, Inf), rate=1.5)
  true_mean <- 1/1.5
  conf_int_90 <- t.test(ars_output)$conf.int
  
  expect_length(ars_output, N)
  expect_true(conf_int_95[1] <= true_mean & conf_int_95[2] >= true_mean)
})


test_that("Output as expected for a normal distribution", {
  N <- 100
  ars_output <- ars(dnorm, N, x0 = c(10,13), mean = 12, sd = 3)
  true_mean <- 12
  true_var <- 9
  
  sample_var <- var(ars_output)
  var_CI_95 <- (N - 1) * sample_var / qchisq(c(0.975, 0.025), N-1)
  mean_CI_95 <- t.test(ars_output)$conf.int
  
  expect_length(ars_output, N)
  expect_true(mean_CI_95[1] <= true_mean & mean_CI_95[2] >= true_mean)
  expect_true(var_CI_95[1] <= true_var & var_CI_95[2] >= true_var)
})


