
test_that("Missing arguments are dealt with", {
  expect_error(ars(), "Missing input arguments")
  expect_error(ars(N=100), "Missing input arguments")
})

test_that("Invalid inputs are dealt with", {
  expect_error(ars(100, 100), "f must be a function")
  expect_error(ars(dnorm, dnorm), "N must be a numeric input") 
  expect_error(ars(dnorm, 100, bounds=2), "length of bounds must be 2")
  expect_error(ars(dnorm, 100, x0 = c(-1, 5), bounds = c(0, 0.1)), "x0 must be inside the bounds")
})

test_that("Inputting log-concave functions should work", {
  expect_error(ars(dgamma, 100, 2, c(0,Inf), shape=2), NA) #gamma(2,1)
  expect_error(ars(dgamma, 100, 2, c(-Inf,Inf), shape=2), "Invalid bounds") 
  expect_error(ars(dgamma, 100, 0.5, c(0,Inf), shape=2), "Invalid x0")
})

test_that("Inputting a log function should work", {
  expect_error(ars(dlnorm, 100, 1, c(0,Inf)), NA)
})

