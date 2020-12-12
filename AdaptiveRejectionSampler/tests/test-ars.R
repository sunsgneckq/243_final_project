
test_that("Missing arguments are dealt with", {
  expect_error(ars(), "Missing input arguments")
  expect_error(ars(N=100), "Missing input arguments")
})

test_that("Invalid inputs are dealt with", {
  expect_error(ars(100, 100), "f must be a function")
  expect_error(ars(dnorm, dnorm), "N must be a numeric input") 
  expect_error(ars(dnorm, 10, bounds=2), "bounds must be length 2")
  expect_error(ars(dnorm, 10, x0 = c(1,2,3)), "x0 must be length 2")
  expect_error(ars(dnorm, 10, x0 = c(-1, 5), bounds = c(0, 0.1)), "x0 must be inside the bounds")
})

test_that("Inputting a log function", {
  h <- function(x) log(dnorm(x))
  expect_error(ars(h, 10), bounds = c(0,100))
})

