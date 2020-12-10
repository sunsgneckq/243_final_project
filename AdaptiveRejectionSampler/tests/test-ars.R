#source("R/args.R")

test_that("Missing arguments are dealt with", {
  expect_error(ars(), "Missing input arguments")
  expect_error(ars(dnorm), "Missing input arguments")
  expect_error(ars(N=100), "Missing input arguments")
})

test_that("Input arguments are valid", {
  expect_error(ars(10,10), "f must be a function")
  expect_error(ars(dnorm,dnorm), "f must be a function") #check is.numeric(N)
})