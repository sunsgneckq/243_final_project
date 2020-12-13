
# tests for the intercept_z_j function

test_that("intercept_z_j works with concave h(x)", {
  dx <- 1E-8
  set.seed(123)
  x <- runif(10,0,2)
  hx <- dnorm(x, log=T)
  dhx <- (dnorm(x+dx, log=T) - dnorm(x, log=T))/dx
  bounds <- c(-Inf, Inf)
  z_2 <- (hx[3] - hx[2] - x[3]*dhx[3] + x[2]*dhx[2]) / (dhx[2] - dhx[3])
  
  expect_length(intercept_z_j(x,hx,dhx,bounds), length(x) + 1)
  expect_equal(intercept_z_j(x, hx, dhx, bounds)[3], z_2)
  expect_true(all(intercept_z_j(x, hx, dhx, bounds) >= 
                    intercept_z_j(x, hx, dhx, bounds)[1]) &
                all(intercept_z_j(x, hx, dhx, bounds) <= 
                      intercept_z_j(x, hx, dhx, bounds)[11])) #intercepts are between bounds
  expect_error(intercept_z_j(x,hx,dhx,c(-Inf, 0)), "x values not within bounds") # don't need if this will never happen
})

