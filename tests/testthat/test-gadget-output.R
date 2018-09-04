context("Gadget output")

library(shinytest)


test_that("Gadget output hasn't changed", {
  expect_pass(testApp("./cases/default/", compareImages = FALSE))
})
