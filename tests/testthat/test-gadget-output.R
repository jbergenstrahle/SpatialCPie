context("Gadget output")

library(shinytest)

test_that("Gadget output hasn't changed", {
    skip_on_bioc()
    expect_pass(testApp("./cases/default/", compareImages = FALSE))
})
