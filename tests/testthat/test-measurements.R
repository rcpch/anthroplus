library(testthat)
library(anthroplus)

test_that("anthroplus_measurements returns median for z = 0", {
  res <- anthroplus_measurements(sex = 1,
                                 age_in_months = 120,
                                 requested_z = 0,
                                 measurement_method = "length",
                                 measurement_precision = 2)
  expect_equal(as.numeric(res[[1]]), 137.78, tolerance = 1e-2)
})

test_that("anthroplus_measurements correct_extreme behaviour", {
  res_no <- anthroplus_measurements(sex = 1,
                                    age_in_months = 120,
                                    requested_z = 4,
                                    measurement_method = "length",
                                    measurement_precision = 2,
                                    correct_extreme = FALSE)

  res_yes <- anthroplus_measurements(sex = 1,
                                     age_in_months = 120,
                                     requested_z = 4,
                                     measurement_method = "length",
                                     measurement_precision = 2,
                                     correct_extreme = TRUE)

  expect_true(is.na(as.numeric(res_no[[1]])))
  expect_equal(as.numeric(res_yes[[1]]), 163.27, tolerance = 1e-2)
})
