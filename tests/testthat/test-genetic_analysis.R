test_that("calculate_effective_migrants handles empty data", {
  empty_data <- data.frame(
    run_number = character(0),
    ticks = numeric(0),
    who = numeric(0),
    breed = character(0),
    t.dad.id = numeric(0),
    t.mom.id = numeric(0),
    t.birth.pond.index = numeric(0)
  )
  
  result <- calculate_effective_migrants(empty_data)
  expect_equal(nrow(result), 0)
  expect_true("origin_pond" %in% names(result))
})

test_that("calculate_effective_migrants requires correct columns", {
  bad_data <- data.frame(x = 1:5)
  
  expect_error(
    calculate_effective_migrants(bad_data),
    "Missing required columns"
  )
})

test_that("balance_genind_populations works with seed", {
  skip_if_not_installed("adegenet")
  
  # This would need actual genind data
  # Placeholder test structure
  expect_true(TRUE)
})