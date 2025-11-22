test_that("apply_breed_na_rules works correctly", {
  # Create sample data
  test_data <- data.frame(
    breed = c("eggs", "larvae", "juveniles", "adults"),
    t.age = c(1, 2, 3, 4),
    t.energy = c(10, 20, 30, 40),
    t.hatch.countdown = c(5, NA, NA, NA),
    t.size.cm = c(1, 2, 3, 4)
  )
  
  result <- apply_breed_na_rules(test_data)
  
  # Test that eggs have NA for t.age and t.energy
  expect_true(is.na(result$t.age[1]))
  expect_true(is.na(result$t.energy[1]))
  
  # Test that larvae don't have hatch.countdown
  expect_true(is.na(result$t.hatch.countdown[2]))
  
  # Test that adults keep their values
  expect_equal(result$t.age[4], 4)
  expect_equal(result$t.energy[4], 40)
})

test_that("apply_breed_na_rules validates input", {
  # Missing breed column
  expect_error(
    apply_breed_na_rules(data.frame(x = 1:5)),
    "must contain a 'breed' column"
  )
  
  # Not a data frame
  expect_error(
    apply_breed_na_rules(c(1, 2, 3)),
    "must be a data.frame"
  )
})

test_that("validate_turtle_data detects issues", {
  # Valid data
  valid_data <- data.frame(
    ticks = 1:5,
    who = 1:5,
    breed = rep("adults", 5),
    xcor = rnorm(5),
    ycor = rnorm(5)
  )
  
  result <- validate_turtle_data(valid_data)
  expect_true(result$valid)
  expect_length(result$errors, 0)
  
  # Missing required column
  invalid_data <- data.frame(x = 1:5)
  result <- validate_turtle_data(invalid_data)
  expect_false(result$valid)
  expect_true(length(result$errors) > 0)
})

test_that("create_analysis_subset works with different methods", {
  test_data <- data.frame(
    ticks = rep(1:10, each = 10),
    who = 1:100,
    run_number = rep("Run 1", 100),
    value = rnorm(100)
  )
  
  # Time subset
  time_subset <- create_analysis_subset(test_data, "time", 5)
  expect_equal(length(unique(time_subset$ticks)), 5)
  
  # Random subset
  random_subset <- create_analysis_subset(test_data, "random", 50, seed = 123)
  expect_equal(nrow(random_subset), 50)
  
  # Systematic subset
  sys_subset <- create_analysis_subset(test_data, "systematic", 50)
  expect_equal(nrow(sys_subset), 50)
})