test_that("assign_turtles_to_ponds validates input", {
  turtle_data <- data.frame(
    who = 1:5,
    xcor = c(1, 2, 3, 4, 5),
    ycor = c(1, 2, 3, 4, 5)
  )
  
  pond_centers <- data.frame(
    pond_id = c(0, 1),
    x = c(1, 5),
    y = c(1, 5)
  )
  
  # Missing coordinates
  expect_error(
    assign_turtles_to_ponds(data.frame(who = 1:5), pond_centers),
    "must contain 'xcor' and 'ycor'"
  )
  
  # Wrong number of columns in pond_centers
  expect_error(
    assign_turtles_to_ponds(turtle_data, data.frame(x = 1)),
    "must have at least 3 columns"
  )
  
  # Not a data frame
  expect_error(
    assign_turtles_to_ponds(c(1, 2, 3), pond_centers),
    "must be a data.frame"
  )
})

test_that("assign_turtles_to_ponds assigns correctly", {
  skip_if_not_installed("FNN")
  
  turtle_data <- data.frame(
    who = 1:4,
    xcor = c(1, 1.2, 5, 5.1),
    ycor = c(1, 1.1, 5, 5.2)
  )
  
  pond_centers <- data.frame(
    pond_id = c(0, 1),
    x = c(1, 5),
    y = c(1, 5)
  )
  
  result <- assign_turtles_to_ponds(turtle_data, pond_centers)
  
  expect_true("corrected.current.pond" %in% names(result))
  expect_equal(result$corrected.current.pond[1], 0)
  expect_equal(result$corrected.current.pond[2], 0)
  expect_equal(result$corrected.current.pond[3], 1)
  expect_equal(result$corrected.current.pond[4], 1)
})

test_that("calculate_turtle_density computes correctly", {
  turtle_data <- data.frame(
    xcor = c(1.1, 1.2, 1.9, 5.1, 5.2),
    ycor = c(1.1, 1.2, 1.8, 5.1, 5.3),
    ticks = c(1, 1, 1, 2, 2)
  )
  
  result <- calculate_turtle_density(turtle_data)
  
  expect_true("mean_density" %in% names(result))
  expect_true("max_density" %in% names(result))
  expect_true(all(result$mean_density >= 0))
})