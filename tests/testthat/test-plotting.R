test_that("plot_allele_frequencies creates ggplot object", {
  # Create sample frequency matrix
  freq_matrix <- matrix(
    c(0.5, 0.3, 0.2, 0.4, 0.4, 0.2),
    nrow = 2,
    ncol = 3,
    dimnames = list(c("Pop1", "Pop2"), c("a1", "a2", "b1"))
  )
  
  p <- plot_allele_frequencies(freq_matrix)
  
  expect_s3_class(p, "ggplot")
})

test_that("plot_allele_frequencies validates input", {
  # Missing column names
  bad_matrix <- matrix(c(0.5, 0.5), nrow = 1)
  expect_error(
    plot_allele_frequencies(bad_matrix),
    "must have column names"
  )
  
  # Not a matrix or data frame
  expect_error(
    plot_allele_frequencies(c(1, 2, 3)),
    "must be a matrix or data frame"
  )
})

test_that("x_axis_angle parameter works", {
  freq_matrix <- matrix(
    c(0.5, 0.5),
    nrow = 1,
    dimnames = list("Pop1", c("a1", "a2"))
  )
  
  # Default 90 degrees
  p1 <- plot_allele_frequencies(freq_matrix)
  expect_s3_class(p1, "ggplot")
  
  # Disabled rotation
  p2 <- plot_allele_frequencies(freq_matrix, x_axis_angle = 0)
  expect_s3_class(p2, "ggplot")
  
  # Custom angle
  p3 <- plot_allele_frequencies(freq_matrix, x_axis_angle = 45)
  expect_s3_class(p3, "ggplot")
})