test_that("wrangle_differentiation_data validates measure type", {
  test_list <- list(
    "Run 1" = list(
      Gst = list(Gst_Est = list(per.locus = c(a = 0.1, b = 0.2)))
    )
  )
  
  result <- wrangle_differentiation_data(test_list, "Gst")
  expect_s3_class(result, "tbl_df")
  expect_true("Gst" %in% names(result))
  
  # Invalid measure type
  expect_error(
    wrangle_differentiation_data(test_list, "Invalid"),
    "Invalid 'measure_type'"
  )
})