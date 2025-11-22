test_that("check_chromosome_structure identifies columns", {
  test_data <- data.frame(
    who = 1:5,
    t.dad.chromosome.1.a = rep("a1", 5),
    t.dad.chromosome.1.b = rep("b1", 5),
    t.mom.chromosome.1.a = rep("a2", 5)
  )
  
  result <- check_chromosome_structure(test_data)
  
  expect_type(result, "list")
  expect_true("expanded" %in% names(result))
  expect_true("loci" %in% names(result))
  expect_true(length(result$loci) > 0)
})

test_that("extract_genotypes handles missing individual", {
  test_data <- data.frame(
    who = 1:3,
    t.dad.chromosome.1.a = c("a1", "a2", "a3")
  )
  
  expect_error(
    extract_genotypes(test_data, individual_id = 999),
    "not found"
  )
})

test_that("extract_genotypes creates proper genotypes", {
  test_data <- data.frame(
    who = c(1),
    t.dad.chromosome.1.a = "a1",
    t.dad.chromosome.2.a = "a2",
    t.mom.chromosome.1.a = "a3",
    t.mom.chromosome.2.a = "a4"
  )
  
  result <- extract_genotypes(test_data, individual_id = 1, loci_names = "a")
  
  expect_type(result, "character")
  expect_equal(names(result), "a")
  expect_true(grepl("/", result["a"]))
})