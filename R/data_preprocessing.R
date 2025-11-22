################################################################################
# R/data_preprocessing.R
################################################################################
# 
# Functions for cleaning, transforming, and preparing ABM simulation data
#
################################################################################

################################################################################
# FUNCTION: apply_breed_na_rules
################################################################################
#' Apply breed-specific NA assignments to turtle data
#'
#' Sets specific column values to NA based on the breed (life stage) of each
#' individual. This is essential for data cleaning where certain metrics are
#' not applicable to particular life stages.
#'
#' @param turtle_data A data.frame or tibble containing simulation data. Must include
#'   a `breed` column with values "eggs", "larvae", "juveniles", "adults" and
#'   the columns to be modified.
#'
#' @return A new data.frame with specified column values set to NA based on
#'   breed-specific rules. The original data.frame remains unchanged.
#'
#' @details The function applies the following NA assignment rules:
#' \describe{
#'   \item{t.age}{NA for "eggs"}
#'   \item{t.age.metamorphosis}{NA for "eggs" and "larvae"}
#'   \item{t.energy}{NA for "eggs"}
#'   \item{t.energy.consumption}{NA for "eggs"}
#'   \item{t.hatch.countdown}{NA for "larvae", "juveniles", and "adults"}
#'   \item{t.max.size}{NA for "eggs" and "larvae"}
#'   \item{t.metamorphic.risk}{NA for "eggs", "juveniles", and "adults"}
#'   \item{t.size.cm}{NA for "eggs"}
#'   \item{t.size.cm.metamorphosis}{NA for "eggs" and "larvae"}
#'   \item{t.stress}{NA for "eggs"}
#' }
#'
#' @examples
#' # Create sample turtle data
#' sample_data <- data.frame(
#'   ticks = 1:10,
#'   who = 1:10,
#'   breed = sample(c("eggs", "larvae", "juveniles", "adults"), 10, replace = TRUE),
#'   t.age = sample(1:10, 10, replace = TRUE),
#'   t.energy = rnorm(10, 50, 10),
#'   t.hatch.countdown = sample(1:5, 10, replace = TRUE)
#' )
#' 
#' # Apply NA rules
#' cleaned_data <- apply_breed_na_rules(sample_data)
#' 
#' # Check that eggs have NA for t.age and t.energy
#' subset(cleaned_data, breed == "eggs", select = c(breed, t.age, t.energy))
#'
#' @importFrom dplyr mutate if_else
#' @export
apply_breed_na_rules <- function(turtle_data) {
  # Input validation
  if (!is.data.frame(turtle_data)) {
    stop("turtle_data must be a data.frame or tibble")
  }
  
  if (!"breed" %in% names(turtle_data)) {
    stop("turtle_data must contain a 'breed' column")
  }
  
  # Check for valid breed values
  valid_breeds <- c("eggs", "larvae", "juveniles", "adults")
  invalid_breeds <- setdiff(unique(turtle_data$breed), valid_breeds)
  if (length(invalid_breeds) > 0) {
    warning("Found unexpected breed values: ", paste(invalid_breeds, collapse = ", "))
  }
  
  # Apply breed-specific NA rules only to columns that exist
  cleaned_data <- turtle_data
  
  # Only apply rules to columns that actually exist in the data
  if ("t.age" %in% names(turtle_data)) {
    cleaned_data$t.age <- ifelse(cleaned_data$breed == "eggs", NA_real_, cleaned_data$t.age)
  }
  
  if ("t.age.metamorphosis" %in% names(turtle_data)) {
    cleaned_data$t.age.metamorphosis <- ifelse(
      cleaned_data$breed %in% c("eggs", "larvae"), 
      NA_real_, 
      cleaned_data$t.age.metamorphosis
    )
  }
  
  if ("t.energy" %in% names(turtle_data)) {
    cleaned_data$t.energy <- ifelse(cleaned_data$breed == "eggs", NA_real_, cleaned_data$t.energy)
  }
  
  if ("t.energy.consumption" %in% names(turtle_data)) {
    cleaned_data$t.energy.consumption <- ifelse(
      cleaned_data$breed == "eggs", 
      NA_real_, 
      cleaned_data$t.energy.consumption
    )
  }
  
  if ("t.hatch.countdown" %in% names(turtle_data)) {
    cleaned_data$t.hatch.countdown <- ifelse(
      cleaned_data$breed %in% c("larvae", "juveniles", "adults"), 
      NA_real_, 
      cleaned_data$t.hatch.countdown
    )
  }
  
  if ("t.metamorphic.risk" %in% names(turtle_data)) {
    cleaned_data$t.metamorphic.risk <- ifelse(
      cleaned_data$breed %in% c("eggs", "juveniles", "adults"), 
      NA_real_, 
      cleaned_data$t.metamorphic.risk
    )
  }
  
  if ("t.size.cm" %in% names(turtle_data)) {
    cleaned_data$t.size.cm <- ifelse(cleaned_data$breed == "eggs", NA_real_, cleaned_data$t.size.cm)
  }
  
  if ("t.size.cm.metamorphosis" %in% names(turtle_data)) {
    cleaned_data$t.size.cm.metamorphosis <- ifelse(
      cleaned_data$breed %in% c("eggs", "larvae"), 
      NA_real_, 
      cleaned_data$t.size.cm.metamorphosis
    )
  }
  
  if ("t.max.size" %in% names(turtle_data)) {
    cleaned_data$t.max.size <- ifelse(
      cleaned_data$breed %in% c("eggs", "larvae"), 
      NA_real_, 
      cleaned_data$t.max.size
    )
  }
  
  if ("t.stress" %in% names(turtle_data)) {
    cleaned_data$t.stress <- ifelse(cleaned_data$breed == "eggs", NA_real_, cleaned_data$t.stress)
  }
  
  return(cleaned_data)
}


################################################################################
# FUNCTION: validate_turtle_data
################################################################################
#' Validate turtle data structure and content
#'
#' Performs comprehensive validation of turtle simulation data, checking for
#' required columns, valid value ranges, and data consistency.
#'
#' @param turtle_data A data.frame containing turtle simulation data.
#' @param required_columns Character vector of required column names. If NULL,
#'   uses a standard set of expected columns.
#' @param strict_validation Whether to apply strict validation rules. Default: FALSE.
#'
#' @return A list containing validation results with components:
#'   \describe{
#'     \item{valid}{Logical indicating if data passed validation}
#'     \item{errors}{Character vector of error messages}
#'     \item{warnings}{Character vector of warning messages}
#'     \item{summary}{List of data summary statistics}
#'   }
#'
#' @examples
#' # Create sample data
#' turtle_data <- data.frame(
#'   ticks = rep(100, 5),
#'   who = 1:5,
#'   breed = c("eggs", "larvae", "juveniles", "adults", "adults"),
#'   t.age = c(NA, 10, 20, 30, 25),
#'   xcor = rnorm(5),
#'   ycor = rnorm(5)
#' )
#' 
#' # Validate the data
#' validation_result <- validate_turtle_data(turtle_data)
#' validation_result$valid
#' validation_result$warnings
#'
#' @export
validate_turtle_data <- function(turtle_data, required_columns = NULL, strict_validation = FALSE) {
  errors <- character(0)
  warnings <- character(0)
  
  # Basic structure checks
  if (!is.data.frame(turtle_data)) {
    errors <- c(errors, "Input must be a data.frame")
    return(list(valid = FALSE, errors = errors, warnings = warnings, summary = NULL))
  }
  
  if (nrow(turtle_data) == 0) {
    errors <- c(errors, "Data frame is empty")
  }
  
  # Define expected columns if not provided
  if (is.null(required_columns)) {
    required_columns <- c("ticks", "who", "breed", "xcor", "ycor")
  }
  
  # Check for required columns
  missing_cols <- setdiff(required_columns, names(turtle_data))
  if (length(missing_cols) > 0) {
    errors <- c(errors, paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Content validation (only if basic structure is valid)
  if (length(errors) == 0) {
    # Check breed values
    if ("breed" %in% names(turtle_data)) {
      valid_breeds <- c("eggs", "larvae", "juveniles", "adults")
      invalid_breeds <- setdiff(unique(turtle_data$breed), valid_breeds)
      if (length(invalid_breeds) > 0) {
        warnings <- c(warnings, paste("Unexpected breed values:", paste(invalid_breeds, collapse = ", ")))
      }
    }
    
    # Check for negative time values
    if ("ticks" %in% names(turtle_data)) {
      if (any(turtle_data$ticks < 0, na.rm = TRUE)) {
        warnings <- c(warnings, "Found negative tick values")
      }
    }
    
    # Check coordinate ranges (if strict validation)
    if (strict_validation) {
      coord_cols <- c("xcor", "ycor")
      for (col in coord_cols) {
        if (col %in% names(turtle_data)) {
          coord_range <- range(turtle_data[[col]], na.rm = TRUE)
          if (coord_range[1] < -1000 || coord_range[2] > 1000) {
            warnings <- c(warnings, paste("Extreme coordinate values in", col, ":", 
                                          round(coord_range[1], 2), "to", round(coord_range[2], 2)))
          }
        }
      }
    }
    
    # Check for duplicate individuals at same time point
    if (all(c("ticks", "who") %in% names(turtle_data))) {
      duplicates <- turtle_data %>%
        dplyr::group_by(ticks, who) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::filter(n > 1)
      
      if (nrow(duplicates) > 0) {
        warnings <- c(warnings, paste("Found", nrow(duplicates), "duplicate individual-time combinations"))
      }
    }
  }
  
  # Create summary statistics
  summary_stats <- list(
    n_rows = nrow(turtle_data),
    n_columns = ncol(turtle_data),
    time_range = if ("ticks" %in% names(turtle_data)) range(turtle_data$ticks, na.rm = TRUE) else NULL,
    n_individuals = if ("who" %in% names(turtle_data)) length(unique(turtle_data$who)) else NULL,
    breed_counts = if ("breed" %in% names(turtle_data)) table(turtle_data$breed) else NULL
  )
  
  # Determine if validation passed
  valid <- length(errors) == 0
  
  return(list(
    valid = valid,
    errors = errors,
    warnings = warnings,
    summary = summary_stats
  ))
}

################################################################################
# FUNCTION: combine_chromosome_columns
################################################################################
#' Combine expanded chromosome columns back to original format
#'
#' Reverses the chromosome expansion process by combining individual allele
#' columns back into bracketed chromosome strings. Useful for data export
#' or compatibility with NetLogo format.
#'
#' @param turtle_data A data.frame with expanded chromosome columns
#'   (e.g., t.dad.chromosome.1.a, t.dad.chromosome.1.b, etc.).
#'
#' @return A data.frame with original chromosome column format
#'   (e.g., t.dad.chromosome.1 = "[a1 b2 c3 d4 e5]").
#'
#' @examples
#' # Create data with expanded chromosome columns
#' expanded_data <- data.frame(
#'   who = 1:3,
#'   t.dad.chromosome.1.a = c("a1", "a2", "a1"),
#'   t.dad.chromosome.1.b = c("b1", "b1", "b2"),
#'   t.dad.chromosome.1.c = c("c1", "c2", "c1"),
#'   t.dad.chromosome.1.d = c("d1", "d1", "d2"),
#'   t.dad.chromosome.1.e = c("e1", "e2", "e1")
#' )
#' 
#' # Combine back to original format
#' combined_data <- combine_chromosome_columns(expanded_data)
#' combined_data$t.dad.chromosome.1
#'
#' @export
combine_chromosome_columns <- function(turtle_data) {
  if (!is.data.frame(turtle_data)) {
    stop("turtle_data must be a data.frame")
  }
  
  # Define chromosome prefixes to process
  chromosome_prefixes <- c("t.dad.chromosome.1", "t.dad.chromosome.2", 
                           "t.mom.chromosome.1", "t.mom.chromosome.2")
  
  # Process each chromosome
  for (prefix in chromosome_prefixes) {
    # Find columns for this chromosome
    chr_cols <- grep(paste0("^", prefix, "\\."), names(turtle_data), value = TRUE)
    
    if (length(chr_cols) > 0) {
      # Create bracketed chromosome string
      chr_matrix <- turtle_data[, chr_cols, drop = FALSE]
      
      turtle_data[[prefix]] <- apply(chr_matrix, 1, function(row) {
        # Remove NA values and combine
        alleles <- row[!is.na(row)]
        if (length(alleles) > 0) {
          paste0("[", paste(alleles, collapse = " "), "]")
        } else {
          NA_character_
        }
      })
      
      # Remove the expanded columns
      turtle_data <- turtle_data[, !names(turtle_data) %in% chr_cols, drop = FALSE]
    }
  }
  
  return(turtle_data)
}

################################################################################
# FUNCTION: create_analysis_subset
################################################################################
#' Create optimized data subsets for analysis
#'
#' Creates smaller, manageable subsets of large turtle datasets for analysis,
#' with options for temporal, replicate-based, or random sampling strategies.
#'
#' @param turtle_data Full turtle dataset.
#' @param subset_method Method for subsetting: "time", "replicate", "random", or "systematic".
#' @param subset_size For "time": number of time points or specific time values.
#'   For "replicate": number of replicates to sample.
#'   For "random"/"systematic": number of total rows.
#' @param seed Random seed for reproducible sampling. Default: NULL.
#'
#' @return A subset of the input data based on the specified method.
#'
#' @examples
#' # Create sample data
#' large_data <- data.frame(
#'   ticks = rep(1:100, each = 50),
#'   who = rep(1:50, 100),
#'   run_number = rep(paste("Run", 1:5), each = 1000),
#'   value = rnorm(5000)
#' )
#' 
#' # Sample 10 evenly spaced time points
#' time_subset <- create_analysis_subset(large_data, "time", 10)
#' length(unique(time_subset$ticks))
#' 
#' # Sample 2 random replicates
#' rep_subset <- create_analysis_subset(large_data, "replicate", 2)
#' length(unique(rep_subset$run_number))
#' 
#' # Random sample of 1000 rows
#' random_subset <- create_analysis_subset(large_data, "random", 1000)
#' nrow(random_subset)
#'
#' @importFrom dplyr filter
#' @export
create_analysis_subset <- function(turtle_data, subset_method, subset_size = NULL, seed = NULL) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Validate inputs
  valid_methods <- c("time", "replicate", "random", "systematic")
  if (!subset_method %in% valid_methods) {
    stop("subset_method must be one of: ", paste(valid_methods, collapse = ", "))
  }
  
  if (subset_method == "time") {
    if ("ticks" %in% names(turtle_data)) {
      unique_times <- sort(unique(turtle_data$ticks))
      
      if (is.null(subset_size)) {
        # Default to 10 evenly spaced time points
        n_times <- min(10, length(unique_times))
        selected_times <- unique_times[round(seq(1, length(unique_times), length.out = n_times))]
      } else if (length(subset_size) == 1 && subset_size <= length(unique_times)) {
        # Evenly spaced time points
        selected_times <- unique_times[round(seq(1, length(unique_times), length.out = subset_size))]
      } else {
        # Specific time values
        selected_times <- intersect(subset_size, unique_times)
        if (length(selected_times) == 0) {
          warning("No specified time points found in data")
          return(turtle_data[FALSE, ])  # Return empty data.frame with same structure
        }
      }
      
      return(turtle_data %>% dplyr::filter(ticks %in% selected_times))
    } else {
      warning("No 'ticks' column found for time-based subsetting")
      return(turtle_data)
    }
    
  } else if (subset_method == "replicate") {
    if ("run_number" %in% names(turtle_data)) {
      unique_runs <- unique(turtle_data$run_number)
      
      if (is.null(subset_size)) {
        subset_size <- min(3, length(unique_runs))
      }
      
      if (subset_size >= length(unique_runs)) {
        return(turtle_data)
      }
      
      selected_runs <- sample(unique_runs, subset_size)
      return(turtle_data %>% dplyr::filter(run_number %in% selected_runs))
    } else {
      warning("No 'run_number' column found for replicate-based subsetting")
      return(turtle_data)
    }
    
  } else if (subset_method == "random") {
    if (is.null(subset_size)) {
      subset_size <- min(100000, nrow(turtle_data))
    }
    
    if (subset_size >= nrow(turtle_data)) {
      return(turtle_data)
    }
    
    sampled_rows <- sample(nrow(turtle_data), subset_size)
    return(turtle_data[sampled_rows, ])
    
  } else if (subset_method == "systematic") {
    if (is.null(subset_size)) {
      subset_size <- min(100000, nrow(turtle_data))
    }
    
    if (subset_size >= nrow(turtle_data)) {
      return(turtle_data)
    }
    
    # Systematic sampling
    interval <- floor(nrow(turtle_data) / subset_size)
    sampled_rows <- seq(1, nrow(turtle_data), by = interval)[1:subset_size]
    return(turtle_data[sampled_rows, ])
  }
}

################################################################################
# FUNCTION: estimate_memory_usage
################################################################################
#' Estimate memory usage for turtle data processing
#'
#' Provides memory usage estimates for turtle data files to help plan
#' data processing strategies for large datasets.
#'
#' @param file_path Path to turtle file or directory containing turtle files.
#' @param pattern File pattern if file_path is a directory. Default: "*turtle".
#' @param sample_rows Number of rows to sample for size estimation. Default: 1000.
#'
#' @return A data.frame with columns: file, rows, size_MB, estimated_memory_MB.
#'
#' @examples
#' \dontrun{
#' # Estimate memory for files in a directory
#' memory_est <- estimate_memory_usage("simulation_output/")
#' print(memory_est)
#' 
#' # Check if total memory requirement is reasonable
#' total_memory_gb <- sum(memory_est$estimated_memory_MB) / 1024
#' cat("Total estimated memory:", round(total_memory_gb, 2), "GB\n")
#' }
#'
#' @importFrom utils object.size
#' @export
estimate_memory_usage <- function(file_path, pattern = "*turtle", sample_rows = 1000) {
  if (file.exists(file_path) && !dir.exists(file_path)) {
    # Single file
    files <- file_path
  } else if (dir.exists(file_path)) {
    # Directory
    files <- list.files(file_path, pattern = pattern, full.names = TRUE)
  } else {
    stop("file_path must be an existing file or directory")
  }
  
  if (length(files) == 0) {
    stop("No files found matching pattern '", pattern, "'")
  }
  
  # Estimate for each file
  estimates <- lapply(files, function(file) {
    tryCatch({
      # Count total lines
      total_lines <- length(readLines(file, warn = FALSE))
      total_rows <- max(0, total_lines - 1)  # Subtract header
      
      if (total_rows == 0) {
        return(data.frame(
          file = basename(file),
          rows = 0,
          size_MB = 0,
          estimated_memory_MB = 0,
          stringsAsFactors = FALSE
        ))
      }
      
      # Sample a subset to estimate memory per row
      sample_size <- min(sample_rows, total_rows)
      sample_data <- read_turtle_space_format(file, sample_size, "head", NULL)
      
      # Calculate memory per row
      memory_per_row <- as.numeric(utils::object.size(sample_data)) / nrow(sample_data)
      
      # Estimate total memory
      estimated_memory <- memory_per_row * total_rows
      
      # File size
      file_size <- file.info(file)$size
      
      data.frame(
        file = basename(file),
        rows = total_rows,
        size_MB = round(file_size / (1024^2), 2),
        estimated_memory_MB = round(estimated_memory / (1024^2), 2),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      warning("Error processing file ", basename(file), ": ", e$message)
      data.frame(
        file = basename(file),
        rows = NA,
        size_MB = NA,
        estimated_memory_MB = NA,
        stringsAsFactors = FALSE
      )
    })
  })
  
  # Combine results
  result <- do.call(rbind, estimates)
  
  # Add summary row
  totals <- data.frame(
    file = "TOTAL",
    rows = sum(result$rows, na.rm = TRUE),
    size_MB = sum(result$size_MB, na.rm = TRUE),
    estimated_memory_MB = sum(result$estimated_memory_MB, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  result <- rbind(result, totals)
  
  return(result)
}