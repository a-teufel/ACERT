################################################################################
# R/data_io.R
################################################################################
# 
# Functions for reading and writing ABM simulation data files
#
################################################################################

################################################################################
# FUNCTION: read_turtles
################################################################################
#' Read turtle files from a directory
#'
#' Reads and combines turtle simulation files from a specified directory with
#' options for sampling, parallel processing, and memory optimization for
#' large datasets. Handles space-separated format with bracketed chromosome data.
#'
#' @param path Directory path containing turtle files. Default: current directory (".").
#' @param pattern File pattern to match using wildcard notation. Default: "*turtle".
#' @param sample_size Number of rows to sample per file. NULL reads all rows. 
#'   Useful for quick analysis of large datasets. Default: NULL.
#' @param sample_method Sampling method to use when \code{sample_size} is specified.
#'   Options: "random" (random sampling), "systematic" (every nth row), 
#'   "head" (first n rows), "tail" (last n rows). Default: "systematic".
#' @param time_points Specific time points (ticks) to keep. NULL keeps all time points.
#'   Useful for analyzing specific simulation stages. Default: NULL.
#' @param chunk_size Number of rows to read at a time for large files. 
#'   Currently unused but reserved for future chunked reading. Default: 100000.
#' @param parallel Whether to use parallel processing for reading multiple files.
#'   Uses number of cores minus 1. Default: FALSE.
#' @param progress Whether to show progress messages during reading. Default: TRUE.
#'
#' @return A data.table with combined turtle data including a \code{run_number} column
#'   that identifies which file each row came from.
#'
#' @details The function automatically:
#' \itemize{
#'   \item Detects file format (space-separated with bracketed chromosomes)
#'   \item Expands chromosome data into individual allele columns
#'   \item Standardizes column types (numeric/character)
#'   \item Extracts run numbers from filenames
#'   \item Combines data from multiple files with consistent structure
#' }
#'
#' Files are expected to follow NetLogo turtle output format with columns for
#' agent properties, coordinates, and genetic data. Chromosome columns in 
#' bracketed format (e.g., "[a1 b2 c1 d2 e1]") are automatically expanded
#' into separate allele columns.
#'
#' @examples
#' \dontrun{
#' # Basic usage: Read all turtle files from current directory
#' turtle_data <- read_turtles()
#' 
#' # Read from specific directory
#' turtle_data <- read_turtles(path = "simulation_output/run1/")
#' 
#' # Read with custom pattern matching
#' turtle_data <- read_turtles(
#'   path = "data/",
#'   pattern = "*_turtles.txt"
#' )
#' 
#' # Sample for quick analysis of large files
#' turtle_sample <- read_turtles(
#'   path = "large_simulation/",
#'   sample_size = 10000,
#'   sample_method = "systematic"  # Every nth row
#' )
#' 
#' # Random sampling
#' turtle_random <- read_turtles(
#'   sample_size = 5000,
#'   sample_method = "random"
#' )
#' 
#' # Read only final time points
#' final_data <- read_turtles(
#'   time_points = c(4000, 5000),
#'   progress = TRUE
#' )
#' 
#' # Read specific time series
#' time_series <- read_turtles(
#'   time_points = seq(0, 5000, by = 500)
#' )
#' 
#' # Use parallel processing for multiple files
#' turtle_data_fast <- read_turtles(
#'   path = "many_runs/",
#'   parallel = TRUE,
#'   progress = TRUE
#' )
#' 
#' # Combine filtering and sampling
#' filtered_sample <- read_turtles(
#'   path = "results/",
#'   time_points = c(1000, 2000, 3000),
#'   sample_size = 20000,
#'   sample_method = "systematic",
#'   parallel = TRUE
#' )
#' 
#' # Read just the beginning of files (quick preview)
#' preview <- read_turtles(
#'   sample_size = 100,
#'   sample_method = "head"
#' )
#' 
#' # Check data structure after reading
#' str(turtle_data)
#' head(turtle_data)
#' table(turtle_data$run_number)
#' }
#'
#' @importFrom data.table fread rbindlist setnames set
#' @importFrom parallel detectCores makeCluster stopCluster clusterEvalQ parLapply
#' @export
read_turtles <- function(path = ".", 
                         pattern = "*turtle",
                         sample_size = NULL,
                         sample_method = "systematic",
                         time_points = NULL,
                         chunk_size = 100000,
                         parallel = FALSE,
                         progress = TRUE) {
  
  # Input validation
  if (!dir.exists(path)) {
    stop("Directory '", path, "' does not exist")
  }
  
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No turtle files found matching pattern '", pattern, "' in directory '", path, "'")
  }
  
  # Process files
  process_file <- function(file, file_index) {
    if (progress && !parallel) {
      message("Reading file ", file_index, " of ", length(files), ": ", basename(file))
    }
    
    dt <- read_turtle_space_format(file, sample_size, sample_method, time_points)
    dt$run_number <- extract_run_number(basename(file))
    return(dt)
  }
  
  # Execute with or without parallel processing
  all_turtle_data <- if (parallel) {
    cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(data.table))
    parallel::parLapply(cl, seq_along(files), function(i) process_file(files[i], i))
  } else {
    lapply(seq_along(files), function(i) process_file(files[i], i))
  }
  
  # Combine all data
  combined_data <- data.table::rbindlist(all_turtle_data, fill = TRUE)
  
  if (progress) {
    message("Total rows loaded: ", format(nrow(combined_data), big.mark = ","))
    message("Memory used: ", format(object.size(combined_data), units = "MB"))
  }
  
  return(combined_data)
}

################################################################################
# FUNCTION: read_turtle_space_format
################################################################################
#' Read turtle file in space-separated format with bracketed chromosomes
#'
#' Internal function for reading turtle files in space-separated format where
#' chromosome data is provided in bracketed arrays that need expansion.
#'
#' @param file Path to turtle file.
#' @param sample_size Number of rows to sample.
#' @param sample_method Sampling method.
#' @param time_points Time points to keep.
#'
#' @return A data.table with expanded chromosome columns.
#'
#' @keywords internal
read_turtle_space_format <- function(file, sample_size = NULL, sample_method = "systematic", time_points = NULL) {
  # Standard header for space-separated format
  standard_header <- c(
    "ticks", "who", "breed", "t.age", "t.age.metamorphosis", "t.birth.pond.index", "t.current.pond.index",
    paste0("t.dad.chromosome.1.", letters[1:5]),
    paste0("t.dad.chromosome.2.", letters[6:10]),
    "t.dad.id", "t.energy", "t.energy.consumption", "t.hatch.countdown", "t.max.size", "t.metamorphic.risk",
    paste0("t.mom.chromosome.1.", letters[1:5]),
    paste0("t.mom.chromosome.2.", letters[6:10]),
    "t.mom.id", "t.sex", "t.size.cm", "t.size.cm.metamorphosis", "t.stress", "ycor", "xcor"
  )
  
  # Read raw data
  dt <- data.table::fread(file, sep = " ", skip = 1, header = FALSE, fill = TRUE, quote = "")
  
  # Adjust header to match actual columns
  actual_cols <- ncol(dt)
  if (actual_cols != length(standard_header)) {
    if (actual_cols < length(standard_header)) {
      standard_header <- standard_header[1:actual_cols]
    } else {
      extra_cols <- actual_cols - length(standard_header)
      standard_header <- c(standard_header, paste0("V", (length(standard_header) + 1):(length(standard_header) + extra_cols)))
    }
  }
  
  # Set column names
  data.table::setnames(dt, standard_header)
  
  # Clean chromosome columns (remove brackets)
  chr_cols <- grep("chromosome", names(dt), value = TRUE)
  for (col in chr_cols) {
    if (col %in% names(dt)) {
      data.table::set(dt, j = col, value = gsub("\\[|\\]", "", dt[[col]]))
    }
  }
  
  # Convert numeric columns
  dt <- standardize_column_types(dt)
  
  # Apply filtering and sampling
  dt <- apply_time_filter(dt, time_points)
  dt <- apply_sampling(dt, sample_size, sample_method, file)
  
  return(dt)
}


################################################################################
# FUNCTION: read_turtles_underscore_format
################################################################################
#' Read turtle files in underscore/hyphen format (newer version)
#'
#' Reads turtle files that use hyphenated column names and underscore-separated
#' chromosome data (e.g., "a1_b2_c1_d2_e1" instead of bracketed format).
#' This format is used in newer NetLogo exports.
#' 
#' @param path Directory path containing turtle files. Default: current directory (".").
#' @param pattern File pattern to match. Default: "*turtle".
#' @param sample_size Number of rows to sample per file. NULL for all rows. Default: NULL.
#' @param sample_method Sampling method: "random", "systematic", "head", "tail".
#'   Default: "systematic".
#' @param time_points Specific time points to keep. NULL for all. Default: NULL.
#' @param chunk_size Reserved for future use. Default: 100000.
#' @param parallel Use parallel processing for multiple files. Default: FALSE.
#' @param progress Show progress messages. Default: TRUE.
#'
#' @return A data.table with harmonized turtle data matching the space-format structure.
#'
#' @details This function handles the alternative turtle file format where:
#' \itemize{
#'   \item Column names use hyphens (e.g., "t-age" instead of "t.age")
#'   \item Chromosomes are underscore-separated strings (e.g., "a1_b2_c1_d2_e1")
#'   \item Output is standardized to match \code{read_turtles()} format
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage for underscore format files
#' turtle_data <- read_turtles_underscore_format()
#' 
#' # From specific directory
#' turtle_data <- read_turtles_underscore_format(
#'   path = "new_format_sims/"
#' )
#' 
#' # With sampling and filtering
#' turtle_sample <- read_turtles_underscore_format(
#'   path = "data/",
#'   sample_size = 15000,
#'   sample_method = "systematic",
#'   time_points = c(1000, 3000, 5000)
#' )
#' 
#' # Parallel processing for speed
#' turtle_data_fast <- read_turtles_underscore_format(
#'   path = "many_files/",
#'   parallel = TRUE
#' )
#' 
#' # Compare formats - columns will be identical after reading
#' space_format <- read_turtles(path = "old_sims/")
#' underscore_format <- read_turtles_underscore_format(path = "new_sims/")
#' 
#' # Both will have same column names
#' identical(names(space_format), names(underscore_format))
#' 
#' # Random sampling for quick preview
#' preview <- read_turtles_underscore_format(
#'   sample_size = 1000,
#'   sample_method = "random",
#'   progress = FALSE
#' )
#' }
#'
#' @export
read_turtles_underscore_format <- function(path = ".", 
                                           pattern = "*turtle",
                                           sample_size = NULL,
                                           sample_method = "systematic",
                                           time_points = NULL,
                                           chunk_size = 100000,
                                           parallel = FALSE,
                                           progress = TRUE) {
  
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No turtle files found in the specified directory")
  
  process_file <- function(file, file_index) {
    if (progress && !parallel) message(paste("Reading file", file_index, "of", length(files), ":", basename(file)))
    dt <- read_turtle_custom_underscore(file, sample_size, sample_method, time_points)
    dt$run_number <- sub("^(\\d+).*", "Run \\1", basename(file))
    return(dt)
  }
  
  all_turtle_data <- if (parallel) {
    cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(data.table))
    parallel::parLapply(cl, seq_along(files), function(i) process_file(files[i], i))
  } else {
    lapply(seq_along(files), function(i) process_file(files[i], i))
  }
  
  combined_data <- data.table::rbindlist(all_turtle_data, fill = TRUE)
  if (progress) {
    message(paste("Total rows loaded:", format(nrow(combined_data), big.mark = ",")))
    message(paste("Memory used:", format(object.size(combined_data), units = "MB")))
  }
  return(combined_data)
}


################################################################################
# FUNCTION: read_turtle_custom_underscore
################################################################################
#' Internal reader for underscore/hyphen turtle format (10 allele columns)
#'
#' @param file Path to file.
#' @param sample_size Number of rows to sample.
#' @param sample_method Sampling method.
#' @param time_points Time points to keep.
#'
#' @return A data.table with expanded chromosome columns.
#' 
#' @keywords internal
read_turtle_custom_underscore <- function(file, sample_size = NULL, sample_method = "systematic", time_points = NULL) {
  dt <- data.table::fread(file, sep = " ", header = TRUE, fill = TRUE, quote = "")
  
  # Rename columns to match expected camel.dot format
  col_rename_map <- list(
    "t-age" = "t.age", "t-age-metamorphosis" = "t.age.metamorphosis",
    "t-birth-pond-index" = "t.birth.pond.index", "t-current-pond-index" = "t.current.pond.index",
    "t-dad-chromosome-1" = "t.dad.chromosome.1", "t-dad-chromosome-2" = "t.dad.chromosome.2",
    "t-dad-id" = "t.dad.id", "t-energy" = "t.energy", "t-energy-consumption" = "t.energy.consumption",
    "t-hatch-countdown" = "t.hatch.countdown", "t-max-size" = "t.max.size", "t-metamorphic-risk" = "t.metamorphic.risk",
    "t-mom-chromosome-1" = "t.mom.chromosome.1", "t-mom-chromosome-2" = "t.mom.chromosome.2",
    "t-mom-id" = "t.mom.id", "t-sex" = "t.sex", "t-size-cm" = "t.size.cm",
    "t-size-cm-metamorphosis" = "t.size.cm.metamorphosis", "t-stress" = "t.stress",
    "ycor" = "ycor", "xcor" = "xcor", "ticks" = "ticks", "who" = "who", "breed" = "breed"
  )
  
  for (old in names(col_rename_map)) {
    if (old %in% names(dt)) {
      data.table::setnames(dt, old, col_rename_map[[old]])
    }
  }
  
  # Expand chromosomes
  expand_chromosome <- function(vec, suffixes) {
    cleaned <- gsub("\\[|\\]", "", vec)
    split <- strsplit(cleaned, "_")
    mat <- do.call(rbind, lapply(split, function(x) {
      x <- x[1:length(suffixes)]
      length(x) <- length(suffixes)
      x
    }))
    colnames(mat) <- suffixes
    return(mat)
  }
  
  # Expand all four chromosomes: 2 per parent, 5 alleles each
  chromosome_expansions <- list(
    "t.dad.chromosome.1" = letters[1:5],
    "t.dad.chromosome.2" = letters[6:10],
    "t.mom.chromosome.1" = letters[1:5],
    "t.mom.chromosome.2" = letters[6:10]
  )
  
  for (chromo in names(chromosome_expansions)) {
    if (chromo %in% names(dt)) {
      suffixes <- chromosome_expansions[[chromo]]
      mat <- expand_chromosome(dt[[chromo]], suffixes)
      for (i in seq_along(suffixes)) {
        dt[[paste0(chromo, ".", suffixes[i])]] <- mat[, i]
      }
      dt[[chromo]] <- NULL
    }
  }
  
  # Cast numeric columns
  numeric_cols <- c("ticks", "t.age", "t.age.metamorphosis", "t.birth.pond.index", "t.current.pond.index",
                    "t.size.cm", "t.energy", "t.energy.consumption", "t.hatch.countdown", "t.max.size",
                    "t.metamorphic.risk", "t.size.cm.metamorphosis", "t.stress", "ycor", "xcor")
  for (col in numeric_cols) {
    if (col %in% names(dt)) {
      suppressWarnings(data.table::set(dt, j = col, value = as.numeric(dt[[col]])))
    }
  }
  
  # Filter by time
  if (!is.null(time_points) && "ticks" %in% names(dt)) {
    dt <- dt[dt$ticks %in% time_points, ]
  }
  
  total_rows <- nrow(dt)
  if (total_rows == 0) {
    warning(sprintf("File '%s' had zero rows after filtering.", basename(file)))
    return(dt)
  }
  
  # Sampling
  if (!is.null(sample_size) && sample_size < total_rows) {
    if (sample_method == "systematic") {
      keep_every <- max(1, floor(total_rows / sample_size))
      indices <- seq(from = 1, to = total_rows, by = keep_every)
      dt <- dt[indices[1:min(length(indices), sample_size)], ]
    } else if (sample_method == "random") {
      dt <- dt[sample.int(total_rows, sample_size), ]
    } else if (sample_method == "head") {
      dt <- dt[1:min(sample_size, total_rows), ]
    } else if (sample_method == "tail") {
      start_row <- max(1, total_rows - sample_size + 1)
      dt <- dt[start_row:total_rows, ]
    } else {
      warning("Unknown sample method; returning all data.")
    }
  }
  
  return(dt)
}


################################################################################
# FUNCTION: read_environments
################################################################################
#' Read environment files from a directory
#'
#' Reads and combines environment simulation files containing patch-level
#' data such as energy, risk, permeability, and regeneration values.
#' Environment files describe the spatial landscape properties at each
#' time point.
#'
#' @param path Directory path containing environment files. Default: current directory (".").
#' @param pattern File pattern to match. Default: "environment".
#' @param sample_size Number of rows to sample per file. NULL for all rows.
#'   Useful for reducing memory usage with large spatial data. Default: NULL.
#' @param time_points Specific time points to keep. NULL for all time points. Default: NULL.
#' @param parallel Whether to use parallel processing. Default: FALSE.
#'
#' @return A data.table with combined environment data including columns for:
#'   \code{ticks}, \code{pxcor}, \code{pycor}, \code{patch_energy}, \code{patch_risk},
#'   \code{permeability}, \code{patch_regeneration}, \code{index}, and identification
#'   columns \code{file_name} and \code{run_number}.
#'
#' @details Environment files contain spatial patch data with properties like:
#' \itemize{
#'   \item Energy levels available in each patch
#'   \item Risk values affecting survival
#'   \item Permeability affecting movement
#'   \item Regeneration rates for resources
#'   \item Patch type index (habitat classification)
#' }
#'
#' @examples
#' \dontrun{
#' # Read all environment files
#' env_data <- read_environments()
#' 
#' # From specific directory
#' env_data <- read_environments(path = "simulation_results/")
#' 
#' # Read specific time points only
#' env_subset <- read_environments(
#'   time_points = c(500, 1000, 1500)
#' )
#' 
#' # Sample patches for large landscapes
#' env_sample <- read_environments(
#'   sample_size = 5000,  # Sample 5000 patches per file
#'   time_points = c(1000, 3000, 5000)
#' )
#' 
#' # Parallel reading for multiple files
#' env_data_fast <- read_environments(
#'   path = "many_runs/",
#'   parallel = TRUE
#' )
#' 
#' # Custom pattern matching
#' env_data <- read_environments(
#'   path = "output/",
#'   pattern = "*_environment.csv"
#' )
#' 
#' # Analyze environment data
#' library(dplyr)
#' 
#' # Mean energy by run and time
#' env_data %>%
#'   group_by(run_number, ticks) %>%
#'   summarise(mean_energy = mean(patch_energy, na.rm = TRUE))
#' 
#' # Risk distribution by patch type
#' env_data %>%
#'   group_by(index) %>%
#'   summarise(
#'     mean_risk = mean(patch_risk),
#'     sd_risk = sd(patch_risk)
#'   )
#' 
#' # Spatial visualization at final time point
#' library(ggplot2)
#' env_final <- read_environments(time_points = 5000)
#' ggplot(env_final, aes(x = pxcor, y = pycor, fill = patch_energy)) +
#'   geom_raster() +
#'   facet_wrap(~run_number)
#' }
#'
#' @importFrom data.table fread rbindlist
#' @importFrom R.utils countLines
#' @export
read_environments <- function(path = ".", 
                              pattern = "environment",
                              sample_size = NULL,
                              time_points = NULL,
                              parallel = FALSE) {
  
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No environment files found matching pattern '", pattern, "' in directory '", path, "'")
  }
  
  process_env_file <- function(file) {
    # Handle sampling if requested
    if (!is.null(sample_size)) {
      total_rows <- R.utils::countLines(file)[1] - 1
      keep_every <- max(1, floor(total_rows / sample_size))
      dt <- data.table::fread(file)
      dt <- dt[seq(1, nrow(dt), by = keep_every)]
    } else {
      dt <- data.table::fread(file)
    }
    
    # Filter by time points
    dt <- apply_time_filter(dt, time_points)
    
    # Convert numeric columns
    numeric_cols <- c("pycor", "pxcor")
    for (col in numeric_cols) {
      if (col %in% names(dt)) {
        dt[[col]] <- as.numeric(dt[[col]])
      }
    }
    
    # Add file identification
    dt$file_name <- basename(file)
    dt$run_number <- extract_run_number(basename(file))
    
    return(dt)
  }
  
  # Process files
  if (parallel) {
    cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(data.table))
    env_data_list <- parallel::parLapply(cl, files, process_env_file)
  } else {
    env_data_list <- lapply(files, process_env_file)
  }
  
  combined_data <- data.table::rbindlist(env_data_list, fill = TRUE)
  return(combined_data)
}

################################################################################
# FUNCTION: read_mutations
################################################################################
#' Read mutation files from a directory
#'
#' Reads and combines mutation event files containing information about
#' genetic changes that occurred during the simulation. Each row represents
#' a mutation event with information about the affected individual, locus,
#' and change in allelic values.
#'
#' @param path Directory path containing mutation files. Default: current directory (".").
#' @param pattern File pattern to match. Default: "mutations".
#' @param sample_size Number of rows to sample per file. NULL for all rows. Default: NULL.
#' @param time_points Specific time points to keep. NULL for all time points. Default: NULL.
#'
#' @return A data.table with combined mutation data including:
#'   \itemize{
#'     \item Original mutation columns from the file
#'     \item \code{file_name}: Source filename
#'     \item \code{run_number}: Extracted run identifier
#'     \item \code{FmI}: Calculated mutation effect size (if data available)
#'   }
#'
#' @details Mutation files typically contain:
#' \itemize{
#'   \item Time (ticks) when mutation occurred
#'   \item Individual ID (turtle)
#'   \item Old and new allele identifiers
#'   \item Old and new additivity values
#' }
#'
#' The function automatically calculates \code{FmI} (fitness mutation impact)
#' as the difference between new and old additivity values when these
#' columns are present.
#'
#' @examples
#' \dontrun{
#' # Read all mutation files
#' mut_data <- read_mutations()
#' 
#' # From specific directory
#' mut_data <- read_mutations(path = "genetic_data/")
#' 
#' # Sample for quick analysis
#' mut_sample <- read_mutations(sample_size = 10000)
#' 
#' # Specific time points
#' mut_subset <- read_mutations(
#'   time_points = c(1000, 2000, 3000, 4000, 5000)
#' )
#' 
#' # Custom pattern
#' mut_data <- read_mutations(
#'   path = "results/",
#'   pattern = "*_mutation_events"
#' )
#' 
#' # Analyze mutation data
#' library(dplyr)
#' 
#' # Mutation rate over time
#' mut_data %>%
#'   group_by(run_number, ticks) %>%
#'   summarise(n_mutations = n())
#' 
#' # Distribution of mutation effects
#' hist(mut_data$FmI, 
#'      main = "Distribution of Mutation Effects",
#'      xlab = "Effect Size (FmI)")
#' 
#' # Mean effect by run
#' mut_data %>%
#'   group_by(run_number) %>%
#'   summarise(
#'     mean_effect = mean(FmI, na.rm = TRUE),
#'     sd_effect = sd(FmI, na.rm = TRUE),
#'     n_mutations = n()
#'   )
#' 
#' # Beneficial vs deleterious mutations
#' mut_data %>%
#'   mutate(type = case_when(
#'     FmI > 0 ~ "Beneficial",
#'     FmI < 0 ~ "Deleterious",
#'     TRUE ~ "Neutral"
#'   )) %>%
#'   group_by(run_number, type) %>%
#'   summarise(count = n())
#' 
#' # Temporal pattern of large-effect mutations
#' mut_data %>%
#'   filter(abs(FmI) > 5) %>%
#'   ggplot(aes(x = ticks, y = FmI, color = run_number)) +
#'   geom_point(alpha = 0.5) +
#'   geom_hline(yintercept = 0, linetype = "dashed") +
#'   labs(title = "Large-Effect Mutations Over Time")
#' }
#'
#' @importFrom data.table fread rbindlist
#' @importFrom R.utils countLines
#' @export
read_mutations <- function(path = ".", 
                           pattern = "mutations",
                           sample_size = NULL,
                           time_points = NULL) {
  
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No mutation files found matching pattern '", pattern, "' in directory '", path, "'")
  }
  
  mut_data_list <- lapply(files, function(file) {
    # Handle sampling if requested
    if (!is.null(sample_size)) {
      total_rows <- R.utils::countLines(file)[1] - 1
      keep_every <- max(1, floor(total_rows / sample_size))
      dt <- data.table::fread(file)
      dt <- dt[seq(1, nrow(dt), by = keep_every)]
    } else {
      dt <- data.table::fread(file)
    }
    
    # Filter by time points
    dt <- apply_time_filter(dt, time_points)
    
    # Add file identification
    dt$file_name <- basename(file)
    dt$run_number <- extract_run_number(basename(file))
    
    return(dt)
  })
  
  combined_data <- data.table::rbindlist(mut_data_list, fill = TRUE)
  
  # Calculate mutation effect size if possible
  if (all(c("new_additivity_value", "old_additivity_value") %in% names(combined_data))) {
    combined_data[, FmI := new_additivity_value - old_additivity_value]
  }
  
  return(combined_data)
}

################################################################################
# FUNCTION: save_abm_data
################################################################################
#' Save processed ABM data to files
#'
#' Saves data in various formats optimized for different use cases.
#' Supports compression for large datasets and cross-platform formats
#' for sharing data with other analysis tools.
#'
#' @param data Data to save. Can be data.frame, data.table, or tibble.
#' @param filename Output filename with appropriate extension (e.g., "data.csv",
#'   "results.rds", "output.parquet").
#' @param format Format to save: "csv", "rds", "txt", "feather", "parquet".
#'   Auto-detected from filename extension if not specified. Default: NULL.
#' @param compress Whether to compress the output. Compression method varies
#'   by format. Default: TRUE.
#'
#' @return Invisibly returns the path to the saved file.
#'
#' @details Format characteristics:
#' \itemize{
#'   \item \strong{csv}: Human-readable, universally compatible, slower for large data
#'   \item \strong{rds}: R-specific, preserves all object attributes, fast
#'   \item \strong{txt}: Tab-separated, human-readable, good for text editors
#'   \item \strong{feather}: Fast, cross-platform (R/Python), preserves types
#'   \item \strong{parquet}: Efficient, cross-platform, good for big data, columnar storage
#' }
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' turtle_data <- read_turtles()
#' 
#' # Save as compressed CSV (most compatible)
#' save_abm_data(turtle_data, "processed_turtles.csv")
#' 
#' # Save as RDS for R-specific use (preserves all attributes)
#' save_abm_data(turtle_data, "processed_turtles.rds", compress = TRUE)
#' 
#' # Save as uncompressed CSV for viewing
#' save_abm_data(turtle_data, "readable_data.csv", compress = FALSE)
#' 
#' # Save as Parquet for cross-platform big data use
#' save_abm_data(turtle_data, "turtles.parquet")
#' 
#' # Save as Feather for fast R/Python interchange
#' save_abm_data(turtle_data, "turtles.feather")
#' 
#' # Explicitly specify format
#' save_abm_data(turtle_data, "output_file", format = "rds")
#' 
#' # Save processed subsets
#' library(dplyr)
#' adults <- turtle_data %>% filter(breed == "adults")
#' save_abm_data(adults, "adults_only.csv")
#' 
#' # Save summary statistics
#' summary_stats <- turtle_data %>%
#'   group_by(run_number, ticks, breed) %>%
#'   summarise(
#'     mean_energy = mean(t.energy),
#'     mean_size = mean(t.size.cm),
#'     count = n()
#'   )
#' save_abm_data(summary_stats, "summary_by_breed.rds")
#' 
#' # Save different formats for different purposes
#' save_abm_data(turtle_data, "for_r_analysis.rds")     # Fast, R-only
#' save_abm_data(turtle_data, "for_python.parquet")     # Cross-platform
#' save_abm_data(turtle_data, "for_excel.csv")          # Universal viewing
#' }
#'
#' @importFrom data.table fwrite
#' @importFrom utils write.table
#' @importFrom R.utils gzip
#' @export
save_abm_data <- function(data, filename, format = NULL, compress = TRUE) {
  # Auto-detect format from filename if not specified
  if (is.null(format)) {
    format <- detect_file_format(filename)
  }
  
  # Validate format
  valid_formats <- c("csv", "rds", "txt", "feather", "parquet")
  if (!format %in% valid_formats) {
    stop("Unsupported format '", format, "'. Use one of: ", paste(valid_formats, collapse = ", "))
  }
  
  # Save based on format
  switch(format,
         csv = {
           data.table::fwrite(data, filename, compress = ifelse(compress, "gzip", "none"))
         },
         rds = {
           saveRDS(data, filename, compress = compress)
         },
         txt = {
           utils::write.table(data, filename, sep = "\t", row.names = FALSE)
           if (compress) {
             R.utils::gzip(filename, overwrite = TRUE)
             filename <- paste0(filename, ".gz")
           }
         },
         feather = {
           if (!requireNamespace("arrow", quietly = TRUE)) {
             stop("Package 'arrow' needed for feather format. Please install it.")
           }
           arrow::write_feather(data, filename, 
                                compression = ifelse(compress, "lz4", "uncompressed"))
         },
         parquet = {
           if (!requireNamespace("arrow", quietly = TRUE)) {
             stop("Package 'arrow' needed for parquet format. Please install it.")
           }
           arrow::write_parquet(data, filename, 
                                compression = ifelse(compress, "snappy", "uncompressed"))
         }
  )
  
  message("Data saved to: ", filename)
  invisible(filename)
}

################################################################################
# FUNCTION: load_abm_data
################################################################################
#' Load processed ABM data efficiently
#'
#' Loads data from various formats with automatic format detection
#' based on file extension. Handles compressed files automatically.
#'
#' @param filename Input filename. Extension determines format unless
#'   \code{format} is explicitly specified.
#' @param format Format of the file: "csv", "rds", "txt", "feather", "parquet".
#'   Auto-detected from extension if NULL. Default: NULL.
#'
#' @return Loaded data as data.table (for most formats) or original object type (for RDS).
#'
#' @details Supports automatic decompression for:
#' \itemize{
#'   \item .csv.gz (gzipped CSV)
#'   \item .txt.gz (gzipped text)
#'   \item Compressed RDS files
#'   \item Compressed feather/parquet (handled by arrow package)
#' }
#'
#' @examples
#' \dontrun{
#' # Auto-detect format and load
#' turtle_data <- load_abm_data("processed_turtles.csv")
#' 
#' # Load compressed CSV
#' turtle_data <- load_abm_data("processed_turtles.csv.gz")
#' 
#' # Load RDS file
#' turtle_data <- load_abm_data("processed_turtles.rds")
#' 
#' # Explicitly specify format
#' turtle_data <- load_abm_data("data_file", format = "rds")
#' 
#' # Load different formats
#' csv_data <- load_abm_data("data.csv")
#' rds_data <- load_abm_data("data.rds")
#' parquet_data <- load_abm_data("data.parquet")
#' feather_data <- load_abm_data("data.feather")
#' 
#' # Load and immediately process
#' library(dplyr)
#' adults <- load_abm_data("turtles.rds") %>%
#'   filter(breed == "adults")
#' 
#' # Load multiple files
#' files <- c("run1_turtles.rds", "run2_turtles.rds", "run3_turtles.rds")
#' all_data <- lapply(files, load_abm_data) %>%
#'   bind_rows()
#' 
#' # Check loaded data
#' data <- load_abm_data("results.csv")
#' str(data)
#' summary(data)
#' 
#' # Load and compare file sizes
#' csv_size <- object.size(load_abm_data("data.csv"))
#' rds_size <- object.size(load_abm_data("data.rds"))
#' parquet_size <- object.size(load_abm_data("data.parquet"))
#' 
#' # Load archived/compressed data
#' old_data <- load_abm_data("archived_results.csv.gz")
#' }
#'
#' @importFrom data.table fread
#' @export
load_abm_data <- function(filename, format = NULL) {
  if (!file.exists(filename)) {
    stop("File '", filename, "' does not exist")
  }
  
  # Auto-detect format if not specified
  if (is.null(format)) {
    format <- detect_file_format(filename)
  }
  
  # Load based on format
  switch(format,
         csv = data.table::fread(filename),
         rds = readRDS(filename),
         txt = data.table::fread(filename, sep = "\t"),
         feather = {
           if (!requireNamespace("arrow", quietly = TRUE)) {
             stop("Package 'arrow' needed for feather format. Please install it.")
           }
           arrow::read_feather(filename)
         },
         parquet = {
           if (!requireNamespace("arrow", quietly = TRUE)) {
             stop("Package 'arrow' needed for parquet format. Please install it.")
           }
           arrow::read_parquet(filename)
         },
         stop("Unsupported format: ", format)
  )
}

################################################################################
# HELPER FUNCTIONS
################################################################################

#' Extract run number from filename
#'
#' Extracts numeric run identifier from filename and formats as "Run X".
#'
#' @param filename Character string filename.
#' @return Character string run number in format "Run X".
#' @keywords internal
extract_run_number <- function(filename) {
  paste("Run", sub("^(\\d+).*", "\\1", filename))
}

#' Detect file format from extension
#'
#' Determines file format based on extension, including compressed files.
#'
#' @param filename Character string filename.
#' @return Character string format identifier.
#' @keywords internal
detect_file_format <- function(filename) {
  if (grepl("\\.csv(\\.gz)?$", filename)) return("csv")
  if (grepl("\\.rds$", filename)) return("rds")
  if (grepl("\\.txt(\\.gz)?$", filename)) return("txt")
  if (grepl("\\.feather$", filename)) return("feather")
  if (grepl("\\.parquet$", filename)) return("parquet")
  stop("Could not determine file format from extension: ", filename)
}

#' Standardize column names
#'
#' Converts column names to consistent format.
#'
#' @param col_names Character vector of column names.
#' @return Character vector of standardized names.
#' @keywords internal
standardize_column_names <- function(col_names) {
  # Convert hyphens to dots for consistency
  col_names <- gsub("-", ".", col_names)
  return(col_names)
}

#' Expand chromosome columns from bracketed format
#'
#' Expands bracketed chromosome arrays into individual allele columns.
#'
#' @param dt Data.table with chromosome columns.
#' @return Data.table with expanded columns.
#' @keywords internal
expand_chromosome_columns <- function(dt) {
  chromosome_cols <- grep("chromosome", names(dt), value = TRUE)
  
  if (length(chromosome_cols) == 0) {
    return(dt)
  }
  
  # Expansion function for individual chromosome column
  expand_single_chromosome <- function(vec, base_name, suffixes) {
    cleaned <- gsub("\\[|\\]", "", vec)
    split_alleles <- strsplit(cleaned, "_")
    
    # Create matrix of alleles
    allele_matrix <- do.call(rbind, lapply(split_alleles, function(x) {
      x <- x[1:length(suffixes)]
      length(x) <- length(suffixes)
      return(x)
    }))
    
    # Add columns to data.table
    for (i in seq_along(suffixes)) {
      col_name <- paste0(base_name, ".", suffixes[i])
      dt[[col_name]] <- allele_matrix[, i]
    }
  }
  
  # Define expansion mappings
  chromosome_mappings <- list(
    "t.dad.chromosome.1" = letters[1:5],
    "t.dad.chromosome.2" = letters[6:10],
    "t.mom.chromosome.1" = letters[1:5],
    "t.mom.chromosome.2" = letters[6:10]
  )
  
  # Expand each chromosome column
  for (chr_col in chromosome_cols) {
    if (chr_col %in% names(chromosome_mappings)) {
      suffixes <- chromosome_mappings[[chr_col]]
      if (chr_col %in% names(dt)) {
        expand_single_chromosome(dt[[chr_col]], chr_col, suffixes)
        # Remove original column
        dt[[chr_col]] <- NULL
      }
    }
  }
  
  return(dt)
}

#' Standardize column types
#'
#' Converts columns to appropriate numeric types.
#'
#' @param dt Data.table to standardize.
#' @return Data.table with proper column types.
#' @keywords internal
standardize_column_types <- function(dt) {
  # Define numeric columns
  numeric_cols <- c(
    "ticks", "t.age", "t.age.metamorphosis", "t.birth.pond.index", "t.current.pond.index",
    "t.size.cm", "t.energy", "t.energy.consumption", "t.hatch.countdown", "t.max.size",
    "t.metamorphic.risk", "t.size.cm.metamorphosis", "t.stress", "ycor", "xcor"
  )
  
  # Convert existing numeric columns
  for (col in numeric_cols) {
    if (col %in% names(dt)) {
      suppressWarnings(data.table::set(dt, j = col, value = as.numeric(dt[[col]])))
    }
  }
  
  return(dt)
}

#' Apply time point filtering
#'
#' Filters data to specified time points.
#'
#' @param dt Data.table to filter.
#' @param time_points Vector of time points to keep.
#' @return Filtered data.table.
#' @keywords internal
apply_time_filter <- function(dt, time_points) {
  if (!is.null(time_points) && "ticks" %in% names(dt)) {
    dt <- dt[dt$ticks %in% time_points, ]
  }
  return(dt)
}

#' Apply sampling to data
#'
#' Samples data according to specified method.
#'
#' @param dt Data.table to sample.
#' @param sample_size Number of rows to sample.
#' @param sample_method Sampling method.
#' @param file_name File name for warnings.
#' @return Sampled data.table.
#' @keywords internal
apply_sampling <- function(dt, sample_size, sample_method, file_name) {
  total_rows <- nrow(dt)
  
  if (total_rows == 0) {
    warning("File '", basename(file_name), "' had zero rows after filtering")
    return(dt)
  }
  
  if (!is.null(sample_size) && sample_size < total_rows) {
    if (sample_method == "systematic") {
      keep_every <- max(1, floor(total_rows / sample_size))
      indices <- seq(from = 1, to = total_rows, by = keep_every)
      if (length(indices) > sample_size) {
        indices <- indices[1:sample_size]
      }
      dt <- dt[indices, ]
    } else if (sample_method == "random") {
      indices <- sample.int(total_rows, sample_size)
      dt <- dt[indices, ]
    } else if (sample_method == "head") {
      dt <- dt[1:min(sample_size, total_rows), ]
    } else if (sample_method == "tail") {
      start_row <- max(1, total_rows - sample_size + 1)
      dt <- dt[start_row:total_rows, ]
    } else {
      warning("Unknown sample method '", sample_method, "'; returning all data")
    }
  }
  
  return(dt)
}


################################################################################
# EXAMPLE DATA GENERATORS (aligned ticks; no external helpers)
################################################################################

#' Generate example MUTATION file with aligned ticks
#'
#' Creates a space-delimited mutations file for testing and demonstration purposes.
#' Generates realistic mutation events distributed across a time schedule.
#'
#' @param file Output file path. Default: "example_mutations.txt".
#' @param n_rows Total mutation events to generate (distributed across ticks). Default: 500.
#' @param max_ticks Maximum simulation time. Default: 5000.
#' @param record_every Recording interval for time points. Default: 1000.
#' @param include_zero If TRUE, include tick 0 as first time point. Default: FALSE.
#' @param seed Random number generator seed for reproducibility. Default: 123.
#'
#' @return Invisibly returns the data.table that was written to file.
#'
#' @details Output columns:
#' \itemize{
#'   \item \code{ticks}: Time point of mutation
#'   \item \code{turtle}: Individual identifier
#'   \item \code{old_allele}: Original allele identifier
#'   \item \code{old_additivity_value}: Original allelic effect
#'   \item \code{new_allele}: Mutated allele identifier
#'   \item \code{new_additivity_value}: New allelic effect
#' }
#'
#' @examples
#' \dontrun{
#' # Generate basic example file
#' generate_example_mutation_file()
#' 
#' # Custom parameters
#' generate_example_mutation_file(
#'   file = "test_mutations.txt",
#'   n_rows = 1000,
#'   max_ticks = 10000,
#'   record_every = 2000
#' )
#' 
#' # Include initial time point
#' generate_example_mutation_file(
#'   file = "mutations_with_t0.txt",
#'   include_zero = TRUE
#' )
#' 
#' # Generate and immediately read
#' generate_example_mutation_file("temp_mut.txt")
#' mut_data <- read_mutations(pattern = "temp_mut")
#' 
#' # Different random seed for variation
#' generate_example_mutation_file(seed = 456)
#' }
#'
#' @importFrom stats rnorm runif
#' @export
generate_example_mutation_file <- function(
    file = "example_mutations.txt",
    n_rows = 500,
    max_ticks = 5000,
    record_every = 1000,
    include_zero = FALSE,
    seed = 123
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
  if (!is.null(seed)) set.seed(seed)
  
  # Build tick schedule locally
  if (max_ticks <= 0L || record_every <= 0L) stop("max_ticks and record_every must be positive.")
  tick_schedule <- seq(from = record_every, to = max_ticks, by = record_every)
  if (include_zero) tick_schedule <- c(0L, tick_schedule)
  tick_schedule <- unique(as.integer(tick_schedule))
  
  # Allocate rows per tick approximately evenly
  k <- length(tick_schedule)
  base <- floor(n_rows / k); rem <- n_rows - base * k
  per_tick <- rep(base, k)
  if (rem > 0) per_tick[seq_len(rem)] <- per_tick[seq_len(rem)] + 1L
  rows_ticks <- unlist(mapply(function(tk, n) rep(tk, n), tick_schedule, per_tick, SIMPLIFY = FALSE))
  n_rows <- length(rows_ticks)
  
  dt <- data.table::data.table(
    ticks  = rows_ticks,
    turtle = paste0("turtle_", sample(10000:999999, n_rows, replace = TRUE))
  )
  
  loci <- letters[1:10]
  old_allele <- paste0(sample(loci, n_rows, TRUE), sample(c("1","2"), n_rows, TRUE))
  old_val <- rnorm(n_rows, mean = 0, sd = 6)
  
  dt[, `:=`(
    old_allele = old_allele,
    old_additivity_value = old_val,
    new_allele = paste0(sub("(.*)", "\\1", old_allele), ".", sub("turtle_", "", turtle)),
    new_additivity_value = old_val + rnorm(.N, mean = 0, sd = 5)
  )]
  
  dt <- dt[, .(ticks, turtle, old_allele, old_additivity_value, new_allele, new_additivity_value)]
  data.table::fwrite(dt, file, sep = " ", quote = FALSE)
  message("Wrote mutations to: ", normalizePath(file))
  invisible(dt)
}

#' Generate example TURTLE file (space-separated, bracketed chromosomes) with aligned ticks
#'
#' Creates an example turtle file matching the space-separated format with
#' bracketed chromosome data. Useful for testing, demonstrations, and vignettes.
#'
#' @param file Output file path. Default: "example_turtles.txt".
#' @param n_rows Total turtle records to generate (distributed across ticks). Default: 1000.
#' @param max_ticks Maximum simulation time. Default: 5000.
#' @param record_every Recording interval. Default: 1000.
#' @param include_zero Include tick 0. Default: FALSE.
#' @param seed Random seed. Default: 123.
#'
#' @return Invisibly returns a data.table of generated data.
#'
#' @details Generates realistic turtle data including:
#' \itemize{
#'   \item Multiple life stages (eggs, larvae, juveniles, adults)
#'   \item Genetic information (chromosomes with 10 alleles total)
#'   \item Spatial coordinates
#'   \item Energy and size attributes
#'   \item Parent IDs
#' }
#'
#' @examples
#' \dontrun{
#' # Generate basic example
#' generate_example_turtle_file()
#' 
#' # Larger dataset
#' generate_example_turtle_file(
#'   file = "large_turtles.txt",
#'   n_rows = 5000,
#'   max_ticks = 10000
#' )
#' 
#' # With initial time point
#' generate_example_turtle_file(
#'   file = "turtles_t0.txt",
#'   include_zero = TRUE
#' )
#' 
#' # Generate and read
#' generate_example_turtle_file("demo.txt")
#' demo_data <- read_turtles(pattern = "demo")
#' 
#' # Check different breeds
#' table(demo_data$breed)
#' }
#'
#' @importFrom stats rnorm runif
#' @export
generate_example_turtle_file <- function(
    file = "example_turtles.txt",
    n_rows = 1000,
    max_ticks = 5000,
    record_every = 1000,
    include_zero = FALSE,
    seed = 123
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
  if (!is.null(seed)) set.seed(seed)
  
  # Tick schedule
  if (max_ticks <= 0L || record_every <= 0L) stop("max_ticks and record_every must be positive.")
  tick_schedule <- seq(from = record_every, to = max_ticks, by = record_every)
  if (include_zero) tick_schedule <- c(0L, tick_schedule)
  tick_schedule <- unique(as.integer(tick_schedule))
  
  # Distribute rows over the schedule
  k <- length(tick_schedule)
  base <- floor(n_rows / k); rem <- n_rows - base * k
  per_tick <- rep(base, k)
  if (rem > 0) per_tick[seq_len(rem)] <- per_tick[seq_len(rem)] + 1L
  ticks_vec <- unlist(mapply(function(tk, n) rep(tk, n), tick_schedule, per_tick, SIMPLIFY = FALSE))
  
  breeds <- c("eggs", "larvae", "juveniles", "adults")
  sexes  <- c("male", "female")
  
  dt <- data.table::data.table(
    ticks = ticks_vec,
    who   = paste0("turtle_", sample(10000:999999, length(ticks_vec), TRUE)),
    breed = sample(breeds, length(ticks_vec), TRUE, prob = c(0.25, 0.15, 0.30, 0.30)),
    t.age = sample(0:220, length(ticks_vec), TRUE)
  )
  
  dt[, `:=`(
    `t-age-metamorphosis` = ifelse(breed %in% c("larvae", "juveniles", "adults"),
                                   sample(10:20, .N, TRUE), NA_integer_),
    `t-birth-pond-index`   = sample(0:2, .N, TRUE),
    `t-current-pond-index` = sample(0:2, .N, TRUE)
  )]
  
  # Inline haplotype strings like "[a1 b2 c2 d1 e1]"
  mk_hap <- function(prefix_letters) {
    paste0("[",
           paste0(paste0(prefix_letters, sample(c("1","2"), length(prefix_letters), TRUE)),
                  collapse = " "),
           "]")
  }
  dad_chr1 <- vapply(seq_len(nrow(dt)), function(i) mk_hap(letters[1:5]), character(1))
  dad_chr2 <- vapply(seq_len(nrow(dt)), function(i) mk_hap(letters[6:10]), character(1))
  mom_chr1 <- vapply(seq_len(nrow(dt)), function(i) mk_hap(letters[1:5]), character(1))
  mom_chr2 <- vapply(seq_len(nrow(dt)), function(i) mk_hap(letters[6:10]), character(1))
  
  dt[, `:=`(
    `t-dad-chromosome-1` = dad_chr1,
    `t-dad-chromosome-2` = dad_chr2,
    `t-dad-id` = paste0("turtle_", sample(10000:999999, .N, TRUE)),
    `t-energy` = runif(.N, 80, 100),
    `t-energy-consumption` = runif(.N, 5, 35),
    `t-hatch-countdown` = ifelse(breed == "eggs", round(runif(.N, 2, 10), 3), NA_real_),
    `t-max-size` = runif(.N, 60, 100),
    `t-metamorphic-risk` = runif(.N, 0, 100),
    `t-mom-chromosome-1` = mom_chr1,
    `t-mom-chromosome-2` = mom_chr2,
    `t-mom-id` = paste0("turtle_", sample(10000:999999, .N, TRUE)),
    `t-sex` = sample(sexes, .N, TRUE),
    `t-size-cm` = ifelse(breed == "eggs", runif(.N, 0, 1), runif(.N, 20, 35)),
    `t-size-cm-metamorphosis` = ifelse(breed %in% c("larvae", "juveniles", "adults"),
                                       runif(.N, 5, 8), NA_real_),
    `t-stress` = runif(.N, 4, 10),
    ycor = runif(.N, 10, 47),
    xcor = runif(.N, 100, 190)
  )]
  
  # Write header line (your reader skips one line), then space-delimited rows
  con <- file(file, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines("ticks who breed t-age t-age-metamorphosis t-birth-pond-index t-current-pond-index t-dad-chromosome-1 t-dad-chromosome-2 t-dad-id t-energy t-energy-consumption t-hatch-countdown t-max-size t-metamorphic-risk t-mom-chromosome-1 t-mom-chromosome-2 t-mom-id t-sex t-size-cm t-size-cm-metamorphosis t-stress ycor xcor", con)
  
  explode <- function(x) {
    y <- gsub("^\\[|\\]$", "", x)
    strsplit(y, " +")[[1]]
  }
  
  for (i in seq_len(nrow(dt))) {
    r <- dt[i]
    tokens <- c(
      r[["ticks"]],
      r[["who"]],
      r[["breed"]],
      r[["t.age"]],
      r[["t-age-metamorphosis"]],
      r[["t-birth-pond-index"]],
      r[["t-current-pond-index"]],
      explode(r[["t-dad-chromosome-1"]]),
      explode(r[["t-dad-chromosome-2"]]),
      r[["t-dad-id"]],
      sprintf("%.10g", r[["t-energy"]]),
      sprintf("%.10g", r[["t-energy-consumption"]]),
      ifelse(is.na(r[["t-hatch-countdown"]]), "NA", sprintf("%.10g", r[["t-hatch-countdown"]])),
      sprintf("%.10g", r[["t-max-size"]]),
      sprintf("%.10g", r[["t-metamorphic-risk"]]),
      explode(r[["t-mom-chromosome-1"]]),
      explode(r[["t-mom-chromosome-2"]]),
      r[["t-mom-id"]],
      r[["t-sex"]],
      sprintf("%.10g", r[["t-size-cm"]]),
      ifelse(is.na(r[["t-size-cm-metamorphosis"]]), "NA", sprintf("%.10g", r[["t-size-cm-metamorphosis"]])),
      sprintf("%.10g", r[["t-stress"]]),
      sprintf("%.10g", r[["ycor"]]),
      sprintf("%.10g", r[["xcor"]])
    )
    writeLines(paste(tokens, collapse = " "), con)
  }
  
  message("Wrote turtles to: ", normalizePath(file))
  invisible(dt)
}

#' Generate example ENVIRONMENT file with aligned ticks
#'
#' Creates an example environment file with spatial patch data across time points.
#'
#' @param file Output file path. Default: "example_environment.txt".
#' @param patches_per_tick Number of patches per time point. Default: 1000.
#' @param max_ticks Maximum simulation time. Default: 5000.
#' @param record_every Recording interval. Default: 1000.
#' @param include_zero Include tick 0. Default: FALSE.
#' @param seed Random seed. Default: 123.
#'
#' @return Invisibly returns a data.table of generated data.
#'
#' @details Generates spatial environment data with realistic patch properties:
#' \itemize{
#'   \item Energy levels
#'   \item Risk values
#'   \item Permeability
#'   \item Regeneration rates
#'   \item Patch type index
#' }
#'
#' @examples
#' \dontrun{
#' # Basic generation
#' generate_example_environment_file()
#' 
#' # Larger landscape
#' generate_example_environment_file(
#'   file = "large_env.txt",
#'   patches_per_tick = 5000
#' )
#' 
#' # More time points
#' generate_example_environment_file(
#'   max_ticks = 10000,
#'   record_every = 500
#' )
#' 
#' # Generate and read
#' generate_example_environment_file("test_env.txt")
#' env <- read_environments(pattern = "test_env")
#' 
#' # Visualize
#' library(ggplot2)
#' ggplot(env, aes(x = pxcor, y = pycor, fill = patch_energy)) +
#'   geom_raster() +
#'   facet_wrap(~ticks)
#' }
#'
#' @importFrom stats rnorm runif
#' @export
generate_example_environment_file <- function(
    file = "example_environment.txt",
    patches_per_tick = 1000,
    max_ticks = 5000,
    record_every = 1000,
    include_zero = FALSE,
    seed = 123
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
  if (!is.null(seed)) set.seed(seed)
  
  if (max_ticks <= 0L || record_every <= 0L) stop("max_ticks and record_every must be positive.")
  tick_schedule <- seq(from = record_every, to = max_ticks, by = record_every)
  if (include_zero) tick_schedule <- c(0L, tick_schedule)
  tick_schedule <- unique(as.integer(tick_schedule))
  
  k <- length(tick_schedule)
  total <- patches_per_tick * k
  
  dt <- data.table::data.table(
    ticks = rep(tick_schedule, each = patches_per_tick),
    pxcor = sample(0:199, total, TRUE),
    pycor = sample(0:99,  total, TRUE),
    patch_energy = runif(total, 130, 175),
    patch_risk   = pmax(0, rnorm(total, mean = 5.5, sd = 1.5)),
    permeability = 100L,
    patch_regeneration = runif(total, 14, 19),
    new  = sample(c("true", "false"), total, TRUE, prob = c(0.1, 0.9)),
    index = sample(c(-1L, 0L, 1L, 2L), total, TRUE, prob = c(0.8, 0.05, 0.1, 0.05)),
    form = sample(c("true", "false"), total, TRUE, prob = c(0.05, 0.95))
  )
  
  dt <- dt[, .(ticks, pxcor, pycor, patch_energy, patch_risk,
               permeability, patch_regeneration, new, index, form)]
  
  data.table::fwrite(dt, file, sep = " ", quote = FALSE)
  message("Wrote environment to: ", normalizePath(file))
  invisible(dt)
}

################################################################################
# CONVENIENCE WRAPPERS (single-run and multi-rep; aligned ticks)
################################################################################

#' Generate a single trio of example files with aligned ticks
#'
#' Convenience wrapper that generates matching turtle, mutation, and environment
#' files for a single simulation run. All files share consistent time points.
#'
#' @param dir Output directory. Created if doesn't exist. Default: current directory.
#' @param prefix File name prefix. Default: "vignette_run".
#' @param max_ticks Maximum simulation time. Default: 5000.
#' @param record_every Recording interval. Default: 1000.
#' @param include_zero Include tick 0. Default: FALSE.
#' @param n_turtles Total turtle records. Default: 1000.
#' @param n_mutations Total mutation records. Default: 500.
#' @param patches_per_tick Patches per time point. Default: 1000.
#' @param seed Random seed. Default: 123.
#'
#' @return Invisibly returns a list with three data.tables: turtles, mutations, environment.
#'
#' @examples
#' \dontrun{
#' # Generate standard example dataset
#' data <- generate_example_abm_triplet()
#' 
#' # Custom output location
#' data <- generate_example_abm_triplet(
#'   dir = "test_data/",
#'   prefix = "demo_sim"
#' )
#' 
#' # Longer simulation
#' data <- generate_example_abm_triplet(
#'   max_ticks = 10000,
#'   record_every = 2000,
#'   n_turtles = 2000
#' )
#' 
#' # Generate and immediately read
#' generate_example_abm_triplet(dir = "temp/", prefix = "test")
#' turtles <- read_turtles(path = "temp/", pattern = "*test*turtle")
#' mutations <- read_mutations(path = "temp/", pattern = "*test*mut")
#' environment <- read_environments(path = "temp/", pattern = "*test*env")
#' }
#'
#' @export
generate_example_abm_triplet <- function(
    dir = ".",
    prefix = "vignette_run",
    max_ticks = 5000,
    record_every = 1000,
    include_zero = FALSE,
    n_turtles = 1000,
    n_mutations = 500,
    patches_per_tick = 1000,
    seed = 123
) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  
  tfile <- file.path(dir, paste0(prefix, "_turtles.txt"))
  mfile <- file.path(dir, paste0(prefix, "_mutations.txt"))
  efile <- file.path(dir, paste0(prefix, "_environment.txt"))
  
  tdt <- generate_example_turtle_file(
    file = tfile,
    n_rows = n_turtles,
    max_ticks = max_ticks,
    record_every = record_every,
    include_zero = include_zero,
    seed = seed
  )
  
  mdt <- generate_example_mutation_file(
    file = mfile,
    n_rows = n_mutations,
    max_ticks = max_ticks,
    record_every = record_every,
    include_zero = include_zero,
    seed = seed + 1
  )
  
  edt <- generate_example_environment_file(
    file = efile,
    patches_per_tick = patches_per_tick,
    max_ticks = max_ticks,
    record_every = record_every,
    include_zero = include_zero,
    seed = seed + 2
  )
  
  invisible(list(turtles = tdt, mutations = mdt, environment = edt))
}

#' Generate N replicate example ABM triplets with aligned ticks
#'
#' Creates multiple sets of example files representing replicate simulation runs.
#' Files are named with numeric prefixes (01_, 02_, etc.) so run numbers are
#' correctly extracted.
#'
#' @param dir Output directory. Default: current directory.
#' @param prefix Base filename after numeric prefix. Default: "acert_demo".
#' @param n_reps Number of replicate sets to generate. Default: 3.
#' @param max_ticks Maximum simulation time. Default: 5000.
#' @param record_every Recording interval. Default: 1000.
#' @param include_zero Include tick 0. Default: FALSE.
#' @param n_turtles Turtle records per replicate. Default: 600.
#' @param mutations_per_tick Mutations per recorded tick per replicate. Default: 120.
#' @param patches_per_tick Patches per recorded tick per replicate. Default: 1200.
#' @param seed Base random seed (each rep gets offset). Default: 123.
#'
#' @return Invisibly returns a data.frame index of all generated file paths.
#'
#' @examples
#' \dontrun{
#' # Generate 3 replicate datasets
#' index <- generate_example_abm_replicates()
#' 
#' # More replicates
#' index <- generate_example_abm_replicates(
#'   dir = "multi_run_test/",
#'   n_reps = 10
#' )
#' 
#' # Custom parameters
#' index <- generate_example_abm_replicates(
#'   dir = "large_study/",
#'   prefix = "experiment_1",
#'   n_reps = 5,
#'   n_turtles = 2000,
#'   max_ticks = 10000
#' )
#' 
#' # Generate and read all
#' generate_example_abm_replicates(dir = "demo/", n_reps = 5)
#' all_turtles <- read_turtles(path = "demo/", pattern = "*turtle")
#' all_mutations <- read_mutations(path = "demo/", pattern = "*mut")
#' all_env <- read_environments(path = "demo/", pattern = "*env")
#' 
#' # Check run numbers
#' table(all_turtles$run_number)
#' 
#' # Use returned index
#' index <- generate_example_abm_replicates(n_reps = 3)
#' print(index)
#' 
#' # Read specific replicates
#' rep1_turtles <- read_turtles(path = dirname(index$turtles[1]))
#' }
#'
#' @export
generate_example_abm_replicates <- function(
    dir = ".",
    prefix = "acert_demo",
    n_reps = 3,
    max_ticks = 5000,
    record_every = 1000,
    include_zero = FALSE,
    n_turtles = 600,
    mutations_per_tick = 120,
    patches_per_tick = 1200,
    seed = 123
) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  out <- vector("list", n_reps)
  
  for (r in seq_len(n_reps)) {
    lead <- sprintf("%02d_%s", r, prefix)
    tfile <- file.path(dir, paste0(lead, "_turtles.txt"))
    mfile <- file.path(dir, paste0(lead, "_mutations.txt"))
    efile <- file.path(dir, paste0(lead, "_environment.txt"))
    
    generate_example_turtle_file(
      file = tfile,
      n_rows = n_turtles,
      max_ticks = max_ticks,
      record_every = record_every,
      include_zero = include_zero,
      seed = seed + r - 1
    )
    
    generate_example_mutation_file(
      file = mfile,
      n_rows = length(seq(if (include_zero) record_every else record_every,
                          max_ticks, by = record_every)) * mutations_per_tick,
      max_ticks = max_ticks,
      record_every = record_every,
      include_zero = include_zero,
      seed = seed + 1000 + r - 1
    )
    
    generate_example_environment_file(
      file = efile,
      patches_per_tick = patches_per_tick,
      max_ticks = max_ticks,
      record_every = record_every,
      include_zero = include_zero,
      seed = seed + 2000 + r - 1
    )
    
    out[[r]] <- data.frame(
      replicate = sprintf("Run %02d", r),
      turtles = normalizePath(tfile),
      mutations = normalizePath(mfile),
      environment = normalizePath(efile),
      stringsAsFactors = FALSE
    )
  }
  
  idx <- do.call(rbind, out)
  invisible(idx)
}