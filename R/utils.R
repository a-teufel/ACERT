################################################################################
# R/utils.R
################################################################################
# 
# Utility functions for ABM analysis that don't fit in other modules
#
################################################################################

################################################################################
# FUNCTION: check_chromosome_structure
################################################################################
#' Check chromosome data structure in turtle data
#'
#' Analyzes the structure of chromosome data to understand how genetic
#' information is organized in the dataset. Useful for debugging and
#' ensuring data integrity.
#'
#' @param turtle_data A data.frame containing turtle data to analyze.
#'
#' @return A list containing information about chromosome column structure:
#'   \describe{
#'     \item{original}{Character vector of original chromosome column names}
#'     \item{expanded}{Character vector of expanded chromosome column names}
#'     \item{loci}{Character vector of detected locus names}
#'   }
#'
#' @details Examines turtle data to identify chromosome columns and determines
#'   whether they are in original bracketed format or expanded individual
#'   allele format. Provides summary information about detected loci and
#'   data completeness.
#'
#' @examples
#' \dontrun{
#' # Check chromosome structure
#' chr_info <- check_chromosome_structure(turtle_data)
#' 
#' # View detected loci
#' chr_info$loci
#' 
#' # Check if data is expanded
#' length(chr_info$expanded) > 0
#' }
#'
#' @importFrom utils head
#' @export
check_chromosome_structure <- function(turtle_data) {
  if (!is.data.frame(turtle_data)) {
    stop("turtle_data must be a data.frame")
  }
  
  # Find all chromosome columns
  chr_cols <- grep("chromosome", names(turtle_data), value = TRUE)
  
  # Separate original vs expanded columns
  original_cols <- chr_cols[!grepl("\\.[a-z]$", chr_cols)]
  expanded_cols <- chr_cols[grepl("\\.[a-z]$", chr_cols)]
  
  # Get unique loci
  loci <- unique(sub(".*\\.", "", expanded_cols))
  
  cat("Chromosome Data Structure:\n")
  cat("-------------------------\n")
  cat("Original chromosome columns:", length(original_cols), "\n")
  if (length(original_cols) > 0) {
    cat("  ", utils::head(original_cols, 4), "...\n")
    # Show sample values
    for (col in utils::head(original_cols, 2)) {
      if (col %in% names(turtle_data)) {
        sample_val <- turtle_data[[col]][!is.na(turtle_data[[col]])][1]
        if (!is.na(sample_val)) {
          cat("  Sample", col, ":", sample_val, "\n")
        }
      }
    }
  }
  
  cat("\nExpanded chromosome columns:", length(expanded_cols), "\n")
  cat("Number of loci detected:", length(loci), "\n")
  cat("Loci names:", paste(sort(loci), collapse = ", "), "\n")
  
  if (length(expanded_cols) > 0) {
    cat("\nExample expanded columns:\n")
    for (col in utils::head(expanded_cols, 6)) {
      sample_val <- turtle_data[[col]][!is.na(turtle_data[[col]])][1]
      if (!is.na(sample_val)) {
        cat("  ", col, "=", sample_val, "\n")
      }
    }
    
    # Check for completeness
    expected_cols <- outer(c("t.dad.chromosome.1", "t.dad.chromosome.2", 
                             "t.mom.chromosome.1", "t.mom.chromosome.2"), 
                           loci, paste, sep = ".")
    missing_cols <- setdiff(as.vector(expected_cols), expanded_cols)
    
    if (length(missing_cols) > 0) {
      cat("\nWarning: Missing expected columns:\n")
      cat("  ", utils::head(missing_cols, 6), "...\n")
    } else {
      cat("\nAll expected chromosome columns are present.\n")
    }
  }
  
  invisible(list(
    original = original_cols,
    expanded = expanded_cols,
    loci = loci
  ))
}

################################################################################
# FUNCTION: extract_genotypes
################################################################################
#' Extract genotypes from expanded chromosome columns for a specific individual
#'
#' Retrieves and formats genetic information for a single individual from
#' turtle data with expanded chromosome columns. Useful for detailed
#' genetic analysis of specific individuals.
#'
#' @param turtle_data A data.frame containing turtle data with expanded chromosome
#'   columns (e.g., t.dad.chromosome.1.a, t.mom.chromosome.2.b, etc.).
#' @param individual_id Identifier of the individual (value from 'who' column).
#' @param loci_names Character vector of locus names to extract. Default: letters[1:10].
#'
#' @return Named character vector of genotypes in "allele1/allele2" format,
#'   where names are locus names.
#'
#' @details Combines alleles from maternal and paternal chromosomes to create
#'   diploid genotypes. Extracts numeric allele values from allele names
#'   (e.g., "a1" becomes "1"). If fewer than 2 alleles are found, creates
#'   homozygous genotypes or returns NA.
#'
#' @examples
#' \dontrun{
#' # Extract genotypes for individual 42
#' genotypes <- extract_genotypes(turtle_data, individual_id = 42)
#' 
#' # View genotypes for first 5 loci
#' genotypes[1:5]
#' 
#' # Extract specific loci only
#' subset_genotypes <- extract_genotypes(
#'   turtle_data, 
#'   individual_id = 42,
#'   loci_names = c("a", "b", "c")
#' )
#' }
#'
#' @export
extract_genotypes <- function(turtle_data, individual_id, loci_names = letters[1:10]) {
  # Input validation
  if (!is.data.frame(turtle_data)) {
    stop("turtle_data must be a data.frame")
  }
  
  if (!"who" %in% names(turtle_data)) {
    stop("turtle_data must contain a 'who' column")
  }
  
  # Find the row for this individual
  ind_row <- turtle_data[turtle_data$who == individual_id, ]
  
  if (nrow(ind_row) == 0) {
    stop("Individual ", individual_id, " not found")
  }
  
  if (nrow(ind_row) > 1) {
    warning("Multiple records for individual ", individual_id, " - using first")
    ind_row <- ind_row[1, ]
  }
  
  # Extract genotypes for each locus
  genotypes <- character(length(loci_names))
  names(genotypes) <- loci_names
  
  for (i in seq_along(loci_names)) {
    locus <- loci_names[i]
    
    # Get alleles from all four chromosomes
    dad1_col <- paste0("t.dad.chromosome.1.", locus)
    dad2_col <- paste0("t.dad.chromosome.2.", locus)
    mom1_col <- paste0("t.mom.chromosome.1.", locus)
    mom2_col <- paste0("t.mom.chromosome.2.", locus)
    
    alleles <- c()
    for (col in c(dad1_col, dad2_col, mom1_col, mom2_col)) {
      if (col %in% names(ind_row) && !is.na(ind_row[[col]])) {
        # Extract numeric part from allele (e.g., "a1" -> "1")
        allele <- ind_row[[col]]
        if (nchar(allele) > 1) {
          allele_value <- substr(allele, 2, nchar(allele))
          alleles <- c(alleles, allele_value)
        }
      }
    }
    
    # Create genotype from first two alleles
    if (length(alleles) >= 2) {
      genotypes[i] <- paste0(alleles[1], "/", alleles[2])
    } else if (length(alleles) == 1) {
      genotypes[i] <- paste0(alleles[1], "/", alleles[1])  # Homozygote
    } else {
      genotypes[i] <- NA_character_
    }
  }
  
  return(genotypes)
}

################################################################################
# FUNCTION: create_time_animation
################################################################################
#' Create animated plots over time (requires gganimate)
#'
#' Creates animated visualizations that show changes over time in ABM
#' simulation data. Requires the gganimate package for rendering.
#'
#' @param data A data.frame containing time-series data.
#' @param plot_function Function that creates a ggplot for a single time point.
#'   Should accept data as first argument.
#' @param time_var Name of the time variable column. Default: "ticks".
#' @param output_file Output filename for animation. Default: "animation.gif".
#' @param fps Frames per second for animation. Default: 10.
#' @param width Width in pixels. Default: 800.
#' @param height Height in pixels. Default: 600.
#'
#' @return Character string path to created animation file.
#'
#' @details Creates animations by applying a plotting function to time-series
#'   data and rendering transitions between time points. Requires gganimate
#'   package to be installed.
#'
#' @examples
#' \dontrun{
#' # Define a plotting function
#' plot_population <- function(data) {
#'   ggplot(data, aes(x = xcor, y = ycor, color = breed)) +
#'     geom_point() +
#'     theme_minimal()
#' }
#' 
#' # Create animation
#' animation_path <- create_time_animation(
#'   turtle_data,
#'   plot_population,
#'   output_file = "population_dynamics.gif"
#' )
#' }
#'
#' @importFrom rlang sym
#' @export
create_time_animation <- function(data, 
                                  plot_function, 
                                  time_var = "ticks",
                                  output_file = "animation.gif", 
                                  fps = 10, 
                                  width = 800, 
                                  height = 600) {
  
  # Check for gganimate package
  if (!requireNamespace("gganimate", quietly = TRUE)) {
    stop("Package 'gganimate' needed for animations. Please install it with: install.packages('gganimate')")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for animations. Please install it.")
  }
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }
  
  if (!time_var %in% names(data)) {
    stop("time_var '", time_var, "' not found in data")
  }
  
  if (!is.function(plot_function)) {
    stop("plot_function must be a function")
  }
  
  # Create base plot
  base_plot <- plot_function(data)
  
  # Add animation
  animated_plot <- base_plot + 
    gganimate::transition_time(!!rlang::sym(time_var)) +
    ggplot2::labs(title = paste("Time: {frame_time}"))
  
  # Render animation
  gganimate::animate(
    animated_plot, 
    fps = fps, 
    width = width, 
    height = height,
    renderer = gganimate::gifski_renderer(output_file)
  )
  
  message("Animation saved to: ", output_file)
  return(output_file)
}