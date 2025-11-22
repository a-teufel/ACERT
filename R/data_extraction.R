################################################################################
# R/data_extraction.R
################################################################################
# 
# Functions for processing and transforming ABM analysis data
#
################################################################################

################################################################################
# FUNCTION: wrangle_differentiation_data
################################################################################
#' Extract and wrangle locus-specific differentiation data
#'
#' @param dm_list Named list where each element represents a replicate run
#' @param measure_type Differentiation measure to extract: "Gst", "Gst_H", or "Jost_D"
#'
#' @return Tibble with columns: File, Locus, and the specified measure
#'
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @importFrom dplyr rename
#' @importFrom rlang .data !!
#' @export
wrangle_differentiation_data <- function(dm_list, measure_type) {
  
  # Input validation
  valid_measures <- c("Gst", "Gst_H", "Jost_D")
  if (!(measure_type %in% valid_measures)) {
    stop("Invalid 'measure_type'. Must be one of: ", paste(valid_measures, collapse = ", "))
  }
  
  # Map measure types to list paths
  list_path_key <- switch(measure_type,
                          "Gst" = "Gst",
                          "Gst_H" = "Gst_Hed",
                          "Jost_D" = "Jost_D")
  
  est_key <- switch(measure_type,
                    "Gst" = "Gst_Est",
                    "Gst_H" = "Gst_H_Est", 
                    "Jost_D" = "Jost_D_Est")
  
  # Extract data from each run
  result <- purrr::map_dfr(names(dm_list), function(run_name) {
    run_data <- dm_list[[run_name]]
    
    # Navigate nested list structure
    per_locus_data <- NULL
    tryCatch({
      per_locus_data <- run_data[[list_path_key]][[est_key]]$per.locus
    }, error = function(e) {
      warning("Could not extract '", measure_type, "' data for run '", run_name, "'")
    })
    
    if (!is.null(per_locus_data) && length(per_locus_data) > 0) {
      tibble::tibble(
        File = run_name,
        Locus = names(per_locus_data),
        Value = as.numeric(per_locus_data)
      )
    } else {
      tibble::tibble(File = character(), Locus = character(), Value = numeric())
    }
  })
  
  # Rename Value column to measure type
  result %>% dplyr::rename(!!measure_type := "Value")
}

################################################################################
# FUNCTION: wrangle_heterozygosity_data
################################################################################
#' Wrangle heterozygosity data for plotting
#'
#' @param gd_list Named list where each element represents a replicate run containing
#'   Obs_Het and Exp_Het matrices with loci as rows and populations as columns
#'
#' @return Tibble with columns: File, Metric_Population, Mean_Value, Error_Value
#'
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr mutate bind_rows group_by summarise select
#' @importFrom tidyr pivot_longer
#' @importFrom stats na.omit sd
#' @importFrom rlang .data
#' @export
wrangle_heterozygosity_data <- function(gd_list) {
  
  # Helper function for Standard Error of the Mean
  calculate_sem <- function(x) {
    x <- stats::na.omit(x)
    if (length(x) < 2) return(NA_real_)
    stats::sd(x) / sqrt(length(x))
  }
  
  # Process each run
  all_runs_data <- purrr::map_dfr(names(gd_list), function(run_name) {
    run_data <- gd_list[[run_name]]
    
    obs_het_matrix <- run_data$Obs_Het
    exp_het_matrix <- run_data$Exp_Het
    
    # Validate matrices exist and have data
    if (is.null(obs_het_matrix) || is.null(exp_het_matrix) ||
        nrow(obs_het_matrix) == 0 || ncol(obs_het_matrix) == 0 ||
        nrow(exp_het_matrix) == 0 || ncol(exp_het_matrix) == 0) {
      warning("Missing or empty heterozygosity data for run '", run_name, "'")
      return(tibble::tibble())
    }
    
    # Convert observed heterozygosity to long format
    obs_het_long <- tibble::as_tibble(obs_het_matrix, rownames = "Locus") %>%
      tidyr::pivot_longer(
        cols = -.data$Locus,
        names_to = "Population", 
        values_to = "Value"
      ) %>%
      dplyr::mutate(Het_Type = "Ho")
    
    # Convert expected heterozygosity to long format
    exp_het_long <- tibble::as_tibble(exp_het_matrix, rownames = "Locus") %>%
      tidyr::pivot_longer(
        cols = -.data$Locus,
        names_to = "Population",
        values_to = "Value"
      ) %>%
      dplyr::mutate(Het_Type = "He")
    
    # Combine and summarize
    combined_het_long <- dplyr::bind_rows(obs_het_long, exp_het_long)
    
    summary_df <- combined_het_long %>%
      dplyr::group_by(.data$Population, .data$Het_Type) %>%
      dplyr::summarise(
        Mean_Value = mean(.data$Value, na.rm = TRUE),
        Error_Value = 2 * calculate_sem(.data$Value),
        .groups = 'drop'
      ) %>%
      dplyr::mutate(File = run_name)
    
    return(summary_df)
  })
  
  # Create ordered factor levels for plotting
  all_runs_data <- dplyr::mutate(all_runs_data, Population = as.character(.data$Population))
  unique_populations <- sort(unique(all_runs_data$Population))
  
  ordered_metric_levels <- character()
  for (pop in unique_populations) {
    ordered_metric_levels <- c(ordered_metric_levels, paste0("Ho", pop), paste0("He", pop))
  }
  
  # Finalize data structure
  all_runs_data %>%
    dplyr::mutate(
      Metric_Population = paste0(.data$Het_Type, .data$Population),
      Metric_Population = factor(.data$Metric_Population, levels = ordered_metric_levels)
    ) %>%
    dplyr::select(.data$File, .data$Metric_Population, .data$Mean_Value, .data$Error_Value)
}