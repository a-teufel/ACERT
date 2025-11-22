################################################################################
# R/plotting.R
################################################################################
# 
# Visualization functions for ABM analysis results
#
################################################################################

################################################################################
# FUNCTION: plot_allele_frequencies
################################################################################
#' Plot allele frequencies by locus and population
#'
#' Creates a stacked bar chart visualization showing allele frequencies across
#' populations, faceted by locus. Useful for visualizing genetic diversity
#' patterns in agent-based model outputs.
#'
#' @param pop_freq_matrix Numeric matrix where rows are populations and columns are alleles.
#'   Column names must be in format 'LocusAllele' (e.g., 'A1', 'A2', 'B1', 'B2').
#' @param min_frequency Minimum frequency threshold for displaying alleles. 
#'   Alleles below this frequency will be filtered out. Default: 0.001.
#' @param population_labels Character vector of population labels. If NULL, 
#'   populations will be labeled 0, 1, 2, etc. Default: NULL.
#' @param locus_labels Named character vector for custom locus labels. 
#'   Names should be locus identifiers, values should be display labels. Default: NULL.
#' @param title Plot title. Default: "Allele Frequencies by Locus and Population".
#' @param show_legend Whether to display the legend showing allele colors. Default: FALSE.
#' @param bar_outline_color Color for bar outlines. Default: "black".
#' @param bar_outline_width Width of bar outline in mm. Default: 0.3.
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#'
#' @return A ggplot object showing stacked bar charts faceted by locus.
#'
#' @examples
#' \dontrun{
#' # Create sample allele frequency matrix
#' # Columns represent alleles at two loci (A and B)
#' freq_matrix <- matrix(
#'   c(0.4, 0.3, 0.2, 0.1,  # Population 0
#'     0.2, 0.5, 0.2, 0.1,  # Population 1
#'     0.1, 0.2, 0.6, 0.1), # Population 2
#'   nrow = 3, byrow = TRUE
#' )
#' colnames(freq_matrix) <- c("A1", "A2", "B1", "B2")
#' 
#' # Basic plot
#' plot_allele_frequencies(freq_matrix)
#' 
#' # With custom population labels and title
#' plot_allele_frequencies(
#'   freq_matrix,
#'   population_labels = c("Pond A", "Pond B", "Pond C"),
#'   title = "Genetic Structure Across Ponds"
#' )
#' 
#' # With custom locus labels and legend
#' plot_allele_frequencies(
#'   freq_matrix,
#'   locus_labels = c("A" = "Microsatellite A", "B" = "Microsatellite B"),
#'   show_legend = TRUE,
#'   x_axis_angle = 0
#' )
#' 
#' # Filter rare alleles
#' plot_allele_frequencies(
#'   freq_matrix,
#'   min_frequency = 0.15,
#'   bar_outline_width = 0.5
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col facet_wrap scale_fill_manual theme_minimal labs theme element_text
#' @importFrom dplyr mutate filter
#' @importFrom tidyr pivot_longer
#' @importFrom stats setNames
#' @importFrom rlang .data
#' @export
plot_allele_frequencies <- function(pop_freq_matrix, 
                                    min_frequency = 0.001,
                                    population_labels = NULL,
                                    locus_labels = NULL,
                                    title = "Allele Frequencies by Locus and Population",
                                    show_legend = FALSE,
                                    bar_outline_color = "black",
                                    bar_outline_width = 0.3,
                                    x_axis_angle = 45) {
  
  # Input validation
  if (!is.matrix(pop_freq_matrix) && !is.data.frame(pop_freq_matrix)) {
    stop("pop_freq_matrix must be a matrix or data frame")
  }
  if (is.null(colnames(pop_freq_matrix))) {
    stop("pop_freq_matrix must have column names in format 'LocusAllele'")
  }
  
  # Convert to data frame and handle missing values
  if (is.matrix(pop_freq_matrix)) {
    pop_freq_matrix <- as.data.frame(pop_freq_matrix)
  }
  pop_freq_matrix[is.na(pop_freq_matrix)] <- 0
  
  # Set population labels
  n_populations <- nrow(pop_freq_matrix)
  if (is.null(population_labels)) {
    population_labels <- as.character(seq_len(n_populations) - 1)
  }
  
  # Transform to long format
  df_long <- pop_freq_matrix %>%
    dplyr::mutate(Population = factor(population_labels, levels = population_labels)) %>%
    tidyr::pivot_longer(-"Population", names_to = "Allele", values_to = "Frequency") %>%
    dplyr::mutate(Locus = substr(.data$Allele, 1, 1)) %>%
    dplyr::filter(.data$Frequency > min_frequency)
  
  # Create locus labels
  unique_loci <- unique(df_long$Locus)
  if (is.null(locus_labels)) {
    locus_labels <- stats::setNames(paste("Locus", unique_loci), unique_loci)
  }
  
  # Define colors
  base_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                   "#0072B2", "#D55E00", "#CC79A7", "#999999")
  n_alleles <- length(unique(df_long$Allele))
  colors <- rep(base_colors, length.out = n_alleles)
  
  # Create plot
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data$Population, y = .data$Frequency, fill = .data$Allele)) +
    ggplot2::geom_col(position = "stack", color = bar_outline_color, linewidth = bar_outline_width) +
    ggplot2::facet_wrap(~.data$Locus, labeller = ggplot2::labeller(Locus = locus_labels)) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Population", y = "Frequency") +
    ggplot2::theme(legend.position = if (show_legend) "right" else "none")
  
  # Add x-axis rotation
  if (x_axis_angle != 0) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1))
  }
  
  return(p)
}

################################################################################
# FUNCTION: plot_allele_frequency_spectra
################################################################################
#' Plot faceted allele frequency spectra
#'
#' Creates faceted bar plots showing the distribution of alleles across
#' frequency bins for multiple simulation runs. This is useful for comparing
#' the allele frequency spectrum across different parameter combinations
#' or replicate runs.
#'
#' @param paf_list Named list where each element is a numeric vector of allele 
#'   frequency bin counts. Each vector should be named with bin labels 
#'   (e.g., "0-0.1", "0.1-0.2", etc.). List names represent run identifiers.
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#'
#' @return A ggplot object showing faceted bar plots of allele frequency spectra.
#'
#' @examples
#' \dontrun{
#' # Create sample allele frequency spectrum data
#' bin_names <- c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
#'                "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
#' 
#' paf_list <- list(
#'   "Run 1" = setNames(c(45, 23, 15, 8, 5, 3, 2, 1, 1, 0), bin_names),
#'   "Run 2" = setNames(c(50, 20, 12, 10, 4, 2, 1, 1, 0, 0), bin_names),
#'   "Run 3" = setNames(c(42, 25, 14, 9, 6, 2, 1, 1, 0, 0), bin_names),
#'   "Run 4" = setNames(c(48, 22, 13, 8, 5, 3, 1, 0, 0, 0), bin_names)
#' )
#' 
#' # Basic plot
#' plot_allele_frequency_spectra(paf_list)
#' 
#' # Without x-axis rotation
#' plot_allele_frequency_spectra(paf_list, x_axis_angle = 0)
#' 
#' # With many runs (will automatically arrange in 5 columns)
#' paf_list_large <- lapply(1:15, function(i) {
#'   setNames(rpois(10, lambda = 20 - i), bin_names)
#' })
#' names(paf_list_large) <- paste("Run", 1:15)
#' plot_allele_frequency_spectra(paf_list_large)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col facet_wrap labs theme_minimal theme element_text element_blank element_rect unit
#' @importFrom dplyr tibble
#' @importFrom purrr map_dfr
#' @export
plot_allele_frequency_spectra <- function(paf_list, x_axis_angle = 45) {
  if (length(paf_list) == 0 || is.null(names(paf_list[[1]]))) {
    stop("Input 'paf_list' must be a named list with named numeric vectors inside.")
  }
  
  # Transform to tidy data frame
  df_afs <- purrr::map_dfr(paf_list, ~{
    dplyr::tibble(
      bin = names(.x),
      count = as.numeric(.x)
    )
  }, .id = "Run")
  
  # Ensure bins are ordered factors
  ordered_bins <- names(paf_list[[1]])
  df_afs$bin <- factor(df_afs$bin, levels = ordered_bins)
  
  # Create base plot
  p <- ggplot2::ggplot(df_afs, ggplot2::aes(x = .data$bin, y = .data$count)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::facet_wrap(~ .data$Run, scales = "free_y", ncol = 5) +
    ggplot2::labs(
      x = "Allele Frequency Bin",
      y = "Number of Alleles",
      title = "Allele Frequency Spectra by Run"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7),
      panel.spacing = ggplot2::unit(0.1, "lines"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
  
  # Add x-axis rotation
  if (x_axis_angle != 0) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1, size = 7))
  } else {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7))
  }
  
  return(p)
}

################################################################################
# FUNCTION: plot_breed_population
################################################################################
#' Plot population counts over time for a specific breed type
#'
#' Creates line plots showing how population sizes change over time for a
#' specific life stage (breed type) across multiple simulation runs. Each
#' run is shown in a separate facet panel.
#'
#' @param tdat Data frame containing turtle (agent) data with required columns: 
#'   \code{ticks} (time steps), \code{breed} (life stage), and \code{run_number} (simulation ID).
#' @param breed_type Breed category to plot. Must be one of: "eggs", "larvae", 
#'   "juveniles", or "adults".
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#'
#' @return A ggplot object showing population counts over time faceted by run.
#'
#' @examples
#' \dontrun{
#' # Create sample turtle data
#' tdat <- data.frame(
#'   ticks = rep(0:100, each = 50),
#'   breed = sample(c("eggs", "larvae", "juveniles", "adults"), 
#'                  5050, replace = TRUE),
#'   run_number = rep(1:5, each = 1010),
#'   xcor = runif(5050, -50, 50),
#'   ycor = runif(5050, -50, 50)
#' )
#' 
#' # Plot adult population dynamics
#' plot_breed_population(tdat, "adults")
#' 
#' # Plot larvae population dynamics without axis rotation
#' plot_breed_population(tdat, "larvae", x_axis_angle = 0)
#' 
#' # Plot juveniles
#' plot_breed_population(tdat, "juveniles")
#' 
#' # Compare all life stages by creating multiple plots
#' library(patchwork)
#' p1 <- plot_breed_population(tdat, "eggs")
#' p2 <- plot_breed_population(tdat, "larvae")
#' p3 <- plot_breed_population(tdat, "juveniles")
#' p4 <- plot_breed_population(tdat, "adults")
#' (p1 + p2) / (p3 + p4)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_y_continuous labs facet_wrap theme_minimal theme element_text element_line
#' @importFrom dplyr filter group_by summarise n mutate
#' @importFrom rlang .data
#' @export
plot_breed_population <- function(tdat, breed_type, x_axis_angle = 45) {
  valid_breeds <- c("eggs", "larvae", "juveniles", "adults")
  if (!breed_type %in% valid_breeds) {
    stop("Invalid 'breed_type'. Must be one of: ", paste(valid_breeds, collapse = ", "))
  }
  
  breed_colors <- c(
    "eggs" = "#1f77b4", "larvae" = "#ff7f0e", 
    "juveniles" = "#2ca02c", "adults" = "#d62728"
  )
  
  plot_data <- tdat %>%
    dplyr::filter(.data$breed == breed_type) %>%
    dplyr::group_by(.data$run_number, .data$ticks, .data$breed) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(breed = factor(.data$breed, levels = breed_type))
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$ticks, y = .data$count, group = 1)) +
    ggplot2::geom_line(linewidth = 1, color = breed_colors[breed_type]) +
    ggplot2::geom_point(size = 2, alpha = 0.8, color = breed_colors[breed_type]) +
    ggplot2::scale_y_continuous(name = paste("Count (", tools::toTitleCase(breed_type), ")")) +
    ggplot2::labs(
      title = paste(tools::toTitleCase(breed_type), "Population Counts Over Time by Run"),
      x = "Ticks (Time Steps)"
    ) +
    ggplot2::facet_wrap(~ .data$run_number, scales = "free_y", ncol = 5) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  
  # Add x-axis rotation
  if (x_axis_angle != 0) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1))
  }
  
  return(p)
}

################################################################################
# FUNCTION: plot_diff_heatmap
################################################################################
#' Plot pairwise differentiation heatmap
#'
#' Creates heatmap visualizations of pairwise genetic differentiation between
#' populations. Supports multiple differentiation measures (Fst, Gst, etc.)
#' and displays results from multiple runs in faceted panels.
#'
#' @param data Named list where each element contains differentiation matrices 
#'   as dist objects. List names should follow the format "Run X Y" where X is
#'   the run number and Y is additional information.
#' @param measure Differentiation measure to plot. Must be one of: "Fst", "Gst", 
#'   "Gst_H", or "Jost_D". This must match a name in the differentiation objects.
#' @param viridis_option Viridis color palette option. Options include "A" (magma),
#'   "B" (inferno), "C" (plasma), "D" (viridis), "E" (cividis). Default: "D".
#' @param viridis_direction Direction of color scale: 1 for normal, -1 for reversed. 
#'   Default: 1.
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#'
#' @return A ggplot object showing differentiation heatmap faceted by run.
#'
#' @examples
#' \dontrun{
#' # Create sample differentiation data
#' # Simulate Fst values between 4 populations
#' create_diff_matrix <- function() {
#'   mat <- matrix(runif(16, 0, 0.3), nrow = 4)
#'   diag(mat) <- 0
#'   mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
#'   rownames(mat) <- colnames(mat) <- 0:3
#'   return(mat)
#' }
#' 
#' diff_data <- list(
#'   "Run 1 100" = list(
#'     Fst = as.dist(create_diff_matrix()),
#'     Gst = as.dist(create_diff_matrix())
#'   ),
#'   "Run 2 100" = list(
#'     Fst = as.dist(create_diff_matrix()),
#'     Gst = as.dist(create_diff_matrix())
#'   ),
#'   "Run 3 100" = list(
#'     Fst = as.dist(create_diff_matrix()),
#'     Gst = as.dist(create_diff_matrix())
#'   )
#' )
#' 
#' # Basic Fst heatmap
#' plot_diff_heatmap(diff_data, measure = "Fst")
#' 
#' # Gst with different color palette
#' plot_diff_heatmap(diff_data, measure = "Gst", viridis_option = "A")
#' 
#' # Reversed color scale
#' plot_diff_heatmap(diff_data, measure = "Fst", viridis_direction = -1)
#' 
#' # Without axis rotation
#' plot_diff_heatmap(diff_data, measure = "Fst", x_axis_angle = 0)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c facet_wrap labs theme_minimal theme element_text label_bquote
#' @importFrom dplyr bind_rows
#' @export
plot_diff_heatmap <- function(data, measure, viridis_option = "D", viridis_direction = 1, x_axis_angle = 45) {
  valid_measures <- c("Fst", "Gst", "Gst_H", "Jost_D")
  if (!measure %in% valid_measures) {
    stop("measure must be one of: ", paste(valid_measures, collapse = ", "))
  }
  
  # Extract run numbers and convert data
  run_names <- names(data)
  run_numbers <- gsub("Run\\s+(\\d+)\\s+\\d+", "\\1", run_names)
  
  plot_data_list <- vector("list", length(data))
  for (i in seq_along(data)) {
    if (!measure %in% names(data[[i]])) {
      stop("Measure '", measure, "' not found in ", run_names[i])
    }
    
    diff_matrix <- as.matrix(data[[i]][[measure]])
    pond_names <- rownames(diff_matrix)
    if (is.null(pond_names)) {
      pond_names <- 0:(nrow(diff_matrix)-1)
    }
    
    pond_pairs <- expand.grid(
      pond1 = pond_names, pond2 = pond_names, stringsAsFactors = FALSE
    )
    
    plot_data_list[[i]] <- data.frame(
      run_number = run_numbers[i],
      pond1 = pond_pairs$pond1,
      pond2 = pond_pairs$pond2,
      diff_value = as.vector(diff_matrix),
      stringsAsFactors = FALSE
    )
  }
  
  plot_data <- dplyr::bind_rows(plot_data_list)
  plot_data$run_number <- factor(plot_data$run_number, 
                                 levels = sort(as.numeric(unique(plot_data$run_number))))
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$pond2, y = .data$pond1, fill = .data$diff_value)) +
    ggplot2::geom_tile(color = "white", size = 0.1) +
    ggplot2::scale_fill_viridis_c(option = viridis_option, direction = viridis_direction, name = measure) +
    ggplot2::facet_wrap(~.data$run_number, labeller = ggplot2::label_bquote(Run~.(run_number))) +
    ggplot2::labs(title = paste("Pairwise", measure, "Differentiation"), x = "Pond", y = "Pond") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  # Add x-axis rotation
  if (x_axis_angle != 0) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1))
  }
  
  return(p)
}

################################################################################
# FUNCTION: plot_differentiation_strip
################################################################################
#' Create genetic differentiation strip plot
#'
#' Creates a strip plot showing differentiation values for multiple loci across
#' different files/replicates. Each locus is shown with a unique color and shape,
#' and mean values are overlaid as horizontal bars.
#'
#' @param data Data frame with columns: \code{File} (replicate identifier), 
#'   \code{Locus} (locus name), and a differentiation measure column.
#' @param y_var_col Column name containing differentiation values (as a string).
#' @param y_label Y-axis label for the plot.
#' @param title Plot title.
#' @param locus_names_for_mapping Character vector of locus names for color/shape mapping.
#'   Should match the unique values in the \code{Locus} column.
#' @param my_hardcoded_colors Character vector of color codes for each locus.
#'   Must be same length as \code{locus_names_for_mapping}.
#' @param my_custom_shapes Numeric vector of shape codes for each locus (use integers 0-25).
#'   Must be same length as \code{locus_names_for_mapping}.
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#'
#' @return A ggplot object showing strip plot with custom colors and shapes.
#'
#' @examples
#' \dontrun{
#' # Create sample differentiation data
#' diff_data <- data.frame(
#'   File = rep(paste("Rep", 1:5), each = 3),
#'   Locus = rep(c("LocusA", "LocusB", "LocusC"), 5),
#'   Fst = runif(15, 0.05, 0.25)
#' )
#' 
#' # Define loci, colors, and shapes
#' loci <- c("LocusA", "LocusB", "LocusC")
#' colors <- c("#E69F00", "#56B4E9", "#009E73")
#' shapes <- c(16, 17, 18)  # circle, triangle, diamond
#' 
#' # Basic strip plot
#' plot_differentiation_strip(
#'   data = diff_data,
#'   y_var_col = "Fst",
#'   y_label = "Fst",
#'   title = "Genetic Differentiation by Locus",
#'   locus_names_for_mapping = loci,
#'   my_hardcoded_colors = colors,
#'   my_custom_shapes = shapes
#' )
#' 
#' # Without x-axis rotation
#' plot_differentiation_strip(
#'   diff_data, "Fst", "Fst", "Differentiation",
#'   loci, colors, shapes, x_axis_angle = 0
#' )
#' 
#' # With more loci and different shapes
#' diff_data_large <- data.frame(
#'   File = rep(paste("Rep", 1:8), each = 5),
#'   Locus = rep(paste("Locus", LETTERS[1:5]), 8),
#'   Jost_D = runif(40, 0.1, 0.4)
#' )
#' 
#' plot_differentiation_strip(
#'   diff_data_large, "Jost_D", "Jost's D", 
#'   "Differentiation Across Replicates",
#'   paste("Locus", LETTERS[1:5]),
#'   rainbow(5),
#'   c(15, 16, 17, 18, 19)
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_jitter geom_point scale_color_manual scale_shape_manual labs theme_minimal theme guide_legend guides element_text
#' @importFrom dplyr group_by summarise
#' @importFrom rlang .data sym
#' @importFrom stats setNames
#' @export
plot_differentiation_strip <- function(data, y_var_col, y_label, title,
                                       locus_names_for_mapping,
                                       my_hardcoded_colors,
                                       my_custom_shapes,
                                       x_axis_angle = 45) {
  
  y_var_sym <- rlang::sym(y_var_col)
  
  mean_diff_data <- data %>%
    dplyr::group_by(.data$File) %>%
    dplyr::summarise(Mean_Diff = mean(!!y_var_sym, na.rm = TRUE), .groups = 'drop')
  
  custom_locus_colors <- stats::setNames(my_hardcoded_colors, locus_names_for_mapping)
  custom_locus_shapes <- stats::setNames(my_custom_shapes, locus_names_for_mapping)
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$File, y = !!y_var_sym)) +
    ggplot2::geom_jitter(ggplot2::aes(color = .data$Locus, shape = .data$Locus),
                         width = 0.2, height = 0, size = 3, alpha = 0.8) +
    ggplot2::geom_point(data = mean_diff_data, ggplot2::aes(y = .data$Mean_Diff),
                        color = "black", size = 6, shape = 95, stroke = 1.5) +
    ggplot2::scale_color_manual(values = custom_locus_colors) +
    ggplot2::scale_shape_manual(values = custom_locus_shapes) +
    ggplot2::labs(title = title, x = "Replicate", y = y_label) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1)),
                    shape = ggplot2::guide_legend())
  
  # Add x-axis rotation
  if (x_axis_angle != 0) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1))
  }
  
  return(p)
}

################################################################################
# FUNCTION: plot_distribution_by_facet
################################################################################
#' Plot distribution with violin, boxplot, and mean
#'
#' Creates faceted violin plots combined with boxplots and mean indicators.
#' Useful for visualizing distributions of metrics across different categories.
#'
#' @param data Data frame containing the data to plot.
#' @param y_var Column name for y-axis variable (as string).
#' @param facet_var Column name for faceting variable (as string).
#' @param title Plot title. Default: NULL (no title).
#' @param y_label Y-axis label. Default: NULL (uses column name).
#' @param fill_color Fill color for violin plots. Default: "lightblue".
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#' 
#' @return A ggplot object showing faceted violin plots with boxplots and means.
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' het_data <- data.frame(
#'   run_number = rep(1:5, each = 20),
#'   heterozygosity = rnorm(100, mean = 0.3, sd = 0.05),
#'   Fst = rnorm(100, mean = 0.15, sd = 0.03)
#' )
#' 
#' # Basic violin plot by run
#' plot_distribution_by_facet(
#'   data = het_data,
#'   y_var = "heterozygosity",
#'   facet_var = "run_number",
#'   title = "Heterozygosity Distribution by Run",
#'   y_label = "Observed Heterozygosity"
#' )
#' 
#' # Plot with custom color
#' plot_distribution_by_facet(
#'   het_data, "Fst", "run_number",
#'   title = "Fst Distribution",
#'   y_label = "Fst Value",
#'   fill_color = "coral"
#' )
#' 
#' # Without axis rotation
#' plot_distribution_by_facet(
#'   het_data, "heterozygosity", "run_number",
#'   y_label = "Het", x_axis_angle = 0
#' )
#' 
#' # With categorical faceting variable
#' het_data$scenario <- rep(c("Low Migration", "High Migration"), 50)
#' plot_distribution_by_facet(
#'   het_data, "heterozygosity", "scenario",
#'   title = "Heterozygosity by Migration Scenario",
#'   fill_color = "steelblue"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes_string geom_violin geom_boxplot stat_summary facet_wrap labs theme_minimal theme element_blank element_text
#' @export
plot_distribution_by_facet <- function(data, y_var, facet_var, title = NULL, y_label = NULL, fill_color = "lightblue",x_axis_angle = 45) {
  if (!all(c(y_var, facet_var) %in% colnames(data))) {
    stop("Specified columns not found in the data frame.")
  }
  
  p<- ggplot2::ggplot(data = data, ggplot2::aes_string(x = "''", y = y_var)) +
    ggplot2::geom_violin(fill = fill_color, alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.8) +
    ggplot2::stat_summary(fun = "mean", geom = "crossbar", width = 0.5, color = "black") +
    ggplot2::facet_wrap(as.formula(paste("~", facet_var))) +
    ggplot2::labs(title = title, y = y_label, x = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  # Add x-axis rotation
  if (x_axis_angle != 0) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1))
  }
  
  return(p)
}

################################################################################
# FUNCTION: plot_heterozygosity_bars
################################################################################
#' Plot mean heterozygosity with error bars
#'
#' Creates grouped bar plots showing mean heterozygosity values with error bars
#' across multiple replicates. Useful for comparing different heterozygosity
#' metrics (observed vs. expected) across populations.
#'
#' @param data Tibble with columns: \code{File} (replicate identifier), 
#'   \code{Metric_Population} (heterozygosity type and population), 
#'   \code{Mean_Value} (mean heterozygosity), and 
#'   \code{Error_Value} (standard error or standard deviation).
#' @param title Plot title.
#' @param y_label Y-axis label.
#' @param fill_colors Optional character vector of colors for bars. If NULL,
#'   uses ColorBrewer Set1 palette. Default: NULL.
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#'
#' @return A ggplot object showing grouped bar plot with error bars.
#'
#' @examples
#' \dontrun{
#' # Create sample heterozygosity data
#' het_summary <- data.frame(
#'   File = rep(paste("Rep", 1:5), each = 4),
#'   Metric_Population = rep(c("Ho_Pop0", "Ho_Pop1", "He_Pop0", "He_Pop1"), 5),
#'   Mean_Value = runif(20, 0.25, 0.45),
#'   Error_Value = runif(20, 0.01, 0.05)
#' )
#' 
#' # Basic plot with default colors
#' plot_heterozygosity_bars(
#'   data = het_summary,
#'   title = "Heterozygosity Across Populations",
#'   y_label = "Heterozygosity"
#' )
#' 
#' # With custom colors
#' custom_colors <- c("Ho_Pop0" = "#1f77b4", "Ho_Pop1" = "#ff7f0e",
#'                    "He_Pop0" = "#2ca02c", "He_Pop1" = "#d62728")
#' plot_heterozygosity_bars(
#'   het_summary, 
#'   "Observed vs Expected Heterozygosity",
#'   "Het Value",
#'   fill_colors = custom_colors
#' )
#' 
#' # Without x-axis rotation
#' plot_heterozygosity_bars(
#'   het_summary,
#'   "Heterozygosity Summary",
#'   "Heterozygosity",
#'   x_axis_angle = 0
#' )
#' 
#' # With many populations
#' het_large <- data.frame(
#'   File = rep(paste("Rep", 1:8), each = 6),
#'   Metric_Population = rep(paste("Ho_Pop", 0:5, sep = ""), 8),
#'   Mean_Value = runif(48, 0.2, 0.5),
#'   Error_Value = runif(48, 0.02, 0.08)
#' )
#' plot_heterozygosity_bars(het_large, "Multi-Population Het", "Value")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_errorbar labs theme_minimal theme position_dodge scale_fill_brewer scale_fill_manual element_text
#' @importFrom rlang .data
#' @export
plot_heterozygosity_bars <- function(data, title, y_label, fill_colors = NULL, x_axis_angle = 45) {
  pd <- ggplot2::position_dodge(width = 0.9)
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$File, y = .data$Mean_Value, fill = .data$Metric_Population)) +
    ggplot2::geom_bar(stat = "identity", position = pd) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$Mean_Value - .data$Error_Value,
                                        ymax = .data$Mean_Value + .data$Error_Value),
                           position = pd, width = 0.25) +
    ggplot2::labs(title = title, x = "Replicate", y = y_label, fill = "Heterozygosity Type and Population") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  if (!is.null(fill_colors)) {
    p <- p + ggplot2::scale_fill_manual(values = fill_colors)
  } else {
    p <- p + ggplot2::scale_fill_brewer(palette = "Set1")
  }
  
  # Add x-axis rotation
  if (x_axis_angle != 0) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1))
  }
  
  return(p)
}

################################################################################
# FUNCTION: plot_migration_heatmap
################################################################################
#' Plot effective migration heatmap
#'
#' Creates heatmap visualizations showing migration patterns between populations.
#' Can display either detailed data for a specific time point or summarized
#' data across time points.
#'
#' @param data Tibble containing migration data. Required columns depend on \code{data_type}.
#'   For "detailed": \code{run_number}, \code{ticks}, \code{origin_pond}, 
#'   \code{destination_pond}, \code{effective_migrants}.
#'   For "summarized": \code{run_number}, \code{origin_pond}, \code{destination_pond}, 
#'   \code{median_effective_migrants}.
#' @param data_type Type of data: "detailed" or "summarized".
#' @param tick Specific time point for detailed data. Required when \code{data_type = "detailed"}.
#'   Ignored for summarized data. Default: NULL.
#' @param viridis_option Viridis color palette option ("A", "B", "C", "D", "E"). Default: "D".
#' @param viridis_direction Direction of color scale: 1 for normal, -1 for reversed. Default: 1.
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#'
#' @return A ggplot object showing migration heatmap faceted by run_number.
#'
#' @examples
#' \dontrun{
#' # Create sample detailed migration data
#' detailed_mig <- expand.grid(
#'   run_number = 1:3,
#'   ticks = c(50, 100, 150),
#'   origin_pond = 0:3,
#'   destination_pond = 0:3
#' )
#' detailed_mig$effective_migrants <- rpois(nrow(detailed_mig), lambda = 5)
#' 
#' # Plot detailed data at specific time point
#' plot_migration_heatmap(
#'   data = detailed_mig,
#'   data_type = "detailed",
#'   tick = 100
#' )
#' 
#' # With different color palette
#' plot_migration_heatmap(
#'   detailed_mig, "detailed", tick = 100,
#'   viridis_option = "A", viridis_direction = -1
#' )
#' 
#' # Create sample summarized migration data
#' summarized_mig <- expand.grid(
#'   run_number = 1:4,
#'   origin_pond = 0:4,
#'   destination_pond = 0:4
#' )
#' summarized_mig$median_effective_migrants <- rpois(nrow(summarized_mig), 8)
#' 
#' # Plot summarized data
#' plot_migration_heatmap(
#'   data = summarized_mig,
#'   data_type = "summarized"
#' )
#' 
#' # Without axis rotation
#' plot_migration_heatmap(
#'   summarized_mig, "summarized",
#'   x_axis_angle = 0
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c facet_wrap labs theme_minimal theme element_text
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @export
plot_migration_heatmap <- function(data, data_type, tick = NULL, viridis_option = "D", viridis_direction = 1, x_axis_angle = 45) {
  if (!data_type %in% c("detailed", "summarized")) {
    stop("data_type must be either 'detailed' or 'summarized'")
  }
  
  if (data_type == "detailed") {
    required_cols <- c("run_number", "ticks", "origin_pond", "destination_pond", "effective_migrants")
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    if (is.null(tick)) {
      stop("tick parameter is required for detailed data")
    }
    
    plot_data <- data %>% dplyr::filter(.data$ticks == tick)
    plot_title <- paste("Effective Migration Patterns at Tick", tick)
    fill_var <- "effective_migrants"
  } else {
    required_cols <- c("run_number", "origin_pond", "destination_pond", "median_effective_migrants")
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    plot_data <- data
    plot_title <- "Effective Migration Patterns (Summarized)"
    fill_var <- "median_effective_migrants"
  }
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$destination_pond, y = .data$origin_pond, fill = .data[[fill_var]])) +
    ggplot2::geom_tile(color = "white", size = 0.1) +
    ggplot2::scale_fill_viridis_c(option = viridis_option, direction = viridis_direction, name = "Effective\nMigrants") +
    ggplot2::facet_wrap(~.data$run_number) +
    ggplot2::labs(title = plot_title, x = "Destination Pond", y = "Origin Pond") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  # Add x-axis rotation
  if (x_axis_angle != 0) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1))
  }
  
  return(p)
}

################################################################################
# FUNCTION: plot_migration_violins
################################################################################
#' Plot migration distributions with faceted violin plots
#'
#' Creates violin plots showing the distribution of migration values between
#' different pond combinations, faceted by destination pond.
#'
#' @param data Data frame containing migration data.
#' @param origin_col Column name for origin category (as string).
#' @param destination_col Column name for destination category (as string).
#' @param value_col Column name with numeric values to plot (as string).
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#'
#' @return A ggplot object showing violin plots faceted by destination.
#'
#' @examples
#' \dontrun{
#' # Create sample migration data
#' migration_data <- expand.grid(
#'   origin_pond = 0:3,
#'   destination_pond = 0:3,
#'   replicate = 1:20
#' )
#' migration_data$migrants <- rpois(nrow(migration_data), lambda = 10)
#' 
#' # Basic violin plot
#' plot_migration_violins(
#'   data = migration_data,
#'   origin_col = "origin_pond",
#'   destination_col = "destination_pond",
#'   value_col = "migrants"
#' )
#' 
#' # Without axis rotation
#' plot_migration_violins(
#'   migration_data,
#'   "origin_pond",
#'   "destination_pond",
#'   "migrants",
#'   x_axis_angle = 0
#' )
#' 
#' # With different variable names
#' migration_data2 <- data.frame(
#'   source = sample(c("A", "B", "C"), 100, replace = TRUE),
#'   target = sample(c("A", "B", "C"), 100, replace = TRUE),
#'   gene_flow = rnorm(100, mean = 15, sd = 5)
#' )
#' 
#' plot_migration_violins(
#'   migration_data2,
#'   origin_col = "source",
#'   destination_col = "target",
#'   value_col = "gene_flow"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_jitter facet_wrap vars labs theme_minimal theme element_text
#' @importFrom rlang sym
#' @export
plot_migration_violins <- function(data, origin_col, destination_col, value_col, x_axis_angle = 45) {
  origin_sym <- rlang::sym(origin_col)
  destination_sym <- rlang::sym(destination_col)
  value_sym <- rlang::sym(value_col)
  
  tryCatch({
    p <- ggplot2::ggplot(data, ggplot2::aes(x = factor(!!origin_sym), y = !!value_sym)) +
      ggplot2::geom_violin(ggplot2::aes(fill = factor(!!destination_sym)), trim = FALSE, scale = "width") +
      ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
      ggplot2::facet_wrap(ggplot2::vars(!!destination_sym), scales = "free_y") +
      ggplot2::labs(
        title = paste("Distribution of", value_col, "by Pond Combination"),
        subtitle = paste("Faceted by", destination_col),
        x = origin_col, y = value_col, fill = destination_col
      ) +
      ggplot2::theme_minimal()
    
    # Add x-axis rotation
    if (x_axis_angle != 0) {
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1))
    }
    
    return(p)
  }, error = function(e) {
    stop("Error during plotting. Check column names and data types. Original error: ", e)
  })
}

################################################################################
# FUNCTION: plot_netlogo_landscape
################################################################################
#' Plot NetLogo turtle positions on landscape
#'
#' Creates a spatial visualization showing agent (turtle) positions overlaid
#' on an environmental landscape with continuous variable coloring.
#'
#' @param tdat Data frame with turtle data. Required columns: \code{xcor}, \code{ycor}, 
#'   \code{breed}, \code{ticks}.
#' @param edat Data frame with environment data. Required columns: \code{pxcor}, \code{pycor}, 
#'   \code{index}, and the continuous variable specified in \code{continuous_var}.
#' @param continuous_var Continuous variable for landscape coloring. Must be one of:
#'   "patch_risk", "permeability", "patch_regeneration", or "patch_energy".
#' @param color_palette_option Color palette option for the landscape. Default: "viridis".
#' @param title Optional plot title. If NULL, generates automatic title. Default: NULL.
#' @param show_legend Whether to display legend. Default: TRUE.
#' @param turtle_alpha Transparency of turtle points (0-1). Default: 0.7.
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#'
#' @return A ggplot object showing landscape with turtle positions.
#'
#' @examples
#' \dontrun{
#' # Create sample turtle data
#' turtles <- data.frame(
#'   xcor = runif(200, -25, 25),
#'   ycor = runif(200, -25, 25),
#'   breed = sample(c("juveniles", "adults"), 200, replace = TRUE),
#'   ticks = rep(c(0, 50, 100), length.out = 200)
#' )
#' 
#' # Create sample environment data
#' env <- expand.grid(
#'   pxcor = -25:25,
#'   pycor = -25:25
#' )
#' env$index <- sample(c(-1, 0, 1, 2), nrow(env), replace = TRUE)
#' env$patch_risk <- runif(nrow(env), 0, 1)
#' env$permeability <- runif(nrow(env), 0.5, 1)
#' env$patch_regeneration <- runif(nrow(env), 0, 0.1)
#' env$patch_energy <- runif(nrow(env), 50, 100)
#' 
#' # Basic landscape plot with risk
#' plot_netlogo_landscape(
#'   tdat = turtles,
#'   edat = env,
#'   continuous_var = "patch_risk"
#' )
#' 
#' # Show energy landscape with custom title
#' plot_netlogo_landscape(
#'   turtles, env,
#'   continuous_var = "patch_energy",
#'   title = "Agent Distribution on Energy Landscape",
#'   turtle_alpha = 0.5
#' )
#' 
#' # Permeability landscape without legend
#' plot_netlogo_landscape(
#'   turtles, env,
#'   continuous_var = "permeability",
#'   color_palette_option = "plasma",
#'   show_legend = FALSE
#' )
#' 
#' # Multiple time points (will automatically facet)
#' plot_netlogo_landscape(
#'   turtles, env,
#'   continuous_var = "patch_regeneration",
#'   x_axis_angle = 0
#' )
#' }
#'
#' @importFrom ggplot2 ggplot geom_raster geom_tile geom_point scale_fill_viridis_c theme_minimal labs theme element_blank facet_wrap
#' @importFrom rlang .data
#' @export
plot_netlogo_landscape <- function(tdat, edat, continuous_var, 
                                   color_palette_option = "viridis", 
                                   title = NULL, show_legend = TRUE, turtle_alpha = 0.7, x_axis_angle = 45) {
  
  valid_vars <- c("patch_risk", "permeability", "patch_regeneration", "patch_energy")
  if (!(continuous_var %in% valid_vars)) {
    stop("Invalid continuous_var. Choose from: ", paste(valid_vars, collapse = ", "))
  }
  
  # Base plot with continuous variable
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = edat[edat$index < 0, ],
                         mapping = ggplot2::aes(x = .data$pxcor, y = .data$pycor, fill = .data[[continuous_var]])) +
    ggplot2::scale_fill_viridis_c(option = color_palette_option, name = continuous_var) +
    ggplot2::geom_tile(data = edat[edat$index >= 0, ],
                       mapping = ggplot2::aes(x = .data$pxcor, y = .data$pycor), fill = "lightblue") +
    ggplot2::geom_point(data = tdat,
                        mapping = ggplot2::aes(x = .data$xcor, y = .data$ycor, color = .data$breed),
                        size = 1.5, alpha = turtle_alpha) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "X-coordinate", y = "Y-coordinate", color = "Breed") +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank())
  
  # Add title
  if (!is.null(title)) {
    p <- p + ggplot2::labs(title = title)
  } else {
    p <- p + ggplot2::labs(title = paste0("Turtle Locations on a ", continuous_var, " Landscape"))
  }
  
  # Add faceting for multiple time points
  if (length(unique(tdat$ticks)) > 1) {
    p <- p + ggplot2::facet_wrap(~ .data$ticks)
  }
  
  # Control legend
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  if (x_axis_angle != 0) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1))
  }
  
  return(p)
}

################################################################################
# FUNCTION: plot_selected_mutations
################################################################################
#' Plot selected mutations over time
#'
#' Visualizes allele frequency trajectories over time for alleles showing
#' significant selection signals. Creates faceted panels showing both observed
#' frequencies and fitted trend lines.
#'
#' @param mutation_results Tibble output from \code{find_selected_mutations()}.
#'   Must contain columns: \code{locus}, \code{allele_name}, \code{population}, 
#'   \code{p_value}, \code{slope}.
#' @param genind_list List of genind objects used in \code{find_selected_mutations()}.
#'   Names should indicate time points.
#' @param target_loci Character vector of locus names used in \code{find_selected_mutations()}.
#' @param max_panels Maximum number of panels to plot. Default: 20.
#' @param rank_by How to rank results: "p_value" (most significant first) or 
#'   "slope" (steepest slopes first). Default: "p_value".
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#'
#' @return A ggplot object showing allele frequency trajectories over time.
#'
#' @examples
#' \dontrun{
#' # Assume you have genind objects across time points
#' library(adegenet)
#' 
#' # Create sample genind data (simplified example)
#' genind_t0 <- adegenet::genind(...)  # Your genind object at time 0
#' genind_t50 <- adegenet::genind(...) # Your genind object at time 50
#' genind_t100 <- adegenet::genind(...)# Your genind object at time 100
#' 
#' genind_list <- list(
#'   "Run 1 0" = genind_t0,
#'   "Run 1 50" = genind_t50,
#'   "Run 1 100" = genind_t100
#' )
#' 
#' # Find selected mutations (hypothetical function)
#' mutation_results <- find_selected_mutations(
#'   genind_list,
#'   target_loci = c("LocusA", "LocusB", "LocusC")
#' )
#' 
#' # Plot selected mutations ranked by p-value
#' plot_selected_mutations(
#'   mutation_results = mutation_results,
#'   genind_list = genind_list,
#'   target_loci = c("LocusA", "LocusB", "LocusC"),
#'   max_panels = 12
#' )
#' 
#' # Rank by slope instead
#' plot_selected_mutations(
#'   mutation_results, genind_list,
#'   target_loci = c("LocusA", "LocusB", "LocusC"),
#'   rank_by = "slope",
#'   max_panels = 15
#' )
#' 
#' # Without axis rotation
#' plot_selected_mutations(
#'   mutation_results, genind_list,
#'   target_loci = c("LocusA", "LocusB", "LocusC"),
#'   max_panels = 20,
#'   x_axis_angle = 0
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth facet_wrap labs theme_minimal theme element_text
#' @importFrom adegenet tab nInd seppop
#' @importFrom dplyr as_tibble bind_rows filter arrange slice_head inner_join mutate desc row_number
#' @importFrom tidyr separate
#' @importFrom rlang .data
#' @importFrom utils installed.packages
#' @export
plot_selected_mutations <- function(mutation_results, genind_list, target_loci,
                                    max_panels = 20, rank_by = "p_value", x_axis_angle = 45) {
  
  if (nrow(mutation_results) == 0) {
    stop("No results to plot. mutation_results is empty.")
  }
  if (!(rank_by %in% c("p_value", "slope"))) {
    stop("rank_by must be either 'p_value' or 'slope'.")
  }
  
  # Limit and rank results
  if (nrow(mutation_results) > max_panels) {
    if (rank_by == "p_value") {
      mutation_results <- mutation_results %>%
        dplyr::arrange(.data$p_value) %>%
        dplyr::slice_head(n = max_panels)
    } else {
      mutation_results <- mutation_results %>%
        dplyr::arrange(desc(.data$slope)) %>%
        dplyr::slice_head(n = max_panels)
    }
  }
  
  # Extract time points
  list_names <- names(genind_list)
  ticks <- as.numeric(gsub(".* (\\d+)$", "\\1", list_names))
  
  # Determine frequency scope
  n_populations <- length(unique(mutation_results$population))
  use_population_facets <- n_populations > 1 && !"global" %in% mutation_results$population
  
  # Recalculate frequencies for plotting
  frequency_data <- list()
  
  for (i in seq_along(genind_list)) {
    genind_obj <- genind_list[[i]]
    current_tick <- ticks[i]
    
    # Filter to target loci
    genind_obj <- genind_obj[loc = target_loci]
    
    if (use_population_facets) {
      # Population-specific frequencies
      pop_list <- adegenet::seppop(genind_obj)
      pop_names <- names(pop_list)
      
      pop_freq_list <- list()
      for (j in seq_along(pop_list)) {
        pop_obj <- pop_list[[j]]
        pop_name <- pop_names[j]
        
        if (adegenet::nInd(pop_obj) > 0) {
          pop_allele_counts <- adegenet::tab(pop_obj, freq = FALSE, NA.method = "mean")
          pop_allele_freqs <- colSums(pop_allele_counts) / (2 * adegenet::nInd(pop_obj))
          
          pop_freq_list[[j]] <- data.frame(
            allele = names(pop_allele_freqs),
            frequency = pop_allele_freqs,
            ticks = current_tick,
            population = pop_name,
            stringsAsFactors = FALSE
          )
        }
      }
      
      freq_df <- dplyr::bind_rows(pop_freq_list)
      
    } else {
      # Global frequencies
      allele_counts <- adegenet::tab(genind_obj, freq = FALSE, NA.method = "mean")
      allele_freqs <- colSums(allele_counts) / (2 * adegenet::nInd(genind_obj))
      
      freq_df <- data.frame(
        allele = names(allele_freqs),
        frequency = allele_freqs,
        ticks = current_tick,
        population = "global",
        stringsAsFactors = FALSE
      )
    }
    
    frequency_data[[i]] <- freq_df
  }
  
  # Combine and process frequency data
  all_freqs <- dplyr::bind_rows(frequency_data) %>%
    dplyr::as_tibble() %>%
    tidyr::separate(.data$allele, into = c("locus", "allele_name"), sep = "\\.", extra = "merge") %>%
    dplyr::filter(.data$frequency > 0)
  
  # Filter to only significant alleles
  mutation_results_ordered <- mutation_results %>%
    dplyr::mutate(original_order = row_number()) %>%
    dplyr::select(.data$locus, .data$allele_name, .data$population, .data$original_order)
  
  plot_data <- all_freqs %>%
    dplyr::inner_join(mutation_results_ordered, by = c("locus", "allele_name", "population")) %>%
    dplyr::mutate(allele_id = paste(.data$locus, .data$allele_name, sep = ".")) %>%
    dplyr::arrange(.data$original_order) %>%
    dplyr::mutate(allele_id = factor(.data$allele_id, levels = unique(.data$allele_id)))
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$ticks, y = .data$frequency)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.5) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 1) +
    ggplot2::facet_wrap(~ .data$allele_id, scales = "free_y") +
    ggplot2::labs(
      x = "Time (ticks)",
      y = "Allele Frequency", 
      title = "Selected Allele Trajectories Over Time"
    ) +
    ggplot2::theme_minimal()
  
  # Add population faceting if needed
  if (use_population_facets) {
    p <- p + ggplot2::facet_wrap(~ .data$allele_id + .data$population, scales = "free_y")
  }
  
  # Add x-axis rotation
  if (x_axis_angle != 0) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1))
  }
  
  return(p)
}

################################################################################
# FUNCTION: visualize_patches
################################################################################
#' Visualize patch properties
#'
#' Creates raster visualizations of environmental patch properties across space,
#' with options for temporal and spatial faceting.
#' 
#' @param env_data Environment data frame. Required columns: \code{pxcor}, \code{pycor}, 
#'   \code{ticks}, \code{run_number}, and the property specified in \code{property}.
#' @param property Property to visualize. Options include: "index", "patch_energy",
#'   "patch_risk", "permeability", "patch_regeneration".
#' @param time_points Specific time points to visualize. If NULL, shows all time points.
#'   Default: NULL.
#' @param facet_by Whether to facet by "time" or "run". Default: "time".
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#' 
#' @return A ggplot object showing spatial distribution of patch properties.
#' 
#' @examples
#' \dontrun{
#' # Create sample environment data
#' env_data <- expand.grid(
#'   pxcor = -25:25,
#'   pycor = -25:25,
#'   ticks = c(0, 50, 100),
#'   run_number = 1:3
#' )
#' env_data$index <- sample(c(-1, 0, 1, 2), nrow(env_data), replace = TRUE)
#' env_data$patch_energy <- runif(nrow(env_data), 50, 150)
#' env_data$patch_risk <- runif(nrow(env_data), 0, 1)
#' env_data$permeability <- runif(nrow(env_data), 0.3, 1)
#' env_data$patch_regeneration <- runif(nrow(env_data), 0, 0.2)
#' 
#' # Visualize patch index (habitat types) across time
#' visualize_patches(
#'   env_data = env_data,
#'   property = "index",
#'   facet_by = "time"
#' )
#' 
#' # Visualize energy at specific time points
#' visualize_patches(
#'   env_data,
#'   property = "patch_energy",
#'   time_points = c(0, 100)
#' )
#' 
#' # Compare risk across runs
#' visualize_patches(
#'   env_data,
#'   property = "patch_risk",
#'   time_points = 50,
#'   facet_by = "run"
#' )
#' 
#' # Permeability without axis rotation
#' visualize_patches(
#'   env_data,
#'   property = "permeability",
#'   facet_by = "time",
#'   x_axis_angle = 0
#' )
#' 
#' # Regeneration patterns
#' visualize_patches(
#'   env_data,
#'   property = "patch_regeneration",
#'   time_points = c(0, 50, 100)
#' )
#' }
#' 
#' @importFrom stats dist
#' @importFrom grDevices colorRampPalette
#' @export
visualize_patches <- function(env_data, property, time_points = NULL, facet_by = "time",x_axis_angle = 45) {
  # Filter time points if specified
  if (!is.null(time_points)) {
    env_data <- env_data %>%
      dplyr::filter(ticks %in% time_points)
  }
  
  # Create base plot
  p <- ggplot2::ggplot(env_data, aes(x = pxcor, y = pycor)) +
    ggplot2::geom_raster(aes(fill = .data[[property]]))
  
  # Add appropriate color scale
  if (property == "index") {
    # Convert numeric index to factor for discrete scale
    env_data[[property]] <- as.factor(env_data[[property]])
    
    # Recreate the plot with factor data
    p <- ggplot2::ggplot(env_data, aes(x = pxcor, y = pycor)) +
      ggplot2::geom_raster(aes(fill = .data[[property]]))
    
    # Create a custom color function for index values
    # Positive numbers (including 0) get shades of blue
    # Negative numbers get shades of green
    
    # Get the range of index values in the data (convert back to numeric for logic)
    index_values <- unique(as.numeric(as.character(env_data[[property]])))
    
    # Separate positive (including 0) and negative values
    pos_values <- index_values[index_values >= 0]
    neg_values <- index_values[index_values < 0]
    
    # Create color mappings
    colors <- c()
    
    # Handle negative values with green shades
    if (length(neg_values) > 0) {
      neg_values <- sort(neg_values)
      green_colors <- colorRampPalette(c("darkgreen", "lightgreen"))(length(neg_values))
      names(green_colors) <- as.character(neg_values)
      colors <- c(colors, green_colors)
    }
    
    # Handle positive values (including 0) with blue shades
    if (length(pos_values) > 0) {
      pos_values <- sort(pos_values)
      blue_colors <- colorRampPalette(c("lightblue", "darkblue"))(length(pos_values))
      names(blue_colors) <- as.character(pos_values)
      colors <- c(colors, blue_colors)
    }
    
    p <- p + ggplot2::scale_fill_manual(
      values = colors,
      name = "Patch type"
    )
  } else if (property == "patch_energy") {
    p <- p + ggplot2::scale_fill_gradient(low = "blue", high = "red", name = "Energy")
  } else if (property == "patch_risk") {
    p <- p + ggplot2::scale_fill_gradient(low = "white", high = "black", name = "Risk")
  } else if (property == "permeability") {
    p <- p + ggplot2::scale_fill_gradient(low = "darkgreen", high = "lightgreen", name = "Permeability")
  } else if (property == "patch_regeneration") {
    p <- p + ggplot2::scale_fill_gradient(low = "blue", high = "yellow", name = "Regeneration")
  }
  
  
  
  # Add faceting
  if (facet_by == "time") {
    p <- p + ggplot2::facet_wrap(~ ticks)
  } else if (facet_by == "run") {
    p <- p + ggplot2::facet_wrap(~ run_number)
  }
  
  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste(property, "across", facet_by),
      x = "X coordinate",
      y = "Y coordinate"
    )
  
  # Add x-axis rotation
  if (x_axis_angle != 0) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1))
  }
  
  
  return(p)
}



################################################################################
# FUNCTION: plot_turtle_property
################################################################################
#' Create turtle property plots
#'
#' Creates boxplots showing the distribution of turtle (agent) properties across
#' replicates, time steps, and optional grouping variables.
#'
#' @param data Analyzed turtle data frame. Must contain columns: \code{ticks},
#'   \code{mean_value}, and any grouping variable specified.
#' @param property Property name for labeling the y-axis.
#' @param grouping Optional grouping variable for coloring (column name as string).
#'   Default: NULL.
#' @param experiment_name Name of the experiment for plot title. Default: "" (no name).
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#' 
#' @return A ggplot object showing boxplots of the property.
#' 
#' @examples
#' \dontrun{
#' # Create sample turtle property data
#' turtle_data <- data.frame(
#'   ticks = rep(c(0, 25, 50, 75, 100), each = 20),
#'   mean_value = rnorm(100, mean = 75, sd = 10),
#'   scenario = rep(c("Low", "High"), 50),
#'   replicate = rep(1:20, 5)
#' )
#' 
#' # Basic plot without grouping
#' plot_turtle_property(
#'   data = turtle_data,
#'   property = "Energy"
#' )
#' 
#' # With grouping by scenario
#' plot_turtle_property(
#'   turtle_data,
#'   property = "Energy",
#'   grouping = "scenario",
#'   experiment_name = "Migration Study"
#' )
#' 
#' # Without axis rotation
#' plot_turtle_property(
#'   turtle_data,
#'   property = "Energy",
#'   grouping = "scenario",
#'   x_axis_angle = 0
#' )
#' 
#' # Multiple properties example
#' turtle_data$breeding_status <- sample(c("Breeding", "Non-breeding"), 
#'                                        100, replace = TRUE)
#' plot_turtle_property(
#'   turtle_data,
#'   property = "Body Condition",
#'   grouping = "breeding_status",
#'   experiment_name = "Breeding Dynamics"
#' )
#' }
#' 
#' @export
plot_turtle_property <- function(data, property, grouping = NULL, experiment_name = "",x_axis_angle = 45) {
  plot_title <- paste("Replicate Mean", property, "by", 
                      ifelse(is.null(grouping), "Time", paste(grouping, "and Time")),
                      ifelse(experiment_name == "", "", paste("in", experiment_name)))
  
  p <- ggplot2::ggplot(data, aes(x = factor(ticks), y = mean_value))
  
  if (!is.null(grouping)) {
    p <- p + ggplot2::geom_boxplot(
      aes(fill = factor(.data[[grouping]]),
          group = interaction(factor(ticks), factor(.data[[grouping]]))),
      position = ggplot2::position_dodge(width = 0.75),
      alpha = 0.5,
      outlier.shape = 16
    )
  } else {
    p <- p + ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = 16)
  }
  
  p <- p +
    ggplot2::labs(
      x = "Time steps (weeks)",
      y = paste("Mean", property),
      title = plot_title,
      fill = grouping
    ) +
    ggplot2::theme_minimal()
  
  # Add x-axis rotation
  if (x_axis_angle != 0) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1))
  }
  
  return(p)
}


################################################################################
# FUNCTION: plot_major_alleles
################################################################################
#' Plot Major Allele Frequencies Over Time
#'
#' Visualizes the output from \code{wrangle_major_alleles()}, showing the frequency
#' of the two major alleles over time with optional trend lines.
#'
#' @param plot_data A tibble generated by \code{wrangle_major_alleles()}.
#'   Required columns: \code{tick}, \code{locus}, \code{allele}, \code{frequency}.
#'   Optional: \code{population} (if analysis was by population).
#' @param trend_line Character string specifying type of trend line: "linear", 
#'   "quadratic", or "loess". Ignored if \code{type = "line"}. Default: "quadratic".
#' @param type Character string: "scatter" (points with optional trend line) or 
#'   "line" (connected points without trend). Default: "scatter".
#' @param panels_per_row Integer specifying number of panels per row in faceted plot.
#'   Default: 4.
#' @param analysis_mode Character string: "global" or "by_population".
#'   Used for labeling the plot title and facets. Default: "global".
#' @param x_axis_angle Angle for x-axis text rotation in degrees. 
#'   Set to 0 to disable rotation. Default: 45.
#'
#' @return A ggplot object showing major allele frequencies over time.
#'
#' @examples
#' \dontrun{
#' # Assume you have wrangle_major_alleles function and genind_list
#' library(adegenet)
#' 
#' # Example: Create sample data structure
#' plot_data <- data.frame(
#'   tick = rep(c(0, 25, 50, 75, 100), each = 8),
#'   locus = rep(rep(c("LocusA", "LocusB"), each = 4), 5),
#'   allele = rep(c("A1", "A2", "B1", "B2"), 10),
#'   frequency = runif(40, 0.2, 0.8)
#' )
#' 
#' # Basic scatter plot with quadratic trend
#' plot_major_alleles(
#'   plot_data,
#'   trend_line = "quadratic",
#'   type = "scatter",
#'   analysis_mode = "global"
#' )
#' 
#' # Line plot connecting points
#' plot_major_alleles(
#'   plot_data,
#'   type = "line",
#'   panels_per_row = 2
#' )
#' 
#' # Scatter with linear trend
#' plot_major_alleles(
#'   plot_data,
#'   trend_line = "linear",
#'   type = "scatter",
#'   panels_per_row = 3
#' )
#' 
#' # LOESS smoothing
#' plot_major_alleles(
#'   plot_data,
#'   trend_line = "loess",
#'   type = "scatter",
#'   x_axis_angle = 0
#' )
#' 
#' # By population analysis
#' plot_data_pop <- data.frame(
#'   tick = rep(c(0, 50, 100), each = 12),
#'   locus = rep(rep(c("LocusA", "LocusB"), each = 6), 3),
#'   allele = rep(c("A1", "A2", "B1", "B2"), 9),
#'   frequency = runif(36, 0.1, 0.9),
#'   population = rep(c("Pop0", "Pop1", "Pop2"), each = 4)
#' )
#' 
#' plot_major_alleles(
#'   plot_data_pop,
#'   trend_line = "quadratic",
#'   analysis_mode = "by_population"
#' )
#' }
#'
#' @import ggplot2
#' @importFrom tools toTitleCase
#' @export
plot_major_alleles <- function(plot_data,
                               trend_line = "quadratic",
                               type = "scatter",
                               panels_per_row = 4,
                               analysis_mode = "global",x_axis_angle = 45) {
  
  # Input validation
  if (!trend_line %in% c("linear", "quadratic", "loess")) {
    stop("trend_line must be 'linear', 'quadratic', or 'loess'")
  }
  
  if (!type %in% c("scatter", "line")) {
    stop("type must be 'scatter' or 'line'")
  }
  
  # Prepare facet labels
  if (analysis_mode == "global") {
    plot_data$facet_label <- paste("Locus:", plot_data$locus)
  } else {
    plot_data$facet_label <- paste("Locus:", plot_data$locus, "\nPop:", plot_data$population)
  }
  
  # Create base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = tick, y = frequency, color = allele)) +
    ggplot2::labs(
      title = paste0(
        "Major Allele Frequencies Over Time (",
        tools::toTitleCase(gsub("_", " ", analysis_mode)),
        ")"
      ),
      x = "Tick Number",
      y = "Allele Frequency",
      color = "Allele"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 9),
      legend.position = "bottom"
    )
  
  # Add geometries based on type
  if (type == "line") {
    p <- p + ggplot2::geom_line(linewidth = 0.8)
  } else if (type == "scatter") {
    p <- p + ggplot2::geom_point(alpha = 0.7, size = 1.5)
    
    # Add trend line
    if (trend_line == "linear") {
      p <- p + ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.8)
    } else if (trend_line == "quadratic") {
      p <- p + ggplot2::geom_smooth(
        method = "lm",
        formula = y ~ poly(x, 2),
        se = FALSE,
        linewidth = 0.8
      )
    } else {
      p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.8)
    }
  }
  
  # Add faceting
  p <- p + ggplot2::facet_wrap(~ facet_label, scales = "fixed", ncol = panels_per_row)
  
  # Add x-axis rotation
  if (x_axis_angle != 0) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_axis_angle, hjust = 1, vjust = 1))
  }
  
  return(p)
}