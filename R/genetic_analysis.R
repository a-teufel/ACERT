################################################################################
# R/genetic_analysis.R
################################################################################
# 
# Functions for genetic analysis including selection detection, allele tracking,
# mutation analysis, and population genetics of ABM simulation data
#
################################################################################

################################################################################
# FUNCTION: find_selected_mutations
################################################################################
#' Identify putatively selected mutations from time-series genetic data
#'
#' Analyzes time-series genetic data to identify alleles showing statistically
#' significant increases in frequency over time, indicating potential positive
#' selection. Uses linear regression to test for positive trends and integrates
#' mutation data for biological context.
#'
#' @param genind_list A named list of genind objects ordered by time. Names should
#'   contain tick numbers (e.g., "Run 8 500", "Run 8 1000").
#' @param mutation_data A data.frame with mutation information including columns:
#'   ticks, turtle, old_additivity_value, new_additivity_value, new_allele, 
#'   run_number, and FmI.
#' @param target_loci Character vector of locus names to analyze.
#' @param frequency_scope Either "global" (pool all individuals) or "population"
#'   (analyze each population separately). Default: "global".
#' @param analyze_all_alleles If TRUE, analyze all alleles. If FALSE, only analyze
#'   alleles absent at the initial time point. Default: FALSE.
#' @param significance_threshold P-value threshold for significance. Default: 0.05.
#'
#' @return A tibble with significant positive trends containing columns for
#'   locus, allele_name, population, slope, p_value, r_squared, std_error,
#'   t_value, new_additivity_value, turtle_origin, FmI, and run_number.
#'
#' @examples
#' \dontrun{
#' # Create genind list from turtle data
#' genind_list <- create_genind_time_series(turtle_data, time_points = c(500, 1000, 1500))
#' 
#' # Find selected mutations
#' selected <- find_selected_mutations(
#'   genind_list = genind_list,
#'   mutation_data = mutation_data,
#'   target_loci = letters[1:5],
#'   frequency_scope = "global"
#' )
#' }
#'
#' @importFrom adegenet tab seppop nInd
#' @importFrom dplyr %>% as_tibble bind_rows filter group_by mutate select distinct inner_join anti_join pull
#' @importFrom tidyr separate nest pivot_longer complete
#' @importFrom purrr map map_dbl
#' @importFrom broom tidy glance
#' @importFrom rlang .data
#' @importFrom stats lm
#' @importFrom stringr str_extract
#' @importFrom tibble tibble
#' @export
find_selected_mutations <- function(genind_list,
                                    mutation_data,
                                    target_loci,
                                    frequency_scope = "global",
                                    analyze_all_alleles = FALSE,
                                    significance_threshold = 0.05) {
  
  # Input validation
  if (!is.list(genind_list) || length(genind_list) == 0) {
    stop("genind_list must be a non-empty list of genind objects")
  }
  
  if (!all(sapply(genind_list, function(x) inherits(x, "genind")))) {
    stop("All elements in genind_list must be genind objects")
  }
  
  if (is.null(names(genind_list))) {
    stop("genind_list must have named elements containing time information")
  }
  
  if (!frequency_scope %in% c("global", "population")) {
    stop("frequency_scope must be either 'global' or 'population'")
  }
  
  if (!is.data.frame(mutation_data)) {
    stop("mutation_data must be a data.frame")
  }
  
  # Extract run and time information
  list_names <- names(genind_list)
  run_number_str <- unique(gsub(" \\d+$", "", list_names))
  ticks <- as.numeric(gsub(".* (\\d+)$", "\\1", list_names))
  
  if (length(run_number_str) != 1) {
    stop("All genind objects must be from the same simulation run")
  }
  
  if (any(is.na(ticks))) {
    stop("Could not extract time information from genind_list names")
  }
  
  message("Analyzing ", length(genind_list), " time points from ", run_number_str)
  
  # Calculate allele frequencies over time
  message("Calculating allele frequencies...")
  frequency_data <- calculate_allele_frequencies_time_series(
    genind_list, target_loci, frequency_scope
  )
  
  if (nrow(frequency_data) == 0) {
    warning("No frequency data calculated")
    return(create_empty_selection_results())
  }
  
  message("Total frequency records: ", nrow(frequency_data))
  
  # Filter alleles based on analysis scope
  if (!analyze_all_alleles) {
    frequency_data <- filter_new_alleles(frequency_data, genind_list[[1]], target_loci)
  }
  
  if (nrow(frequency_data) == 0) {
    warning("No alleles found for analysis after filtering")
    return(create_empty_selection_results())
  }
  
  message("Alleles to analyze: ", 
          length(unique(paste(frequency_data$locus, frequency_data$allele_name))))
  
  # Perform linear regression analysis
  message("Testing for selection signatures...")
  regression_results <- analyze_frequency_trends(frequency_data, significance_threshold)
  
  if (nrow(regression_results) == 0) {
    warning("No significant positive trends found")
    return(create_empty_selection_results())
  }
  
  message("Significant positive trends found: ", nrow(regression_results))
  
  # Cross-reference with mutation data
  message("Cross-referencing with mutation data...")
  final_results <- integrate_mutation_data(regression_results, mutation_data, run_number_str)
  
  message("Final results after joining with mutation data: ", nrow(final_results))
  
  return(final_results)
}
################################################################################
# FUNCTION: wrangle_major_alleles
################################################################################
#' Wrangle Major Allele Frequencies into a Tibble
#'
#' Takes a list of genind objects and wrangles them into a single tidy tibble
#' suitable for plotting. Efficiently calculates frequencies and identifies the
#' two most frequent alleles for each locus at the earliest time point.
#'
#' @param genind_list A list of genind objects. Each element should be named
#'   with a format containing a numerical tick value (e.g., "Run 8 500").
#' @param analysis_mode Character string: "global" or "by_population".
#'   If "global", major alleles are determined across all individuals.
#'   If "by_population", major alleles are determined within each population.
#'   Default: "by_population".
#'
#' @return A tibble with columns for tick, population, locus, allele,
#'   locus_allele, and frequency, ready for plotting.
#'
#' @examples
#' \dontrun{
#' # Wrangle major alleles globally
#' plot_data <- wrangle_major_alleles(
#'   genind_list,
#'   analysis_mode = "global"
#' )
#' 
#' # Wrangle by population
#' plot_data <- wrangle_major_alleles(
#'   genind_list,
#'   analysis_mode = "by_population"
#' )
#' }
#'
#' @importFrom adegenet tab pop nInd indNames
#' @importFrom stringr str_extract str_split_fixed
#' @importFrom dplyr %>% mutate group_by summarise ungroup inner_join filter arrange slice_head select desc
#' @importFrom tidyr crossing replace_na pivot_longer
#' @importFrom tibble tibble rownames_to_column
#' @export
wrangle_major_alleles <- function(genind_list, analysis_mode = "by_population") {
  
  # Input validation
  if (!is.list(genind_list) || length(genind_list) == 0) {
    stop("genind_list must be a non-empty list of genind objects")
  }
  
  if (!all(sapply(genind_list, function(x) inherits(x, "genind")))) {
    stop("All elements must be genind objects")
  }
  
  if (is.null(names(genind_list))) {
    stop("genind_list elements must be named")
  }
  
  if (!analysis_mode %in% c("global", "by_population")) {
    stop("analysis_mode must be 'global' or 'by_population'")
  }
  
  message("Wrangling allele frequency data for ", length(genind_list), " time points...")
  
  # Extract frequencies for each time point
  extract_func <- function(genind_obj, tick_name) {
    tick_number <- as.numeric(stringr::str_extract(tick_name, "\\d+$"))
    if (is.na(tick_number)) {
      warning("Could not extract tick from: ", tick_name)
      return(NULL)
    }
    
    if (analysis_mode == "global") {
      counts_tab <- adegenet::tab(genind_obj, NA.method = "mean")
      total_alleles <- 2 * adegenet::nInd(genind_obj)
      global_freqs <- colSums(counts_tab, na.rm = TRUE) / total_alleles
      
      freq_df <- tibble::tibble(
        locus_allele = names(global_freqs),
        frequency = global_freqs
      ) %>%
        dplyr::mutate(population = "Global")
      
    } else {
      if (is.null(adegenet::pop(genind_obj))) {
        warning("No populations defined for tick: ", tick_number)
        return(NULL)
      }
      
      # Calculate frequencies by population
      pop_counts_matrix <- adegenet::tab(genind_obj, pop = FALSE, NA.method = "mean")
      pop_data <- tibble::tibble(
        population = adegenet::pop(genind_obj),
        genotype = as.character(adegenet::indNames(genind_obj))
      )
      
      freq_df <- as.data.frame(pop_counts_matrix) %>%
        tibble::rownames_to_column(var = "genotype") %>%
        dplyr::inner_join(pop_data, by = "genotype") %>%
        tidyr::pivot_longer(
          cols = -c(genotype, population),
          names_to = "locus_allele",
          values_to = "count"
        ) %>%
        dplyr::group_by(population, locus_allele) %>%
        dplyr::summarise(
          total_alleles = 2 * dplyr::n(),
          total_count = sum(.data$count, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(frequency = total_count / total_alleles) %>%
        dplyr::select(population, locus_allele, frequency)
    }
    
    # Split locus_allele into separate columns
    split_cols <- stringr::str_split_fixed(freq_df$locus_allele, pattern = "\\.", n = 2)
    freq_df %>%
      dplyr::mutate(
        tick = tick_number,
        locus = split_cols[, 1],
        allele = split_cols[, 2]
      )
  }
  
  # Extract frequencies for all time points
  all_freq_data <- mapply(extract_func, genind_list, names(genind_list), SIMPLIFY = FALSE) %>%
    dplyr::bind_rows()
  
  # Identify major alleles at first time point
  first_tick_data <- all_freq_data %>%
    dplyr::filter(tick == min(all_freq_data$tick))
  
  major_alleles <- first_tick_data %>%
    dplyr::group_by(locus, population) %>%
    dplyr::arrange(dplyr::desc(frequency)) %>%
    dplyr::slice_head(n = 2) %>%
    dplyr::ungroup() %>%
    dplyr::select(locus, population, allele, locus_allele)
  
  message("Identified ", nrow(major_alleles), " major alleles to track")
  
  # Create complete time series grid
  all_ticks <- unique(all_freq_data$tick)
  complete_grid <- major_alleles %>%
    tidyr::crossing(tick = all_ticks)
  
  # Join with actual frequency data
  all_freq_data_clean <- all_freq_data %>%
    dplyr::distinct(tick, population, locus_allele, .keep_all = TRUE)
  
  plot_data <- complete_grid %>%
    dplyr::left_join(all_freq_data_clean, by = c("tick", "population", "locus_allele")) %>%
    dplyr::mutate(
      locus = dplyr::if_else(is.na(locus.y), locus.x, locus.y),
      allele = dplyr::if_else(is.na(allele.y), allele.x, allele.y),
      frequency = tidyr::replace_na(frequency, 0)
    ) %>%
    dplyr::select(locus, population, allele, locus_allele, tick, frequency)
  
  return(plot_data)
}

################################################################################
# FUNCTION: calculate_effective_migrants
################################################################################
#' Calculate effective migration between populations
#'
#' Analyzes agent-based model data to quantify the number of effective migrants
#' between populations at each time step.
#'
#' @param turtle_data A data.frame containing simulation data with columns:
#'   run_number, ticks, who, breed, t.dad.id, t.mom.id, t.birth.pond.index.
#' @param summarize_across_time If FALSE, returns detailed migration counts for
#'   each time step. If TRUE, returns median migration rates. Default: FALSE.
#'
#' @return A tibble with migration data.
#'
#' @examples
#' \dontrun{
#' # Calculate detailed migration for each time step
#' migration_detailed <- calculate_effective_migrants(turtle_data)
#' }
#'
#' @importFrom dplyr %>% filter select left_join rename group_by summarise n_distinct arrange
#' @importFrom tidyr pivot_longer complete
#' @importFrom stats median
#' @export
calculate_effective_migrants <- function(turtle_data, summarize_across_time = FALSE) {
  
  # Input validation
  required_cols <- c("run_number", "ticks", "who", "breed", "t.dad.id", "t.mom.id", "t.birth.pond.index")
  missing_cols <- setdiff(required_cols, names(turtle_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  message("Identifying parent-offspring relationships...")
  
  # Identify offspring and their parents
  offspring_data <- turtle_data %>%
    dplyr::filter(breed %in% c("eggs", "larvae")) %>%
    dplyr::select(run_number, ticks, who, t.dad.id, t.mom.id, t.birth.pond.index)
  
  if (nrow(offspring_data) == 0) {
    warning("No offspring found in data")
    return(create_empty_migration_results(summarize_across_time))
  }
  
  # Get parent birth pond information
  parent_birth_ponds <- turtle_data %>%
    dplyr::select(who, parent_birth_pond = t.birth.pond.index) %>%
    dplyr::distinct()
  
  message("Tracing migration events...")
  
  # Create parent-offspring migration events
  migrant_events <- offspring_data %>%
    # Add dad birth pond
    dplyr::left_join(parent_birth_ponds, by = c("t.dad.id" = "who")) %>%
    dplyr::rename(dad_birth_pond = parent_birth_pond) %>%
    # Add mom birth pond
    dplyr::left_join(parent_birth_ponds, by = c("t.mom.id" = "who")) %>%
    dplyr::rename(mom_birth_pond = parent_birth_pond) %>%
    # Reshape to have one row per parent
    tidyr::pivot_longer(
      cols = c(t.dad.id, t.mom.id),
      names_to = "parent_type", 
      values_to = "parent_id"
    ) %>%
    tidyr::pivot_longer(
      cols = c(dad_birth_pond, mom_birth_pond),
      names_to = "pond_type",
      values_to = "parent_origin_pond"
    ) %>%
    # Keep matching parent-pond combinations
    dplyr::filter(
      (parent_type == "t.dad.id" & pond_type == "dad_birth_pond") |
        (parent_type == "t.mom.id" & pond_type == "mom_birth_pond")
    ) %>%
    dplyr::select(-parent_type, -pond_type)
  
  # Calculate effective migrants per time step
  effective_migrants_by_tick <- migrant_events %>%
    dplyr::filter(!is.na(parent_id)) %>%
    dplyr::filter(t.birth.pond.index != parent_origin_pond) %>%
    dplyr::group_by(run_number, ticks, 
                    origin_pond = parent_origin_pond, 
                    destination_pond = t.birth.pond.index) %>%
    dplyr::summarise(effective_migrants = dplyr::n_distinct(parent_id), .groups = "drop")
  
  # Complete missing combinations with zeros
  if (nrow(effective_migrants_by_tick) > 0) {
    effective_migrants_by_tick <- effective_migrants_by_tick %>%
      tidyr::complete(run_number, ticks, origin_pond, destination_pond, 
                      fill = list(effective_migrants = 0))
  }
  
  message("Migration events calculated for ", length(unique(effective_migrants_by_tick$ticks)), " time points")
  
  # Handle summarization option
  if (summarize_across_time) {
    if (nrow(effective_migrants_by_tick) == 0) {
      return(create_empty_migration_results(TRUE))
    }
    
    final_output <- effective_migrants_by_tick %>%
      dplyr::group_by(run_number, origin_pond, destination_pond) %>%
      dplyr::summarise(
        median_effective_migrants = median(effective_migrants, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    final_output <- effective_migrants_by_tick %>%
      dplyr::arrange(run_number, ticks, origin_pond, destination_pond)
  }
  
  return(final_output)
}

################################################################################
# FUNCTION: analyze_mutation_effects
################################################################################
#' Analyze mutation effect sizes and their distribution
#'
#' Provides comprehensive analysis of mutation effects including statistical
#' summaries, temporal patterns, and effect size distributions.
#'
#' @param mutation_data A data.frame containing mutation information.
#' @param group_by Character vector of columns to group analysis by. 
#'   Default: c("run_number", "ticks").
#' @param effect_size_column Name of the column containing effect sizes. Default: "FmI".
#'
#' @return A tibble with mutation effect statistics.
#'
#' @examples
#' \dontrun{
#' # Basic mutation effect analysis
#' mutation_stats <- analyze_mutation_effects(mutation_data)
#' }
#'
#' @importFrom dplyr group_by summarise across all_of n
#' @importFrom stats sd
#' @export
analyze_mutation_effects <- function(mutation_data, 
                                     group_by = c("run_number", "ticks"),
                                     effect_size_column = "FmI") {
  
  # Input validation
  if (!is.data.frame(mutation_data)) {
    stop("mutation_data must be a data.frame")
  }
  
  if (!effect_size_column %in% names(mutation_data)) {
    # Try to calculate FmI if possible
    if (all(c("new_additivity_value", "old_additivity_value") %in% names(mutation_data))) {
      mutation_data$FmI <- mutation_data$new_additivity_value - mutation_data$old_additivity_value
      message("Calculated FmI from additivity values")
    } else {
      stop("Effect size column '", effect_size_column, "' not found and cannot be calculated")
    }
  }
  
  missing_group_cols <- setdiff(group_by, names(mutation_data))
  if (length(missing_group_cols) > 0) {
    stop("Grouping columns not found: ", paste(missing_group_cols, collapse = ", "))
  }
  
  # Remove rows with missing effect sizes
  clean_data <- mutation_data[!is.na(mutation_data[[effect_size_column]]), ]
  
  if (nrow(clean_data) == 0) {
    stop("No non-missing values found in effect size column")
  }
  
  # Calculate statistics
  effect_stats <- clean_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_by))) %>%
    dplyr::summarise(
      n_mutations = dplyr::n(),
      mean_effect = mean(.data[[effect_size_column]], na.rm = TRUE),
      median_effect = median(.data[[effect_size_column]], na.rm = TRUE),
      sd_effect = sd(.data[[effect_size_column]], na.rm = TRUE),
      min_effect = min(.data[[effect_size_column]], na.rm = TRUE),
      max_effect = max(.data[[effect_size_column]], na.rm = TRUE),
      prop_beneficial = mean(.data[[effect_size_column]] > 0, na.rm = TRUE),
      prop_neutral = mean(.data[[effect_size_column]] == 0, na.rm = TRUE),
      prop_deleterious = mean(.data[[effect_size_column]] < 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(effect_stats)
}

################################################################################
# FUNCTION: create_genind_time_series
################################################################################
#' Create time series of genind objects from turtle data
#'
#' Converts turtle simulation data into a time-ordered list of genind objects
#' suitable for temporal genetic analysis.
#'
#' @param turtle_data A data.frame containing turtle genetic data.
#' @param time_points Numeric vector of specific time points to include. Default: NULL.
#' @param loci_names Character vector of locus names to include. Default: letters[1:10].
#' @param min_individuals_per_pop Minimum number of individuals required per
#'   population to include in analysis. Default: 5.
#' @param run_identifier Column name containing run information. Default: "run_number".
#'
#' @return A named list of genind objects.
#'
#' @examples
#' \dontrun{
#' # Create genind time series for specific time points
#' genind_series <- create_genind_time_series(
#'   turtle_data,
#'   time_points = c(500, 1000, 1500, 2000),
#'   loci_names = letters[1:5]
#' )
#' }
#'
#' @importFrom adegenet df2genind
#' @export
create_genind_time_series <- function(turtle_data,
                                      time_points = NULL,
                                      loci_names = letters[1:10],
                                      min_individuals_per_pop = 5,
                                      run_identifier = "run_number") {
  
  # Input validation
  if (!is.data.frame(turtle_data)) {
    stop("turtle_data must be a data.frame")
  }
  
  required_cols <- c("ticks", "t.birth.pond.index", run_identifier)
  missing_cols <- setdiff(required_cols, names(turtle_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Filter to specified time points
  if (!is.null(time_points)) {
    turtle_data <- turtle_data %>%
      dplyr::filter(.data$ticks %in% time_points)
    
    if (nrow(turtle_data) == 0) {
      stop("No data found for specified time points")
    }
  }
  
  # Get unique run-time combinations
  time_combos <- turtle_data %>%
    dplyr::select(dplyr::all_of(c(run_identifier, "ticks"))) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data[[run_identifier]], .data$ticks)
  
  message("Creating genind objects for ", nrow(time_combos), " time points")
  
  # Convert to pegas format
  pegas_data <- convert_to_pegas_format(turtle_data, loci_names)
  
  genind_list <- list()
  
  for (i in 1:nrow(time_combos)) {
    run_id <- time_combos[[run_identifier]][i]
    time_point <- time_combos$ticks[i]
    
    # Filter data for this run-time combination
    time_data <- pegas_data %>%
      dplyr::filter(.data[[run_identifier]] == run_id, .data$ticks == time_point)
    
    if (nrow(time_data) == 0) {
      warning("No data for ", run_id, " at time ", time_point)
      next
    }
    
    # Check population sizes
    pop_sizes <- table(time_data$pop)
    valid_pops <- names(pop_sizes)[pop_sizes >= min_individuals_per_pop]
    
    if (length(valid_pops) == 0) {
      warning("No populations with sufficient individuals at ", run_id, " time ", time_point)
      next
    }
    
    # Filter to valid populations
    time_data_filtered <- time_data %>%
      dplyr::filter(.data$pop %in% valid_pops)
    
    # Create genind object
    tryCatch({
      genind_obj <- adegenet::df2genind(
        time_data_filtered[, loci_names],
        sep = "/",
        pop = as.factor(time_data_filtered$pop),
        ploidy = 2
      )
      
      # Create descriptive name
      list_name <- paste(run_id, time_point)
      genind_list[[list_name]] <- genind_obj
      
    }, error = function(e) {
      warning("Failed to create genind for ", run_id, " at time ", time_point, ": ", e$message)
    })
  }
  
  if (length(genind_list) == 0) {
    stop("No genind objects could be created")
  }
  
  message("Successfully created ", length(genind_list), " genind objects")
  
  return(genind_list)
}


################################################################################
# FUNCTION: convert_to_pegas_format
################################################################################
#' Convert turtle chromosome data to pegas format for analysis
#'
#' Transforms ABM turtle genetic data with expanded chromosome columns into
#' pegas-compatible format for population genetic analysis.
#'
#' @param turtle_data Data.frame containing turtle data with expanded chromosome columns.
#' @param loci_names Character vector of locus names. Default: letters[1:10].
#' @param chunk_size Number of rows to process at once for memory efficiency. Default: 10000.
#'
#' @return Data.frame suitable for conversion to pegas format with genotype
#'   columns, population assignments, and time information.
#'
#' @examples
#' \dontrun{
#' # Convert turtle data to pegas format
#' pegas_data <- convert_to_pegas_format(
#'   turtle_data, 
#'   loci_names = letters[1:5]
#' )
#' }
#'
#' @importFrom stats setNames
#' @export
convert_to_pegas_format <- function(turtle_data, 
                                    loci_names = letters[1:10], 
                                    chunk_size = 10000) {
  
  # Input validation
  if (!is.data.frame(turtle_data)) {
    stop("turtle_data must be a data.frame")
  }
  
  # Check for expanded chromosome columns
  chr_cols <- grep("chromosome.*\\.[a-z]$", names(turtle_data), value = TRUE)
  if (length(chr_cols) == 0) {
    stop("No expanded chromosome columns found. Data must have expanded chromosome format.")
  }
  
  # Validate required columns
  required_cols <- c("who", "t.birth.pond.index", "ticks")
  missing_cols <- setdiff(required_cols, names(turtle_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Get unique loci from column names
  detected_loci <- unique(sub(".*\\.", "", chr_cols))
  actual_loci <- intersect(loci_names, detected_loci)
  
  if (length(actual_loci) == 0) {
    stop("No matching loci found between loci_names and data columns")
  }
  
  message("Processing ", length(actual_loci), " loci: ", paste(actual_loci, collapse = ", "))
  
  # Process data in chunks for memory efficiency
  n_rows <- nrow(turtle_data)
  n_chunks <- ceiling(n_rows / chunk_size)
  
  #result_chunks <- list()
  result_chunks <- vector("list", n_chunks)
  
  for (chunk in 1:n_chunks) {
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, n_rows)
    
    chunk_data <- turtle_data[start_idx:end_idx, ]
    
    # Create genotype matrix for this chunk
    genotype_matrix <- matrix(NA_character_, 
                              nrow = nrow(chunk_data), 
                              ncol = length(actual_loci))
    colnames(genotype_matrix) <- actual_loci
    
    # Process each locus
    for (j in seq_along(actual_loci)) {
      locus <- actual_loci[j]
      
      # Find chromosome columns for this locus
      dad1_col <- paste0("t.dad.chromosome.1.", locus)
      dad2_col <- paste0("t.dad.chromosome.2.", locus)
      mom1_col <- paste0("t.mom.chromosome.1.", locus)
      mom2_col <- paste0("t.mom.chromosome.2.", locus)
      
      # Extract alleles for each individual
      for (i in 1:nrow(chunk_data)) {
        alleles <- character(0)
        
        # Collect alleles from all four chromosomes
        for (col in c(dad1_col, dad2_col, mom1_col, mom2_col)) {
          if (col %in% names(chunk_data) && !is.na(chunk_data[[col]][i])) {
            # Extract numeric part of allele (e.g., "a1" -> "1")
            allele_string <- chunk_data[[col]][i]
            if (nchar(allele_string) > 1) {
              allele_value <- substr(allele_string, 2, nchar(allele_string))
              alleles <- c(alleles, allele_value)
            }
          }
        }
        
        # Create diploid genotype from first two alleles
        if (length(alleles) >= 2) {
          genotype_matrix[i, j] <- paste0(alleles[1], "/", alleles[2])
        } else if (length(alleles) == 1) {
          # Homozygote
          genotype_matrix[i, j] <- paste0(alleles[1], "/", alleles[1])
        }
        # else remains NA
      }
    }
    
    # Combine with metadata
    chunk_result <- data.frame(
      genotype_matrix,
      pop = chunk_data$t.birth.pond.index,
      ticks = chunk_data$ticks,
      individual_id = chunk_data$who,
      run_number = if("run_number" %in% names(chunk_data)) chunk_data$run_number else NA,
      stringsAsFactors = FALSE
    )
    
    result_chunks[[chunk]] <- chunk_result
    
    # Clean up intermediate objects
    rm(chunk_data, genotype_matrix)
    if(chunk %% 10 == 0) gc(verbose = FALSE)
  }
  
  # Combine all chunks
  final_result <- do.call(rbind, result_chunks)
  rownames(final_result) <- final_result$individual_id
  
  return(final_result)
}

################################################################################
# FUNCTION: balance_genind_populations
################################################################################
#' Balance genind Populations
#'
#' This function takes a genind object and samples each population down to a
#' specified target size. Populations with fewer individuals than the target
#' size are left as is.
#'
#' @param x A genind object containing the populations to be balanced.
#' @param target_sample_size An integer specifying the target number of
#'   individuals per population.
#' @param seed An optional integer to set the random seed for reproducible
#'   sampling. Defaults to NULL.
#'
#' @return A new genind object where each population has been sampled down to
#'   the target size or remains its original size if smaller.
#'   
#' @export
balance_genind_populations <- function(x, target_sample_size, seed = NULL) {
  # Check if a seed is provided and set it
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Separate the genind object into a list of populations
  pops_list <- adegenet::seppop(x)
  
  # Initialize a list to hold the sampled populations
  sampled_pops_list <- list()
  
  # Iterate through each population
  for (pop_name in names(pops_list)) {
    current_pop <- pops_list[[pop_name]]
    
    # Check if the population size is greater than the target size
    if (adegenet::nInd(current_pop) > target_sample_size) {
      # Sample down to the target size without replacement
      sampled_indices <- sample(
        1:adegenet::nInd(current_pop),
        size = target_sample_size,
        replace = FALSE
      )
      sampled_pops_list[[pop_name]] <- current_pop[sampled_indices, ]
    } else {
      # If the population is already smaller, keep it as is
      sampled_pops_list[[pop_name]] <- current_pop
    }
  }
  
  # Recombine the sampled populations into a single genind object
  final_genind_obj <- adegenet::repool(sampled_pops_list)
  
  return(final_genind_obj)
}

################################################################################
# FUNCTION: get_effective_alleles_from_freq
################################################################################
#' @title Calculate Effective Number of Alleles and Expected Heterozygosity from Frequencies
#'
#' @description Calculates both the effective number of alleles (\eqn{A_e}) and the
#' expected heterozygosity (\eqn{H_e}) for each locus within each population.
#' \eqn{H_e} is calculated as \eqn{1 - \sum p_i^2}. \eqn{A_e} is calculated as \eqn{1 / (1 - H_e)}.
#'
#' @param pop_freq_list A list of matrices or data frames, typically the `pop.freq`
#' component from the output of `hierfstat::basic.stats`. Each element in the list
#' represents a locus, with rows being alleles and columns being populations,
#' containing allele frequencies.
#'
#' @return A list containing two data frames:
#' \itemize{
#'   \item \code{Ae}: A data frame where rows are locus names and columns are population names.
#'         Each cell contains the effective number of alleles ($A_e$) for that locus
#'         within that specific population.
#'   \item \code{He}: A data frame where rows are locus names and columns are population names.
#'         Each cell contains the expected heterozygosity ($H_e$) for that locus
#'         within that specific population.
#' }
#'
#' @examples
#' # Example using a simplified structure similar to rbsout$pop.freq
#' simple_pop_freq_example <- list(
#'   a = data.frame(
#'     pop1 = c(0.5341, 0.4545, 0.0114), # Sums to 1.0000
#'     pop2 = c(0.4886, 0.5114, 0.0000), # Sums to 1.0000
#'     row.names = c("allele1", "allele2", "allele3")
#'   ),
#'   b = data.frame(
#'     pop1 = c(0.4886, 0.0114, 0.0341, 0.3523, 0.1136), # Sums to 1.0000
#'     pop2 = c(0.6023, 0.0000, 0.0000, 0.3977, 0.0000), # Sums to 1.0000
#'     row.names = c("alleleA", "alleleB", "alleleC", "alleleD", "alleleE")
#'   )
#' )
#'
#' # Use the function
#' results <- get_effective_alleles_from_freq(simple_pop_freq_example)
#'
#' @export
get_effective_alleles_from_freq <- function(pop_freq_list) {
  expected_heterozygosity_matrix <- sapply(pop_freq_list, function(locus_df) {
    1 - colSums(locus_df^2)
  })
  
  effective_alleles_matrix <- 1 / expected_heterozygosity_matrix
  effective_alleles_matrix[is.infinite(effective_alleles_matrix)] <- 1
  
  # Transpose and convert to data frames
  df_He <- as.data.frame(t(expected_heterozygosity_matrix))
  df_Ae <- as.data.frame(t(effective_alleles_matrix))
  
  # Return a list containing both data frames
  list(Ae = df_Ae, He = df_He)
}

################################################################################
# FUNCTION: get_num_observed_alleles_per_locus
################################################################################
#' @title Count Observed Alleles Per Locus Per Population
#'
#' @description Calculates the number of alleles with frequencies greater than zero
#' for each locus, specific to each population.
#' This function expects a list where each element is a matrix or data frame
#' (representing a locus), with rows as alleles and columns as populations.
#'
#' @param pop_freq_list A list of matrices/data frames, typically `rbsout$pop.freq`
#' from `hierfstat::basic.stats`.
#'
#' @return A data frame where rows are locus names and columns are population names.
#' Each cell contains the count of alleles with frequencies > 0 for that locus
#' within that specific population.
#'
#' @examples
#' # Example using a simplified structure similar to rbsout$pop.freq
#' simple_pop_freq_example <- list(
#'   locusA = data.frame(
#'     pop1 = c(0.5, 0.4, 0.1),
#'     pop2 = c(0.6, 0.0, 0.4), # Allele A2 has 0 frequency in pop2
#'     row.names = c("A1", "A2", "A3")
#'   ),
#'   locusB = data.frame(
#'     pop1 = c(0.7, 0.3),
#'     pop2 = c(0.0, 1.0), # Allele B1 has 0 frequency in pop2
#'     row.names = c("B1", "B2")
#'   )
#' )
#'
#' # Use the function
#' get_num_observed_alleles_per_locus(simple_pop_freq_example)
#'
#' @export
get_num_observed_alleles_per_locus <- function(pop_freq_list) {
  allele_counts_matrix <- sapply(pop_freq_list, function(locus_df) {
    colSums(locus_df > 0)
  })
  as.data.frame(t(allele_counts_matrix))
}

################################################################################
# FUNCTION: restore_basic_stats_names
################################################################################
#' Restore Original Names to basic.stats Output
#'
#' Restores original locus and allele names from a genind object to the output
#' of hierfstat's basic.stats() function. The genind2hierfstat() conversion
#' replaces meaningful locus names with generic ones (X1, X2, etc.) and
#' renumbers alleles. This function reverses that process.
#'
#' @param genind_obj A genind object (from adegenet package) that was used
#'   to generate the basic.stats output.
#' @param basic_stats_output Output from hierfstat's basic.stats() function,
#'   typically generated after converting the genind object with genind2hierfstat().
#' @param verbose Logical indicating whether to print mapping information.
#'   Default is FALSE.
#'
#' @return A modified basic.stats output object with original locus and allele
#'   names restored where applicable.
#'
#' @import adegenet
#' @importFrom methods slot
#'
#' @examples
#' \dontrun{
#' library(adegenet)
#' library(hierfstat)
#' 
#' # Load example data
#' data(nancycats)
#' 
#' # Run basic.stats
#' hierfstat_data <- genind2hierfstat(nancycats)
#' bs_output <- basic.stats(hierfstat_data)
#' 
#' # Restore original names
#' bs_restored <- restore_basic_stats_names(nancycats, bs_output)
#' }
#'
#' @export
restore_basic_stats_names <- function(genind_obj, 
                                      basic_stats_output, 
                                      verbose = FALSE) {
  
  # Input validation
  if (!inherits(genind_obj, "genind")) {
    stop("genind_obj must be a genind object from the adegenet package")
  }
  
  if (!is.list(basic_stats_output)) {
    stop("basic_stats_output must be a list (output from hierfstat::basic.stats)")
  }
  
  # Extract original names
  original_loci <- adegenet::locNames(genind_obj)
  original_alleles <- genind_obj@all.names
  
  if (verbose) {
    cat("Original loci:", length(original_loci), "\n")
    cat("Original locus names:", paste(original_loci, collapse = ", "), "\n")
  }
  
  # Create locus mapping (X1 -> original_name)
  hierfstat_loci <- paste0("X", 1:length(original_loci))
  locus_mapping <- setNames(original_loci, hierfstat_loci)
  
  # Create allele mapping
  allele_mapping <- list()
  for (i in 1:length(original_alleles)) {
    locus_name <- names(original_alleles)[i]
    hierfstat_locus <- paste0("X", i)
    
    # hierfstat typically uses 01, 02, 03, etc. for alleles
    n_alleles <- length(original_alleles[[i]])
    hierfstat_alleles <- sprintf("%02d", 1:n_alleles)
    
    allele_mapping[[hierfstat_locus]] <- setNames(
      original_alleles[[i]], 
      hierfstat_alleles
    )
  }
  
  # Function to replace locus names in a character vector
  replace_locus_names <- function(name_vector) {
    if (is.null(name_vector)) return(NULL)
    
    # Replace each hierfstat locus name with original name
    new_names <- name_vector
    for (i in seq_along(locus_mapping)) {
      hierfstat_name <- names(locus_mapping)[i]
      original_name <- locus_mapping[i]
      new_names <- gsub(paste0("^", hierfstat_name, "$"), original_name, new_names)
    }
    
    return(new_names)
  }
  
  # Function to replace allele names in a character vector
  replace_allele_names <- function(name_vector, locus_context = NULL) {
    if (is.null(name_vector)) return(NULL)
    
    new_names <- sapply(name_vector, function(x) {
      # Check if it looks like a locus.allele format
      if (grepl("\\.", x)) {
        parts <- strsplit(x, "\\.")[[1]]
        if (length(parts) == 2) {
          locus_part <- parts[1]
          allele_part <- parts[2]
          
          # Try to map both parts
          if (locus_part %in% names(locus_mapping) && 
              locus_part %in% names(allele_mapping) &&
              allele_part %in% names(allele_mapping[[locus_part]])) {
            
            original_locus <- locus_mapping[locus_part]
            original_allele <- allele_mapping[[locus_part]][allele_part]
            
            return(paste0(original_locus, ".", original_allele))
          }
        }
      }
      
      # Check if we have a locus context and this looks like a simple allele code
      if (!is.null(locus_context) && locus_context %in% names(allele_mapping)) {
        # Format as zero-padded string to match our mapping
        allele_code <- sprintf("%02d", as.numeric(x))
        if (allele_code %in% names(allele_mapping[[locus_context]])) {
          return(allele_mapping[[locus_context]][allele_code])
        }
      }
      
      return(x)  # Return unchanged if no mapping found
    }, USE.NAMES = FALSE)
    
    return(new_names)
  }
  
  # Recursive function to process any object and replace names
  replace_names_recursive <- function(obj, parent_name = NULL) {
    if (is.null(obj)) {
      return(obj)
    }
    
    # Handle vectors
    if (is.vector(obj) && !is.list(obj)) {
      # Replace names of the vector
      if (!is.null(names(obj))) {
        names(obj) <- replace_locus_names(names(obj))
        names(obj) <- replace_allele_names(names(obj), parent_name)
      }
      return(obj)
    }
    
    # Handle matrices
    if (is.matrix(obj)) {
      # Replace rownames
      if (!is.null(rownames(obj))) {
        rownames(obj) <- replace_locus_names(rownames(obj))
        
        # For pop.freq structure, use parent_name as locus context
        if (!is.null(parent_name) && parent_name %in% locus_mapping) {
          # Find the hierfstat locus name for this original locus name
          hierfstat_locus <- names(locus_mapping)[locus_mapping == parent_name]
          if (length(hierfstat_locus) > 0) {
            rownames(obj) <- replace_allele_names(rownames(obj), hierfstat_locus[1])
          }
        } else {
          rownames(obj) <- replace_allele_names(rownames(obj), parent_name)
        }
      }
      # Replace colnames
      if (!is.null(colnames(obj))) {
        colnames(obj) <- replace_locus_names(colnames(obj))
        colnames(obj) <- replace_allele_names(colnames(obj), parent_name)
      }
      return(obj)
    }
    
    # Handle data frames
    if (is.data.frame(obj)) {
      # Replace rownames
      if (!is.null(rownames(obj))) {
        rownames(obj) <- replace_locus_names(rownames(obj))
        
        # For pop.freq structure, use parent_name as locus context
        if (!is.null(parent_name) && parent_name %in% locus_mapping) {
          # Find the hierfstat locus name for this original locus name
          hierfstat_locus <- names(locus_mapping)[locus_mapping == parent_name]
          if (length(hierfstat_locus) > 0) {
            rownames(obj) <- replace_allele_names(rownames(obj), hierfstat_locus[1])
          }
        } else {
          rownames(obj) <- replace_allele_names(rownames(obj), parent_name)
        }
      }
      # Replace colnames
      if (!is.null(colnames(obj))) {
        colnames(obj) <- replace_locus_names(colnames(obj))
        colnames(obj) <- replace_allele_names(colnames(obj), parent_name)
      }
      return(obj)
    }
    
    # Handle arrays
    if (is.array(obj)) {
      # Replace dimnames
      if (!is.null(dimnames(obj))) {
        for (i in seq_along(dimnames(obj))) {
          if (!is.null(dimnames(obj)[[i]])) {
            dimnames(obj)[[i]] <- replace_locus_names(dimnames(obj)[[i]])
            dimnames(obj)[[i]] <- replace_allele_names(dimnames(obj)[[i]], parent_name)
          }
        }
      }
      return(obj)
    }
    
    # Handle lists (recursively)
    if (is.list(obj)) {
      # Replace names of list elements
      if (!is.null(names(obj))) {
        names(obj) <- replace_locus_names(names(obj))
        names(obj) <- replace_allele_names(names(obj), parent_name)
      }
      
      # Recursively process each element, passing the element name as context
      for (i in seq_along(obj)) {
        element_name <- names(obj)[i]
        obj[[i]] <- replace_names_recursive(obj[[i]], element_name)
      }
      return(obj)
    }
    
    # For other object types, return as-is
    return(obj)
  }
  
  # Create a copy of the output to modify
  restored_output <- basic_stats_output
  
  # Process the entire structure recursively
  restored_output <- replace_names_recursive(restored_output)
  
  # Add mapping information as attributes
  attr(restored_output, "locus_mapping") <- locus_mapping
  attr(restored_output, "allele_mapping") <- allele_mapping
  attr(restored_output, "original_genind") <- deparse(substitute(genind_obj))
  
  if (verbose) {
    cat("Name restoration complete!\n")
    cat("Processed components:", paste(names(restored_output), collapse = ", "), "\n")
  }
  
  return(restored_output)
}

################################################################################
# HELPER FUNCTIONS
################################################################################

#' Calculate allele frequencies across time series
#' @param genind_list List of genind objects
#' @param loci_names Vector of locus names or NULL for all
#' @param frequency_scope "global" or "population"
#' @return Tibble with frequency data
#' @keywords internal
calculate_allele_frequencies_time_series <- function(genind_list, loci_names, frequency_scope) {
  frequency_data <- list()
  
  for (i in seq_along(genind_list)) {
    genind_obj <- genind_list[[i]]
    tick_number <- as.numeric(stringr::str_extract(names(genind_list)[i], "\\d+$"))
    
    # Filter to target loci if specified
    if (!is.null(loci_names)) {
      genind_obj <- genind_obj[loc = loci_names]
    }
    
    if (frequency_scope == "global") {
      # Global frequencies
      allele_counts <- adegenet::tab(genind_obj, freq = FALSE, NA.method = "mean")
      total_alleles <- 2 * adegenet::nInd(genind_obj)
      allele_freqs <- colSums(allele_counts, na.rm = TRUE) / total_alleles
      
      freq_df <- tibble::tibble(
        allele = names(allele_freqs),
        frequency = allele_freqs,
        tick = tick_number,
        population = "global"
      )
    } else {
      # Population-specific frequencies
      if (is.null(adegenet::pop(genind_obj))) {
        warning("No population assignments found for tick ", tick_number)
        next
      }
      
      pop_list <- adegenet::seppop(genind_obj)
      pop_freq_list <- list()
      
      for (j in seq_along(pop_list)) {
        pop_obj <- pop_list[[j]]
        pop_name <- names(pop_list)[j]
        
        if (adegenet::nInd(pop_obj) > 0) {
          pop_allele_counts <- adegenet::tab(pop_obj, freq = FALSE, NA.method = "mean")
          pop_total_alleles <- 2 * adegenet::nInd(pop_obj)
          pop_allele_freqs <- colSums(pop_allele_counts, na.rm = TRUE) / pop_total_alleles
          
          pop_freq_list[[j]] <- tibble::tibble(
            allele = names(pop_allele_freqs),
            frequency = pop_allele_freqs,
            tick = tick_number,
            population = pop_name
          )
        }
      }
      
      freq_df <- dplyr::bind_rows(pop_freq_list)
    }
    
    frequency_data[[i]] <- freq_df
  }
  
  # Combine and format
  all_freqs <- dplyr::bind_rows(frequency_data) %>%
    dplyr::filter(frequency > 0) %>%
    tidyr::separate(.data$allele, into = c("locus", "allele_name"), sep = "\\.", extra = "merge") %>%
    dplyr::mutate(locus_allele = paste(.data$locus, .data$allele_name, sep = "."))
  
  return(all_freqs)
}

#' Filter to new alleles not present initially
#' @param frequency_data Frequency data tibble
#' @param initial_genind First genind object
#' @param target_loci Vector of target loci
#' @return Filtered frequency data
#' @keywords internal
filter_new_alleles <- function(frequency_data, initial_genind, target_loci) {
  # Get alleles present at first time point
  initial_genind_filtered <- initial_genind[loc = target_loci]
  initial_alleles <- colnames(adegenet::tab(initial_genind_filtered, freq = FALSE))
  
  # Convert to same format as frequency data
  initial_alleles_df <- tibble::tibble(allele = initial_alleles) %>%
    tidyr::separate(.data$allele, into = c("locus", "allele_name"), sep = "\\.", extra = "merge")
  
  # Filter out initial alleles
  filtered_data <- frequency_data %>%
    dplyr::anti_join(initial_alleles_df, by = c("locus", "allele_name"))
  
  return(filtered_data)
}

#' Analyze frequency trends with linear regression
#' @param frequency_data Frequency data tibble
#' @param significance_threshold P-value threshold
#' @return Regression results tibble
#' @keywords internal
analyze_frequency_trends <- function(frequency_data, significance_threshold) {
  results <- frequency_data %>%
    dplyr::group_by(.data$locus, .data$allele_name, .data$population) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model = purrr::map(.data$data, ~{
        if (nrow(.x) < 3) return(NULL)
        stats::lm(frequency ~ tick, data = .x)
      }),
      slope = purrr::map_dbl(.data$model, ~{
        if (is.null(.x)) return(NA_real_)
        broom::tidy(.x) %>% 
          dplyr::filter(.data$term == "tick") %>% 
          dplyr::pull(.data$estimate)
      }),
      p_value = purrr::map_dbl(.data$model, ~{
        if (is.null(.x)) return(NA_real_)
        broom::tidy(.x) %>% 
          dplyr::filter(.data$term == "tick") %>% 
          dplyr::pull(.data$p.value)
      }),
      std_error = purrr::map_dbl(.data$model, ~{
        if (is.null(.x)) return(NA_real_)
        broom::tidy(.x) %>% 
          dplyr::filter(.data$term == "tick") %>% 
          dplyr::pull(.data$std.error)
      }),
      t_value = purrr::map_dbl(.data$model, ~{
        if (is.null(.x)) return(NA_real_)
        broom::tidy(.x) %>% 
          dplyr::filter(.data$term == "tick") %>% 
          dplyr::pull(.data$statistic)
      }),
      r_squared = purrr::map_dbl(.data$model, ~{
        if (is.null(.x)) return(NA_real_)
        broom::glance(.x) %>% 
          dplyr::pull(.data$r.squared)
      })
    ) %>%
    dplyr::select(-.data$data, -.data$model) %>%
    dplyr::filter(!is.na(.data$slope)) %>%
    dplyr::filter(.data$slope > 0, .data$p_value < significance_threshold)
  
  return(results)
}

#' Integrate regression results with mutation data
#' @param regression_results Results from regression analysis
#' @param mutation_data Original mutation data
#' @param run_number_str Run number string
#' @return Integrated results tibble
#' @keywords internal
integrate_mutation_data <- function(regression_results, mutation_data, run_number_str) {
  # Process mutation data to match genind format
  mutation_processed <- mutation_data %>%
    dplyr::filter(.data$run_number == run_number_str) %>%
    dplyr::mutate(
      locus = substr(.data$new_allele, 1, 1),
      base_allele = substr(.data$new_allele, 2, 2),
      mutation_id = sub(".*\\.", "", .data$new_allele),
      allele_name = paste0(.data$base_allele, "_", .data$mutation_id)
    ) %>%
    dplyr::select(
      .data$locus, .data$allele_name,
      new_additivity_value = .data$new_additivity_value,
      turtle_origin = .data$turtle,
      FmI = .data$FmI
    ) %>%
    dplyr::distinct()
  
  # Join with regression results
  final_results <- regression_results %>%
    dplyr::inner_join(mutation_processed, by = c("locus", "allele_name")) %>%
    dplyr::mutate(run_number = run_number_str) %>%
    dplyr::distinct() %>%
    dplyr::as_tibble()
  
  return(final_results)
}

#' Create empty selection results tibble
#' @return Empty tibble with correct column structure
#' @keywords internal
create_empty_selection_results <- function() {
  tibble::tibble(
    locus = character(0),
    allele_name = character(0),
    population = character(0),
    slope = numeric(0),
    p_value = numeric(0),
    std_error = numeric(0),
    t_value = numeric(0),
    r_squared = numeric(0),
    new_additivity_value = numeric(0),
    turtle_origin = character(0),
    FmI = numeric(0),
    run_number = character(0)
  )
}

#' Create empty migration results tibble
#' @param summarized Whether results should be summarized format
#' @return Empty tibble with correct column structure
#' @keywords internal
create_empty_migration_results <- function(summarized = FALSE) {
  if (summarized) {
    tibble::tibble(
      run_number = character(0),
      origin_pond = numeric(0),
      destination_pond = numeric(0),
      median_effective_migrants = numeric(0)
    )
  } else {
    tibble::tibble(
      run_number = character(0),
      ticks = numeric(0),
      origin_pond = numeric(0),
      destination_pond = numeric(0),
      effective_migrants = numeric(0)
    )
  }
}