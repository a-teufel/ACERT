################################################################################
# R/population_genetics.R
################################################################################
# 
# Functions for population genetic analysis of ABM simulation data
#
################################################################################

################################################################################
# FUNCTION: calculate_popgen_stats
################################################################################
#' Calculate population genetic summary statistics
#'
#' Computes comprehensive population genetic statistics including F-statistics,
#' Gst variants, and Jost's D with optional bootstrap confidence intervals.
#'
#' @param genind_obj A genind object from the adegenet package.
#' @param pop_column Numeric value indicating the column containing population
#'   identifiers for Pegas F-statistics. Default: NULL.
#' @param n_bootstrap Number of bootstrap replicates for confidence intervals.
#'   Default: 1000.
#' @param use_bootstrap Whether to calculate bootstrap confidence intervals.
#'   Default: TRUE.
#'
#' @return A list containing differentiation measures:
#'   \describe{
#'     \item{Fstats}{F-statistics from Pegas}
#'     \item{Gst}{Nei's Gst with optional confidence intervals}
#'     \item{Gst_Hedrick}{Hedrick's standardized Gst with optional CIs}
#'     \item{Jost_D}{Jost's D with optional confidence intervals}
#'   }
#'
#' @details Bootstrap confidence intervals are calculated using the mmod package's
#'   chao_bootstrap function when use_bootstrap = TRUE. Each differentiation
#'   measure includes both point estimates and confidence intervals when available.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(nancycats)
#' 
#' # Calculate all statistics with bootstrap CIs
#' popgen_stats <- calculate_popgen_stats(nancycats, n_bootstrap = 100)
#' 
#' # View Gst results
#' popgen_stats$Gst$Gst_Est
#' popgen_stats$Gst$Gst_CI
#' 
#' # Calculate without bootstrap (faster)
#' quick_stats <- calculate_popgen_stats(nancycats, use_bootstrap = FALSE)
#' }
#'
#' @importFrom pegas genind2loci Fst
#' @importFrom mmod chao_bootstrap Gst_Nei Gst_Hedrick D_Jost summarise_bootstrap
#' @export
calculate_popgen_stats <- function(genind_obj, 
                                   pop_column = NULL, 
                                   n_bootstrap = 1000, 
                                   use_bootstrap = TRUE) {
  
  # Input validation
  if (!inherits(genind_obj, "genind")) {
    stop("genind_obj must be a genind object from the adegenet package")
  }
  
  if (use_bootstrap && n_bootstrap < 10) {
    warning("Very low number of bootstrap replicates may give unreliable confidence intervals")
  }
  
  # Convert to pegas format for F-statistics
  pegas_obj <- pegas::genind2loci(genind_obj)
  
  # Calculate bootstrap confidence intervals if requested
  if (use_bootstrap) {
    tryCatch({
      bootstrap_results <- mmod::chao_bootstrap(genind_obj, nreps = n_bootstrap)
      
      gst_ci <- mmod::summarise_bootstrap(bootstrap_results, mmod::Gst_Nei)
      gst_hedrick_ci <- mmod::summarise_bootstrap(bootstrap_results, mmod::Gst_Hedrick)
      jost_d_ci <- mmod::summarise_bootstrap(bootstrap_results, mmod::D_Jost)
    }, error = function(e) {
      warning("Bootstrap failed: ", e$message, ". Proceeding without confidence intervals.")
      gst_ci <- NULL
      gst_hedrick_ci <- NULL
      jost_d_ci <- NULL
    })
  } else {
    gst_ci <- NULL
    gst_hedrick_ci <- NULL
    jost_d_ci <- NULL
  }
  
  # Calculate point estimates
  fstats <- pegas::Fst(pegas_obj, pop = pop_column)
  gst_est <- mmod::Gst_Nei(genind_obj)
  gst_hedrick_est <- mmod::Gst_Hedrick(genind_obj)
  jost_d_est <- mmod::D_Jost(genind_obj)
  
  # Compile results
  results <- list(
    Fstats = fstats,
    Gst = list(
      Gst_Est = gst_est,
      Gst_CI = gst_ci
    ),
    Gst_Hed = list(
      Gst_H_Est = gst_hedrick_est,
      Gst_H_CI = gst_hedrick_ci
    ),
    Jost_D = list(
      Jost_D_Est = jost_d_est,
      Jost_D_CI = jost_d_ci
    )
  )
  
  return(results)
}

################################################################################
# FUNCTION: calculate_pairwise_differentiation
################################################################################
#' Calculate pairwise differentiation measures between populations
#'
#' Computes multiple pairwise differentiation statistics including Fst, Gst,
#' Hedrick's Gst, and Jost's D for all population pairs.
#'
#' @param genind_obj A genind object with population assignments.
#' @param linearized Whether to linearize differentiation measures for use as
#'   distances. Default: FALSE.
#'
#' @return A list containing pairwise distance matrices:
#'   \describe{
#'     \item{Fst}{Weir & Cockerham's Fst}
#'     \item{Gst}{Nei's Gst}
#'     \item{Gst_H}{Hedrick's standardized Gst}
#'     \item{Jost_D}{Jost's D}
#'   }
#'
#' @details All measures are returned as dist objects suitable for clustering
#'   or other distance-based analyses. Linearized versions are appropriate
#'   for use as genetic distances in phylogenetic or UPGMA analyses.
#'
#' @examples
#' \dontrun{
#' # Calculate pairwise differentiation
#' pairwise_diff <- calculate_pairwise_differentiation(nancycats)
#' 
#' # Extract Fst matrix
#' fst_matrix <- as.matrix(pairwise_diff$Fst)
#' 
#' # Use for clustering
#' fst_tree <- hclust(pairwise_diff$Fst)
#' plot(fst_tree)
#' }
#'
#' @importFrom hierfstat pairwise.WCfst genind2hierfstat
#' @importFrom mmod pairwise_Gst_Nei pairwise_Gst_Hedrick pairwise_D
#' @importFrom stats as.dist
#' @export
calculate_pairwise_differentiation <- function(genind_obj, linearized = FALSE) {
  
  # Input validation
  if (!inherits(genind_obj, "genind")) {
    stop("genind_obj must be a genind object")
  }
  
  if (is.null(genind_obj$pop)) {
    stop("genind_obj must have population assignments (genind_obj$pop)")
  }
  
  n_pops <- length(levels(genind_obj$pop))
  if (n_pops < 2) {
    stop("Need at least 2 populations for pairwise analysis")
  }
  
  # Convert to hierfstat format for Fst calculation
  hierfstat_data <- hierfstat::genind2hierfstat(genind_obj)
  
  # Calculate pairwise measures
  pairwise_results <- list(
    Fst = stats::as.dist(hierfstat::pairwise.WCfst(hierfstat_data)),
    Gst = mmod::pairwise_Gst_Nei(genind_obj, linearized = linearized),
    Gst_H = mmod::pairwise_Gst_Hedrick(genind_obj, linearized = linearized),
    Jost_D = mmod::pairwise_D(genind_obj, linearized = linearized)
  )
  
  return(pairwise_results)
}

################################################################################
# FUNCTION: calculate_genetic_diversity
################################################################################
#' Calculate genetic diversity indices for populations
#'
#' Computes comprehensive genetic diversity measures including allelic richness,
#' heterozygosity, and F-statistics using hierfstat's basic.stats function
#' with restored original names.
#'
#' @param genind_obj A genind object with population assignments.
#' @param restore_names Whether to restore original locus and allele names
#'   from the genind object. Default: TRUE.
#'
#' @return A list containing diversity measures:
#'   \describe{
#'     \item{Num_Alleles}{Number of observed alleles per locus per population}
#'     \item{Num_Eff_Alleles}{Effective number of alleles (1/He)}
#'     \item{Obs_Het}{Observed heterozygosity by locus and population}
#'     \item{Exp_Het}{Expected heterozygosity by locus and population}
#'     \item{Allele_Freqs}{Allele frequencies by population}
#'     \item{Fis}{Inbreeding coefficient by locus and population}
#'   }
#'
#' @examples
#' \dontrun{
#' # Calculate diversity indices
#' diversity <- calculate_genetic_diversity(nancycats)
#' 
#' # View observed heterozygosity
#' diversity$Obs_Het
#' 
#' # Calculate mean diversity across loci
#' mean_ho <- rowMeans(diversity$Obs_Het, na.rm = TRUE)
#' mean_he <- rowMeans(diversity$Exp_Het, na.rm = TRUE)
#' }
#'
#' @importFrom hierfstat basic.stats
#' @export
calculate_genetic_diversity <- function(genind_obj, restore_names = TRUE) {
  
  # Input validation
  if (!inherits(genind_obj, "genind")) {
    stop("genind_obj must be a genind object")
  }
  
  if (is.null(genind_obj$pop)) {
    stop("genind_obj must have population assignments")
  }
  
  # Calculate basic statistics using hierfstat
  basic_stats_result <- hierfstat::basic.stats(genind_obj)
  
  # Restore original names if requested
  if (restore_names) {
    basic_stats_result <- restore_basic_stats_names(genind_obj, basic_stats_result)
  }
  
  # Calculate derived measures
  num_alleles <- get_num_observed_alleles_per_locus(basic_stats_result$pop.freq)
  effective_alleles_result <- get_effective_alleles_from_freq(basic_stats_result$pop.freq)
  
  # Compile diversity indices
  diversity_indices <- list(
    Num_Alleles = num_alleles,
    Num_Eff_Alleles = effective_alleles_result$Ae,
    Obs_Het = basic_stats_result$Ho,
    Exp_Het = effective_alleles_result$He,
    Allele_Freqs = basic_stats_result$pop.freq,
    Fis = basic_stats_result$Fis
  )
  
  return(diversity_indices)
}

################################################################################
# FUNCTION: calculate_allele_frequency_spectrum
################################################################################
#' Calculate allele frequency spectra for populations
#'
#' Computes frequency distributions of alleles across specified frequency bins,
#' useful for detecting demographic changes and selection signatures.
#'
#' @param genind_obj A genind object with population assignments.
#' @param by_population Whether to calculate spectra separately for each
#'   population. Default: FALSE (global spectrum).
#' @param freq_bins Numeric vector defining frequency bin boundaries.
#'   Default: seq(0, 1, by = 0.05).
#' @param min_frequency Minimum allele frequency to include. Default: 0.
#'
#' @return If by_population = FALSE, returns a table of allele counts per
#'   frequency bin. If by_population = TRUE, returns a list with one spectrum
#'   per population.
#'
#' @examples
#' \dontrun{
#' # Global allele frequency spectrum
#' global_afs <- calculate_allele_frequency_spectrum(nancycats)
#' barplot(global_afs, main = "Global Allele Frequency Spectrum")
#' 
#' # Population-specific spectra
#' pop_afs <- calculate_allele_frequency_spectrum(nancycats, by_population = TRUE)
#' 
#' # Custom frequency bins
#' rare_afs <- calculate_allele_frequency_spectrum(
#'   nancycats, 
#'   freq_bins = c(0, 0.01, 0.05, 0.1, 0.5, 1.0)
#' )
#' }
#'
#' @importFrom adegenet makefreq seppop
#' @export
calculate_allele_frequency_spectrum <- function(genind_obj, 
                                                by_population = FALSE,
                                                freq_bins = seq(0, 1, by = 0.05),
                                                min_frequency = 0) {
  
  # Input validation
  if (!inherits(genind_obj, "genind")) {
    stop("genind_obj must be a genind object")
  }
  
  if (min(freq_bins) < 0 || max(freq_bins) > 1) {
    stop("freq_bins must be between 0 and 1")
  }
  
  if (min_frequency < 0 || min_frequency > 1) {
    stop("min_frequency must be between 0 and 1")
  }
  
  if (by_population) {
    if (is.null(genind_obj$pop)) {
      stop("Population assignments required for by_population = TRUE")
    }
    
    # Calculate spectrum for each population
    pop_list <- adegenet::seppop(genind_obj)
    pop_names <- names(pop_list)
    spectra <- list()
    
    for (pop_name in pop_names) {
      pop_obj <- pop_list[[pop_name]]
      
      # Get allele frequencies for this population
      allele_freqs <- adegenet::makefreq(pop_obj)
      mean_freqs <- colMeans(allele_freqs, na.rm = TRUE)
      
      # Filter by minimum frequency
      filtered_freqs <- mean_freqs[mean_freqs >= min_frequency]
      
      # Create frequency bins
      if (length(filtered_freqs) > 0) {
        freq_categories <- cut(filtered_freqs, 
                               breaks = freq_bins,
                               include.lowest = TRUE,
                               right = FALSE)
        spectra[[pop_name]] <- table(freq_categories)
      } else {
        # Empty spectrum
        empty_bins <- cut(numeric(0), breaks = freq_bins, include.lowest = TRUE)
        spectra[[pop_name]] <- table(empty_bins)
      }
    }
    
    return(spectra)
    
  } else {
    # Global spectrum across all populations
    allele_freqs <- adegenet::makefreq(genind_obj)
    mean_freqs <- colMeans(allele_freqs, na.rm = TRUE)
    
    # Filter by minimum frequency
    filtered_freqs <- mean_freqs[mean_freqs >= min_frequency]
    
    # Create frequency bins
    if (length(filtered_freqs) > 0) {
      freq_categories <- cut(filtered_freqs,
                             breaks = freq_bins,
                             include.lowest = TRUE,
                             right = FALSE)
      return(table(freq_categories))
    } else {
      # Return empty spectrum
      empty_bins <- cut(numeric(0), breaks = freq_bins, include.lowest = TRUE)
      return(table(empty_bins))
    }
  }
}

################################################################################
# FUNCTION: perform_amova
################################################################################
#' Perform Analysis of Molecular Variance (AMOVA)
#'
#' Conducts hierarchical analysis of genetic variance using population
#' structure to partition total genetic variation.
#'
#' @param genind_obj A genind object with population assignments.
#' @param hierarchical_structure Data.frame defining hierarchical population
#'   structure with columns for different levels (e.g., region, population).
#'   Default: NULL (uses only population level).
#' @param distance_method Distance method for genetic distances. Options include
#'   "euclidean", "manhattan", "canberra". Default: "euclidean".
#'
#' @return AMOVA results object from poppr package containing variance
#'   components, phi-statistics, and significance tests.
#'
#' @examples
#' \dontrun{
#' # Simple AMOVA with population structure only
#' amova_result <- perform_amova(nancycats)
#' amova_result
#' 
#' # Hierarchical AMOVA with region and population structure
#' strata_df <- data.frame(
#'   region = rep(c("North", "South"), each = 9),
#'   population = levels(nancycats$pop)
#' )
#' 
#' hierarchical_amova <- perform_amova(nancycats, strata_df)
#' }
#'
#' @importFrom poppr as.genclone poppr.amova
#' @importFrom adegenet strata
#' @export
perform_amova <- function(genind_obj, 
                          hierarchical_structure = NULL,
                          distance_method = "euclidean") {
  
  # Input validation
  if (!inherits(genind_obj, "genind")) {
    stop("genind_obj must be a genind object")
  }
  
  if (is.null(genind_obj$pop)) {
    stop("genind_obj must have population assignments")
  }
  
  # Convert to genclone object required by poppr
  genclone_obj <- poppr::as.genclone(genind_obj)
  
  # Set hierarchical structure if provided
  if (!is.null(hierarchical_structure)) {
    if (!is.data.frame(hierarchical_structure)) {
      stop("hierarchical_structure must be a data.frame")
    }
    
    if (nrow(hierarchical_structure) != length(levels(genind_obj$pop))) {
      stop("hierarchical_structure must have one row per population")
    }
    
    # Expand to individual level
    individual_strata <- hierarchical_structure[as.numeric(genind_obj$pop), , drop = FALSE]
    adegenet::strata(genclone_obj) <- individual_strata
    
    # Create formula based on column names
    formula_terms <- paste(names(hierarchical_structure), collapse = "/")
    amova_formula <- as.formula(paste("~", formula_terms))
  } else {
    # Simple population-level AMOVA
    amova_formula <- ~ pop
  }
  
  # Perform AMOVA
  amova_result <- poppr::poppr.amova(genclone_obj, 
                                     amova_formula,
                                     method = distance_method)
  
  return(amova_result)
}

################################################################################
# FUNCTION: perform_dapc
################################################################################
#' Perform Discriminant Analysis of Principal Components (DAPC)
#'
#' Conducts DAPC analysis to identify genetic clusters and assess population
#' differentiation with dimensionality reduction.
#'
#' @param genind_obj A genind object with or without population assignments.
#' @param n_pca Number of principal components to retain. If NULL, uses
#'   optimization to determine optimal number. Default: NULL.
#' @param n_da Number of discriminant functions to retain. If NULL, uses
#'   number of populations minus 1. Default: NULL.
#' @param find_clusters Whether to identify genetic clusters if no population
#'   assignments exist. Default: FALSE.
#' @param max_clusters Maximum number of clusters to test when find_clusters = TRUE.
#'   Default: 20.
#'
#' @return DAPC results object containing:
#'   \describe{
#'     \item{eigenvalues}{Eigenvalues of discriminant functions}
#'     \item{ind.coord}{Individual coordinates on discriminant axes}
#'     \item{grp.coord}{Group centroids on discriminant axes}
#'     \item{posterior}{Posterior assignment probabilities}
#'   }
#'
#' @examples
#' \dontrun{
#' # DAPC with existing population assignments
#' dapc_result <- perform_dapc(nancycats, n_pca = 50, n_da = 2)
#' 
#' # Plot results
#' scatter(dapc_result)
#' 
#' # DAPC with cluster identification
#' dapc_clusters <- perform_dapc(nancycats, find_clusters = TRUE)
#' }
#'
#' @importFrom adegenet dapc find.clusters optim.a.score
#' @export
perform_dapc <- function(genind_obj, 
                         n_pca = NULL, 
                         n_da = NULL,
                         find_clusters = FALSE,
                         max_clusters = 20) {
  
  # Input validation
  if (!inherits(genind_obj, "genind")) {
    stop("genind_obj must be a genind object")
  }
  
  # Determine population assignments
  if (is.null(genind_obj$pop) || find_clusters) {
    if (find_clusters) {
      message("Identifying genetic clusters...")
      cluster_result <- adegenet::find.clusters(genind_obj, 
                                                max.n.clust = max_clusters,
                                                n.pca = 100,
                                                choose.n.clust = FALSE,
                                                criterion = "diffNgroup")
      population_assignments <- cluster_result$grp
      optimal_k <- length(levels(cluster_result$grp))
      message("Optimal number of clusters: ", optimal_k)
    } else {
      stop("No population assignments found. Set find_clusters = TRUE or provide population assignments.")
    }
  } else {
    population_assignments <- genind_obj$pop
  }
  
  # Optimize number of PCA components if not specified
  if (is.null(n_pca)) {
    message("Optimizing number of PCA components...")
    # Initial DAPC to get optimization data
    temp_dapc <- adegenet::dapc(genind_obj, 
                                population_assignments, 
                                n.pca = min(100, nrow(genind_obj@tab) - 1), 
                                n.da = 2)
    
    optim_result <- adegenet::optim.a.score(temp_dapc)
    n_pca <- optim_result$best
    message("Optimal number of PCA components: ", n_pca)
  }
  
  # Set number of discriminant axes if not specified
  if (is.null(n_da)) {
    n_da <- min(length(unique(population_assignments)) - 1, n_pca)
  }
  
  # Perform final DAPC
  dapc_result <- adegenet::dapc(genind_obj, 
                                population_assignments,
                                n.pca = n_pca, 
                                n.da = n_da)
  
  # Add metadata to results
  dapc_result$call$n_pca_used <- n_pca
  dapc_result$call$n_da_used <- n_da
  dapc_result$call$find_clusters <- find_clusters
  
  return(dapc_result)
}


################################################################################
# FUNCTION: process_popgen
################################################################################
#' Process population genetics data with memory optimization
#'
#' Converts turtle data to genind objects and calculates population genetic
#' statistics including Fst values across time points and replicates.
#'
#' @param all_turtle_data List of turtle data frames or single data.table
#' @param loci_names Names of loci to analyze. Default: letters[1:10].
#' @param max_individuals Maximum individuals per population for memory management.
#'   Default: NULL.
#' @param progress Show progress messages. Default: TRUE.
#'
#' @return List containing Fst data and genind objects:
#'   \describe{
#'     \item{combined_fst_data}{Data frame with Fst values by replicate, time, and locus}
#'     \item{genind_data}{List of genind objects for further analysis}
#'   }
#'
#' @examples
#' \dontrun{
#' # Process population genetics for multiple replicates
#' popgen_results <- process_popgen(turtle_data_list, max_individuals = 1000)
#' 
#' # View Fst results
#' head(popgen_results$combined_fst_data)
#' }
#'
#' @importFrom adegenet df2genind
#' @importFrom hierfstat basic.stats
#' @importFrom stats setNames
#' @export
process_popgen <- function(all_turtle_data, 
                           loci_names = letters[1:10],
                           max_individuals = NULL,
                           progress = TRUE) {
  
  # Handle both list and data.table input
  if (is.data.frame(all_turtle_data)) {
    all_turtle_data <- split(all_turtle_data, all_turtle_data$run_number)
  }
  
  combined_fst_data <- data.frame()
  genind_data <- list()
  
  for (i in seq_along(all_turtle_data)) {
    replicate <- all_turtle_data[[i]]
    replicate_num <- unique(replicate[['run_number']])
    
    if (progress) {
      message("Processing population genetics for replicate ", replicate_num)
    }
    
    # Sample if too many individuals
    if (!is.null(max_individuals) && nrow(replicate) > max_individuals) {
      replicate <- replicate[sample(nrow(replicate), max_individuals), ]
    }
    
    # Convert to pegas format - using the convert_to_pegas function
    data_pegas <- convert_to_pegas(replicate, loci_names = loci_names)
    
    # Process each time point
    ticks <- unique(data_pegas[['ticks']])
    for (tick in ticks) {
      data_for_tick <- data_pegas[data_pegas[['ticks']] == tick, ]
      
      # Skip if too few individuals
      if (nrow(data_for_tick) < 10) next
      
      # Get loci columns (all except pop and ticks)
      loci_cols <- setdiff(names(data_for_tick), c("pop", "ticks"))
      loci_data <- data_for_tick[, loci_cols]
      population_info <- data_for_tick$pop
      
      # Create genind object
      tryCatch({
        data_genind <- adegenet::df2genind(
          loci_data, 
          sep = "/",
          pop = as.factor(population_info), 
          ploidy = 2
        )
        
        # Store genind object
        genind_data[[paste(replicate_num, as.character(tick))]] <- data_genind
        
        # Calculate Fst
        fst_stats <- hierfstat::basic.stats(data_genind)
        fst_values <- fst_stats[["perloc"]][["Fst"]]
        
        # Combine results
        combined_fst_data <- rbind(
          combined_fst_data,
          data.frame(
            File = paste(replicate_num, as.character(tick)),
            Tick = tick,
            Locus = loci_names,
            Fst = as.numeric(fst_values),
            stringsAsFactors = FALSE
          )
        )
      }, error = function(e) {
        if (progress) {
          warning("Error processing ", replicate_num, " at time ", tick, ": ", e$message)
        }
      })
    }
    if(i %% 5 == 0) {
      gc(verbose = FALSE)
      message(sprintf("Processed %d/%d replicates", i, length(all_turtle_data)))
    }
  }
  
  return(list(
    combined_fst_data = combined_fst_data,
    genind_data = genind_data
  ))
}