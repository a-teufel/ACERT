################################################################################
# R/analysis_summary.R
################################################################################
# 
# Functions for summarizing analysis results, creating reports, and assessing
# replicate similarity in ABM simulation data
#
################################################################################

################################################################################
# FUNCTION: assess_replicate_similarity
################################################################################
#' Assess similarity of replicate simulation runs
#'
#' Analyzes the similarity between replicate runs using hierarchical clustering
#' and k-means analysis based on summary statistics. Helps identify potential
#' outlier runs or batch effects in simulation studies.
#'
#' @param turtle_data Combined turtle data from all replicates.
#' @param env_data Combined environment data from all replicates. Default: NULL.
#' @param time_points Time points to include in assessment. If NULL, uses all
#'   time points. Default: NULL.
#' @param properties Turtle properties to include in similarity calculation.
#'   Default: c("t.size.cm.metamorphosis", "t.age.metamorphosis").
#'
#' @return A list containing clustering results and distance matrix:
#'   \describe{
#'     \item{replicate_summaries}{Summary statistics per replicate}
#'     \item{distance_matrix}{Distance matrix between replicates}
#'     \item{hclust}{Hierarchical clustering results}
#'     \item{kmeans_results}{K-means clustering for different k values}
#'     \item{silhouette_scores}{Silhouette scores for cluster validation}
#'     \item{optimal_k}{Optimal number of clusters}
#'     \item{optimal_clusters}{Cluster assignments for optimal k}
#'   }
#'
#' @details The function calculates summary statistics for each replicate,
#'   standardizes the data, and performs both hierarchical and k-means clustering.
#'   Silhouette analysis is used to determine the optimal number of clusters.
#'
#' @examples
#' \dontrun{
#' # Assess similarity of turtle replicates
#' similarity_result <- assess_replicate_similarity(
#'   turtle_data,
#'   properties = c("t.size.cm", "t.age", "t.energy")
#' )
#' 
#' # View clustering results
#' print(similarity_result$optimal_k)
#' table(similarity_result$optimal_clusters)
#' 
#' # Include environmental data
#' full_assessment <- assess_replicate_similarity(
#'   turtle_data, 
#'   env_data,
#'   time_points = c(1000, 2000, 3000)
#' )
#' }
#'
#' @importFrom stats dist hclust kmeans var
#' @importFrom cluster silhouette
#' @importFrom dplyr group_by summarize left_join filter
#' @export
assess_replicate_similarity <- function(turtle_data, 
                                        env_data = NULL, 
                                        time_points = NULL,
                                        properties = c("t.size.cm.metamorphosis", 
                                                       "t.age.metamorphosis")) {
  
  # Input validation
  if (!is.data.frame(turtle_data)) {
    stop("turtle_data must be a data.frame")
  }
  
  if (!"run_number" %in% names(turtle_data)) {
    stop("turtle_data must contain a 'run_number' column")
  }
  
  # Filter time points if specified
  if (!is.null(time_points)) {
    if (!"ticks" %in% names(turtle_data)) {
      stop("time_points specified but 'ticks' column not found in turtle_data")
    }
    turtle_data <- turtle_data %>%
      dplyr::filter(ticks %in% time_points)
    
    if (!is.null(env_data) && "ticks" %in% names(env_data)) {
      env_data <- env_data %>%
        dplyr::filter(ticks %in% time_points)
    }
  }
  
  # Calculate summary statistics for each replicate
  replicate_summaries <- list()
  
  # Turtle property summaries
  for (prop in properties) {
    if (prop %in% names(turtle_data)) {
      prop_summary <- turtle_data %>%
        dplyr::group_by(run_number) %>%
        dplyr::summarize(
          !!paste0(prop, "_mean") := mean(.data[[prop]], na.rm = TRUE),
          !!paste0(prop, "_var") := stats::var(.data[[prop]], na.rm = TRUE),
          .groups = "drop"
        )
      
      if (length(replicate_summaries) == 0) {
        replicate_summaries <- prop_summary
      } else {
        replicate_summaries <- dplyr::left_join(
          replicate_summaries, 
          prop_summary, 
          by = "run_number"
        )
      }
    } else {
      warning("Property '", prop, "' not found in turtle_data")
    }
  }
  
  # Add environment summaries if provided
  if (!is.null(env_data)) {
    env_props <- c("patch_energy", "patch_risk", "permeability", "patch_regeneration")
    for (prop in env_props) {
      if (prop %in% names(env_data)) {
        env_summary <- env_data %>%
          dplyr::group_by(run_number) %>%
          dplyr::summarize(
            !!paste0(prop, "_mean") := mean(.data[[prop]], na.rm = TRUE),
            !!paste0(prop, "_var") := stats::var(.data[[prop]], na.rm = TRUE),
            .groups = "drop"
          )
        
        replicate_summaries <- dplyr::left_join(
          replicate_summaries,
          env_summary,
          by = "run_number"
        )
      }
    }
  }
  
  if (nrow(replicate_summaries) < 2) {
    stop("Need at least 2 replicates for similarity analysis")
  }
  
  # Prepare data for clustering
  clustering_data <- replicate_summaries %>%
    dplyr::select(-run_number) %>%
    as.matrix()
  
  # Remove columns with all NA or zero variance
  valid_cols <- apply(clustering_data, 2, function(x) {
    !all(is.na(x)) && stats::var(x, na.rm = TRUE) > 0
  })
  
  if (sum(valid_cols) == 0) {
    stop("No valid columns for clustering analysis")
  }
  
  clustering_data <- clustering_data[, valid_cols, drop = FALSE]
  rownames(clustering_data) <- replicate_summaries$run_number
  
  # Standardize the data
  clustering_data_scaled <- scale(clustering_data)
  
  # Handle any remaining NA or infinite values
  clustering_data_scaled[is.na(clustering_data_scaled) | is.infinite(clustering_data_scaled)] <- 0
  
  # Calculate distance matrix
  dist_matrix <- stats::dist(clustering_data_scaled)
  
  # Perform hierarchical clustering
  hclust_result <- stats::hclust(dist_matrix, method = "ward.D2")
  
  # Perform k-means for different numbers of clusters
  kmeans_results <- list()
  max_k <- min(10, nrow(clustering_data_scaled) - 1)
  
  if (max_k > 1) {
    for (k in 2:max_k) {
      tryCatch({
        kmeans_results[[paste0("k", k)]] <- stats::kmeans(clustering_data_scaled, centers = k, nstart = 25)
      }, error = function(e) {
        warning("K-means failed for k = ", k, ": ", e$message)
      })
    }
  }
  
  # Calculate silhouette scores
  silhouette_scores <- numeric(max_k - 1)
  optimal_k <- 2  # Default
  
  if (length(kmeans_results) > 0) {
    for (i in 2:max_k) {
      k_name <- paste0("k", i)
      if (k_name %in% names(kmeans_results)) {
        tryCatch({
          sil <- cluster::silhouette(kmeans_results[[k_name]]$cluster, dist_matrix)
          silhouette_scores[i-1] <- mean(sil[, 3])
        }, error = function(e) {
          silhouette_scores[i-1] <- NA
        })
      }
    }
    
    # Find optimal number of clusters
    valid_scores <- !is.na(silhouette_scores)
    if (any(valid_scores)) {
      optimal_k <- which.max(silhouette_scores[valid_scores]) + 1
    }
  }
  
  # Get optimal cluster assignments
  optimal_clusters <- NULL
  if (paste0("k", optimal_k) %in% names(kmeans_results)) {
    optimal_clusters <- kmeans_results[[paste0("k", optimal_k)]]$cluster
  }
  
  # Return results
  results <- list(
    replicate_summaries = replicate_summaries,
    distance_matrix = dist_matrix,
    hclust = hclust_result,
    kmeans_results = kmeans_results,
    silhouette_scores = silhouette_scores,
    optimal_k = optimal_k,
    optimal_clusters = optimal_clusters,
    n_replicates = nrow(replicate_summaries),
    properties_analyzed = properties
  )
  
  return(results)
}

################################################################################
# FUNCTION: plot_replicate_assessment
################################################################################
#' Plot replicate similarity assessment
#'
#' @param assessment_results Results from assess_replicate_similarity
#' @param plot_type Type of plot ("dendrogram", "pca", "silhouette")
#' @return A ggplot object
#' @importFrom ggplot2 ggplot aes geom_segment geom_text theme_minimal labs geom_point geom_line geom_vline
#' @importFrom ggdendro dendro_data 
#' @importFrom dplyr select
#' @importFrom stats prcomp as.dendrogram
#' @export
plot_replicate_assessment <- function(assessment_results, plot_type = "dendrogram") {
  if (plot_type == "dendrogram") {
    # Create dendrogram
    dend <- as.dendrogram(assessment_results$hclust)
    dend_data <- ggdendro::dendro_data(dend)
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = dend_data$segments,
                            aes(x = x, y = y, xend = xend, yend = yend)) +
      ggplot2::geom_text(data = dend_data$labels,
                         aes(x = x, y = y, label = label),
                         hjust = 1, angle = 90, size = 3) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Hierarchical Clustering of Replicate Runs",
                    x = "Replicate", y = "Distance")
    
  } else if (plot_type == "pca") {
    # Perform PCA on the summary data
    pca_data <- assessment_results$replicate_summaries %>%
      dplyr::select(-run_number) %>%
      as.matrix()
    
    pca_result <- prcomp(pca_data, scale = TRUE)
    
    # Create PCA plot
    pca_df <- data.frame(
      PC1 = pca_result$x[, 1],
      PC2 = pca_result$x[, 2],
      run_number = assessment_results$replicate_summaries$run_number,
      cluster = factor(assessment_results$optimal_clusters)
    )
    
    var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
    
    p <- ggplot2::ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_text(aes(label = run_number), hjust = 0, vjust = 0, size = 3) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "PCA of Replicate Runs",
        x = paste0("PC1 (", var_explained[1], "% variance)"),
        y = paste0("PC2 (", var_explained[2], "% variance)"),
        color = "Cluster"
      )
    
  } else if (plot_type == "silhouette") {
    # Create silhouette plot
    sil_df <- data.frame(
      k = 2:(length(assessment_results$silhouette_scores) + 1),
      silhouette = assessment_results$silhouette_scores
    )
    
    p <- ggplot2::ggplot(sil_df, aes(x = k, y = silhouette)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_vline(xintercept = assessment_results$optimal_k,
                          linetype = "dashed", color = "red") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Silhouette Scores for Different Numbers of Clusters",
        x = "Number of Clusters",
        y = "Average Silhouette Score"
      )
  }
  
  return(p)
}

################################################################################
# FUNCTION: create_abm_report
################################################################################
#' Create a comprehensive report of ABM analysis results
#'
#' Generates a complete analysis report including plots, summary statistics,
#' and data quality assessments. Optimized for large datasets with sampling
#' options for visualization.
#'
#' @param turtle_data Turtle simulation data.
#' @param env_data Environment data. Default: NULL.
#' @param mutation_data Mutation data. Default: NULL.
#' @param output_dir Directory to save report files. Default: current directory.
#' @param experiment_name Name of the experiment for file naming. Default: "ABM_Analysis".
#' @param sample_for_plots Whether to sample data for plotting to improve performance.
#'   Default: TRUE.
#' @param max_plot_points Maximum number of points to show in plots. Default: 50000.
#'
#' @return Character string path to the created report directory.
#'
#' @details Creates a comprehensive report including:
#'   \itemize{
#'     \item Summary statistics and data quality checks
#'     \item Plots of key turtle properties over time
#'     \item Population genetics analysis (if chromosome data available)
#'     \item Environmental analysis (if env_data provided)
#'     \item Mutation analysis (if mutation_data provided)
#'     \item Replicate similarity assessment
#'   }
#'
#' @examples
#' \dontrun{
#' # Create basic report
#' report_path <- create_abm_report(
#'   turtle_data,
#'   experiment_name = "Pond_Evolution_Study"
#' )
#' 
#' # Full report with all data types
#' full_report <- create_abm_report(
#'   turtle_data = turtle_data,
#'   env_data = env_data,
#'   mutation_data = mutation_data,
#'   output_dir = "final_analysis/",
#'   experiment_name = "Complete_Analysis"
#' )
#' }
#'
#' @importFrom ggplot2 ggsave
#' @importFrom utils object.size
#' @export
create_abm_report <- function(turtle_data, 
                              env_data = NULL, 
                              mutation_data = NULL,
                              output_dir = ".", 
                              experiment_name = "ABM_Analysis",
                              sample_for_plots = TRUE, 
                              max_plot_points = 50000) {
  
  # Input validation
  if (!is.data.frame(turtle_data)) {
    stop("turtle_data must be a data.frame")
  }
  
  # Create output directory structure
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
  }
  
  message("Creating ABM analysis report...")
  message("Output directory: ", output_dir)
  
  # Sample data for plotting if requested
  if (sample_for_plots && nrow(turtle_data) > max_plot_points) {
    message("Sampling data for plotting (", max_plot_points, " points)...")
    plot_data <- turtle_data[sample(nrow(turtle_data), max_plot_points), ]
  } else {
    plot_data <- turtle_data
  }
  
  # Generate summary statistics
  message("Generating summary statistics...")
  data_summary <- summarize_simulation_data(turtle_data, env_data, mutation_data)
  
  # Create plots for key turtle properties
  message("Creating turtle property plots...")
  
  # Check for available columns with both naming conventions
  size_col <- if("t.size.cm.metamorphosis" %in% names(plot_data)) {
    "t.size.cm.metamorphosis" 
  } else if("t_size_cm_metamorphosis" %in% names(plot_data)) {
    "t_size_cm_metamorphosis"
  } else {
    NULL
  }
  
  age_col <- if("t.age.metamorphosis" %in% names(plot_data)) {
    "t.age.metamorphosis"
  } else if("t_age_metamorphosis" %in% names(plot_data)) {
    "t_age_metamorphosis"
  } else {
    NULL
  }
  
  # Create plots for available properties
  if (!is.null(size_col) && "t.birth.pond.index" %in% names(plot_data)) {
    tryCatch({
      size_data <- analyze_turtle_property(plot_data, size_col, "t.birth.pond.index")
      size_plot <- plot_turtle_property(size_data, "Size at Metamorphosis", "t.birth.pond.index", experiment_name)
      ggplot2::ggsave(file.path(plots_dir, "size_at_metamorphosis.png"), size_plot, width = 10, height = 6)
      message("  Created size at metamorphosis plot")
    }, error = function(e) {
      warning("Failed to create size plot: ", e$message)
    })
  }
  
  if (!is.null(age_col) && "t.birth.pond.index" %in% names(plot_data)) {
    tryCatch({
      age_data <- analyze_turtle_property(plot_data, age_col, "t.birth.pond.index")
      age_plot <- plot_turtle_property(age_data, "Age at Metamorphosis", "t.birth.pond.index", experiment_name)
      ggplot2::ggsave(file.path(plots_dir, "age_at_metamorphosis.png"), age_plot, width = 10, height = 6)
      message("  Created age at metamorphosis plot")
    }, error = function(e) {
      warning("Failed to create age plot: ", e$message)
    })
  }
  
  # Population genetics analysis if chromosome data available
  chr_cols <- grep("chromosome.*\\.[a-z]$", names(turtle_data), value = TRUE)
  if (length(chr_cols) > 0) {
    message("Processing population genetics...")
    
    # Sample for population genetics if needed
    if (nrow(turtle_data) > 100000) {
      message("  Sampling for population genetics analysis...")
      popgen_data <- create_analysis_subset(turtle_data, "random", 100000)
    } else {
      popgen_data <- turtle_data
    }
    
    tryCatch({
      # This function may not exist in all versions, so wrap in tryCatch
      popgen_results <- process_popgen(popgen_data, max_individuals = 1000)
      
      if (!is.null(popgen_results$combined_fst_data)) {
        # Create Fst plot
        avg_fst <- popgen_results$combined_fst_data %>%
          dplyr::group_by(Tick, Locus) %>%
          dplyr::summarize(
            Mean_Fst = mean(Fst, na.rm = TRUE),
            SE_Fst = sd(Fst, na.rm = TRUE) / sqrt(dplyr::n()),
            .groups = "drop"
          )
        
        if (nrow(avg_fst) > 0) {
          fst_plot <- ggplot2::ggplot(avg_fst, ggplot2::aes(x = Tick, y = Mean_Fst, color = Locus)) +
            ggplot2::geom_line(size = 1) +
            ggplot2::geom_point(size = 2) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean_Fst - SE_Fst, ymax = Mean_Fst + SE_Fst), width = 0.1) +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = paste("Mean Fst Over Time in", experiment_name),
                          x = "Time (ticks)", y = "Mean Fst +/- SE")
          
          ggplot2::ggsave(file.path(plots_dir, "fst_over_time.png"), fst_plot, width = 10, height = 6)
          message("  Created Fst plot")
        }
      }
    }, error = function(e) {
      message("  Population genetics analysis failed: ", e$message)
    })
  }
  
  # Environment analysis if available
  if (!is.null(env_data)) {
    message("Creating environment analysis...")
    
    tryCatch({
      # Get representative time points
      time_points <- unique(env_data$ticks)
      if (length(time_points) > 6) {
        time_points <- time_points[seq(1, length(time_points), length.out = 6)]
      }
      
      # Create environment visualizations if functions are available
      env_props <- c("patch_energy", "patch_risk", "permeability", "patch_regeneration")
      for (prop in env_props) {
        if (prop %in% names(env_data)) {
          tryCatch({
            # This assumes visualize_patches function exists
            p <- visualize_patches(env_data, prop, time_points = time_points)
            ggplot2::ggsave(file.path(plots_dir, paste0(prop, ".png")), p, width = 12, height = 8)
            message("  Created ", prop, " visualization")
          }, error = function(e) {
            message("  Failed to create ", prop, " plot: ", e$message)
          })
        }
      }
    }, error = function(e) {
      message("Environment analysis failed: ", e$message)
    })
  }
  
  # Mutation analysis if available
  if (!is.null(mutation_data)) {
    message("Creating mutation analysis...")
    
    tryCatch({
      # Sample mutations if too many
      mutation_plot_data <- if (nrow(mutation_data) > 100000) {
        mutation_data[sample(nrow(mutation_data), 100000), ]
      } else {
        mutation_data
      }
      
      if ("FmI" %in% names(mutation_plot_data) && "ticks" %in% names(mutation_plot_data)) {
        mut_plot <- ggplot2::ggplot(mutation_plot_data, ggplot2::aes(x = factor(ticks), y = FmI)) +
          ggplot2::geom_boxplot(ggplot2::aes(fill = factor(ticks)), show.legend = FALSE, outlier.shape = NA) +
          ggplot2::geom_jitter(height = 0, width = 0.1, alpha = 0.05, size = 0.5) +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = paste("Mutation Effect Sizes Over Time in", experiment_name),
                        x = "Time (ticks)", y = "Mutation Effect Size (FmI)")
        
        ggplot2::ggsave(file.path(plots_dir, "mutation_effects.png"), mut_plot, width = 10, height = 6)
        message("  Created mutation effects plot")
      }
    }, error = function(e) {
      message("Mutation analysis failed: ", e$message)
    })
  }
  
  # Assess replicate similarity
  message("Assessing replicate similarity...")
  tryCatch({
    # Sample for replicate assessment
    assess_data <- if (nrow(turtle_data) > 50000) {
      create_analysis_subset(turtle_data, "random", 50000)
    } else {
      turtle_data
    }
    
    assessment <- assess_replicate_similarity(assess_data, env_data)
    
    # Save assessment plots would go here (requires plotting functions)
    message("  Replicate similarity assessment completed")
    
  }, error = function(e) {
    message("Replicate similarity assessment failed: ", e$message)
    assessment <- NULL
  })
  
  # Save summary report
  message("Saving summary report...")
  summary_file <- file.path(output_dir, paste0(experiment_name, "_summary.txt"))
  
  tryCatch({
    sink(summary_file)
    cat("ABM Analysis Summary Report\n")
    cat("===========================\n\n")
    cat("Experiment:", experiment_name, "\n")
    cat("Generated:", Sys.time(), "\n\n")
    
    # Data summary
    cat("Data Summary:\n")
    cat("-------------\n")
    cat("Total turtle records:", format(nrow(turtle_data), big.mark = ","), "\n")
    
    if ("run_number" %in% names(turtle_data)) {
      cat("Number of replicates:", length(unique(turtle_data$run_number)), "\n")
    }
    
    if ("ticks" %in% names(turtle_data)) {
      time_range <- range(turtle_data$ticks, na.rm = TRUE)
      cat("Time range:", paste(time_range, collapse = " - "), "ticks\n")
    }
    
    # Chromosome structure
    cat("Chromosome columns:", length(chr_cols), "expanded chromosome columns found\n")
    if (length(chr_cols) > 0) {
      loci <- unique(sub(".*\\.", "", chr_cols))
      cat("Loci detected:", paste(loci, collapse = ", "), "\n")
    }
    
    if (!is.null(env_data)) {
      cat("Environment records:", format(nrow(env_data), big.mark = ","), "\n")
    }
    
    if (!is.null(mutation_data)) {
      cat("Mutation records:", format(nrow(mutation_data), big.mark = ","), "\n")
    }
    
    # Replicate clustering results
    if (!is.null(assessment)) {
      cat("\nReplicate Analysis:\n")
      cat("-------------------\n")
      cat("Optimal number of clusters:", assessment$optimal_k, "\n")
      if (!is.null(assessment$optimal_clusters)) {
        cluster_dist <- table(assessment$optimal_clusters)
        cat("Cluster distribution:", paste(names(cluster_dist), "=", cluster_dist, collapse = ", "), "\n")
      }
    }
    
    cat("\nFiles created:\n")
    cat("- Plots directory:", plots_dir, "\n")
    cat("- Summary file:", summary_file, "\n")
    
    sink()
    
  }, error = function(e) {
    if (sink.number() > 0) sink()  # Ensure sink is closed
    warning("Failed to create summary file: ", e$message)
  })
  
  message("Report creation completed!")
  message("Report saved to: ", output_dir)
  
  return(output_dir)
}

################################################################################
# FUNCTION: summarize_simulation_data
################################################################################
#' Create summary statistics for simulation data
#'
#' Generates comprehensive summary statistics for turtle, environment, and
#' mutation data to provide an overview of simulation results.
#'
#' @param turtle_data Turtle simulation data.
#' @param env_data Environment data. Default: NULL.
#' @param mutation_data Mutation data. Default: NULL.
#'
#' @return A list containing summary statistics for each data type.
#'
#' @examples
#' \dontrun{
#' # Create data summary
#' summary_stats <- summarize_simulation_data(
#'   turtle_data,
#'   env_data,
#'   mutation_data
#' )
#' 
#' # View turtle summary
#' summary_stats$turtle_summary
#' }
#'
#' @importFrom dplyr summarize across where
#' @export
summarize_simulation_data <- function(turtle_data, env_data = NULL, mutation_data = NULL) {
  
  summary_list <- list()
  
  # Turtle data summary
  if (!is.null(turtle_data)) {
    turtle_summary <- list(
      n_records = nrow(turtle_data),
      n_individuals = length(unique(turtle_data$who)),
      n_replicates = if("run_number" %in% names(turtle_data)) length(unique(turtle_data$run_number)) else 1,
      time_range = if("ticks" %in% names(turtle_data)) range(turtle_data$ticks, na.rm = TRUE) else c(NA, NA),
      breeds = if("breed" %in% names(turtle_data)) table(turtle_data$breed) else "Not available",
      memory_size_mb = round(as.numeric(utils::object.size(turtle_data)) / (1024^2), 2)
    )
    
    # Numeric column summaries
    numeric_cols <- names(turtle_data)[sapply(turtle_data, is.numeric)]
    numeric_cols <- setdiff(numeric_cols, c("who", "ticks", "xcor", "ycor"))
    
    if (length(numeric_cols) > 0) {
      turtle_summary$numeric_summaries <- turtle_data %>%
        dplyr::summarize(dplyr::across(all_of(numeric_cols), 
                                       list(mean = ~mean(.x, na.rm = TRUE),
                                            sd = ~sd(.x, na.rm = TRUE),
                                            min = ~min(.x, na.rm = TRUE),
                                            max = ~max(.x, na.rm = TRUE)),
                                       .names = "{.col}_{.fn}"))
    }
    
    summary_list$turtle_summary <- turtle_summary
  }
  
  # Environment data summary
  if (!is.null(env_data)) {
    env_summary <- list(
      n_records = nrow(env_data),
      n_replicates = if("run_number" %in% names(env_data)) length(unique(env_data$run_number)) else 1,
      time_range = if("ticks" %in% names(env_data)) range(env_data$ticks, na.rm = TRUE) else c(NA, NA),
      spatial_extent = if(all(c("pxcor", "pycor") %in% names(env_data))) {
        list(x_range = range(env_data$pxcor, na.rm = TRUE),
             y_range = range(env_data$pycor, na.rm = TRUE))
      } else {
        "Coordinates not available"
      },
      memory_size_mb = round(as.numeric(utils::object.size(env_data)) / (1024^2), 2)
    )
    
    summary_list$env_summary <- env_summary
  }
  
  # Mutation data summary
  if (!is.null(mutation_data)) {
    mut_summary <- list(
      n_mutations = nrow(mutation_data),
      n_replicates = if("run_number" %in% names(mutation_data)) length(unique(mutation_data$run_number)) else 1,
      time_range = if("ticks" %in% names(mutation_data)) range(mutation_data$ticks, na.rm = TRUE) else c(NA, NA),
      memory_size_mb = round(as.numeric(utils::object.size(mutation_data)) / (1024^2), 2)
    )
    
    # Effect size summary if available
    if ("FmI" %in% names(mutation_data)) {
      mut_summary$effect_size_summary <- list(
        mean_effect = mean(mutation_data$FmI, na.rm = TRUE),
        sd_effect = sd(mutation_data$FmI, na.rm = TRUE),
        min_effect = min(mutation_data$FmI, na.rm = TRUE),
        max_effect = max(mutation_data$FmI, na.rm = TRUE),
        prop_beneficial = mean(mutation_data$FmI > 0, na.rm = TRUE),
        prop_deleterious = mean(mutation_data$FmI < 0, na.rm = TRUE),
        prop_neutral = mean(mutation_data$FmI == 0, na.rm = TRUE)
      )
    }
    
    summary_list$mutation_summary <- mut_summary
  }
  
  return(summary_list)
}