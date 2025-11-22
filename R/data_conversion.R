################################################################################
# R/data_conversion.R
################################################################################
# 
# Functions for converting genetic data between different analysis software formats
#
################################################################################

################################################################################
# FUNCTION: genind_to_geneland
################################################################################
#' Convert genind object to Geneland input files
#'
#' Converts a genind object and associated coordinates to the input format
#' required by Geneland software for spatial population genetics analysis.
#' Creates both genotype (G.txt) and coordinate (XY.txt) files.
#'
#' @param genind_obj An object of class 'genind' from the adegenet package.
#' @param coords_obj A matrix or data.frame with two columns representing
#'   geographic coordinates (e.g., Longitude, Latitude) for each individual.
#' @param output_dir Directory where output files will be saved. 
#'   Default: current working directory.
#' @param ploidy Ploidy level of the organisms. Default: 2.
#' @param missing_val Integer value for missing data in output. Default: 0.
#'
#' @return Invisibly returns NULL. Creates G.txt (genotypes) and XY.txt 
#'   (coordinates) files in the specified directory.
#'
#' @details Alleles are recoded to consecutive integers (1, 2, 3, ...) for each 
#'   locus based on alphabetical order. Missing genotypes are converted to the
#'   specified missing_val value. The output format follows Geneland requirements
#'   with individuals in rows and alleles in columns.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(nancycats)
#' 
#' # Create coordinate data (normally from GPS/sampling locations)
#' coords <- matrix(runif(nrow(nancycats@tab) * 2), ncol = 2)
#' colnames(coords) <- c("Longitude", "Latitude")
#' 
#' # Convert to Geneland format
#' genind_to_geneland(
#'   genind_obj = nancycats,
#'   coords_obj = coords,
#'   output_dir = "geneland_input/"
#' )
#' 
#' # Check created files
#' list.files("geneland_input/")
#' }
#'
#' @importFrom adegenet genind2df nInd locNames
#' @importFrom utils write.table
#' @importFrom stats na.omit setNames
#' @export
genind_to_geneland <- function(genind_obj, coords_obj, output_dir = getwd(), ploidy = 2, missing_val = 0) {
  # Input validation
  if (!inherits(genind_obj, "genind")) {
    stop("Input 'genind_obj' must be an object of class 'genind'.")
  }
  if (!is.matrix(coords_obj) && !is.data.frame(coords_obj)) {
    stop("Input 'coords_obj' must be a matrix or data frame.")
  }
  if (ncol(coords_obj) != 2) {
    stop("'coords_obj' must have exactly 2 columns (e.g., Longitude, Latitude).")
  }
  if (nrow(genind_obj@tab) != nrow(coords_obj)) {
    stop("Number of individuals in 'genind_obj' (", nrow(genind_obj@tab), 
         ") and 'coords_obj' (", nrow(coords_obj), ") do not match.")
  }
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: '", output_dir, "'")
  }
  
  # Process coordinates (XY.txt)
  coords_output_path <- file.path(output_dir, "XY.txt")
  utils::write.table(coords_obj, coords_output_path,
                     row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  message("Coordinates saved to: '", coords_output_path, "'")
  
  # Process genotypes (G.txt)
  num_individuals <- adegenet::nInd(genind_obj)
  num_loci <- length(adegenet::locNames(genind_obj))
  
  # Initialize genotype matrix
  geneland_genotypes <- matrix(NA, nrow = num_individuals, ncol = num_loci * ploidy)
  
  # Process each locus
  for (i in 1:num_loci) {
    locus_name <- adegenet::locNames(genind_obj)[i]
    
    # Extract genotypes as strings
    locus_genotypes_df <- adegenet::genind2df(genind_obj[loc = locus_name], sep = "/", oneColPerAll = FALSE)
    locus_genotypes_vec <- locus_genotypes_df[[locus_name]]
    
    # Get unique alleles for integer mapping
    clean_genotypes_vec <- stats::na.omit(locus_genotypes_vec)
    alleles_split <- unlist(strsplit(clean_genotypes_vec, "/"))
    unique_alleles_for_locus <- sort(unique(alleles_split))
    
    # Create allele to integer mapping
    allele_map <- stats::setNames(1:length(unique_alleles_for_locus), unique_alleles_for_locus)
    
    # Process each individual's genotype
    for (j in 1:num_individuals) {
      genotype_str <- locus_genotypes_vec[j]
      
      if (is.na(genotype_str) || genotype_str == "NA/NA") {
        geneland_genotypes[j, ((i - 1) * ploidy + 1):(i * ploidy)] <- rep(missing_val, ploidy)
      } else {
        alleles_char <- strsplit(genotype_str, "/")[[1]]
        alleles_int <- as.numeric(allele_map[alleles_char])
        geneland_genotypes[j, ((i - 1) * ploidy + 1):(i * ploidy)] <- alleles_int
      }
    }
  }
  
  # Save genotypes
  genotypes_output_path <- file.path(output_dir, "G.txt")
  utils::write.table(geneland_genotypes, genotypes_output_path,
                     row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  message("Genotypes saved to: '", genotypes_output_path, "'")
  
  invisible(NULL)
}

################################################################################
# FUNCTION: genind_to_lea
################################################################################
#' Convert genind object to LEA input format
#'
#' Converts a genind object to the format required by the LEA (Landscape and
#' Ecological Association) R package for population structure and environmental
#' association analysis.
#'
#' @param genind_obj An object of class 'genind' from the adegenet package.
#' @param output_dir Directory where output file will be saved. If NULL,
#'   no file is saved. Default: NULL.
#' @param output_format LEA output format: "geno" (loci in rows) or "lfmm" 
#'   (individuals in rows). Default: "geno".
#'
#' @return Numeric matrix of genotypes coded as 0, 1, or 2. Missing data coded as -9.
#'
#' @details Converts genind to genlight then to LEA format using dartR::gl2geno.
#'   Genotypes are coded as 0, 1, or 2 representing the number of reference alleles.
#'
#' @examples
#' \dontrun{
#' # Convert to LEA format
#' lea_matrix <- genind_to_lea(nancycats, output_format = "geno")
#' 
#' # Save to file
#' genind_to_lea(nancycats, 
#'               output_dir = "lea_analysis/", 
#'               output_format = "lfmm")
#' }
#'
#' @importFrom adegenet as.genlight
#' @importFrom dartR gl2geno
#' @importFrom utils read.table
#' @export
genind_to_lea <- function(genind_obj, output_dir = NULL, output_format = "geno") {
  # Input validation
  if (!inherits(genind_obj, "genind")) {
    stop("Input 'genind_obj' must be an object of class 'genind'.")
  }
  if (!output_format %in% c("geno", "lfmm")) {
    stop("'output_format' must be either 'geno' or 'lfmm'.")
  }
  
  # Convert genind to genlight
  genlight_obj <- adegenet::as.genlight(genind_obj)
  
  # Convert to LEA format using dartR
  temp_geno_file <- tempfile(fileext = ".geno")
  lea_matrix_geno_format <- dartR::gl2geno(
    x = genlight_obj,
    outfile = temp_geno_file,
    outpath = tempdir(),
    verbose = 0
  )
  
  # Read matrix and clean up temp file
  lea_matrix_geno_format <- as.matrix(utils::read.table(temp_geno_file, header = FALSE, sep = "\t", na.strings = "-9"))
  file.remove(temp_geno_file)
  
  # Adjust format if needed
  if (output_format == "lfmm") {
    lea_matrix_output <- t(lea_matrix_geno_format)
  } else {
    lea_matrix_output <- lea_matrix_geno_format
  }
  
  # Clean row/column names
  rownames(lea_matrix_output) <- NULL
  colnames(lea_matrix_output) <- NULL
  
  # Save to file if requested
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      message("Created output directory: '", output_dir, "'")
    }
    
    file_name_base <- deparse(substitute(genind_obj))
    if (length(file_name_base) > 1 || file_name_base == ".") {
      file_name_base <- "converted_genotypes"
    }
    
    output_file_path <- file.path(output_dir, paste0(file_name_base, ".", output_format))
    
    # Use LEA write functions
    if (output_format == "geno") {
      LEA::write.geno(lea_matrix_output, output_file_path)
    } else {
      LEA::write.lfmm(lea_matrix_output, output_file_path)
    }
    message("LEA data saved to: '", output_file_path, "'")
  }
  
  return(lea_matrix_output)
}

################################################################################
# FUNCTION: genind_to_tess3r
################################################################################
#' Convert genind object to tess3r input format
#'
#' Converts a genind object and coordinates to the format required by tess3r
#' for spatial population structure analysis using tensor factorization.
#'
#' @param genind_obj An object of class 'genind' from the adegenet package.
#' @param coords_obj A matrix or data.frame with two columns (X, Y coordinates)
#'   for each individual.
#'
#' @return List with elements:
#'   \describe{
#'     \item{X}{Numeric genotype matrix (individuals in rows, loci in columns)}
#'     \item{coord}{Numeric coordinate matrix (individuals in rows, X and Y in columns)}
#'   }
#'
#' @details Genotypes coded as 0, 1, or 2 (number of reference alleles). 
#'   Missing data will be NA.
#'
#' @examples
#' \dontrun{
#' # Create coordinate data
#' coords <- matrix(runif(nrow(nancycats@tab) * 2), ncol = 2)
#' 
#' # Convert to tess3r format
#' tess3r_data <- genind_to_tess3r(nancycats, coords)
#' 
#' # Use with tess3r package
#' library(tess3r)
#' tess3_result <- tess3(X = tess3r_data$X, coord = tess3r_data$coord, K = 1:5)
#' }
#'
#' @importFrom adegenet as.genlight
#' @export
genind_to_tess3r <- function(genind_obj, coords_obj) {
  # Input validation
  if (!inherits(genind_obj, "genind")) {
    stop("Input 'genind_obj' must be an object of class 'genind'.")
  }
  if (!is.matrix(coords_obj) && !is.data.frame(coords_obj)) {
    stop("Input 'coords_obj' must be a matrix or data frame.")
  }
  if (ncol(coords_obj) != 2) {
    stop("'coords_obj' must have exactly 2 columns (X, Y coordinates).")
  }
  if (nrow(genind_obj@tab) != nrow(coords_obj)) {
    stop("Number of individuals in 'genind_obj' (", nrow(genind_obj@tab), 
         ") and 'coords_obj' (", nrow(coords_obj), ") do not match.")
  }
  
  # Convert genind to genlight and extract genotype matrix
  genlight_obj <- adegenet::as.genlight(genind_obj)
  genotype_matrix_X <- as.matrix(genlight_obj)
  
  # Clean row/column names
  rownames(genotype_matrix_X) <- NULL
  colnames(genotype_matrix_X) <- NULL
  
  # Prepare coordinate matrix
  coordinate_matrix_coord <- as.matrix(coords_obj)
  rownames(coordinate_matrix_coord) <- NULL
  colnames(coordinate_matrix_coord) <- NULL
  
  return(list(X = genotype_matrix_X, coord = coordinate_matrix_coord))
}

################################################################################
# FUNCTION: convert_to_pegas
################################################################################
#' Convert turtle chromosome data to pegas format for analysis
#'
#' Transforms ABM turtle genetic data with expanded chromosome columns into
#' pegas-compatible format for population genetic analysis. Processes data
#' in chunks for memory efficiency with large datasets.
#'
#' @param data Turtle data containing chromosome information with expanded
#'   chromosome columns (e.g., t.dad.chromosome.1.a, t.mom.chromosome.2.b).
#' @param loci_names Names of loci to include in conversion. Default: letters[1:10].
#' @param chunk_size Number of rows to process at once for memory efficiency.
#'   Default: 10000.
#'
#' @return A data.frame suitable for pegas analysis with genotype columns,
#'   population assignments, and metadata.
#'
#' @details Combines alleles from maternal and paternal chromosomes to create
#'   diploid genotypes in "allele1/allele2" format. Extracts numeric allele
#'   values from allele names (e.g., "a1" becomes "1"). Processes data in
#'   chunks to handle large datasets efficiently.
#'
#' @examples
#' \dontrun{
#' # Convert turtle data to pegas format
#' pegas_data <- convert_to_pegas(
#'   turtle_data, 
#'   loci_names = letters[1:5]
#' )
#' 
#' # Convert to pegas loci object
#' library(pegas)
#' loci_obj <- df2loci(pegas_data[, 1:5], 
#'                     pop = pegas_data$pop,
#'                     ploidy = 2, 
#'                     sep = "/")
#' }
#'
#' @importFrom stats setNames
#' @export
convert_to_pegas <- function(data, loci_names = letters[1:10], chunk_size = 10000) {
  # Check if we have the expanded chromosome columns
  chr_cols <- grep("chromosome.*\\.[a-z]$", names(data), value = TRUE)
  
  if (length(chr_cols) == 0) {
    stop("No expanded chromosome columns found. Make sure data was read with read_turtles()")
  }
  
  # Get the unique loci from column names
  loci <- unique(sub(".*\\.", "", chr_cols))
  if (!missing(loci_names)) {
    loci <- loci_names[1:length(loci)]
  }
  
  # Process in chunks for large datasets
  n_rows <- nrow(data)
  n_chunks <- ceiling(n_rows / chunk_size)
  
  gene_data_chunks <- list()
  
  for (chunk in 1:n_chunks) {
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, n_rows)
    
    chunk_data <- data[start_idx:end_idx, ]
    
    gene_data_list <- lapply(1:nrow(chunk_data), function(i) {
      genotypes <- character(length(loci))
      names(genotypes) <- loci
      
      for (j in seq_along(loci)) {
        locus <- loci[j]
        
        # Get alleles from all four chromosomes
        dad1_col <- paste0("t.dad.chromosome.1.", locus)
        dad2_col <- paste0("t.dad.chromosome.2.", locus)
        mom1_col <- paste0("t.mom.chromosome.1.", locus)
        mom2_col <- paste0("t.mom.chromosome.2.", locus)
        
        alleles <- c()
        if (dad1_col %in% names(chunk_data) && !is.na(chunk_data[[dad1_col]][i])) {
          # Extract just the numeric part from the allele (e.g., "a1" -> "1")
          allele <- chunk_data[[dad1_col]][i]
          allele_value <- substr(allele, 2, nchar(allele))
          alleles <- c(alleles, allele_value)
        }
        if (dad2_col %in% names(chunk_data) && !is.na(chunk_data[[dad2_col]][i])) {
          allele <- chunk_data[[dad2_col]][i]
          allele_value <- substr(allele, 2, nchar(allele))
          alleles <- c(alleles, allele_value)
        }
        if (mom1_col %in% names(chunk_data) && !is.na(chunk_data[[mom1_col]][i])) {
          allele <- chunk_data[[mom1_col]][i]
          allele_value <- substr(allele, 2, nchar(allele))
          alleles <- c(alleles, allele_value)
        }
        if (mom2_col %in% names(chunk_data) && !is.na(chunk_data[[mom2_col]][i])) {
          allele <- chunk_data[[mom2_col]][i]
          allele_value <- substr(allele, 2, nchar(allele))
          alleles <- c(alleles, allele_value)
        }
        
        # Create genotype (take first two alleles)
        if (length(alleles) >= 2) {
          genotypes[j] <- paste0(alleles[1], "/", alleles[2])
        } else if (length(alleles) == 1) {
          genotypes[j] <- paste0(alleles[1], "/", alleles[1])  # Homozygote
        } else {
          genotypes[j] <- NA
        }
      }
      
      return(genotypes)
    })
    
    gene_data_chunks[[chunk]] <- gene_data_list
  }
  
  # Combine all chunks
  all_gene_data <- unlist(gene_data_chunks, recursive = FALSE)
  
  # Create genotype matrix
  genotype_data <- do.call(rbind, all_gene_data)
  colnames(genotype_data) <- loci
  rownames(genotype_data) <- data$who
  
  # Create final data frame
  final_data <- data.frame(
    genotype_data,
    pop = data$t.birth.pond.index,
    ticks = data$ticks,
    stringsAsFactors = FALSE
  )
  
  return(final_data)
}