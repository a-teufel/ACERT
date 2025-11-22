################################################################################
# R/spatial_analysis.R
################################################################################
# 
# Functions for spatial and environmental analysis of ABM simulation data
#
################################################################################

################################################################################
# FUNCTION: assign_turtles_to_ponds
################################################################################
#' Assign turtles to closest ponds based on spatial coordinates
#'
#' Assigns each turtle to the pond whose center point is geographically closest
#' using Euclidean distance. This function corrects potential inaccuracies in
#' the original pond assignments by using actual spatial coordinates.
#'
#' @param turtle_data A data.frame containing turtle data with `xcor` and `ycor`
#'   columns representing turtle locations.
#' @param pond_centers A data.frame containing pond center coordinates with exactly
#'   3 columns: pond identifier, x-coordinate, y-coordinate (in that order).
#'
#' @return A new data.frame based on input turtle_data with an additional
#'   `corrected.current.pond` column containing the nearest pond identifier.
#'
#' @details Uses the FNN package's get.knnx function for efficient nearest
#'   neighbor search. The pond_centers data.frame column order is critical:
#'   column 1 must be pond identifiers, columns 2-3 must be x,y coordinates.
#'
#' @examples
#' # Create sample turtle data
#' turtle_sample <- data.frame(
#'   who = 1:6,
#'   xcor = c(1, 1.5, 5, 5.2, 9, 8.5),
#'   ycor = c(1, 1.2, 5, 5.5, 9, 8.3),
#'   t.current.pond.index = c(0, 0, 1, 1, 2, 2)
#' )
#' 
#' # Create pond centers
#' pond_centers <- data.frame(
#'   pond_id = c(0, 1, 2),
#'   center_x = c(1, 5, 9),
#'   center_y = c(1, 5, 9)
#' )
#' 
#' # Assign turtles to nearest ponds
#' corrected_data <- assign_turtles_to_ponds(turtle_sample, pond_centers)
#' 
#' # Compare original vs corrected assignments
#' corrected_data[c("who", "xcor", "ycor", "t.current.pond.index", "corrected.current.pond")]
#'
#' @importFrom dplyr select mutate
#' @importFrom FNN get.knnx
#' @export
assign_turtles_to_ponds <- function(turtle_data, pond_centers) {
  # Input validation
  if (!is.data.frame(turtle_data)) {
    stop("turtle_data must be a data.frame")
  }
  
  if (!all(c("xcor", "ycor") %in% names(turtle_data))) {
    stop("turtle_data must contain 'xcor' and 'ycor' columns")
  }
  
  if (!is.data.frame(pond_centers)) {
    stop("pond_centers must be a data.frame")
  }
  
  if (ncol(pond_centers) < 3) {
    stop("pond_centers must have at least 3 columns: pond_id, x_coord, y_coord")
  }
  
  if (nrow(pond_centers) == 0) {
    stop("pond_centers cannot be empty")
  }
  
  # Extract turtle coordinates
  turtle_coords <- turtle_data %>%
    dplyr::select(xcor, ycor) %>%
    as.matrix()
  
  # Extract pond information
  pond_labels <- pond_centers[[1]]  # First column is pond identifier
  pond_coords <- pond_centers %>%
    dplyr::select(2, 3) %>%  # Columns 2 and 3 are x,y coordinates
    as.matrix()
  
  # Find nearest pond for each turtle
  nn_results <- FNN::get.knnx(
    data = pond_coords, 
    query = turtle_coords, 
    k = 1
  )
  
  # Extract nearest pond labels
  nearest_indices <- nn_results$nn.index[, 1]
  corrected_pond_labels <- pond_labels[nearest_indices]
  
  # Add corrected pond assignment to turtle data
  result_data <- turtle_data %>%
    dplyr::mutate(corrected.current.pond = corrected_pond_labels)
  
  return(result_data)
}

################################################################################
# FUNCTION: calculate_turtle_density
################################################################################
#' Calculate turtle density on spatial patches
#'
#' Computes turtle density per spatial patch using floor coordinates to create
#' a regular grid. Provides both mean and maximum density statistics across
#' time points.
#'
#' @param turtle_data A data.frame with turtle location data including `xcor`,
#'   `ycor`, and `ticks` columns.
#' @param environment_data Optional data.frame with patch information for
#'   merging additional environmental variables. Default: NULL.
#' @param time_points Specific time points to include in density calculation.
#'   If NULL, uses all time points. Default: NULL.
#' @param grid_resolution Grid resolution for density calculation. Default: 1
#'   (uses floor coordinates as-is).
#'
#' @return A data.frame with turtle density per patch containing:
#'   \describe{
#'     \item{xcor_grid, ycor_grid}{Grid coordinates}
#'     \item{mean_density}{Mean number of turtles across time points}
#'     \item{max_density}{Maximum number of turtles at any time point}
#'     \item{total_turtle_time}{Total turtle-time units}
#'     \item{n_time_points}{Number of time points with data}
#'   }
#'
#' @examples
#' # Create sample turtle movement data
#' turtle_movement <- data.frame(
#'   who = rep(1:10, 5),
#'   xcor = runif(50, 0, 10),
#'   ycor = runif(50, 0, 10),
#'   ticks = rep(1:5, each = 10)
#' )
#' 
#' # Calculate density
#' density_result <- calculate_turtle_density(turtle_movement)
#' 
#' # View high-density areas
#' high_density <- subset(density_result, mean_density > 2)
#' print(high_density)
#' 
#' # Calculate density for specific time points
#' recent_density <- calculate_turtle_density(
#'   turtle_movement, 
#'   time_points = c(4, 5)
#' )
#'
#' @importFrom dplyr group_by summarise mutate left_join
#' @importFrom tidyr replace_na
#' @export
calculate_turtle_density <- function(turtle_data, 
                                     environment_data = NULL,
                                     time_points = NULL,
                                     grid_resolution = 1) {
  
  # Input validation
  required_cols <- c("xcor", "ycor", "ticks")
  missing_cols <- setdiff(required_cols, names(turtle_data))
  if (length(missing_cols) > 0) {
    stop("turtle_data missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!is.numeric(grid_resolution) || grid_resolution <= 0) {
    stop("grid_resolution must be a positive number")
  }
  
  # Filter time points if specified
  if (!is.null(time_points)) {
    turtle_data <- turtle_data %>%
      dplyr::filter(ticks %in% time_points)
    
    if (nrow(turtle_data) == 0) {
      stop("No data found for specified time points")
    }
  }
  
  # Create grid coordinates
  density_data <- turtle_data %>%
    dplyr::mutate(
      xcor_grid = floor(xcor / grid_resolution) * grid_resolution,
      ycor_grid = floor(ycor / grid_resolution) * grid_resolution
    ) %>%
    dplyr::group_by(xcor_grid, ycor_grid, ticks) %>%
    dplyr::summarise(turtle_count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(xcor_grid, ycor_grid) %>%
    dplyr::summarise(
      mean_density = mean(turtle_count),
      max_density = max(turtle_count),
      total_turtle_time = sum(turtle_count),
      n_time_points = dplyr::n(),
      .groups = "drop"
    )
  
  # Merge with environment data if provided
  if (!is.null(environment_data)) {
    if (all(c("pxcor", "pycor") %in% names(environment_data))) {
      # Match environment coordinates to grid coordinates
      env_summary <- environment_data %>%
        dplyr::mutate(
          xcor_grid = floor(pxcor / grid_resolution) * grid_resolution,
          ycor_grid = floor(pycor / grid_resolution) * grid_resolution
        ) %>%
        dplyr::group_by(xcor_grid, ycor_grid) %>%
        dplyr::summarise(
          dplyr::across(where(is.numeric), mean, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::select(-pxcor, -pycor)  # Remove original coordinates
      
      density_data <- density_data %>%
        dplyr::left_join(env_summary, by = c("xcor_grid", "ycor_grid"))
    } else {
      warning("environment_data missing 'pxcor' and 'pycor' columns; not merged")
    }
  }
  
  return(density_data)
}


################################################################################
# HELPER FUNCTIONS
################################################################################

#' Assign environmental values to turtle locations
#' @param turtle_data Turtle location data
#' @param environment_data Environmental grid data
#' @param habitat_variables Variables to extract
#' @return Data frame with environmental values
#' @keywords internal
assign_environmental_values <- function(turtle_data, environment_data, habitat_variables) {
  # Simple nearest neighbor assignment using floor coordinates
  turtle_with_env <- turtle_data %>%
    dplyr::mutate(
      pxcor_match = floor(xcor),
      pycor_match = floor(ycor)
    ) %>%
    dplyr::left_join(
      environment_data,
      by = c("pxcor_match" = "pxcor", "pycor_match" = "pycor")
    )
  
  return(turtle_with_env)
}


#' Calculate Hopkins statistic for spatial randomness
#' @param coords Matrix of coordinates
#' @param bounds Study area boundaries
#' @return Hopkins statistic value
#' @keywords internal
calculate_hopkins_statistic <- function(coords, bounds) {
  n <- nrow(coords)
  
  # Sample size for Hopkins test (typically n/10)
  m <- max(5, min(50, floor(n / 10)))
  
  # Generate random points in study area
  random_points <- matrix(c(
    stats::runif(m, bounds[1], bounds[2]),
    stats::runif(m, bounds[3], bounds[4])
  ), ncol = 2)
  
  # Calculate distances from random points to nearest data points
  u_distances <- numeric(m)
  for (i in 1:m) {
    distances_to_data <- sqrt(rowSums((coords - matrix(random_points[i, ], nrow = n, ncol = 2, byrow = TRUE))^2))
    u_distances[i] <- min(distances_to_data)
  }
  
  # Calculate distances from random sample of data points to nearest neighbors
  sample_indices <- sample(1:n, m)
  w_distances <- numeric(m)
  for (i in 1:m) {
    idx <- sample_indices[i]
    distances_to_others <- sqrt(rowSums((coords[-idx, ] - matrix(coords[idx, ], nrow = n-1, ncol = 2, byrow = TRUE))^2))
    w_distances[i] <- min(distances_to_others)
  }
  
  # Calculate Hopkins statistic
  hopkins_stat <- sum(u_distances^2) / (sum(u_distances^2) + sum(w_distances^2))
  
  return(hopkins_stat)
}
