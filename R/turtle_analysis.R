################################################################################
# R/turtle_analysis.R
################################################################################

################################################################################
# FUNCTION: analyze_turtle_property
################################################################################
#' Analyze turtle properties across different groupings
#'
#' Computes summary statistics for a specified turtle property, optionally
#' grouped by categorical variables. Useful for comparing traits across
#' time points, life stages, populations, or experimental treatments.
#'
#' @param data Turtle data frame or tibble. Must contain the specified property
#'   column along with \code{ticks} and \code{run_number} columns.
#' @param property Character string specifying the property to analyze.
#'   Common properties include: "t.energy", "t.size.cm", "t.age", "t.stress".
#' @param grouping Character vector of grouping variable(s). Use NULL for no grouping
#'   beyond run_number and ticks. Common grouping variables: "breed", "t.sex",
#'   "t.birth.pond.index". Default: NULL.
#' @param summarize Logical indicating whether to calculate summary statistics
#'   (mean, SD, SE, n). If FALSE, returns filtered raw data. Default: TRUE.
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item \code{run_number}: Simulation replicate identifier
#'     \item \code{ticks}: Time point
#'     \item Grouping variable(s) if specified
#'     \item \code{mean_value}: Mean of the property
#'     \item \code{sd_value}: Standard deviation
#'     \item \code{n}: Sample size
#'     \item \code{se_value}: Standard error of the mean
#'   }
#'
#' @details This function automatically:
#' \itemize{
#'   \item Removes NA values from the analyzed property
#'   \item Groups by run_number and ticks
#'   \item Adds additional grouping levels if specified
#'   \item Calculates standard errors for plotting confidence intervals
#' }
#'
#' The function is designed to work seamlessly with plotting functions and
#' downstream statistical analyses.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' turtle_data <- read_turtles(path = "simulation_output/")
#' 
#' # Basic analysis: energy levels over time
#' energy_summary <- analyze_turtle_property(
#'   data = turtle_data,
#'   property = "t.energy"
#' )
#' 
#' # View results
#' head(energy_summary)
#' 
#' # Analyze size by life stage
#' size_by_breed <- analyze_turtle_property(
#'   data = turtle_data,
#'   property = "t.size.cm",
#'   grouping = "breed"
#' )
#' 
#' # Compare across multiple groupings
#' energy_by_sex_and_pond <- analyze_turtle_property(
#'   data = turtle_data,
#'   property = "t.energy",
#'   grouping = c("t.sex", "t.current.pond.index")
#' )
#' 
#' # Get raw data instead of summaries
#' raw_stress <- analyze_turtle_property(
#'   data = turtle_data,
#'   property = "t.stress",
#'   summarize = FALSE
#' )
#' 
#' # Analyze age by sex
#' age_by_sex <- analyze_turtle_property(
#'   turtle_data,
#'   property = "t.age",
#'   grouping = "t.sex"
#' )
#' 
#' # Visualize results
#' library(ggplot2)
#' library(dplyr)
#' 
#' # Plot energy over time by life stage
#' energy_by_breed %>%
#'   ggplot(aes(x = ticks, y = mean_value, color = breed)) +
#'   geom_line() +
#'   geom_ribbon(aes(ymin = mean_value - se_value,
#'                   ymax = mean_value + se_value,
#'                   fill = breed), alpha = 0.3) +
#'   facet_wrap(~run_number) +
#'   labs(title = "Energy Levels by Life Stage",
#'        y = "Mean Energy", x = "Time (ticks)")
#' 
#' # Statistical summary across runs
#' energy_summary %>%
#'   group_by(ticks) %>%
#'   summarise(
#'     grand_mean = mean(mean_value),
#'     between_run_sd = sd(mean_value),
#'     n_runs = n()
#'   )
#' 
#' # Compare final vs initial values
#' energy_summary %>%
#'   filter(ticks %in% c(min(ticks), max(ticks))) %>%
#'   group_by(run_number, ticks) %>%
#'   summarise(mean_energy = mean(mean_value))
#' 
#' # Analyze stress levels in adults only
#' adult_stress <- turtle_data %>%
#'   filter(breed == "adults") %>%
#'   analyze_turtle_property(property = "t.stress")
#' 
#' # Multiple properties workflow
#' properties <- c("t.energy", "t.size.cm", "t.stress", "t.age")
#' 
#' results <- lapply(properties, function(prop) {
#'   analyze_turtle_property(turtle_data, property = prop, grouping = "breed")
#' })
#' names(results) <- properties
#' 
#' # Check sample sizes
#' energy_by_breed %>%
#'   group_by(breed) %>%
#'   summarise(
#'     total_measurements = sum(n),
#'     mean_n_per_timepoint = mean(n)
#'   )
#' }
#'
#' @importFrom stats sd var
#' @importFrom dplyr group_by summarize mutate select filter across all_of n
#' @importFrom rlang .data
#' @export
analyze_turtle_property <- function(data, property, grouping = NULL, summarize = TRUE) {
  # Check if property exists in data
  if (!property %in% names(data)) {
    stop(paste("Property", property, "not found in data"))
  }
  
  # Create analysis data
  analysis_data <- data %>%
    dplyr::select(all_of(c(property, grouping, "ticks", "run_number"))) %>%
    dplyr::filter(!is.na(.data[[property]]))
  
  if (summarize) {
    if (is.null(grouping)) {
      result <- analysis_data %>%
        dplyr::group_by(run_number, ticks) %>%
        dplyr::summarize(
          mean_value = mean(.data[[property]], na.rm = TRUE),
          sd_value = sd(.data[[property]], na.rm = TRUE),
          n = n(),
          se_value = sd_value / sqrt(n),
          .groups = "drop"
        )
    } else {
      result <- analysis_data %>%
        dplyr::group_by(across(all_of(c("run_number", "ticks", grouping)))) %>%
        dplyr::summarize(
          mean_value = mean(.data[[property]], na.rm = TRUE),
          sd_value = sd(.data[[property]], na.rm = TRUE),
          n = n(),
          se_value = sd_value / sqrt(n),
          .groups = "drop"
        )
    }
  } else {
    result <- analysis_data
  }
  
  return(result)
}


################################################################################
# FUNCTION: calculate_sex_ratio
################################################################################
#' Calculate sex ratios for turtle populations
#'
#' Computes sex ratios and sex-specific counts across time points and
#' simulation replicates. Useful for analyzing sex determination patterns,
#' demographic changes, and population viability.
#'
#' @param data Turtle data frame or tibble. Must contain a \code{t.sex} column
#'   with values "male" and "female", plus \code{run_number} and \code{ticks}.
#' @param grouping Character vector of additional grouping variable(s) beyond
#'   run_number and ticks. Useful for comparing sex ratios across populations,
#'   life stages, or other categorical factors. Default: NULL.
#'
#' @return A tibble with sex ratio statistics:
#'   \itemize{
#'     \item \code{run_number}: Simulation replicate
#'     \item \code{ticks}: Time point
#'     \item Additional grouping variables if specified
#'     \item \code{n_total}: Total count of individuals
#'     \item \code{n_male}: Count of males
#'     \item \code{n_female}: Count of females
#'     \item \code{prop_male}: Proportion male (0-1)
#'     \item \code{prop_female}: Proportion female (0-1)
#'     \item \code{sex_ratio}: Male to female ratio (males/females)
#'   }
#'
#' @details Sex ratio interpretation:
#' \itemize{
#'   \item sex_ratio = 1.0: Equal males and females (1:1)
#'   \item sex_ratio > 1.0: Male-biased (e.g., 2.0 = 2 males per female)
#'   \item sex_ratio < 1.0: Female-biased (e.g., 0.5 = 1 male per 2 females)
#' }
#'
#' The function handles missing sex data by excluding NA values from calculations.
#'
#' @examples
#' \dontrun{
#' # Load turtle data
#' turtle_data <- read_turtles()
#' 
#' # Basic sex ratio calculation
#' sex_ratios <- calculate_sex_ratio(turtle_data)
#' 
#' # View results
#' head(sex_ratios)
#' summary(sex_ratios$sex_ratio)
#' 
#' # Sex ratio by life stage
#' sex_ratio_by_breed <- calculate_sex_ratio(
#'   data = turtle_data,
#'   grouping = "breed"
#' )
#' 
#' # Sex ratio by population (pond)
#' sex_ratio_by_pond <- calculate_sex_ratio(
#'   turtle_data,
#'   grouping = "t.current.pond.index"
#' )
#' 
#' # Multiple groupings
#' sex_ratio_detailed <- calculate_sex_ratio(
#'   turtle_data,
#'   grouping = c("breed", "t.birth.pond.index")
#' )
#' 
#' # Visualize sex ratio over time
#' library(ggplot2)
#' library(dplyr)
#' 
#' sex_ratios %>%
#'   ggplot(aes(x = ticks, y = sex_ratio)) +
#'   geom_line(aes(color = run_number)) +
#'   geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
#'   labs(title = "Sex Ratio Over Time",
#'        y = "Sex Ratio (M:F)",
#'        x = "Time (ticks)") +
#'   scale_y_continuous(limits = c(0, 2))
#' 
#' # Check for sex ratio bias
#' sex_ratios %>%
#'   mutate(bias = case_when(
#'     sex_ratio > 1.2 ~ "Male-biased",
#'     sex_ratio < 0.8 ~ "Female-biased",
#'     TRUE ~ "Balanced"
#'   )) %>%
#'   group_by(bias) %>%
#'   summarise(n_timepoints = n())
#' 
#' # Compare across life stages
#' sex_ratio_by_breed %>%
#'   ggplot(aes(x = ticks, y = sex_ratio, color = breed)) +
#'   geom_line() +
#'   geom_hline(yintercept = 1, linetype = "dashed") +
#'   facet_wrap(~run_number) +
#'   labs(title = "Sex Ratio by Life Stage")
#' 
#' # Summary statistics
#' sex_ratios %>%
#'   group_by(run_number) %>%
#'   summarise(
#'     mean_sex_ratio = mean(sex_ratio, na.rm = TRUE),
#'     sd_sex_ratio = sd(sex_ratio, na.rm = TRUE),
#'     min_sex_ratio = min(sex_ratio, na.rm = TRUE),
#'     max_sex_ratio = max(sex_ratio, na.rm = TRUE)
#'   )
#' 
#' # Temporal trend analysis
#' sex_ratios %>%
#'   group_by(ticks) %>%
#'   summarise(
#'     mean_ratio = mean(sex_ratio),
#'     se_ratio = sd(sex_ratio) / sqrt(n())
#'   ) %>%
#'   ggplot(aes(x = ticks, y = mean_ratio)) +
#'   geom_line() +
#'   geom_ribbon(aes(ymin = mean_ratio - se_ratio,
#'                   ymax = mean_ratio + se_ratio),
#'               alpha = 0.3) +
#'   geom_hline(yintercept = 1, linetype = "dashed", color = "red")
#' 
#' # Population size and sex composition
#' sex_ratio_by_pond %>%
#'   filter(ticks == max(ticks)) %>%
#'   ggplot(aes(x = factor(t.current.pond.index))) +
#'   geom_bar(aes(y = n_male, fill = "Male"), 
#'            stat = "identity", position = "dodge") +
#'   geom_bar(aes(y = n_female, fill = "Female"), 
#'            stat = "identity", position = "dodge") +
#'   facet_wrap(~run_number) +
#'   labs(title = "Final Population Sex Composition by Pond",
#'        x = "Pond", y = "Count")
#' 
#' # Adults only sex ratio
#' turtle_data %>%
#'   filter(breed == "adults") %>%
#'   calculate_sex_ratio() %>%
#'   summary()
#' }
#'
#' @export
calculate_sex_ratio <- function(data, grouping = NULL) {
  if (!"t.sex" %in% names(data)) {
    stop("Sex information (t.sex) not found in data")
  }
  
  group_vars <- c("run_number", "ticks")
  if (!is.null(grouping)) {
    group_vars <- c(group_vars, grouping)
  }
  
  sex_summary <- data %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::summarize(
      n_total = n(),
      n_male = sum(t.sex == "male", na.rm = TRUE),
      n_female = sum(t.sex == "female", na.rm = TRUE),
      prop_male = n_male / n_total,
      prop_female = n_female / n_total,
      sex_ratio = n_male / n_female,
      .groups = "drop"
    )
  
  return(sex_summary)
}

################################################################################
# FUNCTION: summarize_turtle_data
################################################################################
#' Generate comprehensive summary statistics for turtle data
#'
#' Computes mean, standard deviation, minimum, and maximum for all numeric
#' turtle properties. Provides a high-level overview of simulation outputs
#' with flexible grouping options.
#'
#' @param turtle_data Turtle data frame or tibble containing numeric properties
#'   to summarize.
#' @param by_replicate Logical indicating whether to group by run_number.
#'   If FALSE, aggregates across all replicates. Default: TRUE.
#' @param by_time Logical indicating whether to group by ticks (time points).
#'   If FALSE, aggregates across all time points. Default: TRUE.
#'
#' @return A data frame with summary statistics for each numeric column:
#'   \itemize{
#'     \item Grouping variables (run_number and/or ticks if enabled)
#'     \item \code{n_turtles}: Count of individuals
#'     \item \code{[property]_mean}: Mean value for each property
#'     \item \code{[property]_sd}: Standard deviation
#'     \item \code{[property]_min}: Minimum value
#'     \item \code{[property]_max}: Maximum value
#'   }
#'
#' @details The function automatically:
#' \itemize{
#'   \item Identifies all numeric columns (excluding ID and coordinate columns)
#'   \item Handles NA values in calculations
#'   \item Returns NA for statistics when all values are missing
#'   \item Generates descriptive column names with suffixes (_mean, _sd, etc.)
#' }
#'
#' Excluded from summaries: who (ID), ticks (if used as grouping), xcor, ycor
#' (spatial coordinates).
#'
#' @examples
#' \dontrun{
#' # Load data
#' turtle_data <- read_turtles()
#' 
#' # Complete summary by replicate and time
#' full_summary <- summarize_turtle_data(turtle_data)
#' 
#' # View structure
#' head(full_summary)
#' names(full_summary)
#' 
#' # Overall summary (no grouping)
#' overall_summary <- summarize_turtle_data(
#'   turtle_data,
#'   by_replicate = FALSE,
#'   by_time = FALSE
#' )
#' print(overall_summary)
#' 
#' # By time only (average across replicates)
#' time_summary <- summarize_turtle_data(
#'   turtle_data,
#'   by_replicate = FALSE,
#'   by_time = TRUE
#' )
#' 
#' # By replicate only (average across time)
#' replicate_summary <- summarize_turtle_data(
#'   turtle_data,
#'   by_replicate = TRUE,
#'   by_time = FALSE
#' )
#' 
#' # Examine specific properties
#' library(dplyr)
#' 
#' full_summary %>%
#'   select(run_number, ticks, n_turtles, 
#'          t.energy_mean, t.energy_sd,
#'          t.size.cm_mean, t.size.cm_sd) %>%
#'   head(10)
#' 
#' # Track population size over time
#' full_summary %>%
#'   ggplot(aes(x = ticks, y = n_turtles, color = run_number)) +
#'   geom_line() +
#'   labs(title = "Population Size Over Time",
#'        y = "Number of Individuals")
#' 
#' # Compare energy levels across runs
#' full_summary %>%
#'   ggplot(aes(x = ticks, y = t.energy_mean, color = run_number)) +
#'   geom_line() +
#'   geom_ribbon(aes(ymin = t.energy_mean - t.energy_sd,
#'                   ymax = t.energy_mean + t.energy_sd,
#'                   fill = run_number), alpha = 0.2) +
#'   labs(title = "Mean Energy ± SD")
#' 
#' # Summary table for final time point
#' full_summary %>%
#'   filter(ticks == max(ticks)) %>%
#'   select(run_number, n_turtles, ends_with("_mean")) %>%
#'   arrange(run_number)
#' 
#' # Range of values across simulation
#' time_summary %>%
#'   select(ticks, 
#'          t.age_min, t.age_max,
#'          t.size.cm_min, t.size.cm_max,
#'          t.energy_min, t.energy_max)
#' 
#' # Coefficient of variation over time
#' full_summary %>%
#'   mutate(
#'     energy_cv = t.energy_sd / t.energy_mean,
#'     size_cv = t.size.cm_sd / t.size.cm_mean
#'   ) %>%
#'   select(run_number, ticks, energy_cv, size_cv) %>%
#'   tidyr::pivot_longer(cols = ends_with("_cv"),
#'                       names_to = "property",
#'                       values_to = "cv") %>%
#'   ggplot(aes(x = ticks, y = cv, color = property)) +
#'   geom_line() +
#'   facet_wrap(~run_number) +
#'   labs(title = "Coefficient of Variation Over Time")
#' 
#' # Export summary for reporting
#' write.csv(overall_summary, "simulation_summary.csv", row.names = FALSE)
#' 
#' # Compare first and last time points
#' full_summary %>%
#'   filter(ticks %in% c(min(ticks), max(ticks))) %>%
#'   select(run_number, ticks, t.energy_mean, t.size.cm_mean, t.stress_mean) %>%
#'   tidyr::pivot_wider(names_from = ticks, values_from = ends_with("_mean"))
#' 
#' # Life stage specific summaries
#' turtle_data %>%
#'   filter(breed == "adults") %>%
#'   summarize_turtle_data() %>%
#'   head()
#' 
#' # Check for trends
#' time_summary %>%
#'   select(ticks, ends_with("_mean")) %>%
#'   tidyr::pivot_longer(-ticks, names_to = "property", values_to = "mean") %>%
#'   ggplot(aes(x = ticks, y = mean)) +
#'   geom_line() +
#'   facet_wrap(~property, scales = "free_y") +
#'   labs(title = "All Properties Over Time")
#' }
#'
#' @export
summarize_turtle_data <- function(turtle_data, by_replicate = TRUE, by_time = TRUE) {
  numeric_cols <- names(turtle_data)[sapply(turtle_data, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, c("who", "ticks", "xcor", "ycor"))
  
  group_vars <- c()
  if (by_replicate) group_vars <- c(group_vars, "run_number")
  if (by_time) group_vars <- c(group_vars, "ticks")
  
  summarize_fun <- function(x, fun) {
    if (all(is.na(x))) {
      return(NA_real_)
    } else {
      return(fun(x, na.rm = TRUE))
    }
  }
  
  if (length(group_vars) == 0) {
    summary_data <- turtle_data %>%
      dplyr::summarise(
        n_turtles = dplyr::n(),
        dplyr::across(
          all_of(numeric_cols),
          list(
            mean = ~mean(.x, na.rm = TRUE),
            sd = ~sd(.x, na.rm = TRUE),
            min = ~summarize_fun(.x, min),
            max = ~summarize_fun(.x, max)
          ),
          .names = "{.col}_{.fn}"
        )
      )
  } else {
    summary_data <- turtle_data %>%
      dplyr::group_by(across(all_of(group_vars))) %>%
      dplyr::summarise(
        n_turtles = dplyr::n(),
        dplyr::across(
          all_of(numeric_cols),
          list(
            mean = ~mean(.x, na.rm = TRUE),
            sd = ~sd(.x, na.rm = TRUE),
            min = ~summarize_fun(.x, min),
            max = ~summarize_fun(.x, max)
          ),
          .names = "{.col}_{.fn}"
        ),
        .groups = "drop"
      )
  }
  
  return(summary_data)
}

################################################################################
# FUNCTION: compare_turtle_groups
################################################################################
#' Compare turtle properties between groups with statistical tests
#'
#' Performs statistical hypothesis tests to compare a turtle property across
#' different groups. Supports t-tests, Wilcoxon tests, ANOVA, and Kruskal-Wallis
#' tests with automatic test selection based on number of groups.
#'
#' @param turtle_data Turtle data frame or tibble.
#' @param property Character string specifying the numeric property to compare
#'   (e.g., "t.energy", "t.size.cm", "t.stress").
#' @param group_var Character string specifying the grouping variable
#'   (e.g., "breed", "t.sex", "t.birth.pond.index").
#' @param test_type Character string specifying the statistical test. Options:
#'   \itemize{
#'     \item "auto": Automatic selection (t-test for 2 groups, ANOVA for 3+)
#'     \item "t.test": Student's t-test (requires 2 groups)
#'     \item "wilcox": Wilcoxon rank-sum test (requires 2 groups)
#'     \item "anova": One-way ANOVA (for 2+ groups)
#'     \item "kruskal": Kruskal-Wallis test (for 2+ groups)
#'   }
#'   Default: "auto".
#' @param time_point Numeric value specifying a single time point to analyze.
#'   If NULL, uses all available time points (not recommended for most tests).
#'   Default: NULL.
#'
#' @return Statistical test result object from base R stats functions:
#'   \itemize{
#'     \item \code{htest} object for t.test, wilcox.test, kruskal.test
#'     \item \code{aov} object for ANOVA
#'   }
#'   Use \code{summary()} on the result for detailed output.
#'
#' @details Test selection guide:
#' \itemize{
#'   \item \strong{t-test}: Two groups, normally distributed data
#'   \item \strong{Wilcoxon}: Two groups, non-normal or ordinal data
#'   \item \strong{ANOVA}: Three or more groups, normally distributed
#'   \item \strong{Kruskal-Wallis}: Three or more groups, non-normal data
#' }
#'
#' The function automatically removes NA values from both the property and
#' grouping variable before analysis.
#'
#' @examples
#' \dontrun{
#' # Load data
#' turtle_data <- read_turtles()
#' 
#' # Compare energy between sexes at final time point
#' final_tick <- max(turtle_data$ticks)
#' 
#' energy_sex_test <- compare_turtle_groups(
#'   turtle_data = turtle_data,
#'   property = "t.energy",
#'   group_var = "t.sex",
#'   test_type = "auto",
#'   time_point = final_tick
#' )
#' print(energy_sex_test)
#' 
#' # Wilcoxon test for non-normal data
#' size_sex_test <- compare_turtle_groups(
#'   turtle_data,
#'   property = "t.size.cm",
#'   group_var = "t.sex",
#'   test_type = "wilcox",
#'   time_point = final_tick
#' )
#' 
#' # Compare across life stages (3+ groups)
#' energy_breed_test <- compare_turtle_groups(
#'   turtle_data,
#'   property = "t.energy",
#'   group_var = "breed",
#'   test_type = "anova",
#'   time_point = 1000
#' )
#' summary(energy_breed_test)
#' 
#' # Non-parametric alternative for life stages
#' stress_breed_test <- compare_turtle_groups(
#'   turtle_data,
#'   property = "t.stress",
#'   group_var = "breed",
#'   test_type = "kruskal",
#'   time_point = 1000
#' )
#' 
#' # Compare across ponds
#' energy_pond_test <- compare_turtle_groups(
#'   turtle_data,
#'   property = "t.energy",
#'   group_var = "t.current.pond.index",
#'   time_point = max(turtle_data$ticks)
#' )
#' 
#' # Workflow: Test at multiple time points
#' library(dplyr)
#' 
#' time_points <- seq(0, max(turtle_data$ticks), by = 1000)
#' 
#' test_results <- lapply(time_points, function(tp) {
#'   test <- compare_turtle_groups(
#'     turtle_data,
#'     property = "t.energy",
#'     group_var = "t.sex",
#'     test_type = "t.test",
#'     time_point = tp
#'   )
#'   data.frame(
#'     time_point = tp,
#'     p_value = test$p.value,
#'     statistic = test$statistic,
#'     estimate_diff = diff(test$estimate)
#'   )
#' })
#' 
#' results_df <- do.call(rbind, test_results)
#' 
#' # Plot p-values over time
#' library(ggplot2)
#' ggplot(results_df, aes(x = time_point, y = p_value)) +
#'   geom_line() +
#'   geom_point() +
#'   geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
#'   labs(title = "Sexual Dimorphism in Energy Over Time",
#'        y = "P-value", x = "Time (ticks)")
#' 
#' # ANOVA with post-hoc tests
#' aov_result <- compare_turtle_groups(
#'   turtle_data,
#'   property = "t.size.cm",
#'   group_var = "breed",
#'   test_type = "anova",
#'   time_point = 2000
#' )
#' 
#' # Post-hoc pairwise comparisons
#' TukeyHSD(aov_result)
#' 
#' # Check assumptions before parametric tests
#' library(car)
#' 
#' # Subset data for testing
#' test_data <- turtle_data %>%
#'   filter(ticks == final_tick, !is.na(t.energy), !is.na(t.sex))
#' 
#' # Test normality
#' shapiro.test(test_data$t.energy[test_data$t.sex == "male"])
#' shapiro.test(test_data$t.energy[test_data$t.sex == "female"])
#' 
#' # Test homogeneity of variance
#' leveneTest(t.energy ~ t.sex, data = test_data)
#' 
#' # Extract specific results
#' t_result <- compare_turtle_groups(
#'   turtle_data, "t.energy", "t.sex",
#'   test_type = "t.test", time_point = final_tick
#' )
#' 
#' cat("Mean difference:", diff(t_result$estimate), "\n")
#' cat("95% CI:", t_result$conf.int, "\n")
#' cat("P-value:", t_result$p.value, "\n")
#' 
#' # Compare multiple properties
#' properties <- c("t.energy", "t.size.cm", "t.stress", "t.age")
#' 
#' multi_test <- lapply(properties, function(prop) {
#'   test <- compare_turtle_groups(
#'     turtle_data, prop, "t.sex",
#'     test_type = "t.test", time_point = final_tick
#'   )
#'   data.frame(
#'     property = prop,
#'     p_value = test$p.value,
#'     significant = test$p.value < 0.05
#'   )
#' })
#' 
#' do.call(rbind, multi_test)
#' 
#' # Life stage comparisons
#' breeds <- c("eggs", "larvae", "juveniles", "adults")
#' turtle_data_subset <- turtle_data %>%
#'   filter(breed %in% breeds, ticks == 3000)
#' 
#' kw_result <- compare_turtle_groups(
#'   turtle_data_subset,
#'   property = "t.stress",
#'   group_var = "breed",
#'   test_type = "kruskal",
#'   time_point = 3000
#' )
#' print(kw_result)
#' }
#'
#' @importFrom stats t.test wilcox.test as.formula aov kruskal.test
#' @export
compare_turtle_groups <- function(turtle_data, property, group_var, 
                                  test_type = "auto", time_point = NULL) {
  # Filter to specific time point if requested
  if (!is.null(time_point)) {
    turtle_data <- turtle_data[turtle_data$ticks == time_point, ]
  }
  
  # Remove missing values
  test_data <- turtle_data[!is.na(turtle_data[[property]]) & !is.na(turtle_data[[group_var]]), ]
  
  # Auto-detect test type if needed
  if (test_type == "auto") {
    n_groups <- length(unique(test_data[[group_var]]))
    if (n_groups == 2) {
      test_type <- "t.test"
    } else {
      test_type <- "anova"
    }
  }
  
  # Perform test
  if (test_type == "t.test") {
    groups <- split(test_data[[property]], test_data[[group_var]])
    if (length(groups) != 2) {
      stop("t.test requires exactly 2 groups")
    }
    result <- t.test(groups[[1]], groups[[2]])
    
  } else if (test_type == "wilcox") {
    groups <- split(test_data[[property]], test_data[[group_var]])
    if (length(groups) != 2) {
      stop("wilcox.test requires exactly 2 groups")
    }
    result <- wilcox.test(groups[[1]], groups[[2]])
    
  } else if (test_type == "anova") {
    formula <- as.formula(paste(property, "~", group_var))
    result <- aov(formula, data = test_data)
    
  } else if (test_type == "kruskal") {
    formula <- as.formula(paste(property, "~", group_var))
    result <- kruskal.test(formula, data = test_data)
  }
  
  return(result)
}