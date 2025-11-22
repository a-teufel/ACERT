
################################################################################
# README.md
################################################################################
# ACERT: Agent-Based Model of Complex Life Cycle Evolution - R Tools

## Overview

ACERT is an R package designed for analyzing outputs from agent-based models of complex life cycle evolution. The package provides tools for:

- Data wrangling and import of ABM output files
- Turtle property analysis across different groupings
- Population genetics calculations and visualizations
- Patch/environment analysis and visualization
- Assessment of replicate run similarity

## Installation

### From GitHub

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install ACERT from GitHub
devtools::install_github("a-teufel/ACERT")
```

### Building from Source

1. Clone the repository:
```bash
git clone https://github.com/a-teufel/ACERT.git
cd ACERT
```

2. Build the package:
```r
# In R, from the parent directory of the package
devtools::build("ACERT")
devtools::install("ACERT")
```

## Quick Start

```r
library(ACERT)

# Read turtle data
turtle_data <- read_turtles("path/to/turtle/files")

# Read environment data
env_data <- read_environments("path/to/environment/files")

# Analyze turtle properties
size_analysis <- analyze_turtle_property(turtle_data, 
                                       property = "t.size.cm.metamorphosis",
                                       grouping = "t.birth.pond.index")

# Create plots
size_plot <- plot_turtle_property(size_analysis, 
                                property = "Size at Metamorphosis",
                                grouping = "t.birth.pond.index",
                                experiment_name = "Experiment 1")

# Generate comprehensive report
create_abm_report(turtle_data, env_data, 
                 output_dir = "results",
                 experiment_name = "Experiment 1")
```

## Main Functions

### Data Import
- `read_turtles()`: Import turtle data files
- `read_environments()`: Import environment data files
- `read_mutations()`: Import mutation data files
- `convert_to_pegas()`: Convert chromosome data to PEGAS format

### Turtle Analysis
- `analyze_turtle_property()`: Analyze any turtle property across groupings
- `plot_turtle_property()`: Create visualizations of turtle properties
- `calculate_sex_ratio()`: Calculate sex ratios for populations

### Population Genetics
- `process_popgen()`: Process population genetics data
- `calculate_popgen_stats()`: Calculate summary statistics
- `calculate_pairwise_fst()`: Calculate pairwise Fst values
- `perform_dapc()`: Discriminant Analysis of Principal Components
- `perform_amova()`: Analysis of Molecular Variance

### Patch Analysis
- `visualize_patches()`: Create patch property visualizations
- `calculate_patch_distances()`: Calculate distances between patches
- `calculate_landscape_complexity()`: Compute landscape metrics
- `calculate_turtle_density()`: Calculate turtle density on patches

### Replicate Assessment
- `assess_replicate_similarity()`: Assess similarity between replicate runs
- `plot_replicate_assessment()`: Visualize replicate clustering results

### Utilities
- `create_abm_report()`: Generate comprehensive analysis report
- `save_abm_data()`: Save processed data in various formats

## Dependencies

The package requires:
- R (>= 4.0.0)
- tidyverse packages (dplyr, ggplot2, tidyr)
- Population genetics packages (adegenet, pegas, hierfstat, PopGenReport)
- Additional visualization and analysis packages

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This package is licensed under GPL-3.

## Citation

If you use ACERT in your research, please cite:

```
Lopez, J., Page, R., and Teufel, A. I. (2025). ACERT: Agent-Based Model of Complex Life Cycle Evolution – R Tools. R package version 0.1.0. https://github.com/a-teufel/ACERT

```

## Contact

For questions or issues, please open an issue on GitHub or contact [ateufel@tamusa.edu]
