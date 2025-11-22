# Global variables from NSE (non-standard evaluation)
# This file is used to avoid R CMD check NOTEs about undefined global variables
utils::globalVariables(c(
  # Your existing variables...
  "run_number", "ticks", "sd_value", "mean_value", 
  "t.sex", "t_sex", "n_male", "n_total", 
  "n_female", "xcor", "ycor","xcor_grid", "ycor_grid", "turtle_count", "pxcor", "pycor", "mean_density", 
  "max_density", "Tick", "Locus", "Fst", "Mean_Fst", "SE_Fst", "FmI", 
  "new_additivity_value", "old_additivity_value",
  "t.birth.pond.index", "t_birth_pond_index",
  "t.current.pond.index", "t_current_pond_index",
  "t.size.cm.metamorphosis", "t_size_cm_metamorphosis",
  "t.age.metamorphosis", "t_age_metamorphosis",
  "t.size.cm", "t_size_cm",
  "t.age", "t_age",
  
  # ADD THESE - missing from your R CMD CHECK:
  "breed", "t.energy", "t.energy.consumption", "t.hatch.countdown", 
  "t.max.size", "t.metamorphic.risk", "t.stress", "who", "t.dad.id", 
  "t.mom.id", "parent_birth_pond", "dad_birth_pond", "mom_birth_pond",
  "parent_type", "pond_type", "parent_id", "parent_origin_pond", 
  "origin_pond", "destination_pond", "effective_migrants", "xcor_floor", 
  "ycor_floor", "..cols", "locus", "allele_name", "Population", "Allele", 
  "Frequency", "bin", "pond2", "pond1", "diff_value", "silhouette", 
  "tick", "frequency", "allele", "genotype", "population", "locus_allele", 
  "total_count", "locus.y", "locus.x", "allele.y", "allele.x",
  
  # needed for the example generators
  "patch_energy", "patch_risk", "permeability", "patch_regeneration",
  "new", "index", "form",
  "turtle", "new_allele", ".N",
  
  # ggplot2 aesthetics
  "x", "y", "xend", "yend", "label", "PC1", "PC2", "cluster", "k",
  
  # data.table special symbols
  ":=", "."
))