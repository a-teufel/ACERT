---
title: 'ACERT: An R Package for Analyzing Agent-Based Model Output of Complex Life Cycle Evolution'
tags:
  - R
  - agent-based modeling
  - population genetics
  - life history evolution
  - NetLogo
authors:
  - name: Jesus Lopez
    affiliation: 1
  - name: Robert B. Page
    affiliation: 2
  - name: Ashley I. Teufel
    affiliation: 1
    corresponding: true
    email: ateufel@tamusa.edu
affiliations:
  - name: Texas A&M University - San Antonio, San Antonio, TX, United States 
    index: 1
  - name: Northwestern State University, Natchitoches, LA, United States
    index: 2
date: 26 January 2026
bibliography: paper.bib
---

# Summary

ACERT (Agent-based model of Complex life cycle Evolution - R Tools) is an R package for processing, analyzing, and visualizing output from agent-based models (ABMs) that simulate biological complex life cycles and evolution. The package provides a complete analytical pipeline: from importing large, multi-replicate simulation datasets, through population genetic analysis and spatial landscape characterization, to report generation. By leveraging R's large, preexisting, library of packages related to spatial and genetic analyses, ACERT enables accessible and intuitive output data processing and visualization while supporting exploration of evolutionary, population, and life history dynamics.

# Statement of Need

Many organisms exhibit complex life cycles (CLC) during which they undergo profound, temporally compressed developmental remodeling, known as metamorphosis. Many amphibian species exemplify this pattern: to reach the terrestrial phase, larvae must balance maximizing growth during the larval phase with metamorphosing quickly enough to escape the deteriorating aquatic habitat. Metamorphic timing is generally attributed to complex environmental, physiological, and genetic interactions, making it a challenging, but fascinating, phenomenon to model within an agent-based framework.

ABMs simulating complex life cycles produce large datasets with environmental, physiological, and genetic data on agents across simulated time. A typical simulation tracks thousands of individual agents (called "turtles" in NetLogo terminology) across multiple generations, recording physiological attributes, spatial information, and complete multi-locus genotypes. Simultaneously, the model records environmental conditions across a landscape of discrete habitat patches. When replicated across parameter combinations, researchers face datasets comprising millions of observations distributed across dozens of output files.

Efficiently and comprehensively analyzing such data is important to assessing the evolutionary, population, and life history dynamics that these models produce, as well as how model output varies across different ecological and organismal scenarios. Genetic data must be converted to formats compatible with population genetics software. Individual-level attributes must be summarized across biologically meaningful groupings. Landscape heterogeneity must be quantified and visualized. Variation among replicate runs must be assessed. These tasks are difficult to perform using base R or the standard tidyverse suite of packages alone [@Wickham2019].

ACERT is an R package to analyze the output of CLCs/ABM, an agent-based model that simulates complex life cycles and evolution. This CLCs/ABM, model includes environmental agents called patches and organismal agents called turtles, each with dynamic variables that are tracked and output as the model runs. ACERT allows for accessible and intuitive output data processing and visualization, along with a large preexisting and well-established library of packages related to spatial and genetic analyses. Our package includes functions for importing and cleaning data, moving agent property analysis, population genetic processing, stationary agent visualization, replicate assessment, and miscellaneous other tasks that are encountered when examining model output.

# State of the Field

Several R packages provide population genetics functionality that ACERT builds upon, including `adegenet` [@Jombart2008], `pegas` [@Paradis2010], `hierfstat` [@Goudet2004], and `poppr` [@Kamvar2014]. These packages offer analytical methods but are designed for empirical genetic datasets, not ABM output. Researchers using simulation models must write custom code to parse model output, extract and reformat genetic data, and prepare it for analysis.

General-purpose ABM analysis tools exist but lack domain-specific functionality for evolutionary and population genetic applications. The `nlrx` package [@Salecker2019] facilitates NetLogo experiment design and output collection but does not provide downstream analytical capabilities. The `RNetLogo` package [@Thiele2012] enables R-NetLogo communication but similarly focuses on model control rather than output analysis.

ACERT fills this gap by providing a bridge between ABM output and R's population genetics infrastructure. Rather than duplicating existing analytical methods, ACERT wraps established packages (`adegenet`, `hierfstat`, `mmod`, `poppr`) within functions designed specifically for simulation output, handling data extraction, format conversion, and integration.

![Data processing pipeline in ACERT. Rounded rectangles represent functions, parallelograms represent data objects, and rectangles represent analytical outputs. The three import functions (`read_environments`, `read_turtles`, `read_mutations`) read NetLogo output files and produce corresponding data objects (edat, tdat, mdat). Turtle data (tdat) feeds into multiple analytical pathways: direct outputs for phenotype visualization, sex ratio calculation, and migration analysis; the `process_popgen` function for population genetic processing; and shared outputs with environmental data (edat) for spatial analyses. The `process_popgen` function produces population genetic outputs including allele frequency analyses and genetic subdivision metrics. Mutation data (mdat) contributes to analyses of putatively favored mutations and genetic diversity, the latter also receiving input from `process_popgen` \label{fig:1}](ACERT_pipeline_v3.png)

# Software Design

ACERT's data visualization and statistical analysis features are built to process and extract information from text files containing three different types of information that are recorded by the CLCs/ABM as it runs: (1) environmental data—'edat', (2) organismal data—'tdat', and (3) mutation data—'mdat' (\autoref{fig:1}). This structure mirrors the conceptual organization of agent-based models, where stationary environmental agents (patches) interact with mobile organismal agents (turtles) whose genotypes change through mutation and selection.

The package follows a modular pipeline architecture. The first step in conducting an analysis with ACERT is to import data files, which is done with functions utilizing the `data.table` package to efficiently read in turtle, environment, and mutation datasets while converting genetic data to `adegenet` genind objects or other common population genetics formats [@Jombart2008; @Paradis2010]. Turtle attribute analyses are performed using functions to analyze turtle properties across groupings, visualize turtle properties, and calculate sex ratios. Analysis of turtle population genetics is done via functions that calculate summary statistics, measures of differentiation (e.g., Jost's D, Hedrick's G'st), discriminant analysis of principal components, and analysis of molecular variance [@Adamack2014; @Chessel2004; @Excoffier1992; @Goudet2004; @Jost2008; @Jost2018; @Kamvar2014; @Meirmans2011]. ACERT also provides functions for visualizing and characterizing patches. For example, using ACERT it is straightforward to render raster plots of patch properties and the locations of moving agents on the virtual landscape. Additional functionality includes the ability to assess the similarity between replicate runs, and general utilities such as creating a comprehensive analysis report and saving the analyzed data in various formats.

ACERT was designed to be accessible and straightforward to use by people with a general understanding of base R and the tidyverse, making common analytical workflows straightforward rather than exposing maximum configurability.

# Research Impact Statement

ACERT is an essential companion software to the CLCs/ABM, enabling analyses that would be difficult to perform using base R or the standard tidyverse suite of packages. The package is optimized for handling large datasets to carry out basic exploratory analyses or for conducting more granular analyses of the outcomes of specific replicate runs. As such, ACERT will enable students, wildlife managers, ecologists, and evolutionary biologists to process and interpret output from highly customizable CLCs/ABM software packages.

# AI Usage Disclosure

Generative AI tools (including ChatGPT, Claude, and Gemini) were used during the development of this software for assistance with code writing and debugging. AI tools were also used to assist with editing and drafting portions of this manuscript. All AI-generated content was reviewed, tested, and verified by the authors.

# Acknowledgements

The authors received no specific funding for this work.

# References
