# DCI-after-aSAH-analysis
Code for analysis of DCI complication occurrence after aSAH
# Code for analysis of DCI complication occurrence after aSAH

Associated with manuscript: *[I will fill it later when the title of manuscript will be known, sorry]*

## Overview

This repository contains the R codebase used for the bioinformatic and statistical analyses investigating the occurrence of delayed cerebral ischemia (DCI) complication after aSAH. 

## Workflow & Order of Execution

The analysis relies on a foundational primary script (`00_main_scirpt`) which must be run first in the first place to prepare the workspace, filter data, and initialize the environment. 

Other analysis scripts can be run in any order:

### 1. Core Data Preparation (Run First)
* **`00_main_scirpt.R`** – Primary script handling data import, initial setup, and generation of the base R environment required for all downstream modules.

### 2. Analytical Modules (Run in Any Order)
* **`differences_between_dynamics_of_miRNA_changes_scirpt.R`** – Evaluates temporal changes and dynamics of microRNA expression profiles.
* **`differences_between_groups_day2_and_day6.R`** – Performs comparative statistical analyses between patient groups on Day 2 and Day 6.
* **`heatmap_script.R`** – Generates hierarchical clustering heatmaps to visualize expression patterns.
* **`machine_learning_script.R`** – Implements machine learning models and feature selection.
* **`NormiRazor_script.R`** – Evaluates reference gene stability and normalizers across small RNA-seq datasets.
* **`PCA_script.R`** – Performs Principal Component Analysis for exploratory quality control.
* **`relation_between_gender_and_miRNA_script.R`** – Analyzes potential sex-specific associations with microRNA expression.
* **`top_20_extreme_miRNA_script.R`** – Identifies and visualizes the top 20 most extreme microRNAs.
* **`UMAP_script.R`** – Conducts UMAP dimensionality reduction for sample clustering.
* **`WFNS_and_mFischer_scales_day2_and_day6_script.R`** – Integrates clinical severity scales (WFNS and mFisher) across Day 2 and Day 6 timepoints.

## Data Availability

No data is deposited in this GitHub. All required data available as per the manuscript section (Data availability).
