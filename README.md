# Synergy Score Calculation & Figure Generation (R)

R workflow to calculate drug-combination synergy scores (ZIP/HSA/Bliss/Loewe) using **synergyfinder** and generate manuscript-ready plots (Figure 3A–C and Supplementary Figure 4A–B).

---

## Overview

This repository contains an R script that:

- Loads pre-wrangled, DMSO-normalized drug combination viability data from two experimental batches
- Computes synergy scores per replicate using `synergyfinder`
- Saves ZIP synergy heatmaps for a curated set of combinations
- Aggregates synergy scores across replicates
- Identifies “most synergistic” concentration pairs from an external annotation file
- Generates:
  - A grid of bar plots for selected synergistic concentration pairs (per chemo × cell line)
  - Summary bar plots of mean ZIP synergy across replicates (grouped by drug class)

---

## Repository Structure

Expected directory layout (relative to the project root):

- `data/`
  - `all_cellLines_drug_combo_wrangle_replicate_DMSO_normalized_batch1.RData`
  - `2024-01-08_all_cellLines_drug_combo_wrangle_replicate_DMSO_normalized_batch2.RData`
  - `final_synergistic_combinations_2.txt`
- `scripts/`
  - (place the R script here, if desired)
- `plots/`
  - (auto-created; output PDFs are written here)

---

## Requirements

### R Packages

The script uses:

- `synergyfinder (v3.6.3)`
- `tidyverse (v2.0.0)`
- `patchwork (v1.3.2)`
- `RColorBrewer (v1.1-3)`

Install (example):

```r
install.packages(c("tidyverse", "patchwork", "RColorBrewer"))
# synergyfinder is typically installed from Bioconductor or GitHub depending on your setup
