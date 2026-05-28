# Dengue-phylogeography-South-East-Asia

This repository contains the R script, input files and selected outputs used for the genomic and phylogeographic analysis of a simulated dengue virus outbreak in South-East Asia.

## Repository structure

- `Projet2_Genomic_Analyses.R`: main R script used for data preparation, exploratory analyses, Rt estimation and export of files for external phylogenetic software.
- `data/sequences/`: input metadata and FASTA alignment.
- `data/shapefiles/`: shapefiles used for map visualisation.
- `outputs/`: selected analysis outputs and files generated for external software.

## Main analyses

The workflow includes:

- preliminary assessment of genomic sampling representativeness;
- Rt estimation from reported case data with EpiEstim;
- recombination screening using the Phi test in SplitsTree;
- maximum-likelihood phylogenetic inference using IQ-TREE;
- temporal signal assessment using TempEst;
- Bayesian time-calibrated phylogenetic and discrete phylogeographic inference using BEAST;
- visualisation of inferred spatial spread using Kepler.gl.

## Dynamic visualisation

A short preview of the Kepler.gl dynamic visualisation is shown below:

https://github.com/user-attachments/assets/c3464a74-91ad-465f-b2e1-13bfe6239986

The interactive Kepler.gl visualisation is available in:

- `outputs/kepler.gl.html`
- `outputs/kepler.gl.json`

## Key outputs

Selected outputs include:

- `outputs/Fig1_prephylogeographic_sampling_assessment.pdf`
- `outputs/Fig4_Rt_estimation.pdf`
- `outputs/Rt_estimates.csv`
- `outputs/country_summary.csv`
- `outputs/recombination_phi_test_result.csv`
- `outputs/DENV_ML.iqtree`
- `outputs/DENV_ML.treefile`
- `outputs/Binome_10.log`
- `outputs/MCCtree.svg`
- `outputs/MCCtree.output.csv`

## How to run the R script

Run the script from the repository root folder, which must contain `data/` and `outputs/`:

```r
source("Projet2_Genomic_Analyses.R")
