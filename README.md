# Spatially Resolved Molecular Signatures of Lewy Body Dementia and the Impact of APOE4

This repository contains the computational methods used in the study: "Spatially resolved molecular signatures of Lewy body dementia and the Impact of APOE4" by Jin et al. This study investigates the spatial transcriptomic signatures of Lewy body dementia (LBD) with a focus on the impact of APOE4 genotype. Our analysis examines the molecular differences between APOE3 and APOE4 carriers in LBD patients.


## Code overview

- **processAndIntegrate.r**: QC, filtering, normalization and integration of spatial data using Seurat and Harmony.
- **spacexrDeconvolution.r**: Cell type deconvolution using SpaceXR.
- **spacexr_functions.r**: Custom extensions for SpaceXR.
- **figures.R**: Visualization scripts for all manuscript figures.
- **dotplot_IPA.R**: Visualization of Ingenuity Pathway Analysis results.
- **upsetPlots.r**: UpSet plots for differentially expressed gene intersections.
- **stats.R**: Statistical analysis of LB annotations across layers, genotypes, and disease states.

## Data Analysis Workflow

The analysis follows these key steps:
1. Processing and integration of spatial transcriptomics data
2. Cell type deconvolution using snRNA-seq reference data
3. Annotation of spatial spots for LB status and cortical layers
4. Differential expression analysis between APOE genotypes
5. Pathway analysis and visualization

## Contact

Na Zhao (zhao.na@mayo.edu), Mayo Clinic, Department of Neuroscience for more information.
