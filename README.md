# Thin melanoma spatial multi-omics

This project examines how melanoma-cell states and their tissue environments relate to metastatic outcome. Spatial profiling of thin primary melanomas identified seven metaprograms organized along two inferred trajectories, with further evaluation in independent primary-melanoma cohorts.

## Data modalities

- **Xenium:** targeted single-cell spatial RNA measurements.
- **MILAN:** multiplex immunofluorescence protein measurements, matched to Xenium cells from the same section.
- **Visium CytAssist:** broader spatial gene-expression profiles from corresponding sections.
- **H&E:** image-inferred expression for independent histological validation.
- **Bulk transcriptomics:** external primary-melanoma cohorts for survival analyses.

## Repository structure

Notebooks are grouped by modality and numbered in execution order. Each contains its inputs, outputs, analyses and associated figures.

| Script | Purpose |
|---|---|
| [Xenium/01](Xenium/01_Xenium_QC_normalization_level1_clustering.Rmd) | Import, quality control, normalization and broad-lineage annotation. |
| [Xenium/02](Xenium/02_Xenium_metaprogrammes_trajectories.Rmd) | Identify melanoma metaprograms and infer diffusion-map/Slingshot trajectories. |
| [Xenium/03](Xenium/03_Xenium_cell_subtypes.Rmd) | Annotate non-melanoma cell subtypes. |
| [Xenium/04](Xenium/04_Xenium_combined_annotations.Rmd) | Combine annotations and project reassigned melanoma cells. |
| [Xenium/05](Xenium/05_Xenium_spatial_composition.Rmd) | Analyze cell composition, acquisition concordance, neighborhoods and epidermal depth. |
| [Xenium/06](Xenium/06_Xenium_markers_signatures.Rmd) | Identify markers and compare public signatures and the Kunz cohort. |
| [MILAN/01](MILAN/01_MILAN_matching_protein_profiles.Rmd) | Match Xenium–MILAN cells and summarize protein profiles and lineage concordance. |
| [MILAN/02](MILAN/02_MILAN_RNA_protein_spatial.Rmd) | Compare RNA and protein measurements across space and trajectories. |
| [Visium/01](Visium/01_Visium_preprocessing_trajectories.Rmd) | Preprocess/integrate expression and project metaprograms and trajectories. |
| [Visium/02](Visium/02_Visium_gene_pathway_analysis.Rmd) | Compare genes and Hallmark pathways and model pseudotime trends. |
| [Visium/04](Visium/04_Visium_registration_validation.Rmd) | Evaluate registered Xenium pseudospots, imputation and melanoma-spot selection. |
| [Public_cohorts/01](Public_cohorts/01_Public_cohorts_survival.Rmd) | Test nevus-like scores against survival in Leeds and GSE98394. |
| [virtualHE/01](virtualHE/01_generate_predictions.ipynb) | Predict expression from H&E and select melanoma-enriched spots. |
| [virtualHE/02](virtualHE/02_prepare_expression.ipynb) | Combine predicted expression for metaprogram scoring. |
| [virtualHE/03](virtualHE/03_metaprogrammes_figures.Rmd) | Calculate metaprogram fractions, outcome associations and validation figures. |
| [Registration/01](Registration/01_Registration_coordinates.Rmd) | Summarize prepared registered coordinates; does not perform image registration. |
| [Xenium/helpers.R](Xenium/helpers.R), [Visium/helpers.R](Visium/helpers.R) | Functions shared between notebooks. |

Run Xenium 01–06 first, then MILAN 01–02 and Visium 01, 02 and 04. Visium 03 was incorporated into 01–02. Public-cohort analyses use the prepared Xenium signatures. Run virtual H&E notebooks 01–03 in order; their paths are set in `virtualHE/config.json`.

## Running the code

Data are supplied separately. Set `data_root <- 'path_to_input'`, `output_root <- 'path_to_output'` and, where present, `code_root <- 'path_to_code'` in each notebook's setup chunk. All remaining paths are built automatically. Install the R packages loaded in the notebooks; Visium imputation also uses Python through reticulate.

Syntax and selected small-data tests have passed; full execution of the revised pipeline on the study data remains pending.
