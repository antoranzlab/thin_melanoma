# Thin melanoma spatial multi-omics

Analysis notebooks for the thin melanoma study, covering Xenium spatial transcriptomics, MILAN protein measurements, Visium CytAssist and external bulk-transcriptomic survival cohorts.

## Contents and execution order

- **Xenium:** run notebooks 01–06 in order for quality control, broad-lineage annotation, melanoma metaprograms and trajectories, cell subtypes, combined annotations, spatial composition and signature comparisons.
- **MILAN:** run 01, then 02 after Xenium 04. These notebooks match registered cells and compare RNA and protein profiles.
- **Visium:** run 01, then 02 and 04 after Xenium 04. These notebooks cover preprocessing, metaprogram and trajectory projection, gene/pathway analyses and technical validation. The former notebook 03 was integrated into 01–02; the remaining filenames retain their numbering.
- **Public_cohorts:** survival analysis of the nevus-like signature in Leeds and GSE98394.
- **Registration:** summarizes prepared registered coordinates. Raw-image registration with VALIS is not included in this snapshot.

Each notebook lists its inputs and outputs. Shared functions are in `Xenium/helpers.R` and `Visium/helpers.R`.

## Data and paths

This repository contains code only. Raw data, patient metadata, prepared inputs, analysis objects and generated figures are not included. Execution requires the corresponding inputs; the notebooks are not self-contained examples.

Set the root folders in the first setup chunk of each notebook:

```r
data_root <- 'path_to_input'
output_root <- 'path_to_output'
code_root <- 'path_to_code'
```

Replace the placeholders with your input-data folder, output folder and downloaded repository folder. `code_root` is included where helper imports or the knitting working directory require it. Use the same roots across notebooks. Modality subfolders and cross-modality paths are constructed automatically with `file.path()`; generated outputs are separate from the source code.

The input folder should contain `Xenium`, `MILAN`, `Visium`, `Public_cohorts` and `Registration` subfolders, following the inputs described in each notebook. The code folder contains the corresponding notebook/helper folders. Outputs are written beneath the output root in modality subfolders.

Submitted input data use the final descriptive filenames expected by the notebooks; no filename migration is required.

Rebuild the main analysis objects with these notebooks before using their downstream readers. Visium 01 saves a final `visium_analysis.rds` and retains costly SCT and MAGIC checkpoints; Xenium 02 retains its diffusion-map and Slingshot models within `melanoma_analysis.rds`.

## Software

Use R with R Markdown/knitr and the packages loaded in each notebook. Core dependencies include Seurat, tidyverse, Matrix, destiny, slingshot, SingleR, GSVA, tradeSeq, mgcv, limma, speckle and survival. Some steps additionally require the Python environment accessed through reticulate. Prepared cell2location weights are inputs rather than newly fitted by these notebooks. A pinned software environment is not included in this snapshot.

## Verification status

The curated R code has passed syntax checks and selected small-data tests during reorganization. The full revised pipeline has not yet been run end to end on the study data. Code-generated figure layouts may differ from the manually edited manuscript figures.
