# H&E-based metaprogram analysis

DeepSpot-M predicts gene expression from H&E images. Melanoma-enriched spots are scored for seven melanoma metaprograms, and their patient-level fractions are compared between M+ and M− tumors. The cohort contains 60 slides from 58 patients.

## Run order

Set `workspace` to your input root and `output` to your output root in `config.json`. Both use the same roots as the other analysis folders; virtual H&E inputs and outputs are under `virtualHE/`. Model and decoder paths are relative to the virtual H&E input folder. Run the Python notebooks from this code folder, with cells in order.

1. `01_generate_predictions.ipynb`: predict expression for the 257-gene panel and calculate the melanoma-marker score.
2. `02_prepare_expression.ipynb`: combine expression from selected melanoma spots into a matrix for scoring.
3. `03_metaprogrammes_figures.Rmd`: calculate metaprogram scores, spot assignments, patient fractions, outcome tests and Figure 5A–E panels.

## Data

The code `inputs/` folder contains the gene panel, aliases and metaprogram controls. Slide and patient metadata are supplied separately under the data root at `virtualHE/inputs/`. Paths inside the slide and expression manifests should be relative to the virtual H&E input folder. Image tokens, predictions, H&E images and model weights are supplied separately. Set the six decoder filenames in `config.json` to the files provided with your data. Outputs are written under the output root at `virtualHE/`.

## Analysis

The melanoma-marker score is the mean predicted expression of ten genes, with a threshold of 0.9042399525642395. Metaprogram scoring uses the 257-gene background, ten expression bins and five control genes per signature gene. Each spot is assigned to its highest-scoring metaprogram.

Spot counts are summed across slides from the same patient. Fractions use all melanoma-enriched spots as the denominator; the MES/nevus-like ratio uses their assigned counts without a pseudocount. Outcome comparisons use two-sided Wilcoxon tests with BH correction across the three measures. ROC confidence intervals use 5,000 outcome-stratified bootstrap resamples.

## Software

Python: DeepSpot-M, PyTorch, NumPy, pandas and h5py; Pillow and zarr for image encoding. R: Seurat, Matrix, jsonlite, ggplot2, dplyr, tidyr, png, scales and knitr.
