# Tretow_et_al_2023

This repository comprises all R-scripts and Rmarkdown files used to generate the results and
plots in Tretow et al. 2023.

The analyses are divided into 8 sections. The sections and the R scripts and Rmarkdown files
within the sections must be run in the order given.

The rendered HTML-Reports can be accessed via [the project's GitHub Page](https://mchicken1988.github.io/Tretow_et_al_2023/).

## 1_Quality_assessment

| File type        | File                    |
| ---------------- | ----------------------- |
| ![](png/Rmd.png) | Quality_assessment.Rmd  |

## 2_MAJIQ_AS_analysis

| File type        | File                            |
| ---------------- | ------------------------------- |
| ![](png/R.png)   | MAJIQ_Modulizer_output_to_rds.R |
| ![](png/Rmd.png) | MAJIQ_AS_analysis.Rmd           |

## 3_Dose_response_curve_fitting

| File typ         | File                            |
| ---------------- | ------------------------------- |
| ![](png/Rmd.png) | Dose_response_curve_fitting.Rmd |

## 4_iCLIP_analysis

| File type        | File                                      |
| ---------------- | ----------------------------------------- |
| ![](png/R.png)   | Binding_site_definition.R                 | 
| ![](png/R.png)   | Create_CE_minigenes.R                     |
| ![](png/R.png)   | Calculate_TPMs_in_Controls.R              |
| ![](png/Rmd.png) | Genomic_location_of_binding_sites.Rmd     |
| ![](png/Rmd.png) | Sequence_composition_of_binding_sites.Rmd |
| ![](png/Rmd.png) | RNAmaps_of_regulated_events.Rmd           |
		
## 5_rG4_analysis

| File type        | File                                                |
| ---------------- | --------------------------------------------------- |
| ![](png/R.png)   | rG4_prediction.R                                    |
| ![](png/Rmd.png) | Dependence_of_binding_and_cooperativity_on_rG4s.Rmd |
	
## 6_RTstop profiling

| File type    | File                                                |
| ---------------- | ----------------------------------------------- |
| ![](png/R.png)   | rG4_prediction_on_in_vitro_library_constructs.R |
| ![](png/Rmd.png) | Summary_of_the_in_vitro_library.Rmd             |
| ![](png/Rmd.png) | Comparison_RTstop_profiling_and_predictions.Rmd |

## 7_ivCLIP

| File type    | File                           |
| ---------------- | -------------------------- |
| ![](png/Rmd.png) | in_vitro_CLIP_analysis.Rmd |


## 8_TCGA_analysis

| File type        | File                                                |
| ---------------- | --------------------------------------------------- |
| ![](png/R.png)   | Calculate_expression_BRCA.R                         |
| ![](png/R.png)   | Create_GRanges_object_with_junctions_of_interest.R  |
| ![](png/R.png)   | Calculate_PSIs_BRCA.R                               |
| ![](png/R.png)   | Calculate_adjusted_Rand_Indices.R                   |
| ![](png/Rmd.png) | Correlation_HNRNPH_levels_and_splicing_events.Rmd   |
| ![](png/Rmd.png) | Clustering_BRCA_samples.Rmd                         |
