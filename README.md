# PDTO-growth-modeling

This repository holds raw images and codes used to produce numerical and graphical results reported in the manuscript “Understanding patient-derived tumor organoid growth through an integrated imaging and mathematical modeling framework”. The codes are written in MATLAB. 

## Image and spreadsheet data

Raw images used for NN image analysis and spreadsheet data containing the results of the analysis are provided in separate folders. The files titled "Count and measure" include information on each organoid in each well, whereas files titled "Tracking" include information on tracked organoids. The column "Object ID" can be used to cross-reference between the files. The column "Track ID" identifies the tracked organoids.

## Reproduce results of model fitting for US-UP-UK data

To reproduce the results of model fitting for the US-UP-UK datasets reported in the manuscript, do the following:
1.	Open the data file “data_us_up_uk.mat”.
2.	Run the code “comparison_growth_models_us_up_uk.m”.
The code “comparison_growth_models_us_up_uk” fits the nine growth models shown in Tables 2 and 3 of the manuscript to the organoid time series data. The fitting procedure is described in Section 2.7 of the Materials and Methods section of the manuscript. The matrices means_bic and means_norm contain the numerical results displayed in Tables 2 and 3 of the manuscript.

## Reproduce figures for US-UP-UK data

To reproduce all figures reporting results of model fitting for the US-UP-UK datasets in the manuscript, do the following:
1.	Open the data file “us_up_uk_analysis_results.mat”.
2.	Run the code “figure_codes_published.m” to reproduce figures in the main text of the manuscript and run the code “figure_codes_supplement_published.m” to reproduce figures in the supplementary text of the manuscript.

## Reproduce results of model fitting for US-GFP data

To reproduce the results of model fitting for the US-GFP dataset reported in the manuscript, do the following:
1.	Open the data file “data_usgfp.mat”.
2.	Run the code “comparison_growth_models_usgfp.m”.

## Reproduce figures for US-GFP data

To reproduce all figures reporting results of model fitting for the US-GFP in the manuscript, do the following:
1.	Open the data file "usgfp_analysis_results.mat”.
2.	Run the code “figure_codes_usgfp.m” to reproduce US-GFP figures in the supplementary text of the manuscript.

## Fit growth models to new data

Two codes are provided to fit the growth models of the paper to a new dataset.

In the first code, "comparison_growth_models_new_data_with_tau.m", the user is expected to supply a data matrix where each row corresponds to time series data for a single organoid. The user is also expected to provide an initial organoid size, which is assumed to be the same for all organoids. In this code, all models include the starting time parameter tau (see the accompanying paper). Further instructions are given in comments within the code.

In the second code, "comparison_growth_models_new_data_without_tau.m", the starting time parameter tau is not included. The user supplies a data matrix where each row corresponds to time series data for a single organoid, and the data in the first column corresponds to the initial size of each organoid. Further instructions are given in comments within the code.
