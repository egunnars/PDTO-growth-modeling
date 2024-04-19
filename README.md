# PDTO-growth-modeling

This repository holds codes used to produce numerical and graphical results reported in the manuscript “Understanding patient-derived tumor organoid growth through an integrated imaging and mathematical modeling framework”. The codes are written in MATLAB. 

## Reproduce results of model fitting for US-UP-UK data

To reproduce the results of model fitting for the US-UP-UK datasets reported in the manuscript, do the following:
1.	Open the data file “data_us_up_uk.mat”.
2.	Run the code “comparison_growth_models_us_up_uk.m”.
The code “comparison_growth_models_us_up_uk” fits the nine growth models shown in Tables 2 and 3 of the manuscript to the organoid time series data. The fitting procedure is described in Section 2.7 of the Materials and Methods section of the manuscript. The matrices means_bic and means_norm contain the numerical results displayed in Tables 2 and 3 of the manuscript.

## Reproduce figures for US-UP-UK data

To reproduce all figures reporting results of model fitting for the US-UP-UK datasets in the manuscript, do the following:
1.	Open the data file “us_up_uk_analysis_results.mat”, which contains the results of the model fitting described above.
2.	Run the codes “figure_codes_published.m” to reproduce figures in the main text of the manuscript and run the code “figure_codes_supplement_published.m” to reproduce figures in the supplementary text of the manuscript.
