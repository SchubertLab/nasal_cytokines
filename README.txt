


###################
###### data #######
###################

Data not included in repository. Please contact the corresponding author.

### raw_data ###

# nasal_cytokines_bl_data2020-07-01.csv: This file is the raw data concerning raw abundances of cytokines as well as clinical variables

# meta.csv: This file is the same as "nasal_cytokines_bl_data2020-07-01.csv". However, this file is used in the downstream analysis, and "nasal_cytokines_bl_data2020-07-01.csv" is just used in the pre-processing step.

# plate_correction: This file allocates the total concentration per cytokine per plate to proceed with the normalisation step.

# ERA-BL-codebook: This file has an explanation of the code names used per each clinical variable as well as a brief explanation of the type of measurement used in each case.

### derivative data ###

# 20220803_plate_normalized20220803_data_donwshift_imputed.csv: This file possesses the normalised + downshift imputation values of the 16/20 cytokines analysed through the paper.

# meta2.csv: This file compiles normalised data + downshift imputation and the clinical data and is the main file used in the downstream analysis.

###################
####### src #######
###################

### preprocessing ###

# 01_normalization.R: this file filters, plate normalizes and down-shift imputes cytokine abundances

### downstream_analysis ###

# 02_asthma_severity_analysis.R: This file contains the code to perform all the regression analyses between lung function variables and asthma severity grades.

# 03_differntial_abundance_analysis.R: This file performes differential abundant analysis between asthma severities and cytokine abundances via Kruskal-Wallis with Wilcoxon rank-sum test as post-hoc analysis. In addition, this file calculates the log2 fold-change per cytokine per comparison between asthma severities.

# 04_multiple_regression_analysis.R: This file  performs multiple regression analysis between cytokines and clinical parameters.

# 05_eosinophilic_asthma_classification_analysis.R:  This file contains the code to calculate the F1 score between the different thresholds regarding peripheral blood eosinophilic counts and the reference (sputum), including a bootstrapping (100 times) to establish the F1 score distribution. Additionally, this file calculates the Wilcoxon rank-sum test between T2 high vs T2 low response according to the selected blood eosinophilic count based on the F1 scores from the previous step; this process was repeated in the bootstrapping process to add a sort of stability selection.

# 06_forrest_plot.py: This script creates the forrest plots of the manuscript based on the provided supplementary tables.

