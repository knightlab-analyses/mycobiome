## Mycobiome, microbiome, and immune cell cooccurrence across cancer types

This repository holds the analysis scripts for the multiomic driven comparisons across cancer types.

### Data

##### TCGA data with features intersected with WIS

TCGA data summarized at the species and genus level for fungi and bacteria intersected with the WIS data at those levels.

* genus_intersected_with_WIS

### Code

This directory contains all the code used to generate cooccurrence analysis in TCGA and WIS data. The notebooks are numbered in the order they would need to be run to reproduce the results. Additionally, the `env` directory contains the exact conda environment file that is needed to run the notebooks. 

* 1.0-genus-MMvec-run.ipynb

All MMvec runs for the TCGA data at the genus level intersected with the WIS data.

* 2.0-mmvec-analysis-genus.ipynb

Analysis of the MMvec data inlcuding all plotting and determination of sub-types.

* 2.1-survival-analysis-on-log-ratios.ipynb

Analysis of survival based on above or below medians of the log-ratios between sub-types. 

### Results

* **mmvec-results-genus**

This directory contains the MMvec results for TCGA by center and WIS between bacteria, immune, and fungi pairwise. All tables are generated in notebooks 1.0 and 5.0 

* **tables**

All tables used in the manuscript and raw data used in plots. 

* fungal_bacteria_cooccurrence 
* fungal_immune_cooccurrence 
* log_ratio_comparisons
* log_ratio_survival_analysis

* **figures**

All sub panel plots used in the manuscript. 

* fungal_bacteria_cooccurrence 
* fungal_immune_cooccurrence 
* log_ratio_comparisons
* log_ratio_survival_analysis
