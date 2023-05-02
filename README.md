# ComputationalGenomic_FinalProject

Repo Structure:
- Data:
  - CG11128: gene expression data. Columns: ID (with additional "HS" next to ID) and expression
  - subset_genotype: subset of the genotype table
  - genotype_datasource: resource to download the full genotype data
- association.py: Python code to run association analysis to obtain the effect size estimate and p-value; input needed: subset_genotype and CG11128 (sample files in Data)
- Result:
  - summarystat.csv: sample output for running association.py. this is intended to be used for later PGS calculation


Running Instructions:
- association.py:
  - libraries used: pandas, statsmodels.api, scipy.stats
  - download input data files as described and replace the path to those files in the script. 
- PRScorrelation.py:
  - libraries used: numpy, pandas, scipy.stats, sklearn.metrics, random, warnings
  - download input data files as described and replace the path to those files in the script. 
