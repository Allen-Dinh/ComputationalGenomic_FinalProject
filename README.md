# ComputationalGenomic_FinalProject

Repo Structure:
- Data:
  - CG11128: gene expression data. Columns: ID (with additional "HS" next to ID) and expression
  - subset_genotype: subset of the genotype table
  - genotype_datasource: resource to download the full genotype data
- association.py: Python code to run association analysis to obtain the effect size estimate and p-value; input needed: subset_genotype and CG11128 (sample files in Data)
