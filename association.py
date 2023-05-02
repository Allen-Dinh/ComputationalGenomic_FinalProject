import pandas as pd
import statsmodels.api as sm
import scipy.stats as stats

geno = pd.read_csv("<path>/subset_genotype.csv")
pheno = pd.read_table("<path>/CG11128.txt",header=None)

# for the pheno table, rename columns, modify the ID columns, and convert it to numeric type
column_names = ["ID","expression"]
pheno.columns = column_names
pheno['ID'] = pheno['ID'].str.replace('.HS', '')
pheno['ID'] = pd.to_numeric(pheno['ID'])

# Two functions to modify probability

#Code used to modify probability (approach 2 in the report)

def modify_probability(row):
    # Extract the haplotype columns
    haplotype_cols = row[3:]

    # Check if any value is greater than 0.99
    if (haplotype_cols > 0.99).any():
        # Set values greater than 0.99 to 1, and others to 0
        haplotype_cols[haplotype_cols >= 0.99] = 1.0
        haplotype_cols[haplotype_cols < 0.99] = 0.0
    else:
        # Replace values less than 0.01 with 0
        haplotype_cols = haplotype_cols.where(haplotype_cols < 0.01, 0.0)
        haplotype_cols /= haplotype_cols.sum()

    return row

#Code used to modify/binarize genotype (approach 3 in the report)

def binary_genotype(row):
    # Extract the haplotype columns
    haplotype_cols = row[1:9]

    # Find the index of the maximum value
    max_idx = haplotype_cols.idxmax()

    # Set the maximum value to 1, and others to 0
    haplotype_cols[max_idx] = 1.0
    haplotype_cols[haplotype_cols.index != max_idx] = 0.0

    # Ensure row sums up to 1
    haplotype_cols /= haplotype_cols.sum()

    return row

#group by chr and pos and perform linear regression

grouped = geno.groupby(['chr', 'pos'])
haplotype_results = pd.DataFrame(columns=['beta', 'p-val', 'R-squared'])
for (chr_value, pos_value), group in grouped:
    # Extract the genotypes for the current chr+pos combination
    snp = group.iloc[:, 2:]
    merged = pd.merge(snp, pheno, on="ID", how='inner')
    merged = merged.apply(binary_genotype, axis=1)

    # Perform likelihood ratio test
    reduced_model = sm.OLS(merged["expression"],
                           sm.add_constant(merged.iloc[:, :0]))  # Reduced model with only 'Intercept'
    full_model = sm.OLS(merged["expression"],
                        sm.add_constant(merged.iloc[:, 1:9]))  # Full model with all eight variables
    reduced_results = reduced_model.fit()
    full_results = full_model.fit()
    likelihood_ratio = 2 * (full_results.llf - reduced_results.llf)
    df_diff = full_results.df_model - reduced_results.df_model  # Difference in degrees of freedom
    p_value = 1 - stats.chi2.cdf(likelihood_ratio, df=df_diff)  # Aggregate p-value

    coef = full_results.params.tolist()[1:]
    r_squared = full_results.rsquared

    new_row = {'beta': coef, 'p-val': p_value, 'R-squared': r_squared}
    haplotype_results = haplotype_results.append(new_row, ignore_index=True)
    # Print the extracted genotypes
    print(f'chr: {chr_value}, pos: {pos_value}')
    print('---')

haplotype_results.to_csv('<path>/summarystat.csv', index=False)

print("*****************Association Testing Completed*****************")
print(haplotype_results)
