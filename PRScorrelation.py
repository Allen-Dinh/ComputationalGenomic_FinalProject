
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import scipy.stats
from sklearn.metrics import r2_score
import random
import pandas as pd

# Load summary statistics, genotype and expression data

summarystat = pd.read_csv("<path>/summarystat.txt",index_col=0)
geno = pd.read_table("<path>/geno.txt",header=None)
expression = pd.read_table("<path>/HS_CG11128.txt",header=None)

#Dataframe modification

expression.columns = ['ID','expression']
expression['ID'] = expression['ID'].str.replace('.HS','')
expression['ID'] = expression['ID'].astype(int)
expression = expression.set_index('ID')

last_col_index = geno.columns[-1]
geno = geno.drop(last_col_index, axis=1)
geno.columns = ['chr','pos','ID','haplotype_1','haplotype_2','haplotype_3','haplotype_4','haplotype_5','haplotype_6','haplotype_7','haplotype_8']


# Binarize the HMM haplotype at genomic coordinate. As different parts of the chromosome might come from different haplotypes, to reduce LD, only positions where haplotype changed are recorded.

bigeno = pd.DataFrame()
for id in geno['ID'].unique():
  print('filter for sample',id)
  sample = geno[geno['ID']== id]
  info = sample.iloc[:,0:3]
  haplotype = sample.iloc[:,3:]
  haplotype = haplotype.apply(lambda x: x == x.max(), axis=1).astype(int) 
  result = pd.merge(info, haplotype, left_index=True, right_index=True)
  for chr in result['chr'].unique():
    print('loading chr', chr)
    chrom = result[result['chr'] == chr]
    for col in chrom.columns[3:]:
      chrom['hapid'] = chrom[col].diff()
      chr_filtered = chrom[chrom['hapid'] != 0]
      chr_filtered = chr_filtered[chr_filtered[col] ==1]
      chr_filtered = chr_filtered.iloc[:,:-1]
      bigeno= bigeno.append(chr_filtered)
  bigeno.sort_values(by='pos')


# Calculate effect size based on the haplotype

# Modify summary statistics table


betasum = summarystat[['chr', 'pos', 'beta']]
betasum[['haplotype_1','haplotype_2','haplotype_3','haplotype_4','haplotype_5','haplotype_6','haplotype_7','haplotype_8']]= betasum['beta'].str.split(',',expand=True)
betasum = betasum.drop(['beta'], axis=1)
betasum['haplotype_1'] = betasum['haplotype_1'].str.replace('[', '')
betasum['haplotype_8'] = betasum['haplotype_8'].str.replace(']', '')

effectSize = pd.DataFrame()
for id in geno['ID'].unique():
  print('get',id,'data')
  sample = bigeno[bigeno['ID'] == id]
  for chr in bigeno['chr'].unique():
    print('for chr', chr)
    haplpbeta = betasum[betasum['chr'] == chr]
    chrsample = sample[sample['chr'] == chr]
    chrsample = chrsample.reset_index()
    chrsample = chrsample.drop(columns=['index'])
    betaforsample = haplpbeta[haplpbeta['pos'].isin(chrsample['pos'])]
    prs = chrsample.iloc[:,3:].reset_index() * betaforsample.iloc[:,2:].reset_index()
    prs = prs.iloc[:,1:]
    prs = prs.replace('',np.nan)
    prs = prs.astype(float)
    prs['avg'] = prs.apply(lambda row: np.nanmean(row), axis=1)
    prs = prs.iloc[:,-1:]
    mergedPRS = pd.merge(chrsample.iloc[:,:3], prs, left_index=True, right_index=True)
    mergedPRS = mergedPRS.sort_values(by='pos')
    effectSize= effectSize.append(mergedPRS)

effectSize = effectSize.rename(columns={"avg": "effectSize"})


# Sum up all effect sizes of each sample to calculate PRS.

PRS = effectSize.groupby('ID').sum()
PRS.columns =[ 'pos', 'PRS']
PRS = PRS.iloc[:,1]


# Perform filtertion and calculate correlation of PRS with expression

# merge effect size with summary statistics to obtain p-value at each position


effectSizeP = pd.DataFrame()
for chr in effectSize['chr'].unique():
  p = pd.merge(effectSize[effectSize['chr']== chr], summarystat[summarystat['chr'] == chr], on='pos').iloc[:,[0,1,2,3,5]]
  effectSizeP= effectSizeP.append(p)


# used 10 p-value criteria to find best fitting p-value, using 80% of the data

sumPRSpval = pd.DataFrame(columns=['p-val', 'pearsonr','r^2'])
twtdictpval={}
pval = [0.00005,0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05,]
for p in pval:
  effectSizePFiltered = effectSizeP[effectSizeP['p-val'] <= float(p)] 
  filteredPRS = effectSizePFiltered.groupby('ID').sum()
  filteredPRS = filteredPRS.rename(columns={'effectSize': 'p-val-'+str(p)+'effectSize'})
  filteredPRS = filteredPRS.iloc[:,1]
  filteredPRS = pd.merge(filteredPRS,expression,left_index=True, right_index=True )
  etid = random.sample(list(filteredPRS.index.unique()), k=int(len(filteredPRS.index.unique())*0.8))
  ttid = list(set(filteredPRS.index.unique()) - set(etid))
  et = filteredPRS[filteredPRS.index.isin(etid)]
  tt = filteredPRS[filteredPRS.index.isin(ttid)]
  list_row = [p, scipy.stats.pearsonr(et.iloc[:,0],et.iloc[:,1] )[0],r2_score(et.iloc[:,0],et.iloc[:,1] )]
  sumPRSpval.loc[len(sumPRSpval)] = list_row
  name = 'twt'+ str(p)
  print(name)
  print('for p=',p,':',scipy.stats.pearsonr(et.iloc[:,0],et.iloc[:,1] ) )
  print('for p=',p,'R^2=',r2_score(et.iloc[:,0],et.iloc[:,1] ))
  twtdictpval[name] = tt


# Obtain the optimal P value from the 80% data and perform on the rest 20% of the data

optPval = float(sumPRSpval.iloc[(sumPRSpval['r^2']-1).abs().argsort()[:1]]['p-val'])

twpPRSpval = pd.DataFrame(columns=['p-val', 'pearsonr','r^2'])
list_row = [optPval, scipy.stats.pearsonr(twtdictpval['twt' + str(optPval)].iloc[:,0],twtdictpval['twt' + str(optPval)].iloc[:,1] )[0],r2_score(twtdictpval['twt' + str(optPval)].iloc[:,0],twtdictpval['twt' + str(optPval)].iloc[:,1])]
twpPRSpval .loc[len(twpPRSpval)] = list_row
print(scipy.stats.pearsonr(twtdictpval['twt' + str(optPval)].iloc[:,0],twtdictpval['twt' + str(optPval)].iloc[:,1] ))
print(r2_score(twtdictpval['twt' + str(optPval)].iloc[:,0],twtdictpval['twt' + str(optPval)].iloc[:,1] ))
twpPRSpval

# Use 10 top percentage of data ordered by p-value

sumPRStop = pd.DataFrame(columns=['top-percent', 'pearsonr','r^2'])
twtdicttop={}
effectSizesorted = effectSizeP.sort_values(by='p-val')
percentage = [0.0005,0.001,0.005,0.01,0.025,0.05,0.1,0.2,0.5,0.8]
for perc in percentage:
  effectSizeTop = effectSizesorted.head(int(len(effectSizesorted)*perc))
  filteredPRS = effectSizeTop.groupby('ID').sum()
  filteredPRS = filteredPRS.rename(columns={'effectSize': 'top-'+str(perc)+'effectSize'})
  filteredPRS = filteredPRS.iloc[:,1]
  filteredPRS = pd.merge(filteredPRS,expression,left_index=True, right_index=True )
  etid = random.sample(list(filteredPRS.index.unique()), k=int(len(filteredPRS)*0.8))
  ttid = list(set(filteredPRS.index.unique()) - set(etid))
  et = filteredPRS[filteredPRS.index.isin(etid)]
  tt = filteredPRS[filteredPRS.index.isin(ttid)]
  name = 'top-'+ str(perc)
  list_row = [perc, scipy.stats.pearsonr(et.iloc[:,0],et.iloc[:,1] )[0],r2_score(et.iloc[:,0],et.iloc[:,1] )]
  sumPRStop .loc[len(sumPRStop )] = list_row

  print(name)
  print('for top',perc,':',scipy.stats.pearsonr(et.iloc[:,0],et.iloc[:,1] ) )
  print('for top',perc,'R^2=',r2_score(et.iloc[:,0],et.iloc[:,1] ))
  twtdicttop[name] = tt



# Obtain the optimal top percentage from the 80% data and perform on the rest 20% of the data

optTop = float(sumPRStop.iloc[(sumPRStop['r^2']-1).abs().argsort()[:1]]['top-percent'])

twpPRStop = pd.DataFrame(columns=['top-percent', 'pearsonr','r^2'])
list_row = [optTop, scipy.stats.pearsonr(twtdicttop['top-'+str(optTop)].iloc[:,0],twtdicttop['top-'+str(optTop)].iloc[:,1] )[0],r2_score(twtdicttop['top-'+str(optTop)].iloc[:,0],twtdicttop['top-'+str(optTop)].iloc[:,1])]
twpPRStop .loc[len(twpPRStop )] = list_row

print(scipy.stats.pearsonr(twtdicttop['top-'+str(optTop)].iloc[:,0],twtdicttop['top-'+str(optTop)].iloc[:,1] ))
print(r2_score(twtdicttop['top-'+str(optTop)].iloc[:,0],twtdicttop['top-'+str(optTop)].iloc[:,1] ))



