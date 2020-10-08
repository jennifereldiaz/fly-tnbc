#Aug 2015
#compress CNV data for BRCA dataset for multiple isoforms of genes with different gene coordinates
#genes lost at this stage are those which appear in known common CNVs-removed in the 'no_cnv' files

import pandas as pd
import csv

#read CNV data
print('processing CNVs_genes file...')    
#this will handle the file incorrectly if nothing is named in col 0, row 0. add "TCGA ID" manually if needed.
CNV_genes = pd.read_csv('../BRCA_CNVs_genes_foldchange.csv',index_col=0,dtype={'A1BG':str},header=0,skiprows=[1]) #,skipfooter=1,
CNV_genes.index.rename('TCGA_ID',inplace=True)
print('FILE READ IN:')
print(CNV_genes.head())

#read in candidates file for TN
cands = pd.read_csv('CNV_candidates_compressed.csv',header=0)
print(cands.Symbol.unique().shape[0], 'unique candidate genes') #12621
print(cands.Synonym.unique().shape[0], 'unique candidate genes by synonym') #8920

#check genes that actually have CNV data
sub1 = cands[cands.Symbol.isin(CNV_genes.columns.values)]
rest1 = cands[~(cands.Symbol.isin(CNV_genes.columns.values))]
sub1['CNV data exists'] = 'yes'
rest1['CNV data exists'] = 'no'
print(sub1.Symbol.unique().shape[0], 'candidate genes with CNV data') #12621

#remove those IDs marked as having duplicate CNV data
dupmask = (CNV_genes.iloc[:,0] == 'duplicates')
CNV_genes = CNV_genes[~dupmask]

#convert first column back to floats
CNV_genes.iloc[:,0] = CNV_genes.iloc[:,0].astype(float)

data = CNV_genes.iloc[:,0:-2]

#drop genes that are in known common CNVs
data = data.dropna(axis=1,how='all')
#a few genes get dropped here: ['AFG3L1' 'C16orf55' 'C16orf7' 'C1orf70' 'C8ORFK29' 'CACNA1B' 'CHMP2A'\n 'CN5H6.4' 'DIP2C' 'FAM138E' 'FAM157B' 'HEATR7A' 'HGC6.3' 'LOC646627'\n 'LOC90834' 'METRNL' 'MGC2752' 'MZF1' 'OR4F15' 'OR4F17' 'OR4F4' 'OR4F6'\n 'RPS5' 'SLC27A5' 'TARSL2' 'TM2D3' 'TRIM28' 'TUBB8' 'TUBBP5' 'UBE2M'\n 'WASH3P' 'ZBTB45' 'ZMYND11' 'ZNF132' 'ZNF324' 'ZNF324B' 'ZNF446' 'ZNF584'\n 'ZNF837' 'hsa-mir-1250' 'hsa-mir-1302-10' 'hsa-mir-1302-11' 'hsa-mir-200a'\n 'hsa-mir-3065' 'hsa-mir-3118-1' 'hsa-mir-3118-3' 'hsa-mir-3186'\n 'hsa-mir-338' 'hsa-mir-429' 'hsa-mir-571' 'hsa-mir-657']

print('DATA:')
print(data.head())
print(data.dtypes)

samples = pd.DataFrame(CNV_genes.iloc[:,-2:])
#,columns=['tumor','germ']
data = data.transpose()
data = data.reset_index()
data['index']=[i.split('.')[0] for i in data['index']]

#take mean of isoforms for genes w/ multiple isoforms
data = data.groupby('index').mean().transpose()

#select the genes that are relevant to TN
data = data[list(sub1[sub1.Symbol.isin(data.columns.values)].Symbol.unique())]

#write processed file
CNV_genes = pd.concat([data, samples],axis=1)
print('FILE TO WRITE:')
print(CNV_genes.head())
CNV_genes.to_csv('BRCA_CNVs_genes_foldchange_processed.csv')
print('processed file')

#get genes 
sub2 = sub1[sub1.Symbol.isin(CNV_genes.columns.values)]
rest2 = sub1[~(sub1.Symbol.isin(CNV_genes.columns.values))]
sub2['Gene in known common CNV'] = 'no'
rest2['Gene in known common CNV'] = 'yes'
rest1['Gene in known common CNV'] = ''
cands2 = pd.concat([rest1,sub2,rest2])
print(sub2.Symbol.unique().shape[0], 'unique genes with data in atypical CNVs') #12029
    
cands2.to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_step3.csv',index=False)

#genes lost at this stage are those which don't have data (badly named)
    # or appear in known common CNVs-removed in the 'no_cnv' files
#12029 unique genes with data in atypical CNVs
