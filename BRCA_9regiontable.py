##get matrix of tumors vs cnvs for TN tumors
##nov 13 2015
import pandas as pd
import numpy as np

data = pd.read_csv(
        'BRCA_allGISTIC2.0andISARgenes_folchange_compressed_final_dbupdate_nopseudogenes.csv',
                   header=0,index_col='Symbol')

regions = data[['TCGA basal amp cytoband','TCGA basal del cytoband','TCGA total amp cytoband',
               'TCGA total del cytoband','ISARpeak','CNV type']]
regions.ISARpeak = regions.ISARpeak.astype(str)
repeats = regions[(regions['TCGA basal amp cytoband']==regions['TCGA basal del cytoband'])|
                  (regions['TCGA basal del cytoband']==regions['TCGA total amp cytoband'])|
                  (regions['TCGA basal amp cytoband']==regions['TCGA total amp cytoband'])|
                  (regions['TCGA basal amp cytoband']==regions['TCGA total del cytoband'])|
                  (regions['TCGA basal del cytoband']==regions['TCGA total del cytoband'])|
                  (regions['TCGA total amp cytoband']==regions['TCGA total del cytoband'])]
#regions that show up as both amp and del: 2p25.3,12p13.2,16q12.1
#less than 40 are repeated

amps = pd.Series(regions[['TCGA basal amp cytoband','TCGA total amp cytoband','ISARpeak']].values.ravel()).dropna().unique()
amps = amps[amps != 'nan']
dels = pd.Series(regions[['TCGA basal del cytoband','TCGA total del cytoband']].values.ravel()).dropna().unique()
dels = dels[dels != 'nan']

tumors = pd.read_csv('BRCA_CNVs_TN_filtered2.csv',header=0,index_col='TCGA_ID')

tumorbands = pd.DataFrame(index=tumors.index)


for region in amps:
    genes = pd.Series(regions[(regions['TCGA basal amp cytoband']==region)|
                    (regions['TCGA basal del cytoband']==region)|
                    (regions['TCGA total amp cytoband']==region)|
                    (regions['TCGA total del cytoband']==region)|
                    (regions['ISARpeak']==region)].index.values).dropna()
    genes = genes[genes.isin(tumors.columns)]
    little = tumors[list(genes)]
    #print(region)
    tumorbands[region+'-amp']=little.mean(axis=1)
    #print(tumorbands[region].head())

for region in dels:
    genes = pd.Series(regions[(regions['TCGA basal amp cytoband']==region)|
                    (regions['TCGA basal del cytoband']==region)|
                    (regions['TCGA total amp cytoband']==region)|
                    (regions['TCGA total del cytoband']==region)|
                    (regions['ISARpeak']==region)].index.values).dropna()
    genes = genes[genes.isin(tumors.columns)]
    little = tumors[list(genes)]
    #print(region)
    tumorbands[region+'-del']=little.mean(axis=1)

tumorbands.to_csv('BRCA_TN_CNV_region_values_matrix_update.csv')
