#assign peak information to each gene
#from supp. table 5 of 2012 nature TCGA breast cancer paper
#july 9 2014

import pandas as pd
import numpy as np

basalamp = pd.read_csv('Supplementary Table 5 revised basal amp.txt',sep = '\t',header=0)
totalamp = pd.read_csv('Supplementary Table 5 revised total amp.txt',sep = '\t',header=0)
basaldel = pd.read_csv('Supplementary Table 5 revised basal del.txt',sep = '\t',header=0) 
totaldel = pd.read_csv('Supplementary Table 5 revised total del.txt',sep = '\t',header=0)

def geneinfo(df,genestart):
    df.dropna(axis =1, subset = [2], inplace = True)
    cytoser = pd.Series()
    chromser = pd.Series()
    widestart = pd.Series()
    wideend = pd.Series()
    genesinseg = pd.Series()

    for cytoband in df.columns.values.tolist()[1:]:
        cytoser = pd.concat([cytoser, pd.Series(cytoband,df[cytoband][df[cytoband].notnull()][genestart:])])
        chromser = pd.concat([chromser,pd.Series(df.ix[2,cytoband].split(':')[0].replace('chr',''),df[cytoband][df[cytoband].notnull()][genestart:])])
        widestart = pd.concat([widestart,pd.Series(df.ix[2,cytoband].split(':')[1].split('-')[0],df[cytoband][df[cytoband].notnull()][genestart:])])
        wideend = pd.concat([wideend,pd.Series(df.ix[2,cytoband].split(':')[1].split('-')[1],df[cytoband][df[cytoband].notnull()][genestart:])])
        genesinseg = pd.concat([genesinseg,pd.Series(df[df[cytoband].notnull()].shape[0]-4,df[cytoband][df[cytoband].notnull()][genestart:])])
        
    genesdf = pd.DataFrame({'cytoband':cytoser,'chromosome':chromser,'wide peak start':widestart,'wide peak end':wideend,'genes in peak':genesinseg},
                           columns=['cytoband','chromosome','wide peak start','wide peak end','genes in peak'])
    genesdf['small seg score'] = 1/genesdf['genes in peak']
    return genesdf

basalgenesamp = geneinfo(basalamp,4)
basalgenesdel = geneinfo(basaldel,4)

totalgenesamp = geneinfo(totalamp,5)
totalgenesdel = geneinfo(totaldel,4)

##print('FILES TO WRITE:')                        
##print(basalgenesamp.head())
##print(basalgenesdel.head())
##print(totalgenesamp.head())
##print(totalgenesdel.head())

basalgenesamp.to_csv('BRCA_basalampsetCNVs.csv')
basalgenesdel.to_csv('BRCA_basaldelsetCNVs.csv')
totalgenesamp.to_csv('BRCA_totalampsetCNVs.csv')
totalgenesdel.to_csv('BRCA_totaldelsetCNVs.csv')



isar = pd.read_csv('SupplementaryTables2-3.txt',header=0,sep='\t') #I checked and this is effectively identical to the published supplementary tables
isar.dropna(how='all',inplace=True)
isar.dropna(how='all',axis=1,inplace=True)
isar.set_index('Gene_symbol',inplace=True)
isar = isar[['Peak','Chromosome','PeakSubtype','HeliosScore']]
isar.columns = ['ISARpeak','chromosome','PeakSubtype','HeliosScore']
isar = isar[~(isar['PeakSubtype']=='Luminal')]
peaks = isar['ISARpeak'].drop_duplicates()
isarnew = pd.DataFrame()
for peak in peaks:
    df = isar[isar['ISARpeak']==peak]
    df['genes in peak'] = df.shape[0]
    isarnew = pd.concat([isarnew,df])

isarnew.to_csv('BRCA_isarsetCNVs.csv')



allpeaks = pd.concat([basalgenesamp[['chromosome']],basalgenesdel[['chromosome']],
                      totalgenesdel[['chromosome']],totalgenesamp[['chromosome']],
                      isarnew[['chromosome']]]) #
allpeaks['gene'] = allpeaks.index
allpeaks['chromosome'] = allpeaks['chromosome'].astype(int)
allpeaks.drop_duplicates(inplace=True)

basal = basalgenesamp[['cytoband','genes in peak']].merge(basalgenesdel[['cytoband',
                                                                         'genes in peak']],
                                                          how='outer',left_index=True,
                                                          right_index=True)
basal.columns = ['TCGA basal amp cytoband','TCGA basal amp genes in peak','TCGA basal del cytoband',
                 'TCGA basal del genes in peak']
total = totalgenesamp[['cytoband','genes in peak']].merge(totalgenesdel[['cytoband',
                                                                         'genes in peak']],
                                                          how = 'outer',left_index=True,
                                                          right_index=True)
total.columns = ['TCGA total amp cytoband','TCGA total amp genes in peak','TCGA total del cytoband',
                 'TCGA total del genes in peak']

TCGA = basal.merge(total,how = 'outer',left_index=True,right_index=True)
allpeaks = allpeaks[['chromosome']].merge(TCGA,how='outer',left_index=True,right_index=True)
isarnew.columns = ['ISARpeak','chromosome','ISARPeakSubtype','HeliosScore','ISAR genes in peak']
allpeaks = allpeaks.merge(isarnew[['ISARpeak','ISARPeakSubtype','HeliosScore','ISAR genes in peak']],
                          how='outer',left_index=True,right_index=True)
allpeaks['Symbol'] = allpeaks.index
allpeaks.drop_duplicates(inplace=True)
##allpeaks.sort(columns=['chromosome','ISARpeak','HeliosScore','TCGA basal amp cytoband',
##                       'TCGA basal del cytoband','TCGA total amp cytoband','TCGA total del cytoband'],
##              ascending=[True,True,False,True,True,True,True,],inplace=True)

cands = pd.read_csv('CNV_candidates_compressed.csv',header=0)
cands.columns = ['Symbol(HUGO)','Symbol','chromosome','from','to']
#cands.set_index('Symbol',drop=False,inplace=True)

ls = list(range(1,24)) + ['X','x']
ls = [str(x) for x in ls]
cands = cands[cands.chromosome.isin(ls)]

cands.chromosome.fillna('',inplace=True)
allpeaks.chromosome.fillna('',inplace=True)
def unique(row):
    if row.chromosome != '':
        if row.chromosome != 'X':
            return str(row.Symbol)+'.'+str(int(row.chromosome))
        else: return str(row.Symbol)+'.'+str(row.chromosome)
    else: return str(row.Symbol)

cands['unID'] = cands.apply(unique,axis=1)
allpeaks['unID'] = allpeaks.apply(unique,axis=1)

cands.set_index('unID',inplace=True)
cands.columns = ['Symbol','Synonym','chrom_candfile','from','to']

allpeaks.set_index('unID',inplace=True)

allpeaks = allpeaks.merge(cands[['Symbol','Synonym','chrom_candfile','from','to']],how='outer',left_index=True,right_index=True)
#allpeaks['chrom_candfile'] = allpeaks['chrom_candfile'].astype(float)
allpeaks['from'] = allpeaks['from'].astype(float)
allpeaks.rename(columns={'Symbol_x':'Symbol'},inplace=True)
allpeaks.set_index('Symbol',inplace=True)
##allpeaks.to_csv('BRCA_allGISTIC2.0andISARgenes.csv')
##allpeaks.drop_duplicates(subset=['chromosome','TCGA basal amp cytoband','TCGA basal amp genes in peak','TCGA basal del cytoband',
##                 'TCGA basal del genes in peak','TCGA total amp cytoband','TCGA total amp genes in peak','TCGA total del cytoband',
##                 'TCGA total del genes in peak','ISARpeak','ISARPeakSubtype','HeliosScore','ISAR genes in peak','chrom_candfile'],
##                         inplace=True)
step6 = pd.read_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_step6.csv',header=0)
step6.set_index('Unnamed: 0',inplace=True)
step7 = step6.merge(allpeaks[['TCGA basal amp cytoband', 'TCGA basal amp genes in peak', 'TCGA basal del cytoband',
                              'TCGA basal del genes in peak', 'TCGA total amp cytoband', 'TCGA total amp genes in peak',
                              'TCGA total del cytoband', 'TCGA total del genes in peak', 'ISARpeak', 'ISARPeakSubtype',
                              'HeliosScore', 'ISAR genes in peak']],how='outer',left_index=True,right_index=True)
##allpeaks = allpeaks.merge(step6[['CNV data exists','Gene in known common CNV',
##                      'Copy number skewed in TN','Has RNA data','RNA differentially expressed',
##                      'CNV type','p-value for RNA t-test','Result','CN-RNA relationship']],
##               how='outer',left_index=True,right_index=True)
step7.rename(columns={'Result':'t-test Result'},inplace=True)
step7 = step7[['Synonym', 'Chromosome', 'from', 'to', 'TCGA basal amp cytoband',
                'TCGA basal amp genes in peak','TCGA basal del cytoband',
                'TCGA basal del genes in peak', 'TCGA total amp cytoband',
                'TCGA total amp genes in peak', 'TCGA total del cytoband',
                'TCGA total del genes in peak', 'ISARpeak', 'ISARPeakSubtype',
                'HeliosScore', 'ISAR genes in peak','CNV data exists',
                'Gene in known common CNV', 'Copy number skewed in TN',
                'Has RNA data','RNA differentially expressed',
                'CNV type','CN cutoff','Percent Altered',
                'p-value for RNA t-test','t-test Result', 'CN-RNA relationship']]
step7[['CNV type','TCGA basal amp cytoband','TCGA basal del cytoband','TCGA total amp cytoband',
       'TCGA total del cytoband']] = step7[['CNV type','TCGA basal amp cytoband',
                                 'TCGA basal del cytoband','TCGA total amp cytoband',
       'TCGA total del cytoband']].replace(np.nan,'')

step7['CNV type consistent with dataset'] = np.where(step7['CNV type']=='','',
        np.where((step7['CNV type']=='amp')&((step7['TCGA basal amp cytoband']!='')|(
                step7['TCGA total amp cytoband']!='')|(~np.isnan(step7.ISARpeak))),
        'yes',np.where((step7['CNV type']=='del')&((step7['TCGA basal del cytoband']!='')|(
                step7['TCGA total del cytoband']!='')),
        'yes','no')))
        
step7.sort_values(['Chromosome','from'],inplace=True)
step7.dropna(how='all',axis=0,inplace=True)
#allpeaks.dropna(subset=['chromosome'],inplace=True)
step7.to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_step7.csv')
