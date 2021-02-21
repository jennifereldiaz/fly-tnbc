##examining genes w/in copy number variants in TCGA breast cancer dataset

#the following tier system was used in previous versions
#final system for group assignment is as described in the manuscript


##have p-values for pairwise comparison with RNA expression
##percentage of TN tumors that have each gene amplified (>1.2) or deleted(<0.8)
##and Helios scores from Dana Pe'er
##also whether the gene is synthetic lethal with Myc and whether it came up in another paper.

##then stratify into tiers:
#tier 1: ISAR genes with Helios scores in the top3 or top quartile of that segment for
        #small segments (=<12) and large segments(>12) respectively
#tier 2: ISAR genes with significant CN-RNA p-value and frequency of
        #alteration in the top 3
        #or top quartile for small (=<12) and large segments(>12) respectively
        #and in aure et al or Myc SL or breast or basal essential or in a pan-cancer analysis
#tier 3:rest of ISAR genes with significant CN-RNA p-value and frequency of
        #alteration in the top 3
        #or top quartile for small (=<12) and large segments(>12) respectively
#tier 4:GISTIC genes with significant CN-RNA p-value and frequency of
        #alteration in the top 3
        #or top quartile for small (=<12) and large segments(>12) respectively
        #and in aure et al or Myc SL or breast or basal essential or in a pan-cancer analysis
#tier 5: rest of the GISTIC genes with significant CN-RNA p-value and
        #frequency of alteration in the top 3
        #or top quartile for small (=<12) and large segments(>12) respectively
#tier 6: rest of the ISAR genes with significant p-value for CN-RNA  
#tier 7: rest of the GISTIC genes with p-value significant after multiple hypothesis correction


#First version Nov 21 2014
#Last updated 9/29/20

import pandas as pd
import numpy as np
import math

cands = pd.read_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_step7.csv',header=0) #12906 genes
#cands.rename(columns={'Symbol_x': 'Symbol'},inplace=True)
cands.set_index('Unnamed: 0',inplace=True)
cands = cands[cands.index.astype(str)!='nan']

print('has CNV data',cands[cands['CNV data exists']=='yes'].index.unique().shape) #12621 genes

##it's important to note that a lot of genes simply couldn't be put through all tests
##because I couldn't find the RNA data
print('no RNA data',cands[cands['Has RNA data']=='no'].index.unique().shape) #5606

##start with genes with skewed copy number in TN
#tn = cands[~(np.isnan(cands['Copy number skewed in TN']))]
tn = cands[cands['Copy number skewed in TN']=='yes'] #9708
print('TNBC set',tn.index.dropna().unique().shape)
#tn.set_index('Symbol',inplace=True)

##can check more of these steps with a line like this 
cands[(cands['Gene in known common CNV']=='no')&(
        cands['CNV data exists']=='yes')&(
                cands['Copy number skewed in TN']=='yes')&(
                        cands['Has RNA data']=='yes')&(
                                cands['RNA differentially expressed']=='no')].index.unique().shape

##get amp or del for all skewed genes
parsed = pd.read_csv('BRCA_TN_CNVs_foldchange_parsed_counts.csv')
parsed.set_index('gene',inplace=True)
parsed = parsed.T
amp = parsed[parsed['percent up of altered']>50]
dl = parsed[parsed['percent down of altered']>50]
amp['CNV type'] ='amp'
amp['Percent altered'] = amp['percent up of total']
dl['CNV type'] = 'del'
dl['Percent altered'] = dl['percent down of total']
parsed = pd.concat([amp,dl])
parsed['Percent altered'] = parsed['Percent altered']*100

tn.drop('CNV type',axis=1,inplace=True)
tn = tn.merge(parsed[['Percent altered','CNV type']],how='outer',left_index=True,right_index=True)

##eliminate genes who are del and whose CN-RNA relationship is '-' (doesn't make any biological sense)
##genes that are amp and '-' could be due to chromatin remodeling/epigenetic silencing
amp = tn[tn['CNV type']=='amp']
dl = tn[tn['CNV type']=='del']
dl = dl[(dl['CN-RNA relationship']=='+')|(dl['t-test Result']=='Non significant')]
tn = pd.concat([amp,dl])
print('after removing - dels',tn.index.dropna().unique().shape)
print('TN set without RNA data',tn[tn['Has RNA data']=='no'].index.dropna().unique().shape)
print('TN set with RNA diff exp',tn[tn['RNA differentially expressed']=='yes'].index.dropna().unique().shape)##6694


#check a bunch of other databases - this time i'm putting it in a separate dataframe
db = cands[['Synonym', 'Chromosome', 'from', 'to']].copy()
#genes from pan-cancer analyses (mutations and CNVs)
#cosmic
cosmic = pd.read_excel('../COSMIC_cancer_gene_census.xls',sheetname='List')
cosmic = cosmic[['Symbol','Mutation Type']]
cosmic.set_index('Symbol',inplace=True)
cosmic.columns = ['Cosmic']

db = db.merge(cosmic,how='left',right_index=True,left_index=True)

#pan-cancer
pan = pd.read_excel('../Pan_Cancer.xlsx',sheetname='Sheet1')
pan = pan[(pan.Type=='METHYLATION')|(pan.Type=='MUTATION')]
pan = pan[['Altered Locus/Gene','Type']]
pan.columns = ['Symbol','Pan-cancer']
pan.set_index('Symbol',inplace=True)

db = db.merge(pan,how='left',right_index=True,left_index=True)

#vogelstein
vogmut = pd.read_excel('../Vogelstein-cancer-genes.xlsx',sheetname='Table S2A',
                       skiprows=1)
vogmut['Type'] = 'Mutation'
vogmut = vogmut[['Gene Symbol','Type']]
vogcnv = pd.read_excel('../Vogelstein-cancer-genes.xlsx',sheetname='Table S2B',
                       skiprows=1)
vogcnv['Type'] = 'CNA'
vogcnv = vogcnv[['Gene Symbol','Type']]
voggerm = pd.read_excel('../Vogelstein-cancer-genes.xlsx',sheetname='Table S4',
                        skiprows=1)
voggerm['Type'] = 'Germline'
voggerm = voggerm[['Gene Symbol','Type']]

vog = pd.concat([vogmut,vogcnv,voggerm])
vog.columns = ['Symbol','Vogelstein']
vog.set_index('Symbol',inplace=True)

db = db.merge(vog,how='left',right_index=True,left_index=True)

#civic
civ = pd.read_csv('../nightly-GeneSummaries_CiVICdb_160329.tsv',sep='\t',header=0)
civ = civ[['name','description']]
civ.set_index('name',inplace=True)
civ.columns= ['CiVIC']

db = db.merge(civ,how='left',right_index=True,left_index=True)

#fix dtypes
db[['Cosmic','Pan-cancer','Vogelstein','CiVIC']] = db[[
    'Cosmic','Pan-cancer','Vogelstein','CiVIC']].astype(str)

#kessler RNAi Myc synthetic lethal screen
kessler = pd.read_csv('kesslerS1.txt',sep = ' ',header=0)
kessler.dropna(how='all',inplace=True)
kessler.dropna(axis=1,how='all',inplace=True)
kessler.drop(['v2hs','median.pair.diffs'],axis=1,inplace=True)
kessler.set_index('symbol',inplace=True)
kessler['MYC SL (Kessler et al)'] = 'yes'

db = db.merge(kessler,how='left',right_index=True,left_index=True)

#aure computational study
aure = pd.read_csv('genes_from_Aure2013.txt',header=0)
aure.dropna(how='all',inplace=True)
aure.dropna(axis=1,how='all',inplace=True)
aure['Gene in Aure et al 2013'] = 'yes'
aure.set_index('Gene',inplace=True)

db = db.merge(aure,how='left',right_index=True,left_index=True)

db1 = db.copy() #this is what I'll use to differentiate tier 4 and 5

#ben neel essential genes in cell lines 
#originally I got these before they published. I only sent him genes in 'tn', not the whole brca set.
#neel = pd.read_csv('BRCA_essential_genes_Neel.csv', header=0)
#basalspec = neel[['Basal']].dropna().set_index('Basal')
#breastspec = neel[['Breast-Specific']].dropna().set_index('Breast-Specific')
#basalspec['Basal-specific essential (Neel)'] = 'yes'
#breastspec['Breast-specific essential (Neel)'] ='yes'
#
#db = db.merge(basalspec,how='left',right_index=True,left_index=True)
#db = db.merge(breastspec, how='left',right_index=True,left_index=True)

#now I'm updating to check everything they published.
for c in [1+i*4 for i in range(0,7)]:
    db['Marcotte_'+pd.read_excel('../Marcotte_Neel_BRCA/mmc4.xlsx',sheetname='S3A',
                                 skiprows=2,nrows=1,usecols=[c]).loc[0,'Unnamed: '+str(c)]+
       '_essential'] = np.where(db.index.isin(
            pd.read_excel('../Marcotte_Neel_BRCA/mmc4.xlsx',sheetname='S3A',
                          skiprows=4,usecols=list(range(c+2)),index_col='Gene').index.dropna()),
                            'yes',np.nan)
for c in [1+i*4 for i in range(0,4)]:
    db['Marcotte_Neve_'+pd.read_excel('../Marcotte_Neel_BRCA/mmc4.xlsx',sheetname='S3B',
                                 skiprows=2,nrows=1,usecols=[c]).loc[0,'Unnamed: '+str(c)]+
       '_essential'] = np.where(db.index.isin(
            pd.read_excel('../Marcotte_Neel_BRCA/mmc4.xlsx',sheetname='S3B',
                          skiprows=4,usecols=list(range(c+2)),index_col='Gene').index.dropna()),
                            'yes',np.nan)
for c in [1+i*4 for i in range(0,6)]:
    db['Marcotte_Lehmann_'+pd.read_excel('../Marcotte_Neel_BRCA/mmc4.xlsx',sheetname='S3C',
                                 skiprows=2,nrows=1,usecols=[c]).loc[0,'Unnamed: '+str(c)]+
       '_essential'] = np.where(db.index.isin(
            pd.read_excel('../Marcotte_Neel_BRCA/mmc4.xlsx',sheetname='S3C',
                          skiprows=4,usecols=list(range(c+2)),index_col='Gene').index.dropna()),
                            'yes',np.nan)
for c in [1+i*4 for i in range(0,10)]:
    db['Marcotte_Curtis_'+pd.read_excel('../Marcotte_Neel_BRCA/mmc4.xlsx',sheetname='S3D',
                                 skiprows=2,nrows=1,usecols=[c]).loc[0,'Unnamed: '+str(c)]+
       '_essential'] = np.where(db.index.isin(
            pd.read_excel('../Marcotte_Neel_BRCA/mmc4.xlsx',sheetname='S3D',
                          skiprows=4,usecols=list(range(c+2)),index_col='Gene').index.dropna()),
                            'yes',np.nan)
for c in [1+i*4 for i in range(0,7)]:
    db['Marcotte_RPPA_'+pd.read_excel('../Marcotte_Neel_BRCA/mmc4.xlsx',sheetname='S3E',
                                 skiprows=2,nrows=1,usecols=[c]).loc[0,'Unnamed: '+str(c)]+
       '_essential'] = np.where(db.index.isin(
            pd.read_excel('../Marcotte_Neel_BRCA/mmc4.xlsx',sheetname='S3E',
                          skiprows=4,usecols=list(range(c+2)),index_col='Gene').index.dropna()),
                            'yes',np.nan)
db['Marcotte_basal_essential'] = np.where(db.index.isin(
            pd.read_excel('../Marcotte_Neel_BRCA/mmc4.xlsx',sheetname='S3F',
                          skiprows=3,usecols=list(range(c+2)),index_col='Gene').index.dropna()),
                            'yes',np.nan)

for c in [1+i*3 for i in range(0,45)]:
    db['Marcotte_METABRIC_'+pd.read_excel('../Marcotte_Neel_BRCA/mmc5.xlsx',sheetname='S4A',
                                 skiprows=1,nrows=1,usecols=[c]).loc[0,'Unnamed: '+str(c)]+
       '_essential'] = np.where(db.index.isin(
            pd.read_excel('../Marcotte_Neel_BRCA/mmc5.xlsx',sheetname='S4A',
                          skiprows=3,usecols=list(range(c+2)),index_col='Gene').index.dropna()),
                            'yes',np.nan)
for c in [1+i*3 for i in range(0,45)]:
    db['Marcotte_ISAR_'+pd.read_excel('../Marcotte_Neel_BRCA/mmc5.xlsx',sheetname='S4B',
                                 skiprows=1,nrows=1,usecols=[c]).loc[0,'Unnamed: '+str(c)]+
       '_essential'] = np.where(db.index.isin(
            pd.read_excel('../Marcotte_Neel_BRCA/mmc5.xlsx',sheetname='S4B',
                          skiprows=3,usecols=list(range(c+2)),index_col='Gene').index.dropna()),
                            'yes',np.nan)

db['Marcotte_decreased_expression_essential'] = np.where(db.index.isin(
            pd.read_excel('../Marcotte_Neel_BRCA/mmc7.xlsx',sheetname='S6A',
                          skiprows=2,index_col='Gene').index.dropna()),
                            'yes',np.nan)
db['Marcotte_increased_expression_essential'] = np.where(db.index.isin(
            pd.read_excel('../Marcotte_Neel_BRCA/mmc7.xlsx',sheetname='S6B',
                          skiprows=3,index_col='Gene').index.dropna()),
                            'yes',np.nan)
db['Marcotte_heterozygous_copy_loss_essential'] = np.where(db.index.isin(
            pd.read_excel('../Marcotte_Neel_BRCA/mmc7.xlsx',sheetname='S6C',
                          skiprows=2,index_col='Gene').index.dropna()),
                            'yes',np.nan)
    
#TN specific functional databases
db['Patel_TN_gene_addictions'] = np.where(db.index.isin(
            pd.read_excel(
                    'patel_gene_addictions/41467_2018_3283_MOESM6_ESM.xlsx',
                    sheetname='Supplementary Data 4')[
                    'Integrated Gene List (n=37)'].dropna().tolist()),
                            'yes',np.nan)

db['Koedoot_TN_Hs578T_migration_drivers'] = np.where(db.index.isin(
            pd.read_excel(
                    'koedoot_migration_drivers/41467_2019_11020_MOESM18_ESM.xlsx',
                    sheetname='HS hits',skiprows=2)[
                    'GeneSymbol'].dropna().tolist()),
                            'yes',np.nan)  
db['Koedoot_TN_MDA231_migration_drivers'] = np.where(db.index.isin(
            pd.read_excel(
                    'koedoot_migration_drivers/41467_2019_11020_MOESM18_ESM.xlsx',
                    sheetname='MDA hits',skiprows=2)[
                    'GeneSymbol'].dropna().tolist()),
                            'yes',np.nan)  

db['Miao_sleeping_beauty_BrWSB'] = np.where(db.index.isin(
            pd.read_excel(
                    'miao_sleeping_beauty/41467_2020_16936_MOESM6_ESM_fixed.xlsx',
                    sheetname='Sheet1',skiprows=2)[
                    'Candidate genes in BrWSB group'].str.upper().dropna().tolist()),
                            'yes',np.nan) 
db['Miao_sleeping_beauty_BrMSB'] = np.where(db.index.isin(
            pd.read_excel(
                    'miao_sleeping_beauty/41467_2020_16936_MOESM7_ESM.xlsx',
                    sheetname='Sheet1',skiprows=2)[
                    'Candidate genes in BrMSB group'].str.upper().dropna().tolist()),
                            'yes',np.nan) 

db[(db[db.columns[28:]]=='yes').any(axis=1)].shape 
#export this and use to validate genes from screen
#db.Chromosome=db.Chromosome.astype(int)
#db['from'] = db['from'].astype(int)
#next time exporting this file, replace column names:
#Cosmic (Forbes et al)	Pan-cancer (Ciriello et al)	Parsons et al	CiVIC (Griffith et al)	MYC SL (Kessler et al)	Aure et al
db.drop_duplicates().fillna('').replace('nan','').sort_values(['Chromosome',
                  'from']).to_csv(
    'BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_TNgenes_databases.csv')
#,index=False


#reset index
#tn.Chromosome=tn.Chromosome.astype(int)
#tn['from'] = tn['from'].astype(int)
#db1.Chromosome=db1.Chromosome.astype(int)
#db1['from'] = db1['from'].astype(int)

db1 = db1[db1.index.isin(tn.index)]
db1.sort_values(['Chromosome','from'],inplace=True)
db1 = db1.fillna('nan')[['Cosmic', 'Pan-cancer',
       'Vogelstein', 'CiVIC', 'MYC SL (Kessler et al)',
       'Gene in Aure et al 2013']]

tn.sort_values(['Chromosome','from'],inplace=True)
#db1.reindex(tn.index)#,inplace=True
tn.reset_index(drop=False,inplace=True)
tn.rename(columns ={'index':'Symbol'},inplace=True)
#db1.reset_index(drop=False,inplace=True)
#db1.rename(columns ={'index':'Symbol'},inplace=True)

def segtest(df,column):
    newdf = pd.DataFrame(columns = list(df.columns.values))
    #df = df[df['CNV type consistent with dataset']=='yes']
    for cytoband in df['cytoband'].unique():
        cytodf = df[df['cytoband'] == cytoband]
        if len(cytodf[column].unique()) > 1:
            #print(cytodf.columns)
            cytodf.reset_index(drop=True,inplace=True)
            cytodf.sort_values([column],axis=0,ascending=False,inplace=True)
#I'm going to take the top genes in each list, based on the assumption that no more than
            #3-4 driver genes are in a segment - don't know if that's true
            if cytodf.shape[0] > 12:
                cutoff = cytodf.ix[int(math.floor(cytodf.shape[0]/4)),column] #take top quartile
                last = cytodf[cytodf[column]==cutoff].drop_duplicates(
                    subset=column,keep='last').index.values[0]
                cytodf = cytodf[:last]
            elif cytodf.shape[0] > 3:
                cutoff = cytodf.ix[3,column]
                last = cytodf[cytodf[column]==cutoff].drop_duplicates(
                    subset=column,keep='last').index.values[0]
                cytodf = cytodf[:last]
        newdf = pd.concat([newdf,cytodf])
        #newdf['Gene'] = newdf.index
    #newdf = newdf[['Gene','CNV type','CN cutoff','CN-RNA relationship','Result','p-value for RNA t-test',
      #             column,'N','chromosome','cytoband',
      #                'genes in peak','Seg p-value','small seg score','wide peak start','wide peak end']]
    #newdf['small seg score'] = 1/newdf['genes in peak']
    #newdf.sort(['small seg score','Seg p-value'],ascending = [False,True],inplace=True)
    #newdf.drop_duplicates(cols=['Gene'],inplace=True)
    return newdf

#stratify genes picked up by isar
isar = tn[~np.isnan(tn.ISARpeak)]
isar.rename(columns={'ISARpeak':'cytoband'},inplace=True)

tier1 = segtest(isar,'HeliosScore')
tier1['tier1'] = 1 #128 genes #37 regions
tier1.set_index('Symbol',inplace=True)
tn.set_index('Symbol',inplace=True)
tn = tn.merge(tier1[['tier1']],how='outer',left_index=True,right_index=True)


tn.reset_index(drop=False,inplace=True)



isar = isar[isar['t-test Result']=='Significant']
tier2 = segtest(isar,'Percent altered')
tier3 = tier2[tier2.Symbol.isin(db1[(db1=='nan').all(axis=1)].index)]
tier2 = tier2[tier2.Symbol.isin(db1[(db1!='nan').any(axis=1)].index)]
tier2['tier2'] = 2 #35 regions #139 genes
tier3['tier2'] = 3
tier2 = pd.concat([tier2,tier3])
tier2.sort_values('tier2').drop_duplicates(
        [col for col in tier2.columns if 'tier2' not in col])
tier2.set_index('Symbol',inplace=True)
tn.set_index('Symbol',inplace=True)
tn = tn.merge(tier2[['tier2']],how='outer',left_index=True,right_index=True).drop_duplicates()
tn.reset_index(drop=False,inplace=True)

#stratify genes picked up by gistic but not isar
gistic = tn[np.isnan(tn['ISARpeak'])|((tn['CNV type']=='del')&(
        tn['TCGA basal del cytoband'].astype(str)!='nan')|(
                tn['TCGA total del cytoband'].astype(str)!='nan'))]
##also eliminiate genes not differentially expressed
gistic = gistic[gistic['RNA differentially expressed']=='yes']
gistic = gistic[gistic['t-test Result']=='Significant']
#and consistent cnv type
gistic = gistic[gistic['CNV type consistent with dataset']=='yes']

basalamp = gistic.rename(columns={'TCGA basal amp cytoband':'cytoband'})
basalamp = segtest(basalamp,'Percent altered')
basalamp.rename(columns={'cytoband':'TCGA basal amp cytoband'},inplace=True)

basaldel = gistic.rename(columns={'TCGA basal del cytoband':'cytoband'})
basaldel = segtest(basaldel,'Percent altered')
basaldel.rename(columns={'cytoband':'TCGA basal del cytoband'},inplace=True)

totalamp = gistic.rename(columns={'TCGA total amp cytoband':'cytoband'})
totalamp = segtest(totalamp,'Percent altered')
totalamp.rename(columns={'cytoband':'TCGA total amp cytoband'},inplace=True)

totaldel = gistic.rename(columns={'TCGA total del cytoband':'cytoband'})
totaldel = segtest(totaldel,'Percent altered')
totaldel.rename(columns={'cytoband':'TCGA total del cytoband'},inplace=True)


tier4 = pd.concat([basalamp,basaldel,totalamp,totaldel]).drop_duplicates()
#tier1.5--small segments
tier1point5 = tier4[(tier4['TCGA basal amp genes in peak']<10)|(
    tier4['TCGA basal del genes in peak']<10)|
                    (tier4['TCGA total amp genes in peak']<10)|(
                        tier4['TCGA total del genes in peak']<10)]
#tier5 = tier4[tier4.Symbol.isin(db1[(db1=='nan').all(axis=1)].index)&~(
#        tier4.Symbol.isin(tier1point5.Symbol))]
tier4 = tier4[tier4.Symbol.isin(db1[(db1!='nan').any(axis=1)].index)&~(
        tier4.Symbol.isin(tier1point5.Symbol))]
#gistic.set_index('Symbol',inplace=True)
tier5 = gistic[~gistic.Symbol.isin(tier1point5.index)&~gistic.Symbol.isin(tier4.index)]

tier5['tier4'] = 5
tier4['tier4'] = 4
tier1point5['tier4'] = 1.5 #16 genes #18 regions
tier4 = pd.concat([tier1point5,tier4,tier5]) #tiers 4-5: 1211 genes 68 regions
tier4.sort_values('tier4').drop_duplicates(
        [col for col in tier4.columns if 'tier4' not in col])
tier4.set_index('Symbol',inplace=True)
tn.set_index('Symbol',inplace=True)
tn = tn.merge(tier4[['tier4']],how='outer',left_index=True,right_index=True).drop_duplicates()
tn.reset_index(drop=False,inplace=True)

tn['tier'] = tn[['tier1','tier2','tier4']].min(axis=1)
#tn.drop(['tier1','tier2','tier4'],axis=1,inplace=True)
tn.drop_duplicates(inplace=True)
tn.set_index('Symbol',inplace=True)

#add tier information to full file
cands = cands.merge(tn[['tier']],how='outer',
                    left_index=True,right_index=True)

###tier X: Any genes NOT picked up by the TN-specific analysis - skip this for now
#othertypes = cands[cands['Copy number skewed in TN']=='no']
#TN = cands[cands['Copy number skewed in TN']=='yes']
#known = othertypes[~(othertypes.Cosmic=='nan')|~(othertypes['Pan-cancer']=='nan')
#                              |~(othertypes['Vogelstein']=='nan')|~(othertypes['CiVIC']=='nan')]
#known.tier = 10
#rest = othertypes[(othertypes.Cosmic=='nan')&(othertypes['Pan-cancer']=='nan')
#                              &(othertypes['Vogelstein']=='nan')&(othertypes['CiVIC']=='nan')]
#cands = pd.concat([TN,known,rest])
cands.sort_values(['Chromosome','from'],inplace=True)

#I selected out 'tier 0' manually based on the Helios scores because tier 1 was too much.
#also selected out some of group 2 manually by looking at the different databases.
#import these
tier0 = pd.concat([pd.read_excel('../../../CNV_screen/final_group1_linelevel.xlsx',
                      sheetname='Group 1I',index_col='Unnamed: 0'),
    pd.read_excel('../../../CNV_screen/final_group1_linelevel.xlsx',
                      sheetname='Group 1G',index_col='Unnamed: 0')])
amb = pd.concat([pd.read_excel('../../../CNV_screen/final_group2_linelevel.xlsx',
                                 sheetname='Group 2G',index_col='Unnamed: 0'),
    pd.read_excel('../../../CNV_screen/final_group2_linelevel.xlsx',
                                 sheetname='Group 2I',index_col='Unnamed: 0')])

#raise SystemExit(0)
##reset index and export
cands['group'] = np.where(cands.index.isin(tier0.index),'1I',np.where(
        cands.index.isin(amb.index)&(
                                 cands['CNV type consistent with dataset']=='yes'),
                        '2IG',np.where(cands.index.isin(amb.index)&(
                                 cands['CNV type consistent with dataset']=='no'),
                        '2I',np.where(cands.tier==1.5,'1G',
                         np.where(cands.tier.isin([1,2])&(
                                 cands['CNV type consistent with dataset']=='yes'),
                                    '3I',np.where(
                                 (cands.tier==4)&(
                                 cands['CNV type consistent with dataset']=='yes'),
                                         '3G',
                                         np.where(~np.isnan(cands.ISARpeak),'4I',
                                                  np.where(cands.tier==5,'4G',''))))))))

#fly orthologs
orth = pd.read_csv('New human-fy orthologs library +Diopt.csv',header=0)
#orth.set_index('Human gene',inplace=True)
homs = cands[cands.index.isin(orth['Human gene'])]
nohoms = cands[~cands.index.isin(orth['Human gene'])]
homs['Has fly ortholog'] = 'yes'
nohoms['Has fly ortholog'] = 'no'
cands = pd.concat([homs,nohoms])
cands.sort_values(['Chromosome','from'],inplace=True)


cands.reset_index(drop=False,inplace=True)
cands.rename(columns={'index':'Symbol'},inplace=True)
cands.sort_values(['Chromosome','from'],inplace=True)
cands.drop_duplicates().fillna('').to_csv(
    'BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_dbupdate.csv',index=False)

cands[~cands['Symbol'].str.contains('RN7')&~cands['Symbol'].str.contains(
    'RNU')].drop_duplicates().fillna('').to_csv(
    'BRCA_allGISTIC2.0andISARgenes_folchange_compressed_final_dbupdate_nopseudogenes.csv',index=False)

orth.set_index('Human gene',inplace=True)
homs = homs.merge(orth,how = 'left',left_index=True,right_index=True)
homs.sort_values(['Chromosome','from'],inplace=True)
homs.drop_duplicates(inplace=True)

#ok, now here is a replicate of the 6694 genes that remain after 'selection' step:
cands[(cands['RNA differentially expressed']=='yes')&(
        cands['Copy number skewed in TN']=='yes')&~(
                (cands['CN-RNA relationship']=='-')&(
                        cands['t-test Result']=='Significant')&(
                                cands['CNV type']=='del'))]['Symbol'].unique().shape
##here is the equivalent selection from genes with homologs:
homs[(homs['RNA differentially expressed']=='yes')&(
        homs['Copy number skewed in TN']=='yes')&~(
                (homs['CN-RNA relationship']=='-')&(
                        homs['t-test Result']=='Significant')&(
                                homs['CNV type']=='del'))].index.unique().shape
#4455
#so this is how many genes I have that passed selection and have orthologs.

#new breakdown by group:
#group 1:
homs[homs.group=='1I'].ISARpeak.unique().shape
homs[homs.group=='1I'].index.unique().shape
#count 1G regions manually
homs[homs.group=='1G'].index.unique().shape
#group 2:
homs[homs.group=='2I'].ISARpeak.unique().shape
homs[homs.group=='2I'].index.unique().shape
homs[homs.group=='2IG'].ISARpeak.unique().shape
homs[homs.group=='2IG'].index.unique().shape
# number genes tested, 2IG:
amb[(amb['CNV type consistent with dataset']=='yes')&~(
        amb['Lethality result'].astype(str)=='nan')].index.unique().shape
# number genes tested, 2I:
amb[(amb['CNV type consistent with dataset']=='no')&~(
        amb['Lethality result'].astype(str)=='nan')].index.unique().shape
# number lines tested, 2IG:
amb[(amb['CNV type consistent with dataset']=='yes')&~(
        amb['Lethality result'].astype(str)=='nan')]['Stock tested'].unique().shape
# number lines tested, 2I:
amb[(amb['CNV type consistent with dataset']=='no')&~(
        amb['Lethality result'].astype(str)=='nan')]['Stock tested'].unique().shape
#group 3:
homs[homs.group=='3I'].ISARpeak.unique().shape
homs[homs.group=='3I'].index.unique().shape
len(set(homs[homs.group=='3G'][['TCGA basal amp cytoband',
       'TCGA basal del cytoband',
       'TCGA total amp cytoband',
       'TCGA total del cytoband']].values.ravel())-set([np.nan]))
homs[homs.group=='3G'].index.unique().shape
#group 4:
homs[homs.group=='4I'].ISARpeak.unique().shape
homs[homs.group=='4I'].index.unique().shape
len(set(homs[homs.group=='4G'][['TCGA basal amp cytoband',
       'TCGA basal del cytoband',
       'TCGA total amp cytoband',
       'TCGA total del cytoband']].values.ravel())-set([np.nan]))
homs[homs.group=='4G'].index.unique().shape

#check for duplicates among the groups
groups = homs[['tier','group']]
groups['Symbol'] = groups.index
groups.reset_index(drop=True,inplace=True)
uniq = groups.drop_duplicates()
dups = uniq[uniq.duplicated('Symbol')]
uniq[uniq.Symbol.isin(dups.Symbol)]

#fill in driver group numbers
driv = pd.read_excel(
        '../../../CNV_screen/BRCA_allGISTIC2.0andISARgenes_compressed_final_summary_drivers.xlsx',
                     sheetname='drivers',
                     index_col=0)
drivers = cands[cands['Symbol'].isin(driv.index)]
drivers.Chromosome = drivers.Chromosome.astype(int)
print(drivers.sort_values(['Chromosome',
            'from'])[['Symbol','group']].drop_duplicates())
print(drivers.sort_values(['Chromosome',
            'from'])[['Symbol','group']].drop_duplicates()[['group']].to_csv(index=False))

#driver neel and tn databases
for driver in drivers['Symbol'].unique():
    print(driver)
    print(', '.join([i.replace('_',' ') for i in db.loc[driver][db.loc[driver].astype(
            str)!='nan'].index.tolist() if not i in ['Synonym', 
            'Chromosome', 'from', 'to']]))
    
    

##old breakdown by tier:
#one = homs[homs.tier==1]
#one[(one['CNV type']=='del')&((
#        one['TCGA basal del cytoband'].astype(str)!='nan')|(
#                one['TCGA total del cytoband'].astype(str)!='nan'))].to_csv(
#                    'BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_dbupdate_tier1_ambiguousdeletions.csv')#12hum,16fly, 5reg
#one[(one['CNV type']=='del')&((
#        one['TCGA basal del cytoband'].astype(str)=='nan')&(
#                one['TCGA total del cytoband'].astype(str)=='nan'))].to_csv(
#                    'BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_dbupdate_tier1_paradoxicaldeletions.csv')#119hum,316fly ,27reg
#one[(one['CNV type']!='del')].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_dbupdate_tier1_amps.csv')#89hum,225fly,30reg
#two = homs[homs.tier==2]
#two[two['CNV type']=='amp'].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_dbupdate_tier2.csv')#10hum,20fly,9reg
#homs[(homs.tier==3)&(
#    homs['CNV type']=='amp')].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_dbupdate_tier3.csv')#32hum,63fly,15reg
#homs[homs.tier==4].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_dbupdate_tier4.csv')#191hum,475fly,51reg
#homs[homs.tier==5].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_dbupdate_tier5.csv')#807hum,2302fly,61reg
#homs[homs.tier==1.5].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_dbupdate_tier1.5.csv')#16hum,37fly,16reg
##amps: 6hum,11fly,5regions  #dels: 10hum,26fly,10regions
##homs[homs.tier==10].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_dbupdate_tierX.csv')#65hum,162fly

##here's the new breakdown for the paper
#high priority: tiers 0, 1.5
#intermediate: tiers 1, 2, 4
#low priority -- intersect the remainder of genes with the selection of genes with homologs above:
remainder = homs[(homs.tier==3)|(homs.tier>4)|np.isnan(homs.tier)]
remainder.index.dropna().unique().shape
background = homs[(homs['RNA differentially expressed']=='yes')&(
        homs['Copy number skewed in TN']=='yes')&~(
                (homs['CN-RNA relationship']=='-')&(
                        homs['t-test Result']=='Significant')&(
                                homs['CNV type']=='del'))]
background.index.unique().shape
#human genes --> len(set(remainder.index.dropna())&set(background.index.dropna())) #4138
remainder[remainder.index.isin(background.index)].index.dropna().unique().shape
r = remainder[remainder.index.isin(background.index)]
#low priority regions:
len(set(r['TCGA basal amp cytoband'].dropna())|set(
        r['TCGA basal del cytoband'].dropna())|set(
                r['TCGA total amp cytoband'].dropna())|set(
                        r['TCGA total del cytoband'].dropna())|set(
                                r.ISARpeak.dropna()))
##125 regions (ISAR and GISTIC)

#old:
#everything after tier 5: 4344 humn genes, 16684 fly, in 149 regions

###select 5 random genes from tier 5 for negative controls (ish)
##import random
##df = homs[homs.tier==5].reset_index(drop=False)
##for i in range(5):
##    print(df.ix[random.randrange(df.shape[0]),'index'])

#selected genes: ##didn't end up using these because I had some tier 5 and lower genes on hand
##RDH14
##ARHGEF3
##KIAA1841
##NT5C3L
##VPS33A
