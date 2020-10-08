##examining genes w/in copy number variants in TCGA breast cancer dataset
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

#Nov 21 2014

import pandas as pd
import numpy as np
import math

cands = pd.read_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_step7.csv',header=0) #12906 genes
#cands.rename(columns={'Symbol_x': 'Symbol'},inplace=True)
cands.set_index('Unnamed: 0',inplace=True)

print('has CNV data',cands[cands['CNV data exists']=='yes'].index.unique().shape) #12621 genes

##it's important to note that a lot of genes simply couldn't be put through all tests
##because I couldn't find the RNA data
print('no RNA data',cands[cands['Has RNA data']=='no'].index.unique().shape) #5606

#genes from pan-cancer analyses (mutations and CNVs)
#cosmic
cosmic = pd.read_excel('../COSMIC_cancer_gene_census.xls',sheetname='List')
cosmic = cosmic[['Symbol','Mutation Type']]
cosmic.set_index('Symbol',inplace=True)
cosmic.columns = ['Cosmic']

cands = cands.merge(cosmic,how='left',right_index=True,left_index=True)

#pan-cancer
pan = pd.read_excel('../Pan_Cancer.xlsx',sheetname='Sheet1')
pan = pan[(pan.Type=='METHYLATION')|(pan.Type=='MUTATION')]
pan = pan[['Altered Locus/Gene','Type']]
pan.columns = ['Symbol','Pan-cancer']
pan.set_index('Symbol',inplace=True)

cands = cands.merge(pan,how='left',right_index=True,left_index=True)

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

cands = cands.merge(vog,how='left',right_index=True,left_index=True)

#civic
civ = pd.read_csv('../nightly-GeneSummaries_CiVICdb_160329.tsv',sep='\t',header=0)
civ = civ[['name','description']]
civ.set_index('name',inplace=True)
civ.columns= ['CiVIC']

cands = cands.merge(civ,how='left',right_index=True,left_index=True)

#fix dtypes
cands[['Cosmic','Pan-cancer','Vogelstein','CiVIC']] = cands[[
    'Cosmic','Pan-cancer','Vogelstein','CiVIC']].astype(str)

##start with genes with skewed copy number in TN
#tn = cands[~(np.isnan(cands['Copy number skewed in TN']))]
tn = cands[cands['Copy number skewed in TN']=='yes'] #9708
print('TNBC set',tn.index.dropna().unique().shape)
#tn.set_index('Symbol',inplace=True)

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


#kessler RNAi Myc synthetic lethal screen
kessler = pd.read_csv('kesslerS1.txt',sep = ' ',header=0)
kessler.dropna(how='all',inplace=True)
kessler.dropna(axis=1,how='all',inplace=True)
kessler.drop(['v2hs','median.pair.diffs'],axis=1,inplace=True)
kessler.set_index('symbol',inplace=True)
kessler['MYC SL (Kessler et al)'] = 'yes'

tn = tn.merge(kessler,how='left',right_index=True,left_index=True)

#aure computational study
aure = pd.read_csv('genes_from_Aure2013.txt',header=0)
aure.dropna(how='all',inplace=True)
aure.dropna(axis=1,how='all',inplace=True)
aure['Gene in Aure et al 2013'] = 'yes'
aure.set_index('Gene',inplace=True)

tn = tn.merge(aure,how='left',right_index=True,left_index=True)

#ben neel essential genes in cell lines 
#originally I got these before they published. I only sent him genes in 'tn', not the whole brca set.
#now I am updating to check everything they published.
neel = pd.read_csv('BRCA_essential_genes_Neel.csv', header=0)
basalspec = neel[['Basal']].dropna().set_index('Basal')
breastspec = neel[['Breast-Specific']].dropna().set_index('Breast-Specific')
basalspec['Basal-specific essential (Neel)'] = 'yes'
breastspec['Breast-specific essential (Neel)'] ='yes'

tn = tn.merge(basalspec,how='left',right_index=True,left_index=True)
tn = tn.merge(breastspec, how='left',right_index=True,left_index=True)

#reset index
tn.reset_index(drop=False,inplace=True)
tn.rename(columns ={'index':'Symbol'},inplace=True)


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
tier3 = tier2[~(tier2['MYC SL (Kessler et al)']=='yes')&~(tier2['Gene in Aure et al 2013']=='yes')
              &~(tier2['Basal-specific essential (Neel)']=='yes')&~(
                  tier2['Breast-specific essential (Neel)']=='yes')
              &(tier2.Cosmic=='nan')&(tier2['Pan-cancer']=='nan')&(
                  tier2.Vogelstein=='nan')&(tier2.CiVIC=='nan')]
tier2 = tier2[(tier2['MYC SL (Kessler et al)']=='yes')|(tier2['Gene in Aure et al 2013']=='yes')
              |(tier2['Basal-specific essential (Neel)']=='yes')|(
                  tier2['Breast-specific essential (Neel)']=='yes')
              |~(tier2.Cosmic=='nan')|~(tier2['Pan-cancer']=='nan')|~(
                  tier2.Vogelstein=='nan')|~(tier2.CiVIC=='nan')]
tier2['tier2'] = 2 #35 regions #139 genes
tier3['tier2'] = 3
tier2 = pd.concat([tier2,tier3])
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
tier5 = tier4[~(tier4['MYC SL (Kessler et al)']=='yes')&~(tier4['Gene in Aure et al 2013']=='yes')
              &~(tier4['Basal-specific essential (Neel)']=='yes')&~(
                  tier4['Breast-specific essential (Neel)']=='yes')
              &(tier4.Cosmic=='nan')&(tier4['Pan-cancer']=='nan')&(tier4.Vogelstein=='nan')]
tier4 = tier4[(tier4['MYC SL (Kessler et al)']=='yes')|(tier4['Gene in Aure et al 2013']=='yes')
              |(tier4['Basal-specific essential (Neel)']=='yes')|(
                  tier4['Breast-specific essential (Neel)']=='yes')
              |~(tier4.Cosmic=='nan')|~(tier4['Pan-cancer']=='nan')|~(tier4.Vogelstein=='nan')]

tier5['tier4'] = 5
tier4['tier4'] = 4
tier1point5['tier4'] = 1.5 #16 genes #18 regions
tier4 = pd.concat([tier1point5,tier4,tier5]) #tiers 4-5: 1211 genes 68 regions
tier4.set_index('Symbol',inplace=True)
tn.set_index('Symbol',inplace=True)
tn = tn.merge(tier4[['tier4']],how='outer',left_index=True,right_index=True).drop_duplicates()
tn.reset_index(drop=False,inplace=True)

tn['tier'] = tn[['tier1','tier2','tier4']].min(axis=1)
#tn.drop(['tier1','tier2','tier4'],axis=1,inplace=True)
tn.drop_duplicates(inplace=True)
tn.set_index('Symbol',inplace=True)

#add tier and other information to full file
cands = cands.merge(tn[[
    'MYC SL (Kessler et al)','Gene in Aure et al 2013','Breast-specific essential (Neel)',
                        'Basal-specific essential (Neel)','tier']],how='outer',
                    left_index=True,right_index=True)

##tier X: Any genes NOT picked up by the TN-specific analysis
othertypes = cands[cands['Copy number skewed in TN']=='no']
TN = cands[cands['Copy number skewed in TN']=='yes']
known = othertypes[~(othertypes.Cosmic=='nan')|~(othertypes['Pan-cancer']=='nan')
                              |~(othertypes['Vogelstein']=='nan')|~(othertypes['CiVIC']=='nan')]
known.tier = 10
rest = othertypes[(othertypes.Cosmic=='nan')&(othertypes['Pan-cancer']=='nan')
                              &(othertypes['Vogelstein']=='nan')&(othertypes['CiVIC']=='nan')]
cands = pd.concat([TN,known,rest])
cands.sort_values(['Chromosome','from'],inplace=True)

#fly orthologs
orth = pd.read_csv('New human-fy orthologs library +Diopt.csv',header=0)
#orth.set_index('Human gene',inplace=True)
homs = cands[cands.index.isin(orth['Human gene'])]
nohoms = cands[~cands.index.isin(orth['Human gene'])]
homs['Has fly ortholog'] = 'yes'
nohoms['Has fly ortholog'] = 'no'
cands = pd.concat([homs,nohoms])
cands.sort_values(['Chromosome','from'],inplace=True)

##resent index and export
cands.reset_index(drop=False,inplace=True)
cands.sort_values(['Chromosome','from'],inplace=True)
cands.drop_duplicates().fillna('').to_csv(
    'BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final.csv',index=False)
cands[~cands['index'].str.contains('RN7')&~cands['index'].str.contains(
    'RNU')].drop_duplicates().fillna('').to_csv(
    'BRCA_allGISTIC2.0andISARgenes_folchange_compressed_final_nopseudogenes.csv',index=False)

orth.set_index('Human gene',inplace=True)
homs = homs.merge(orth,how = 'left',left_index=True,right_index=True)
homs.sort_values(['Chromosome','from'],inplace=True)
homs.drop_duplicates(inplace=True)

#ok, now here is a replicate of the 6694 genes that remain after 'selection' step:
#cands[(cands['RNA differentially expressed']=='yes')&(
#        cands['Copy number skewed in TN']=='yes')&~(
#                (cands['CN-RNA relationship']=='-')&(
#                        cands['t-test Result']=='Significant')&(
#                                cands['CNV type']=='del'))]['index'].unique().shape
##here is the equivalent selection from genes with homologs:
#homs[(homs['RNA differentially expressed']=='yes')&(
#        homs['Copy number skewed in TN']=='yes')&~(
#                (homs['CN-RNA relationship']=='-')&(
#                        homs['t-test Result']=='Significant')&(
#                                homs['CNV type']=='del'))].index.unique().shape
#4455
#so this is how many genes I have that passed selection and have orthologs.

one = homs[homs.tier==1]
one[(one['CNV type']=='del')&((
        one['TCGA basal del cytoband'].astype(str)!='nan')|(
                one['TCGA total del cytoband'].astype(str)!='nan'))].to_csv(
                    'BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_tier1_ambiguousdeletions.csv')#12hum,16fly, 5reg
one[(one['CNV type']=='del')&((
        one['TCGA basal del cytoband'].astype(str)=='nan')&(
                one['TCGA total del cytoband'].astype(str)=='nan'))].to_csv(
                    'BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_tier1_paradoxicaldeletions.csv')#119hum,316fly ,27reg
one[(one['CNV type']!='del')].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_tier1_amps.csv')#89hum,225fly,30reg
two = homs[homs.tier==2]
two[two['CNV type']=='amp'].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_tier2.csv')#10hum,20fly,9reg
homs[(homs.tier==3)&(
    homs['CNV type']=='amp')].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_tier3.csv')#32hum,63fly,15reg
homs[homs.tier==4].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_tier4.csv')#191hum,475fly,51reg
homs[homs.tier==5].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_tier5.csv')#807hum,2302fly,61reg
homs[homs.tier==1.5].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_tier1.5.csv')#16hum,37fly,16reg
#amps: 6hum,11fly,5regions  #dels: 10hum,26fly,10regions
homs[homs.tier==10].to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_final_tierX.csv')#65hum,162fly

##here's the new breakdown for the paper
#high priority: tiers 0, 1.5
#intermediate: tiers 1, 2, 4
#low priority -- intersect the remainder of genes with the selection of genes with homologs above:
#remainder = homs[(homs.tier==3)|(homs.tier>4)|np.isnan(homs.tier)]
#remainder.index.dropna().unique().shape
#background = homs[(homs['RNA differentially expressed']=='yes')&(
#        homs['Copy number skewed in TN']=='yes')&~(
#                (homs['CN-RNA relationship']=='-')&(
#                        homs['t-test Result']=='Significant')&(
#                                homs['CNV type']=='del'))]
#background.index.unique().shape
#human genes --> len(set(remainder.index.dropna())&set(background.index.dropna())) #4031
#remainder[remainder.index.isin(background.index)].index.dropna().unique().shape
r = remainder[remainder.index.isin(background.index)]
#low priority regions:
#len(set(r['TCGA basal amp cytoband'].dropna())|set(
#        r['TCGA basal del cytoband'].dropna())|set(
#                r['TCGA total amp cytoband'].dropna())|set(
#                        r['TCGA total del cytoband'].dropna())|set(
#                                r.ISARpeak.dropna()))
##123 regions (ISAR and GISTIC)

#old:
#everything after tier 5: 4344 humn genes, 16684 fly, in 149 regions

###select 5 random genes from tier 5 for negative controls (ish)
##import random
##df = homs[homs.tier==5].reset_index(drop=False)
##for i in range(5):
##    print(df.ix[random.randrange(df.shape[0]),'index'])

#selected genes:
##RDH14
##ARHGEF3
##KIAA1841
##NT5C3L
##VPS33A
