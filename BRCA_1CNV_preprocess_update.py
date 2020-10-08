#june 2014
#concatenate lists of deletions, amplifications from 2012 BRCA TCGA nature paper and Pe'er unpublished data
#supplementary table 5
#identify hg19 coordinates for each gene]==
import pandas as pd
import numpy as np

print('getting CNV candidates...')
basalamp = pd.read_csv('Supplementary Table 5 revised basal amp.txt',sep = '\t') #,names = list(range(0,20))
basalamp = basalamp[4:].drop('cytoband',axis=1)
totalamp = pd.read_csv('Supplementary Table 5 revised total amp.txt',sep = '\t')
totalamp = totalamp[5:].drop('cytoband',axis=1)
basaldel = pd.read_csv('Supplementary Table 5 revised basal del.txt',sep = '\t') #,names = list(range(0,20))
basaldel = basaldel[4:].drop('cytoband',axis=1)
totaldel = pd.read_csv('Supplementary Table 5 revised total del.txt',sep = '\t')
totaldel = totaldel[4:].drop('cytoband',axis=1)

#TCGA GISTIC 2.0 data
CNVs = pd.concat([basalamp, totalamp, basaldel, totaldel], axis=1)
CNVgenes = CNVs.values.ravel()
#concatenates TCGA data


#pe'erdata
peer = pd.read_csv('SupplementaryTables2-3.txt',sep='\t',header=0)

#exclude peaks from only luminal samples
peer = peer[~(peer['PeakSubtype']=='Luminal')]

genes = np.array(peer['Gene_symbol'])

CNVgenes = list(CNVgenes) + list(genes)

CNVgenes = pd.Series(CNVgenes)


CNVgenes = CNVgenes[~CNVgenes.isnull()]


CNVgenes.drop_duplicates(inplace = True)

#table = pd.read_csv('all_data_by_genes.txt',sep='\t',header=0,index_col=0)

#tab = pd.read_csv('all_data_by_genes.txt',sep='\t',header=0)


print('getting gene info...')
#identify genes with multiple isoforms
genetable19 = pd.read_csv('hgTables19_refseq.txt',sep = '\t')
genetab = genetable19[['name','chrom','txStart','txEnd','name2']]

#genetable18 = pd.read_csv('hgTables18_refseq.txt',sep = '\t')
#genetable38 = pd.read_csv('hgTables38_refseq.txt',sep = '\t')

other = pd.read_csv('hgTables19_otherrefseq.txt',sep = '\t')

genetab.chrom = genetab.chrom.str.split('r').str[1]
genetab = genetab.reindex(columns = ['name2','chrom','txStart','txEnd'])
genetab.columns = ['Symbol','chrom','txStart','txEnd']

ensem = pd.read_csv('Homo_sapiens.GRCh37.75.gtf.gtf',skiprows=5,sep='\t',header=None,dtype=str)
ensem = ensem[ensem[2]=='gene']
genes = pd.DataFrame()
genes['Symbol'] = ensem[8].str.split('gene_name').apply(lambda x: x[1])
genes['Symbol'] = genes['Symbol'].str.split(';').apply(lambda x: x[0])
genes['Symbol'] = genes['Symbol'].str.split('"').apply(lambda x: x[1])
genes['chrom'] = ensem[0]
genes['txStart'] = ensem[3]
genes['txEnd'] = ensem[4]
genes.drop_duplicates(inplace=True)

genes = pd.concat([genetab,genes])
genes.drop_duplicates(inplace=True)

mirna = pd.read_csv('hsa.gff3_miRNAs.txt',sep='\t',skiprows=13,names = range(9))
new = pd.DataFrame()
new['Symbol'] = mirna[8].str.split('Name=').apply(lambda x: x[1])
new['Symbol'] = new['Symbol'].str.split(';').apply(lambda x: x[0])
#new['Symbol'] = new['Symbol'].str.split('"').apply(lambda x: x[1])
new['chrom'] = mirna[0].str.split('r').apply(lambda x: x[1])
new['txStart'] = mirna[3]
new['txEnd'] = mirna[4]
new.drop_duplicates(inplace=True)
##copy = new
##copy['Symbol'] = copy['Symbol'].str.strip('a')
##new = pd.concat([new,copy])
new.drop_duplicates(inplace=True)

genes = pd.concat([new,genes])
genes.columns = ['Symbol','Chromosome','from','to']
#missed gene
miss = pd.read_csv('BRCA_CNV_candidates_missinggenecoord_filled.csv',header=0)
miss.drop('Synonym',axis=1,inplace=True)
add = pd.read_csv('AceView.ncbi_37.gene2chromosome2coordinates.txt.txt',header=0,sep='\t')
add.columns = ['Symbol','Chromosome','from','to']
#add.set_index('Symbol',inplace=True)

genes = pd.concat([genes,add,miss])

genes.drop_duplicates(inplace=True)


#gene synonyms
#gene synonyms
print('getting gene synonyms...')
syns = pd.read_csv('HUGOsynonyms.txt',sep='\t',header=0,dtype=str)
misyns = pd.read_csv('HUGOmirnasymbols.txt',sep='\t',header=0,dtype=str)
pseudosyns = pd.read_csv('HUGOpseudogenesymbols.txt',sep='\t',header=0,dtype=str)
syns = pd.concat([syns,misyns,pseudosyns])
syns = syns[['HGNC ID','Approved Symbol','Approved Name','Previous Symbols','Synonyms']]

ncbi = pd.read_csv('Homo_sapiens.gene_info',sep='\t',nrows=1, names = range(50))
ncbi.dropna(how='all',inplace=True)
ncbi.dropna(how='all',inplace=True,axis=1)

columns = ncbi.ix[0,0].split('(')[0].split()[1:]
ncbi = pd.read_csv('Homo_sapiens.gene_info',sep='\t',names=columns,skiprows=1)
ncbi = ncbi[['dbXrefs','GeneID','Symbol','Synonyms','description','Symbol_from_nomenclature_authority',
             'Full_name_from_nomenclature_authority','Other_designations']]
ncbi.GeneID = ncbi.GeneID.astype(str)
IDs = pd.concat([ncbi.Symbol,ncbi.dbXrefs],axis=1)
IDs.set_index('Symbol',inplace=True)
IDs = pd.DataFrame(IDs.dbXrefs.str.split('|').tolist(),index=IDs.index)
IDs = IDs.unstack().dropna()
IDs.reset_index(level=0,drop=True,inplace=True)
hg = IDs[IDs.str.contains('HGNC:')]
hg = pd.DataFrame(hg)
ensem = IDs[IDs.str.contains('Ensembl:')]
ensem = pd.DataFrame(ensem)
IDs = hg.merge(ensem, how='outer',left_index=True,right_index=True)
IDs.columns = ['HGNC ID','Ensembl ID']
ncbi.set_index('Symbol',drop=False,inplace=True)
ncbi = ncbi.merge(IDs,how='outer',left_index=True,right_index=True)

def combine(row):
    ls = []
    if type(row['Approved Symbol']) is str:
        ls += row['Approved Symbol'].split(',')
    if type(row['Approved Name']) is str:
        ls += row['Approved Name'].split(',')
    if type(row['Synonyms']) is str:
        ls += row['Synonyms'].split(',')
    if type(row['Previous Symbols']) is str:
        ls += row['Previous Symbols'].split(',')
    ls = [x.strip() for x in ls]
    ls = list(set(ls + [x.upper() for x in ls] + [x.lower() for x in ls]))
    return ls

syns['All Symbols'] = syns.apply(combine, axis=1)
syns = syns[['Approved Symbol','All Symbols','HGNC ID']]
syns.columns = ['Symbol','All Symbols','HGNC ID']
syns.set_index('Symbol',inplace=True)

def ncombine(row):
    ls = [row.GeneID,row.Full_name_from_nomenclature_authority,row.Symbol,
          row.Symbol_from_nomenclature_authority]
    if row.Synonyms != '-':
        ls += row.Synonyms.split('|')
    if row.Other_designations != '-':
        ls += row.Other_designations.split('|')
    ls = [str(x).strip() for x in ls]
    ls = list(set(ls + [x.upper() for x in ls] + [x.lower() for x in ls]))
    if '-' in ls: 
        ls.remove('-')
    return ls

ncbi['All Symbols'] = ncbi.apply(ncombine, axis=1)
ncbi = ncbi[['Symbol','All Symbols','HGNC ID']]
ncbi['HGNC ID'] = ncbi['HGNC ID'].str[5:]
ncbi.set_index('Symbol',inplace=True)

lost = syns[~(syns.index.isin(ncbi.index))]
missing = ncbi[~(ncbi.index.isin(syns.index))]

def sym(row):
    if type(row['All Symbols']) is list:
        return row['All Symbols']
    else:
        return [row.name]

syns['All Symbols'] = syns.apply(sym,axis=1)
syns = pd.concat([syns['All Symbols'],pd.DataFrame(syns['All Symbols'].tolist(),index=syns.index)],axis=1)
syns['Symbol'] = syns.index
syns = syns.set_index('All Symbols').unstack().dropna()
syns.reset_index(level=0,drop=True,inplace=True)
syns = pd.DataFrame(syns.index,index=syns)
syns.columns = ['All Symbols Hugo']

ncbi['All Symbols'] = ncbi.apply(sym,axis=1)
ncbi = pd.concat([ncbi['All Symbols'],pd.DataFrame(ncbi['All Symbols'].tolist(),index=ncbi.index)],axis=1)
ncbi['Symbol'] = ncbi.index
ncbi = ncbi.set_index('All Symbols').unstack().dropna()
ncbi.reset_index(level=0,drop=True,inplace=True)
ncbi = pd.DataFrame(ncbi.index,index=ncbi)
ncbi.columns = ['All Symbols NCBI']

syns = syns.merge(ncbi,how='outer',right_index=True,left_index=True)

print('adding missed CNV genes...')
missed = pd.read_csv('BRCA_missingCNVgene_synonyms.csv',header=0)
missed.set_index('name',inplace=True)

syns = syns.merge(missed,how='outer',left_index=True, right_index=True)

lastmissed = pd.read_csv('../BRCA_CNV_candidates_additionalmissingsymbols.csv',header=0)
lastmissed.Synonyms = lastmissed.Synonyms.str.split(',')
lastmissed=lastmissed[['CNVgene','Synonyms']]
lastmissed.set_index('CNVgene',inplace=True)

syns = syns.merge(lastmissed,how='outer',left_index=True,right_index=True)

def lcombine(row):
    ls = [row.name]
    if type(row['Symbol']) is str:
        if row['Symbol'] != 'nan':
            ls.append(row['Symbol'])
    if type(row['Synonyms']) is list:
        ls += row['Synonyms']
    if type(row['All Symbols Hugo']) is list:
        ls += row['All Symbols Hugo']
    if type(row['All Symbols NCBI']) is list:
        ls += row['All Symbols NCBI']
    return list(set(ls))
            
syns['All Symbols'] = syns.apply(lcombine,axis=1)

def ls(row):
    if type(row['All Symbols']) is float:
        return []
    else:
        return row['All Symbols']

##
##def lsplus(row):
##    if row['All Symbols'] is list:
##        ls = row['All Symbols']
##        if row['Symbol'] is str:
##            ls.append(row['Symbol'],row.index)
##    else:
##        ls = [row['Symbol'],row.index]
##    return ls
##

#syns['All Symbols'] = syns.apply(lsplus,axis=1)

print('making dataframe...')
syns.reset_index(inplace=True)

#manually added synonyms
man = pd.read_csv('BRCA_CNV_candidates_missing_from_broadandgenecoord_manuallyfilled.csv',
                  header=0,usecols = [0,'All Symbols'])
man = man[['0','All Symbols']]
man.columns = ['Symbol','All Symbols']
man['All Symbols'] = man['All Symbols'].str.strip('[|]').str.split(',')
man['All Symbols'] = man['All Symbols'].apply(lambda y: [x.strip("'") for x in y])

syns = syns[['index','All Symbols']]
syns.columns = ['Symbol','All Symbols']
syns = pd.concat([syns,man])
syns.set_index('Symbol',inplace=True)

new = pd.DataFrame(syns['All Symbols'].tolist(),index=syns.index)
print('unstacking...')
new = pd.DataFrame(new.unstack().dropna().reset_index(level=0,drop=True))
syns = new.merge(syns[['All Symbols']],how='outer',left_index=True,right_index=True)
print('assigning synonyms to CNVgenes...')
have = syns[syns[0].isin(CNVgenes)]
lost = CNVgenes[~(CNVgenes.isin(syns[0]))]
lost = pd.DataFrame(lost,columns=[0])
lost['All Symbols'] = lost[0].apply(lambda x: [x])
CNVsyns = pd.concat([have,lost])



have = CNVsyns[CNVsyns[0].isin(genes.Symbol)]
lost = CNVsyns[~(CNVsyns[0].isin(genes.Symbol))]
lost.reset_index(drop=True,inplace=True)
print('making dataframe...')
df = pd.DataFrame(lost['All Symbols'].tolist(),index=lost.index)
print('unstacking...')
df = df.unstack().dropna()
df = pd.DataFrame(df)
print('resetting index...')
df.reset_index(level=0,drop=True,inplace=True)
df = df.merge(lost[[0]],how='outer',left_index=True,right_index=True)
print('assigning gene info...')
have.set_index(0,inplace=True)
genes.set_index('Symbol',inplace=True)
have = have.merge(genes,how='left',left_index=True,right_index=True)
have = have[['Chromosome','from','to']]
have['Synonym'] = have.index
print('assigning synonyms to all missing genes...')
df.set_index('0_x',inplace=True)
found = df.merge(genes,how='inner',left_index=True,right_index=True)
found.columns = ['Synonym','Chromosome','from','to']
good = pd.concat([have,found])
good['Symbol'] = good.index
good = good[['Symbol','Synonym','Chromosome','from','to']]
good['from'] = good['from'].astype(int)
good['to'] = good['to'].astype(int)
good.drop_duplicates(inplace=True)

stillmissing_genes = CNVgenes[~(CNVgenes.isin(good.index))&~(CNVgenes.isin(good.Synonym))]
stillmissing_genes = CNVsyns[CNVsyns[0].isin(stillmissing_genes)]
stillmissing_genes.drop_duplicates(subset=0).to_csv('BRCA_CNV_candidates_missinggenecoord.csv')

ls = list(range(1,24)) + ['X','x']
ls = [str(x) for x in ls]
good = good[good.Chromosome.isin(ls)]
good.to_csv('CNV_candidates.csv',index=False)

short = good.drop_duplicates(subset = ['Symbol','Synonym','Chromosome'])
short.to_csv('CNV_candidates_compressed.csv',index=False)

##code not using!
#getting gene coordinates
####get data from broad file = problem is I can't find a lot of details about how this was generated
##good = CNVgenes[CNVgenes.isin(tab['Gene Symbol'])]
##missing = CNVgenes[~(CNVgenes.isin(tab['Gene Symbol']))]
##missing = CNVsyns[CNVsyns[0].isin(missing)]
##missing.reset_index(drop=True,inplace=True)
##df = pd.DataFrame(missing['All Symbols'].tolist(),index=missing.index)
##df = df.unstack().dropna()
##df.reset_index(level=0,drop=True,inplace=True)
##df = pd.DataFrame(df)
##df = df.merge(missing[[0]],how='outer',left_index=True,right_index=True)
##found = df[df['0_x'].isin(tab['Gene Symbol'])]
##found.columns = ['Symbol','Symbol in candidate list']
##good = pd.DataFrame(good)
##good.columns = ['Symbol']
##good['Symbol in candidate list'] = good['Symbol']
##CNVlist = pd.concat([good,found])
##CNVlist.drop_duplicates(inplace=True)
##
##stillmissing_tab =CNVgenes[~(CNVgenes.isin(tab['Gene Symbol']))&~(CNVgenes.isin(found['0_y']))]
##stillmissing_tab = CNVsyns[CNVsyns[0].isin(stillmissing_tab)]
##
##
##both_missed = stillmissing_genes[stillmissing_genes[0].isin(stillmissing_tab[0])]
##both_missed.to_csv('BRCA_CNV_candidates_missing_from_broadandgenecoord.csv')
##
##broad_missed_only = stillmissing_tab[~(stillmissing_tab[0].isin(both_missed[0]))]
##broad_missed_genes = good[good.Synonym.isin(broad_missed_only[0])]
##broad_missed_genes['Symbol'] = broad_missed_genes.index
##broad_missed_genes.drop_duplicates(inplace=True)
##broad_missed_genes = broad_missed_genes[['Synonym','chrom','txEnd','txStart']]
##broad_missed_genes.to_csv('BRCA_CNV_missingfrombroad_genecoord.csv')
##br_missed = broad_missed_genes[['Synonym']]
##br_missed.columns = ['Symbol in candidate list']
##br_missed['Symbol'] = br_missed.index
##
##final_w_syns = pd.concat([CNVlist,br_missed])
##final_w_syns.to_csv('BRCA_CNV_candidates_w_syns.csv',index=False)

####
####print('assigning gene info to all synonyms...')
####genes.set_index('Symbol',inplace=True)
####
####allinfo = genes.merge(syns,how='left',left_index=True,right_index=True)
####
####def ls(row):
####    if type(row['All Symbols']) is float:
####        return []
####    else:
####        return row['All Symbols']
####
####allinfo['All Symbols'] = allinfo.apply(ls,axis=1)
####allinfo['Symbol'] = allinfo.index
####new = pd.concat([allinfo[['Symbol']],pd.DataFrame(allinfo['All Symbols'].tolist(),
####                                                                   index=allinfo.index)],axis=1)
####
####new = new.unstack().dropna()
####new.reset_index(level=0,drop=True,inplace=True)
####new = pd.concat([pd.Series(new.index,index=new.index),new],axis=1)
####new.drop_duplicates(inplace=True)
####allinfo = genes.merge(new,how='outer',left_index=True,right_index=True)
####allinfo.set_index(1,inplace=True)
####
####df = pd.DataFrame(CNVgenes,index=CNVgenes)
####
####CNVcands = df.merge(allinfo,how='left',left_index=True,right_index=True)
####CNVcands.dropna(how='all',axis=1,inplace=True)
#####CNVcands.dropna(how='all',inplace=True)
####CNVcands['chrom'] = CNVcands['chrom'].astype(str)
####CNVcands['txStart'] = CNVcands['txStart'].astype(str)
####CNVcands.sort(['chrom','txStart'],inplace=True)
####
#####some of these genes could not be found among all the synonyms I looked for.
#####about 300 will have to be found manually.
####
#####remove isoforms not mapped to the main assembly
####chroms = list(range(1,23)) + ['X']
####chroms = [str(x) for x in chroms]
####good = CNVcands[CNVcands.chrom.isin(chroms)]
####bad = CNVcands[~(CNVcands.chrom.isin(chroms))]
####keep = bad[~(bad.index.isin(good.index))]
####CNVcands = pd.concat([good,keep])
####CNVcands.to_csv('CNV_candidates.csv')
##
