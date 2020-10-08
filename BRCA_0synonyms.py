#Get as large a list as possible of gene synonyms. TCGA lists of GISTIC 2.0
#output, raw CNV data, and RNA data do not share exactly the same gene names
#use synonym list to coordinate.
#11.6.14

import pandas as pd
import itertools

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

def lcombine(row):
    ls = []
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

print('applying synonym list to each synonym...')
syns = pd.concat([syns,lost[['All Symbols']],missing[['All Symbols']]])
print('making dataframe...')
syns = pd.DataFrame(syns['All Symbols'].tolist(),index=syns['All Symbols'])
print('unstacking...')
syns = syns.unstack().dropna()
print('resetting index...')
syns.reset_index(level=0,drop=True,inplace=True)
#print('writing file...')
#syns.to_csv('human_gene_synonym_list.csv')
print('getting missing genes...')
missing = pd.read_csv('BRCA_allGISTIC2.0andISARgenes_compressed_missing.csv',
                      header=0,index_col=0,dtype=str)

##this code is incoporated into the next step-don't need to run
