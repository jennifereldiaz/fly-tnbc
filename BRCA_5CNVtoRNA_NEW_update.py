#june 2014
#dget RNA data for candidate CNV genes


import csv
import math
import numpy as np
import scipy
from scipy import stats
import matplotlib.pyplot as plt
import math
import itertools
from itertools import zip_longest
import pandas as pd
import timeit


#function to transpose
def transpose(mylist):
    return [list(i) for i in zip(*mylist)]

#function for significant digits
from math import log10, floor
def round_to_2(x):
    digits = -int(floor(log10(x))-1)
    digit_str = '.' + str(digits) + 'f'
    return float(format(x, digit_str))

#function for testing if a string is a number
def isnumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#get filtered gene list
with open('BRCA_CNVs_foldchange_all_filtered.csv', 'r') as filtered:
    filtered = csv.reader(filtered)
    filtered_genelist = next(filtered)
    genelist = list(filtered_genelist)[1:]

#assign synonyms
cands = pd.read_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_step4.csv',header=0)
sym = cands[(cands.Symbol.isin(genelist))|(cands.Synonym.isin(genelist))]
sym = sym[['Symbol','Synonym']]
sym.drop_duplicates(inplace=True)


#assign gene IDs
ncbi = pd.read_csv('Homo_sapiens.gene_info',sep='\t',nrows=1, names = range(50))
ncbi.dropna(how='all',inplace=True)
ncbi.dropna(how='all',inplace=True,axis=1)
columns = ncbi.ix[0,0].split('(')[0].split()[1:]
ncbi = pd.read_csv('Homo_sapiens.gene_info',sep='\t',names=columns,skiprows=1)
ncbi = ncbi[['GeneID','Symbol','Synonyms','Symbol_from_nomenclature_authority','Other_designations']]
ncbi.GeneID = ncbi.GeneID.astype(str)

def ncombine(row):
    s = row.GeneID+'|'+row.Symbol+'|'+row.Symbol_from_nomenclature_authority
    if row.Synonyms != '-':
        s += '|'+row.Synonyms
    if row.Other_designations != '-':
        s += '|'+row.Other_designations
    return s

#ncbi['All Symbols'] = ncbi.apply(ncombine, axis=1)
#ncbi = ncbi[['GeneID','All Symbols']]

def Symb(row,x):
    ls = row['All Symbols'].split('|')
    if x in ls:
        return x

def ID(row):
    if (row.name/100).is_integer():
        print(row.name)
    if '?' not in row.Symbol:
        if ncbi[ncbi['Symbol']==row.Symbol].shape[0] ==1:
            frame = ncbi[ncbi['Symbol']==row.Symbol]
            return frame.reset_index()['GeneID'][0]
        elif ncbi[ncbi['Symbol_from_nomenclature_authority']==row.Symbol].shape[0] ==1:
            frame = ncbi[ncbi['Symbol_from_nomenclature_authority']==row.Symbol]
            return frame.reset_index()['GeneID'][0]
        elif ncbi[ncbi['Synonyms'].str.contains(row.Symbol)].shape[0] ==1:
            frame = ncbi[ncbi['Synonyms'].str.contains(row.Symbol)]
            ls = frame.reset_index()['Synonyms'][0].split('|')
            if row.Symbol in ls:
                return frame.reset_index()['GeneID'][0]
        elif ncbi[ncbi['Other_designations'].str.contains(row.Symbol)].shape[0] ==1:
            frame = ncbi[ncbi['Other_designations'].str.contains(row.Symbol)]
            ls = frame.reset_index()['Other_designations'][0].split('|')
            if row.Symbol in ls:
                return frame.reset_index()['GeneID'][0]

print('assigning geneIDs to synonyms for',sym.shape[0],'lines')
sym['ID'] = sym.apply(ID,axis=1)
ls = list(set(sym.values.ravel()))[1:]
##ls.remove('176')
##ls.remove('746')
##ls.remove('1749')
##ls.remove('10992')
##ls.remove('9410')

##function for dealing with gene synonyms
def Syn(row):
    if sym[sym.Symbol==row.gene].shape[0] > 0:
        return row.gene
    elif sym[sym.Symbol==row.gene_id].shape[0] > 0:
        return row.gene_id
    elif sym[sym.Synonym==row.gene].shape[0] ==1:
        frame = sym[sym.Synonym==row.gene]
        return str(frame.reset_index()['Symbol'][0])
    elif sym[sym.Synonym==row.gene_id].shape[0] ==1:
        frame = sym[sym.Synonym==row.gene_id]
        return str(frame.reset_index()['Symbol'][0])
    elif sym[sym.ID==row.id].shape[0] ==1:
        frame = sym[sym.ID==row.id]
        return frame.reset_index()['Symbol'][0]
    else: return ''


#generate summary files of RNAseq data of interest for all tumors by original barcode
#def getRNApandas():
print('getting RNASeq values for gene list with pandas...')
with open('../BRCA_pathology_2014.csv', 'r') as path_file:
    path = list(csv.reader(path_file))
    
filemap = pd.read_csv('../../RNASEQ/FILE_SAMPLE_MAP.txt','\t',header=0,dtype=str)
filemap['tumor-normal'] = filemap['barcode(s)'].str[13:14]
filemap['barcode(s)'] = filemap['barcode(s)'].str[:12]
filemap = filemap[filemap['tumor-normal']=='0']
filemap = filemap[filemap['filename'].str.contains('unc.edu')]
filemap = filemap[filemap['filename'].str.contains('genes.normalized_results')]

df = pd.DataFrame(index=ls,dtype=float)
with open('../../RNASEQ/FILE_SAMPLE_MAP.txt', 'r') as file_map:
    file_map = list(csv.reader(file_map, delimiter='\t'))
for tumorID in path[1:]:
    if float((len(path[1:]) - path[1:].index(tumorID))/50).is_integer(): 
        print(str(len(path[1:]) - path[1:].index(tumorID)) + ' ' + 'tumors left')
    if len(tumorID) > 0:
        tumor_row = filemap[filemap['barcode(s)']==tumorID[0]]
        if tumor_row.shape[0] == 1:
            file_name = '../../RNASEQ/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/' + tumor_row['filename'].max()
            if tumor_row['barcode(s)'].max() in df.columns.values:
                df.drop(tumor_row['barcode(s)'].max(),axis=1,inplace=True)
            else:
                sample_RNA = pd.read_csv(file_name,header=0,sep='\t')
                sample_RNA['gene'] = sample_RNA['gene_id'].str.split('|').apply(lambda x: x[0])
                sample_RNA['id'] = sample_RNA['gene_id'].str.split('|').apply(lambda x: x[1])
                sample_RNA = sample_RNA[(sample_RNA['gene'].isin(ls))|(sample_RNA['gene_id'].isin(ls))|(sample_RNA['id'].isin(ls))]
                sample_RNA.drop_duplicates(subset='gene_id',inplace=True)
                sample_RNA['symbol'] = sample_RNA.apply(Syn,axis=1)
                sample_RNA.symbol = sample_RNA.symbol.astype(str)
                sample_RNA.set_index('symbol',inplace=True)
                sample_RNA = sample_RNA[sample_RNA.index.isin(genelist)]
                #print(sample_RNA.head())
                #for now drop gene names that are duplicated at this point--there are very few
                ser = pd.Series(sample_RNA.index)
                sample_RNA.drop(ser[ser.duplicated()].tolist(),axis=0,inplace=True)
                df[tumorID[0]] = sample_RNA['normalized_count']
                #print(df.head())
                
df.dropna(how='all',inplace=True)
df.dropna(how='all',axis=1,inplace=True)
df = df.T
df.insert(0,'Complete TCGA ID',df.index.values)
df.to_csv('BRCA_RNA_candidates.csv',index=False)
print('created RNA file')


#convert to z-scores
def zscores():
    with open('BRCA_RNA_candidates.csv', 'r') as RNA:
        RNA = csv.reader(RNA)
        RNA = list(RNA)
    RNA_tr = transpose(RNA)
    z_list_tr = []
    z_list_tr.append(RNA_tr[0])
    for cand in range(1,len(RNA[0])):
        #print(RNA[0][cand])
        RNA_list = []
        for i in RNA_tr[cand]:
            if isnumber(i):
                RNA_list.append(float(i))
        normal = scipy.stats.normaltest(RNA_list)
        z_array = scipy.stats.zscore(RNA_list)
        z_list_cand = list(z_array)
        z_list_cand.insert(0, RNA[0][cand])
        z_list_tr.append(z_list_cand)
    z_list = transpose(z_list_tr)
    with open('BRCA_RNA_z_scores.csv','w+') as z_scores:
        z_scores = csv.writer(z_scores)
        for line in z_list:
            z_scores.writerow(line)
    print('created z-scores file')

###RUN
#getRNApandas()
zscores()

cands = pd.read_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_step4.csv',header=0)

with open('BRCA_RNA_z_scores.csv','r') as z:
    z = csv.reader(z)
    z = next(z)
    z = list(z)[1:]

sub = cands[cands.Symbol.isin(z)]
rest = cands[~(cands.Symbol.isin(z))]
sub['Has RNA data'] = 'yes'
rest['Has RNA data'] = 'no'

cands = pd.concat([sub,rest])
cands.to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_step5.csv',index=False)


######    RUN AND TIME    #########
##getpandas = timeit.Timer(stmt='getRNApandas()',setup = 'from __main__ import getRNApandas')
##print('time to get values:',getpandas.timeit(number=1))
##
##scores = timeit.Timer(stmt = 'zscores()',setup='from __main__ import zscores')
##print('time to make scores:',scores.timeit(number=1))







##this function is slow. don't use!
#generate summary files of RNAseq data of interest for all tumors by original barcode
def getRNA():
    print('getting RNASeq values for gene list...')
    with open('BRCA_pathology.csv', 'r') as path_file:
        path = list(csv.reader(path_file))
    with open('BRCA_RNA.csv', 'w+') as RNA_cand:
        RNA_cand = csv.writer(RNA_cand)
        RNA_cand.writerow(['Complete TCGA ID'] + genelist)
        with open('RNASeqV2/FILE_SAMPLE_MAP.txt', 'r') as file_map:
            file_map = list(csv.reader(file_map, delimiter='\t'))
        sample_RNA_list = []
        for tumorID in path[1:]:
            if len(sample_RNA_list) > 1 and not 'duplicate' in sample_RNA_list:
                RNA_cand.writerow(sample_RNA_list)
            if float((len(path[1:]) - path[1:].index(tumorID))/50).is_integer(): 
                print(str(len(path[1:]) - path[1:].index(tumorID)) + ' ' + 'tumors left')
            sample_RNA_list = []
            if len(tumorID) > 0:
                sample_RNA_list = [tumorID[0]]
                ##if any(str(tumorID[0] +'-0') in row1[1] and 'unc.edu' in row1[0] for row1 in file_map if len(row1) >0):
                for gene in genelist:
                    #print(gene)
                    for tumor_row in file_map:
                        if len(tumor_row) > 0 and len(tumor_row[0]) > 10 and tumorID[0] in tumor_row[1]: #selects barcode
                                    if int(tumor_row[1][13]) == 0:#selects the tumor sample for that barcode
                                        if 'unc.edu' in tumor_row[0] and 'genes.normalized_results' in tumor_row[0]:
                                            file_name = 'RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/' + tumor_row[0]
                                            #print(tumor_row[0])
                                            with open(file_name, 'r') as sample_RNA:
                                                sample_RNA = list(csv.reader(sample_RNA, delimiter='\t'))
                                            for RNA in sample_RNA:
                                                if gene == RNA[0].split('|')[0]:
                                                    #print(gene)
                                                    if RNA[1] != '0':
                                                        gene_quant = float(RNA[1])
                                                        if len(sample_RNA_list) == genelist.index(gene) + 1:
                                                            if 'duplicate' not in sample_RNA_list:
                                                                sample_RNA_list.append(gene_quant)
                                                        elif len(sample_RNA_list) > genelist.index(gene) + 1:            
                                                            sample_RNA_list = [tumorID[0], 'duplicate']
                                                        else:
                                                            sample_RNA_list.append('N/A')
                                                    else:
                                                        sample_RNA_list.append('N/A')
                                                                                                
                                                        
                                    else:
                                        #print(tumor_row)
                                        if int(tumor_row[1][13]) != 0 and int(tumor_row[1][13]) != 1:
                                            print('problem barcode')
                                            print(tumor_row[1])
    print('created RNA file')

