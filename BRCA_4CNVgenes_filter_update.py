#june 2014
#determine genes in copy number variants for TNBC
#genes lost at this stage are not relevant to triple negative


import csv
import math
import numpy as np
import scipy
from scipy import stats
from scipy import misc
import matplotlib.pyplot as plt
import math
import itertools
from itertools import zip_longest
import pandas as pd

#in order to create a candidate CNV file for a large number of genes,
#I need to automatically pull out the genomic coordinates for build hg19 for each gene

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

#function for finding binomial probability density with equal probabilities
def BinProb(n,k):
    p = 0.5
    return misc.comb(n,k)*(p**k)*((1-p)**(n-k))

#set amp and del cutoffs
ampcutoff = 2.**0.3 
delcutoff = 2.**-0.3 

#CNV data for candidate genes and pathological information by TCGA barcode
with open('BRCA_CNVs_genes_foldchange_processed.csv', 'r') as CNV_cand:
        cand = csv.reader(CNV_cand)
        cand_genes = next(cand)
        cand_genes = list(cand_genes)[1:-2]
print('Initial Gene List:')
print(len(cand_genes), 'genes') #12029
path = pd.read_csv('../BRCA_pathology_2014.csv',header=0,names=['TCGA_ID','ER_Status','PR_Status','HER2_Status'])
path.set_index('TCGA_ID',drop=False,inplace=True)
#print(path.head())
CNVs = pd.read_csv('BRCA_CNVs_genes_foldchange_processed.csv')
CNVs.set_index('TCGA_ID',drop=True,inplace=True)
allinfo = pd.concat([path,CNVs],axis=1,join='inner')
TNmask = (allinfo.ER_Status == 'Negative') & (allinfo.PR_Status == 'Negative') & (allinfo.HER2_Status == 'Negative')
allTN = allinfo[TNmask]
#print(allTN.head())
allTN.to_csv('BRCA_TN_CNVs_foldchange.csv',index=False)

    
##parse desired cutoffs
with open('BRCA_TN_CNVs_foldchange_parsed.csv', 'w+') as parsed_file:
    parsed = csv.writer(parsed_file)
    with open('BRCA_TN_CNVs_foldchange.csv', 'r') as TN_file:
        TN = csv.reader(TN_file)
        TN = list(TN)
        parsed.writerow(TN[0])
        for row in TN[1:]:
            parsed_row = []
            for i in row:
                if isnumber(i):
                    if float(i) > ampcutoff:
                        parsed_row.append(i)
                    elif float(i) < delcutoff:
                        parsed_row.append(i)
                    else: parsed_row.append('')
                else: parsed_row.append(i)
            parsed.writerow(parsed_row)
print('parsed >1.2 or <0.8')
                       
#add counts. using a p<0.05 cutoff
genes = cand_genes
with open('BRCA_TN_CNVs_foldchange_parsed_counts.csv', 'w+') as parsed_counts_file:
    parsed_counts = csv.writer(parsed_counts_file)
    with open('BRCA_TN_CNVs_foldchange_parsed.csv', 'r') as parsed_file:
        parsed = csv.reader(parsed_file)
        parsed = list(parsed)
        dict_ctr = {}
        dict_upperc = {}
        dict_downperc = {}
        dict_upperctot = {}
        dict_downperctot = {}
        dict_ampdel = {} #dictionary to call amplification or deletion for each gene - need this to determine cutoffs for pairwise comparisons
        first_row = ['gene']
        second_row = ['number altered']
        third_row = ['percent up of altered']
        fourth_row = ['percent down of altered']
        fifth_row = ['percent up of total']
        sixth_row =['percent down of total']
        for gene in genes: 
            ctr = 0
            upctr = 0
            downctr = 0
            for row in parsed:
                if len(row) >= len(genes) + 4:
                    if isnumber(row[genes.index(gene)+4]):
                        ctr += 1
                        if float(row[genes.index(gene)+4]) > ampcutoff:
                                upctr += 1
                        if float(row[genes.index(gene)+4]) < delcutoff:
                                downctr += 1
            if upctr + downctr != ctr:
                print(gene)
                print(upctr, downctr, ctr)
            if ctr>0:
                dict_ctr[gene] = ctr
                dict_upperc[gene] = upctr/ctr*100
                dict_upperctot[gene] = upctr/len(parsed[1:])
                dict_downperc[gene] = downctr/ctr*100
                dict_downperctot[gene] = downctr/len(parsed[1:])
                if BinProb(ctr,upctr) < 0.05:
                        first_row.append(gene)
                        second_row.append(dict_ctr[gene])
                        third_row.append(dict_upperc[gene])
                        fourth_row.append(dict_downperc[gene])
                        fifth_row.append(dict_upperctot[gene])
                        sixth_row.append(dict_downperctot[gene])
                        if upctr > downctr:
                            dict_ampdel[gene] = 'amp'
                        elif downctr > upctr:
                            dict_ampdel[gene] = 'del'
##                else:
##                    print(gene,BinProb(ctr,upctr))
            else: print('no altered samples for', gene)
        parsed_counts.writerow(first_row)
        parsed_counts.writerow(second_row)
        parsed_counts.writerow(third_row)
        parsed_counts.writerow(fourth_row)
        parsed_counts.writerow(fifth_row)
        parsed_counts.writerow(sixth_row)
print('percents in BRCA_TN_CNVs_foldchange_parsed_counts')

##filter by binomial probability
with open('BRCA_TN_CNVs_foldchange_parsed_counts.csv', 'r') as parsed_counts_file:
    counts = csv.reader(parsed_counts_file)
    counts = list(counts)
filtered_genelist = []    #list of filtered genes
for gene in genes:
        if gene in counts[0][1:]:
                filtered_genelist.append(gene)
print('Filtered Gene List:')            
print(len(filtered_genelist)) #9694

TNfiltered = allTN[filtered_genelist]
TNfiltered.to_csv('BRCA_CNVs_foldchange_TN_filtered.csv')

allfiltered = allinfo[filtered_genelist]
allfiltered.to_csv('BRCA_CNVs_foldchange_all_filtered.csv')

cands = pd.read_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_step3.csv',header=0)
sub = cands[cands.Symbol.isin(filtered_genelist)]
rest = cands[~(cands.Symbol.isin(filtered_genelist))]
sub['Copy number skewed in TN'] = 'yes'
rest['Copy number skewed in TN'] = 'no'
cands2 = pd.concat([sub,rest])
cands2.to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_step4.csv',index=False)
print('created filtered CNV files')
