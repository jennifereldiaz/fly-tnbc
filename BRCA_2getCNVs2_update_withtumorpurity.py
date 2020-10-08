#june 2014
#determine most common copy number variants in a set of breast cancer patient
##updated jul 2015 to get lists for both TNBC and ER+


import csv
import math
import numpy as np
#import scipy
#from scipy import stats
#import matplotlib.pyplot as plt
import math
import itertools
from itertools import zip_longest
import pandas as pd


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

#function to pull out a copy number    
def get_copy_number(sample_list, gene_list):
#sample_CNV_file is the file for that sample ID. gene is a list from the candidate file
        copy_number = None
        for row in sample_list:
            if str.isdigit(row[1]):
                if str.isdigit(gene_list[2]):
                    if int(row[1]) == int(gene_list[2]):
                        start_CNV = int(row[2])
                        stop_CNV = int(row[3])
                        start_gene = int(gene_list[3])
                        stop_gene = int(gene_list[4])
                        seg_tuple_CNV = (start_CNV, stop_CNV)
                        seg_tuple_gene = (start_gene, stop_gene) #this figures out if the known CNV is entirely contained in the CNV in the sample 
                        overlap = min(seg_tuple_CNV[1], seg_tuple_gene[1]) - max(seg_tuple_CNV[0], seg_tuple_gene[0])
                        seg_gene = seg_tuple_gene[1] - seg_tuple_gene[0]
                        if overlap >= seg_gene:
                            copy_number = math.pow(2,float(row[5]))
                            return copy_number
            elif row[1] =='X':
                if gene_list[2]=='X':
                    if row[1] == gene_list[2]:
                        start_CNV = int(row[2])
                        stop_CNV = int(row[3])
                        start_gene = int(gene_list[3])
                        stop_gene = int(gene_list[4])
                        seg_tuple_CNV = (start_CNV, stop_CNV)
                        seg_tuple_gene = (start_gene, stop_gene) #this figures out if the known CNV is entirely contained in the CNV in the sample 
                        overlap = min(seg_tuple_CNV[1], seg_tuple_gene[1]) - max(seg_tuple_CNV[0], seg_tuple_gene[0])
                        seg_gene = seg_tuple_gene[1] - seg_tuple_gene[0]
                        if overlap >= seg_gene:
                            copy_number = math.pow(2,float(row[5]))
                            return copy_number
        
                                        

#list of candidate genes
cand_genes = []

tncand = pd.read_csv('new_TN_CNV_run/CNV_candidates.csv',header=0)
ercand = pd.read_csv('ER+_CNV_run/CNV_candidates.csv',header=0)

cand = pd.concat([tncand,ercand])
cand.drop_duplicates(inplace=True)
cand.reset_index(drop=True,inplace=True)

##cand = cand.ix[:4]

cand_genes= cand.Symbol.unique().tolist()

print('Initial Gene List:')
print(len(cand_genes))

##import consensus tumor purity estimate
pur = pd.read_excel('purity_aran2015-s2.xlsx',sheetname='Supp Data 1',skiprows = 3,header=0,index_col='Sample ID')

#generate summary files of candidate amplifications and deletions for all tumors by original barcode, adjusted for tumor purity from aran 2015

with open('BRCA_CNVs_genes_foldchange.csv', 'w+') as CNV_file:
    CNVs = csv.writer(CNV_file)
    with open('BRCA_pathology_2014.csv', 'r') as path_file:
            path = csv.reader(path_file)
            path = list(path)
    secondrow = ['TCGA_ID',]
    firstrow = ['HUGO-->',]
    firstrow += cand.Symbol.tolist()
    secondrow += cand.Synonym.tolist()
    secondrow.append('tumor')
    secondrow.append('normal')
    firstrow.append('tumor')
    firstrow.append('normal')
    CNVs.writerow(firstrow)
    CNVs.writerow(secondrow)
    with open('CNV_data/FILE_SAMPLE_MAP.txt', 'r') as file_map:
        file_map = csv.reader(file_map, delimiter='\t')
        file_map = list(file_map)
        for item in file_map:
                if len(item) <1:
                        file_map.remove(item)
        for tumorID in path:
            print(tumorID[0])
            if len(tumorID[0]) > 0:
                primary = tumorID[0] + '-01' #note that all tumors in my final file are primary, and none were excluded on the basis of having a
                ##met sample.
                met = tumorID[0] + '-06'
                blood = tumorID[0] + '-10'
                solid = tumorID[0] + '-11'
                sample_CNVs_list = []
                germline = ''
                tumor = ''
                if any(primary in row[1] for row in file_map):
                    purityid = primary+'A'
                    if purityid in pur.index.tolist():
                        purity = float(pur.ix[purityid].CPE)
                    else:
                        purity = np.nan
                    #print(primary)
                    for row in file_map:
                        if len(row) > 0:
                            if primary in row[1]:#adds all the CNVs for the tumor sample
                                    if 'nocnv_hg19' in row[0]:
                                            tumor = 'primary'
                                            if sample_CNVs_list == []:
                                                    sample_CNVs_list.append(tumorID[0])
                                            file_name = 'CNV_data/Level_3/' + row[0]
                                            with open(file_name, 'r') as sample_CNVs:
                                                sample_CNVs = csv.reader(sample_CNVs, delimiter='\t')
                                                sample_CNVs = list(sample_CNVs)
                                                for i,row in cand.iterrows():
                                                    if len(sample_CNVs_list) <= i:
                                                        sample_CNVs_list.append('')
                                                    line = row.tolist()
                                                    copy_number_primary = get_copy_number(sample_CNVs, line)
                                                    if copy_number_primary != None:
                                                            if any(blood in row[1] for row in file_map):
                                                                for row2 in file_map:
                                                                    if len(row2) > 0:
                                                                        if blood in row2[1]:#subtracts the CNVs from the germline
                                                                                if 'nocnv_hg19' in row2[0]:
                                                                                    file_name2 = 'CNV_data/Level_3/' + row2[0]
                                                                                    with open(file_name2, 'r') as sample_CNVs_germ:
                                                                                        sample_CNVs_germ = csv.reader(sample_CNVs_germ, delimiter='\t')
                                                                                        sample_CNVs_germ = list(sample_CNVs_germ)
                                                                                        copy_number_germline = get_copy_number(sample_CNVs_germ, line)
                                                                                        if copy_number_germline != None:
                                                                                            if np.isnan(purity): #formula to correct for tumor purity: observedCN = trueCN*purity + germlineCN*(1-purity)
                                                                                                new_copy_number = copy_number_primary/copy_number_germline #fold change over germline
                                                                                            else:
                                                                                                new_copy_number = (copy_number_primary+(purity-1)*copy_number_germline)/(purity*copy_number_germline)
                                                                                            sample_CNVs_list.append(new_copy_number)
                                                                                        else:
                                                                                            if np.isnan(purity):
                                                                                                sample_CNVs_list.append(copy_number_primary)
                                                                                            else:
                                                                                                new_copy_number = (copy_number_primary+(purity-1))/purity
                                                                                                sample_CNVs_list.append(new_copy_number)
                                                                                        germline = 'blood'
                                                            elif any(solid in row[1] for row in file_map):
                                                                for row2 in file_map:
                                                                        if len(row2) > 0:
                                                                            if solid in row2[1]:
                                                                                if 'nocnv_hg19' in row2[0]:
                                                                                        file_name2 = 'CNV_data/Level_3/' + row2[0]
                                                                                        with open(file_name2, 'r') as sample_CNVs_germ_solid:
                                                                                            sample_CNVs_germ_solid = csv.reader(sample_CNVs_germ_solid, delimiter='\t')
                                                                                            sample_CNVs_germ_solid = list(sample_CNVs_germ_solid)
                                                                                            copy_number_germline = get_copy_number(sample_CNVs_germ_solid, line)
                                                                                            if copy_number_germline != None:
                                                                                                if np.isnan(purity): #formula to correct for tumor purity
                                                                                                    new_copy_number = copy_number_primary/copy_number_germline
                                                                                                else:
                                                                                                    new_copy_number = (copy_number_primary+(purity-1)*copy_number_germline)/(purity*copy_number_germline)
                                                                                                sample_CNVs_list.append(new_copy_number)
                                                                                            else:
                                                                                                if np.isnan(purity):
                                                                                                    sample_CNVs_list.append(copy_number_primary)
                                                                                                else:
                                                                                                    new_copy_number = (copy_number_primary+(purity-1))/purity
                                                                                                    sample_CNVs_list.append(new_copy_number)
                                                                                            germline = 'solid'
                                                            else: sample_CNVs_list.append(copy_number_primary)
                else:
                        if any(met in row[1] for row in file_map):
                            purityid = met+'A'
                            if purityid in pur.index.tolist():
                                purity = float(pur.ix[purityid].CPE)
                            else:
                                purity = np.nan
                            #print(met)
                            for row in file_map:
                                if len(row) > 0:
                                    if met in row[1]:#adds all the CNVs for the tumor sample
                                            if 'nocnv_hg19' in row[0]:
                                                    tumor = 'met'
                                                    if sample_CNVs_list == []:
                                                            sample_CNVs_list.append(tumorID[0])
                                                    file_name = 'CNV_data/Level_3/' + row[0]
                                                    with open(file_name, 'r') as sample_CNVs:
                                                        sample_CNVs = csv.reader(sample_CNVs, delimiter='\t')
                                                        sample_CNVs = list(sample_CNVs)
                                                    for i,row in cand.iterrows():
                                                        if len(sample_CNVs_list) <= i:
                                                            sample_CNVs_list.append('')
                                                        line = row.tolist()
                                                        copy_number_met = get_copy_number(sample_CNVs, line)
                                                        if copy_number_met != None:
                                                                if any(blood in row[1] for row in file_map):
                                                                    for row2 in file_map:
                                                                        if len(row2) > 0:
                                                                            if blood in row2[1]:#subtracts the CNVs from the germline
                                                                                    if 'nocnv_hg19' in row2[0]:
                                                                                        file_name2 = 'CNV_data/Level_3/' + row2[0]
                                                                                        with open(file_name2, 'r') as sample_CNVs_germ:
                                                                                            sample_CNVs_germ = csv.reader(sample_CNVs_germ, delimiter='\t')
                                                                                            sample_CNVs_germ = list(sample_CNVs_germ)
                                                                                            copy_number_germline = get_copy_number(sample_CNVs_germ, line)
                                                                                            if copy_number_germline != None:
                                                                                                if np.isnan(purity): #formula to correct for tumor purity
                                                                                                    new_copy_number = copy_number_met/copy_number_germline
                                                                                                else:
                                                                                                    new_copy_number = (copy_number_met+(purity-1)*copy_number_germline)/(purity*copy_number_germline)
                                                                                                sample_CNVs_list.append(new_copy_number)
                                                                                            else:
                                                                                                if purity != np.nan:
                                                                                                    sample_CNVs_list.append(copy_number_met)
                                                                                                else:
                                                                                                    new_copy_number = (copy_number_met+(purity-1))/purity
                                                                                                    sample_CNVs_list.append(new_copy_number)
                                                                                            germline = 'blood'
                                                                elif any(solid in row[1] for row in file_map):
                                                                    for row2 in file_map:
                                                                            if len(row2) > 0:
                                                                                if solid in row2[1]:
                                                                                    if 'nocnv_hg19' in row2[0]:
                                                                                            file_name2 = 'CNV_data/Level_3/' + row2[0]
                                                                                            with open(file_name2, 'r') as sample_CNVs_germ_solid:
                                                                                                sample_CNVs_germ_solid = csv.reader(sample_CNVs_germ_solid, delimiter='\t')
                                                                                                sample_CNVs_germ_solid = list(sample_CNVs_germ_solid)
                                                                                                copy_number_germline = get_copy_number(sample_CNVs_germ, line)
                                                                                                if copy_number_germline != None:
                                                                                                    if np.isnan(purity): #formula to correct for tumor purity
                                                                                                        new_copy_number = copy_number_met/copy_number_germline
                                                                                                    else:
                                                                                                        new_copy_number = (copy_number_met+(purity-1)*copy_number_germline)/(purity*copy_number_germline)
                                                                                                    sample_CNVs_list.append(new_copy_number)
                                                                                                else:
                                                                                                    if np.isnan(purity):
                                                                                                        sample_CNVs_list.append(copy_number_met)
                                                                                                    else:
                                                                                                        new_copy_number = (copy_number_met+(purity-1))/purity
                                                                                                        sample_CNVs_list.append(new_copy_number)
                                                                                                germline = 'solid'
                                                                else: sample_CNVs_list.append(copy_number_met)
                #print(sample_CNVs_list)                                       
                if len(sample_CNVs_list) < len(cand)+1:
                    sample_CNVs_list.append('')
                if len(sample_CNVs_list) > len(cand)+1:
                        del sample_CNVs_list[1:]
                        sample_CNVs_list.append('duplicates')
                if len(sample_CNVs_list) == len(cand)+1:
                    sample_CNVs_list.append(tumor)
                    sample_CNVs_list.append(germline)
                if len(cand)+1 < len(sample_CNVs_list) < len(cand) + 3:
                        sample_CNVs_list.append('')
                if len(sample_CNVs_list) > 1:
                        CNVs.writerow(sample_CNVs_list)




