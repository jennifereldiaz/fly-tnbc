#june 2014
#determine most common copy number variants in a set of breast cancer patients


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

#get filtered gene list
with open('BRCA_CNVs_foldchange_TN_filtered.csv', 'r') as filtered:
    filtered = csv.reader(filtered)
    filtered_genelist = next(filtered)
    filtered_genelist = list(filtered_genelist)[1:]
print('Initial Gene List:')
print(len(filtered_genelist),'genes') #9694


#get amps and del
dict_ampdel = {}
with open('BRCA_TN_CNVs_foldchange_parsed_counts.csv', 'r') as parsed_counts_file:
    parsed_counts = csv.reader(parsed_counts_file)
    parsed_counts = list(parsed_counts)
for gene in filtered_genelist:
    i = parsed_counts[0].index(gene)
    if float(parsed_counts[2][i]) > 50:
        dict_ampdel[gene] = 'amp'
    elif float(parsed_counts[3][i]) > 50:
        dict_ampdel[gene] = 'del'
    else:
        print('PROBLEM GENE:')
        print(gene, 'up:', parsed_counts[2][i], 'down:', parsed_counts[3][i])


#remove genes not differentially expressed (i.e. where the value is 0 in most samples).
#borrowing this step from Akavia et al, 2010, from dana peer's lab
RNA = pd.read_csv('BRCA_RNA_candidates.csv',header=0,index_col=0)
#print('RNA BEFORE FILTERING:')
#print(RNA.head())
RNA = RNA.T
RNA['StDev'] = RNA.std(axis=1)
RNA = RNA[RNA['StDev']>0.25]
RNA.drop('StDev',axis=1)
RNA = RNA.T
#print('RNA AFTER FILTERING:')
#print(RNA.head())
RNA.to_csv('BRCA_RNA_candidates_filtered.csv') 

CNV = pd.read_csv('BRCA_CNVs_foldchange_all_filtered.csv',header=0,index_col=0)
CNV = CNV[list(RNA.columns.values)]
CNV.to_csv('BRCA_CNVs_foldchange_all_filtered2.csv')

TNcounts = pd.read_csv('BRCA_TN_CNVs_foldchange_parsed_counts.csv',header=0,index_col=0)
TNcounts = TNcounts[list(RNA.columns.values)]
TNcounts.to_csv('BRCA_TN_CNVs_foldchange_parsed_counts2.csv')

TNsum = pd.read_csv('BRCA_CNVs_foldchange_TN_filtered.csv',header=0,index_col=0)
TNsum = TNsum[list(RNA.columns.values)]
TNsum.to_csv('BRCA_CNVs_foldchange_TN_filtered2.csv')
filtered_genelist = list(TNsum.columns.values)
print('Filtered Gene List:')
print(len(filtered_genelist),'genes') #6709
    
#convert to z-scores
with open('BRCA_RNA_candidates_filtered.csv', 'r') as RNA:
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


#import RNA file and CNV file as list and concatenate:
with open('BRCA_CNVs_foldchange_all_filtered2.csv', 'r') as CNVs:
    CNVs = csv.reader(CNVs)
    CNVs = list(CNVs)
with open('BRCA_RNA_z_scores.csv', 'r') as RNA:
    RNA = csv.reader(RNA)
    RNA = list(RNA)
with open('BRCA_CNVs_foldchange_and_RNA.csv', 'w+') as comb:
    comb = csv.writer(comb)
    firstrow = ['Complete TCGA ID']
    for CNV_name in CNVs[0][1:len(CNVs[0])]:
        if CNV_name != '':
            CNV_header = CNV_name + '-CNV'
            firstrow.append(CNV_header)
 #   firstrow.append('tumor')
 #   firstrow.append('normal')
    for RNA_name in RNA[0][1:len(RNA[0])]:
        if RNA_name != '':
            RNA_header = RNA_name + '-RNA'
            firstrow.append(RNA_header)
    comb.writerow(firstrow)
    for sample in CNVs[1:]:
        for ID in RNA[1:]:
            if sample[0] == ID[0]:
                sample_list = [sample[0]]
                for i in sample[1:]:
                    sample_list.append(i)
                for i in ID[1:]:
                    sample_list.append(i)
                comb.writerow(sample_list)

##thenconcatenate the CNV and RNA files and do pairwise t-tests.
##set cutoff automatically by amp or del
##then output a list of genes after filtering them by p-value
datadf = pd.read_csv('BRCA_CNVs_foldchange_and_RNA.csv', header=0)
print('COMBINED FILE:')
#print(datadf.head(n=10))
print(datadf.shape)

final_genelist = []
with open('BRCA_CNV_foldchange_siggenes.csv','w+') as sig:
    sig = csv.writer(sig)
    sig.writerow(['Gene','CNV type','CN cutoff','Percent Altered','p-value for RNA t-test','Result','CN-RNA relationship'])
    equal = 0
    unequal = 0
    non_equal = 0
    non_unequal = 0
    upnormal = False
    downnormal = False
    for gene in filtered_genelist:
      if float((len(filtered_genelist) - filtered_genelist.index(gene))/50).is_integer(): 
          print(str(len(filtered_genelist) - filtered_genelist.index(gene)) + ' ' + 'genes left')
      #print(gene)
      CNV_header = gene + '-CNV' 
      RNA_header = gene + '-RNA'
      testdf = datadf[[CNV_header, RNA_header]]
      #print(testdf.head())
      #print(testdf.shape)
      testdf.dropna(inplace=True)
      testdf.columns = ['CNV', 'RNA']
      #print(testdf.head())
      #print(testdf.shape)
      nodup = testdf.RNA #checking to see that there is more than one value in RNA
      nodup.drop_duplicates(inplace=True)
      if nodup.shape[0] > 1:
          if dict_ampdel[gene] == 'amp': #test amplifications. here I will ONLY use the 1.2 cutoff.
              testdf = testdf[testdf.CNV > 2.**-0.3] #remove deletions
              upmask = (testdf.CNV > 2.**0.3) 
              upRNA = testdf[upmask].RNA
              upmean = upRNA.mean()
              upmedian = upRNA.median()
              cutoff = 2.**0.3
              downRNA = testdf[~upmask].RNA
              downmean = downRNA.mean()
              downmedian = downRNA.median()
              perc = TNcounts[gene].ix['percent up of total']*100
          elif dict_ampdel[gene] == 'del': #test deletions
              testdf = testdf[testdf.CNV < 2.**0.3] #remove amplifications
              upmask = (testdf.CNV > 2.**-0.3)
              upRNA = testdf[upmask].RNA
              upmean = upRNA.mean()
              upmedian = upRNA.median()
              downRNA = testdf[~upmask].RNA
              downmean = downRNA.mean()
              downmedian = downRNA.median()
              cutoff = 2.**-0.3
              perc = TNcounts[gene].ix['percent down of total']*100
          if scipy.stats.normaltest(upRNA)[1] > 0.05 or len(upRNA) >=30:
              upnormal = True
          if scipy.stats.normaltest(downRNA)[1] > 0.05 or len(downRNA) >= 30: #using the central limit theorem to say if the sample is large enough in approximates normal.
              downnormal = True
          if upnormal and downnormal: #will use one-sided t tests here, because I only want those cases where upRNA > downRNA, not the other way around.
              if scipy.stats.bartlett(upRNA,downRNA)[1] > 0.05: 
                  p_value = scipy.stats.ttest_ind(upRNA, downRNA)[1]/2
                  equal += 1
              else:
                  p_value = scipy.stats.ttest_ind(upRNA, downRNA,equal_var=False)[1]/2 #Welsch t-test
                  unequal += 1
              if upmean > downmean:
                  relationship = '+'
              else: relationship = '-'
          else:
              if scipy.stats.levene(upRNA,downRNA)[1] > 0.05:
                  p_value = scipy.stats.mannwhitneyu(upRNA, downRNA)[1]/2 #non-parametric test
                  non_equal += 1
              else:
                  p_value = scipy.stats.mannwhitneyu(upRNA, downRNA)[1]/2 ##using the Mann-Whitney U test here. Not robust for unequal variances. could try transform
                  non_unequal += 1
              if upmedian > downmedian:  #can't consider means for a nonnormal samples
                  relationship = '+'
              else: relationship = '-'
              
          #print(p_value)
          if p_value < 0.05/len(filtered_genelist):
              final_genelist.append(gene)
              stat = 'Significant'
          else: stat = 'Non significant'
          sig.writerow([gene,dict_ampdel[gene],cutoff,perc,p_value,stat,relationship])

print('Normal, equal variance:',equal) #
print('Normal, unequal variance:', unequal)#
print('Nonnormal, equal variance:',non_equal) #note that all samples were normal in this run 
print('Nonnormal, unequal variance:',non_unequal) 
if non_unequal > 0:
    print('WARNING: Nonnormal, unequal variance samples were identified. Modify the script to show the distribution of these samples and attempt to transform the samples to achieve equal variances or normality.')
#that is, either they actually were normally distributed, or the sample size was > 30. So it is actually just extra conservative to use only the Welsch test in
#this case, as I suspect the CONEXIC people will have done. I have used either a regular t-test or a welsch t-test where appropriate.
print(len(final_genelist), 'significant genes')
print('created significant genes file')

#read in the sig genes file and output counts
siggenes = pd.read_csv('BRCA_CNV_foldchange_siggenes.csv',header=0)
sigonly = siggenes[siggenes['Result']=='Significant']
#print('Number of significant genes:', sigonly.shape[0])
posonly = sigonly[sigonly['CN-RNA relationship']=='+']
print('Significant genes with + relationship:', posonly.shape[0]) #

siggenes.set_index('Gene',inplace=True)
siggenes = siggenes[['CNV type','CN cutoff','Percent Altered','p-value for RNA t-test','Result','CN-RNA relationship']]
cands = pd.read_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_step5.csv',header=0,index_col='Symbol')
sub1 = cands[cands['Has RNA data']=='yes']
rest1 = cands[~(cands['Has RNA data']=='yes')]
sub2 = sub1[sub1.index.isin(filtered_genelist)]
rest2 = sub1[~(sub1.index.isin(filtered_genelist))]
sub2['RNA differentially expressed'] = 'yes'
rest2['RNA differentially expressed'] = 'no'
rest1['RNA differentially expressed'] = ''
cands = pd.concat([sub2,rest2,rest1])
cands = cands.merge(siggenes,how='outer',left_index=True, right_index=True)
cands.to_csv('BRCA_allGISTIC2.0andISARgenes_foldchange_compressed_step6.csv',index=True)
