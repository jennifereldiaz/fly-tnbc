#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 00:35:57 2020

@author: jenniferlongdiaz
"""
import pandas as pd


muts = pd.read_csv('MutSigCV2015.sig_genes.txt',sep='\t')
muts.set_index('gene',inplace=True)

hits1 = pd.concat([pd.read_excel('../../CNV_screen/final_group1_linelevel.xlsx',
                                 sheetname='Group 1I'),
                   pd.read_excel('../../CNV_screen/final_group1_linelevel.xlsx',
                                 sheetname='Group 1G')])
hits2 = pd.concat([pd.read_excel('../../CNV_screen/final_group2_linelevel.xlsx',
                                 sheetname='Group 2I'),
                   pd.read_excel('../../CNV_screen/final_group2_linelevel.xlsx',
                                 sheetname='Group 2G')])

hits = pd.concat([hits1,hits2])

hits = hits[hits['Validation result'].astype(str).str.contains('\+')] #just the validations (hits.HITS>0)|
hits = hits[['Unnamed: 0','Chromosome','from','to','CNV type', #'HITS',
             'ISARpeak','Stock tested',
             'Type of allele','Validation result','tier']].drop_duplicates()
hits['Unnamed: 0'] = hits['Unnamed: 0'].astype(str)

hits.set_index('Unnamed: 0',inplace=True)

hits = hits.merge(muts,how='left',left_index=True,right_index=True)

hits.sort_values(['Chromosome','from'],inplace=True)

hits.to_csv('driver_mutsig.csv')

hits[hits.q<1]
