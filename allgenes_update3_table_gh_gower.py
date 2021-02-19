#dec 29 2015
#updated feb 17 2016
##add oncodrive,CIVIC
#rerun for fusions
#start looking at transcription factors
#run human RNA seq data through Enrichr for GO and KEGG
#add whether gene has fly homolog for all genes
#updated feb 26 2017
#adjustments to include any validated stock for each gene in screen
#april 23 2016
#adding in reclassifications of the TNBCtypes

#updated june 4 2020 for table:
# include all regions, not just those flagged as TNBC specific in my analysis

import pandas as pd
import numpy as np
import csv
from statsmodels.sandbox.stats.multicomp import multipletests as mult
import seaborn as sns
from matplotlib import pyplot as plt
import scipy as sc
import gower



ampcutoff = 2**0.3
delcutoff = 2**-0.3

#reverse the mutation data back to categorical
dic = {'Missense_Mutation':2.0, 'Splice_Site':0.3, 'Frame_Shift_Del':0.1,
       'In_Frame_Del':1.8, 'Nonsense_Mutation':0.0, 'In_Frame_Ins':1.9,
       'Frame_Shift_Ins':0.2, 'RNA':1.6, 'Nonstop_Mutation':1.7}
dicr= {v: k for k, v in dic.items()}
dicr[1]='None'

#import genomic data
print('importing genomic data...')
cnv = pd.read_csv('../CNVruns/new_TN_CNV_run/BRCA_TN_CNV_region_values_matrix_update.csv',
                  index_col='TCGA_ID') #CNVruns/

mut = pd.read_csv('../Somatic_mutations/BRCA_somatic_muts_matrix.csv', ##mutsig and the databases below are included in this data
                  index_col = 'Unnamed: 0')
#mut = mut.replace(dicr)
mut.columns = [c+'-som' for c in mut.columns]

germ = pd.read_csv('../germline_mutations/BRCA_germline_muts_matrix.csv',index_col='Unnamed: 0')
#germ = germ.replace(dicr)
germ.columns = [c+'-germ' for c in germ.columns]

fusion = pd.read_excel('BRCA_TN_selected_fusion_proteins.xlsx',sheetname = 'Sheet1')


#
print('compiling genomic data...')
#put genomic data together
both = cnv.merge(mut,how='inner',left_index=True,right_index=True)
three = both.merge(germ,how='left',left_index=True,right_index=True).fillna(1.0) #all patients I'm interested in have germline data, but most don't have mutations in this table. fill the rest with normal 1.0
    
#TNBCtypes
mayer = pd.read_excel('source_papers/CCR-13-0583tab1.xlsx',sheetname='Sheet1',skiprows=3)
mayer.set_index('TCGA_ID',drop=False,inplace=True)
rm= dict()
for col in mayer.columns:
    rm[col] = col+'_Mayer'
mayer.rename(columns=rm,inplace=True)
lehmann = pd.read_excel('source_papers/journal.pone.0157368.s008.XLSX',sheetname='TCGA',skiprows=1)
lehmann = lehmann[lehmann.TCGA_SAMPLE.str.contains('-01')]
lehmann.index = lehmann.BARCODE.str[-4:]
rl = dict()
for col in lehmann.columns:
    rl[col] = col+'_Lehmann'
lehmann.rename(columns=rl,inplace=True)
bareche = pd.read_excel('source_papers/bareche_eTable2.xlsx',sheetname='Sheet1',skiprows=2)
bareche = bareche[bareche.Patient.str.contains('TCGA')&bareche.Patient.str.contains('-01')]
bareche.index= bareche.Patient.str[-7:-3]
rb = dict()
for col in bareche.columns:
    rb[col] = col+'_Bareche'
bareche.rename(columns=rb,inplace=True)

types = pd.concat([mayer[['TCGA_ID_Mayer','TNBC_Subtype_Mayer']],
                   lehmann[['BARCODE_Lehmann','TNBC_Lehmann','PAM50_Lehmann',
                            'PAM50lite_Lehmann','TNBCtype_Lehmann','TNBCtype_4_Lehmann']],
                   bareche[['Patient_Bareche','PAM50_Bareche','TNBCtype_Bareche',
                            'TNBCtype Reclassified_Bareche']]],
                  axis=1,join='outer')
#types.TCGA_ID=types.TCGA_ID.astype(str)
#types.Patient=types.Patient.astype(str)

types.dropna(how='all',inplace=True)
types.to_csv('BRCA_TNBCsubtypes_updated.txt',sep='\t')

#put the tnbctypes in with genomic data
three['cut'] = three.index.str.slice(start=8)
three['TCGA_ID'] = three.index
three.set_index('cut',inplace=True)

four = three.merge(types[['TNBC_Subtype_Mayer', 'TNBCtype_Lehmann',
                          'TNBCtype_4_Lehmann','TNBCtype_Bareche', 'TNBCtype Reclassified_Bareche']],how='inner',left_index=True,
                   right_index=True)



#survival
surv = pd.read_excel('../survival/1-s2.0-S0092867418302290-mmc1.xlsx',sheet_name='TCGA-CDR',index_col=1)
surv=surv[surv.type=='BRCA'][['OS.time','DSS.time','DFI.time','PFI.time']]
surv.index =surv.index.str[-4:]
five = four.merge(surv[['OS.time','DSS.time','DFI.time','PFI.time']],how='left',left_index=True,right_index=True)
#for col in ['OS','DSS','DFI','PFI']:
#    five[col+'.gray']=five[col+'.time']/five[col+'.time'].max()
#    five[col+'.rgb']='['+five[col+'.gray'].astype(str)+','+five[col+'.gray'].astype(str)+ \
#                      ','+five[col+'.gray'].astype(str)+']'
#five.drop(['OS.time','DSS.time','DFI.time','PFI.time','OS.gray','DSS.gray','DFI.gray','PFI.gray'],axis=1,inplace=True)

table = five.set_index('TCGA_ID')
table = table.loc[:, (table != 1.0).any(axis=0)]

nosub = table.drop(['TNBC_Subtype_Mayer', 'TNBCtype_Lehmann',
            'TNBCtype_4_Lehmann','TNBCtype_Bareche', 'TNBCtype Reclassified_Bareche',
            'OS.time','DSS.time','DFI.time','PFI.time'],
           axis=1)

lut = dict(zip(table.index,plt.cm.gray_r((table['OS.time']-table['OS.time'].min())/(
        table['OS.time'].max()-table['OS.time'].min()))))

nosubg = nosub.copy()
mcols = [col for col in nosub.columns if any([col in ls for ls in [mut.columns,germ.columns]])]
nosubg[mcols] = nosubg[mcols].replace(dicr)

g = gower.gower_matrix(nosubg)
v = sc.spatial.distance.squareform(g)
l = sc.cluster.hierarchy.linkage(v)

sns.set(font_scale=0.5)

cm = sns.clustermap(nosub,row_linkage=l,
                    col_cluster=False,cmap='RdBu_r',vmin=0,vmax=2,
                    yticklabels=True,xticklabels=False,#'bwr' #'seismic' ,cbar_pos='left'
               row_colors=nosub.index.map(lut),figsize=(11,7))
              # cbar_pos=(0.035,0.05,0.025,0.6))

fig = cm.fig
#fig.savefig('allgenes_figure.pdf')
               #cbar_kws={"ticks":None})


#hm = cm.ax_heatmap.get_position()
#plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=6)
#cm.ax_heatmap.set_position([hm.x0, hm.y0, hm.width*0.25, hm.height])
#col = cm.ax_col_dendrogram.get_position()
#cm.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.25, col.height*0.5])

#table.to_csv('BRCA_all_aberrations_formatlab.txt',sep='\t')
#table.drop(['TNBC_Subtype_Mayer', 'TNBCtype_Lehmann',
#            'TNBCtype_4_Lehmann','TNBCtype_Bareche', 'TNBCtype Reclassified_Bareche',
#            'OS.rgb','DSS.rgb','DFI.rgb','PFI.rgb'],
#           axis=1).to_csv('BRCA_all_aberrations_nosubtype_formatlab.txt',sep='\t')
#table[['OS.rgb']].to_csv('BRCA_survival_labels_OS_formatlab.txt',sep='\t',header=False)
#table[['DSS.rgb']].to_csv('BRCA_survival_labels_DSS_formatlab.txt',sep='\t',header=False)
#table[['PFI.rgb']].to_csv('BRCA_survival_labels_PFI_formatlab.txt',sep='\t',header=False)
#table[['DFI.rgb']].to_csv('BRCA_survival_labels_DFI_formatlab.txt',sep='\t',header=False)


colors = {'BL1':'[0 0.5 0.5]','UNC':'[0.5 1 0]','M':'y','IM':'[1 0.5 0]','MSL':'m','LAR':'r','BL2':'k',np.nan:'w','nan':'w','HER2':'w','ER':'w','UNS':'[0.5 1 0]'}
#BL1 teal, UNC green, M yellow, IM orange, MSL magenta, LAR red, BL2 black
types = table[['TNBC_Subtype_Mayer','TNBCtype_Lehmann','TNBCtype_4_Lehmann','TNBCtype_Bareche',
               'TNBCtype Reclassified_Bareche']]
types.replace(colors,inplace=True)
#types.TNBC_Subtype_Mayer.to_csv('BRCA_subtype_colors_formatlab_Mayer.txt',sep='\t',header=False)
#types.TNBCtype_Lehmann.to_csv('BRCA_subtype_colors_formatlab_Lehmann.txt',sep='\t',header=False)
#types.TNBCtype_4_Lehmann.to_csv('BRCA_subtype_colors_formatlab_Lehmann4.txt',sep='\t',header=False)
#types.TNBCtype_Bareche.to_csv('BRCA_subtype_colors_formatlab_Bareche.txt',sep='\t',header=False)
#types['TNBCtype Reclassified_Bareche'].to_csv('BRCA_subtype_colors_formatlab_Bareche5.txt',sep='\t',header=False)
#types.to_csv('BRCA_subtype_colors_formatlab.txt',sep='\t',header=False)

