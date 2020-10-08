#nov 9 2015

import pandas as pd
import numpy as np

###get list of drivers
###from walsh 2012 paper (mc king)
##drivers = pd.read_csv('walsh2010PNAStable2.txt',sep='\t')

#open maf
maf = pd.read_excel('TCGA_germline-variants_analysis-12-14-11-BERG.Couch.KH_JELD.xlsx',
                    sheetname='Mutation carriers')
maf.dropna(subset=['hugo_symbol'],inplace=True)
maf.rename(columns={'bcr_patient_barcode':'TCGA_ID'},inplace=True)

###remove silent mutations
##allmaf = pd.read_excel('TCGA_germline-variants_analysis-12-14-11-BERG.Couch.KH_JELD.xlsx',
##                    sheetname='Annotated_variants')
##allmaf = allmaf[allmaf.variant_classification != 'Silent']
##
###genes in the walsh paper, not in tcga
##missed = allmaf[allmaf.hugo_symbol.isin(drivers[~drivers.Gene.isin(maf.hugo_symbol)].Gene)]
##missed = missed[missed.revised==3]
###the walsh mutations for these genes don't appear in the TCGA data. So I don't see much of a reason to use the walsh data.


#make a table
dic = {'Missense_Mutation':2.0, 'Splice_Site':0.3, 'Frame_Shift_Del':0.1,
       'In_Frame_Del':1.8, 'Nonsense_Mutation':0, 'In_Frame_Ins':1.9,
       'Frame_Shift_Ins':0.2, 'RNA':1.6, 'Nonstop_Mutation':1.7}
maf['Variant_Value'] = maf['variant_classification'].map(dic)
table = maf.pivot_table(index='TCGA_ID',columns='hugo_symbol',
                         values='Variant_Value',
                         aggfunc=lambda x: np.sum(x)+10*(len(x)-1))
table.fillna(1.0,inplace=True)
missed = pd.Series(maf.TCGA_ID.unique().tolist())
missed = missed[~missed.isin(table.index)]
empty = pd.DataFrame(data=0,index=missed.tolist(),columns=table.columns)
table = pd.concat([table,empty])
table.to_csv('BRCA_germline_muts_matrix.csv')
