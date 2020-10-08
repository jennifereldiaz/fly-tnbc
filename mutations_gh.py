#nov 9 2015

import pandas as pd
import numpy as np

#get list of drivers
#mutsigCV output (from the online tool)
mutsig = pd.read_csv('MutSigCV2015.sig_genes.txt',sep='\t')
mutsig = mutsig[mutsig.q < 0.1]
cosmic = pd.read_excel('COSMIC_cancer_gene_census.xls',sheetname = 'List')
pan = pd.read_excel('Pan_Cancer.xlsx',sheetname='Sheet1')
pan = pan[pan.Type=='MUTATION']
vogel = pd.read_excel('Vogelstein-cancer-genes.xlsx',sheetname='Table S2A',
                      skiprows=1)
civ = pd.read_csv('../nightly-GeneSummaries_CiVICdb_160329.tsv',sep='\t',header=0)
drivers = pd.concat([mutsig.gene,cosmic.Symbol,pan['Altered Locus/Gene'],
                     vogel['Gene Symbol'],civ.name,pd.Series(['TERT'])]).reset_index(drop=True)
#this is WXS so TERT promoter mutations won't show up. Maybe check expression data
#,'TRIO'

#open maf
maf = pd.read_csv('MUTATIONS/Somatic_Mutations/WUSM__IlluminaGA_DNASeq/Level_2/genome.wustl.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf',
                  sep='\t')
maf = maf[['Hugo_Symbol', 'Chrom','Ncbi_Build',
        'Start_Position', 'End_Position', 'Strand', 'Variant_Classification',
       'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1',
       'Tumor_Seq_Allele2', 'Dbsnp_Rs', 'Dbsnp_Val_Status',
       'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode',
       'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2',
       'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2',
       'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2',
        'Mutation_Status']]

maf['TCGA_ID'] = maf.Tumor_Sample_Barcode.str.slice(stop=12)

#remove silent mutations
maf = maf[maf.Variant_Classification != 'Silent']

#select driver genes
dmaf = maf[maf.Hugo_Symbol.isin(drivers)]

#make a table
dic = {'Missense_Mutation':2.0, 'Splice_Site':0.3, 'Frame_Shift_Del':0.1,
       'In_Frame_Del':1.8, 'Nonsense_Mutation':0.0, 'In_Frame_Ins':1.9,
       'Frame_Shift_Ins':0.2, 'RNA':1.6, 'Nonstop_Mutation':1.7}
dmaf['Variant_Value'] = dmaf['Variant_Classification'].map(dic)
table = dmaf.pivot_table(index='TCGA_ID',columns='Hugo_Symbol',
                         values='Variant_Value',
                         aggfunc=lambda x: np.sum(x)+10*(len(x)-1))
table.fillna(1.0,inplace=True)
missed = pd.Series(maf.TCGA_ID.unique().tolist())
missed = missed[~missed.isin(table.index)]
empty = pd.DataFrame(data=1.0,index=missed.tolist(),columns=table.columns)
table = pd.concat([table,empty])
table.to_csv('BRCA_somatic_muts_matrix.csv')
