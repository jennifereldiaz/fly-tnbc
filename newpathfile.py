#create new pathology file to differentiate tumors based on receptor IHC
#oct 7 2014

import pandas as pd

path = pd.read_csv('nationwidechildrens.org_clinical_patient_brca.txt',sep='\t',header=0,
                   skiprows=1)
path = path.dropna(how='all').dropna(how='all',axis=1)

path = path[['bcr_patient_barcode','breast_carcinoma_estrogen_receptor_status',
             'breast_carcinoma_progesterone_receptor_status',
             'lab_proc_her2_neu_immunohistochemistry_receptor_status']]
path.columns= ['TCGA ID','ER Status','PR Status','HER2 Status']

path = path.ix[1:]

path.to_csv('BRCA_pathology_2014.csv',index=False)
