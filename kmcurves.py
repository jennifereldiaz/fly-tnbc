#june 13 2018
#construction kaplan meier curves based on CNV genes from Liu, et al 2018
import pandas as pd
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
from lifelines.utils import median_survival_times
from matplotlib  import pyplot as plt
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests as mult

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'Arial','size': 18}),#'serif':['Times']})
font = {'family':'Arial','size': 16}

#consider pulling these in programmatically
#amps = ['AKIRIN1',
#'HEYL',
#'PPCS',
#'TM2D1',
#'INADL',
#'PVRL4',
#'GRHL1',
#'KLF11',
#'PRKCI',
#'TBL1XR1',
#'PIK3CA',
#'TRIO',
#'E2F3',
#'SOX4',
#'PTP4A1',
#'MYB',
#'EGFR',
#'AUTS2',
#'UBE3C',
#'DNAJB6',
#'ANKRD46',
#'GRHL2',
#'TRPS1',
#'MYC',
#'CCNE1',
#'ADAM15']
#
#dels=['GAB1',
#'MSRA',
#'STXBP1',
#'PTEN',
#'BCL9L',
#'RILPL2',
#'KATNAL1',
#'RB1',
#'ARF6']

hits = pd.concat([pd.read_excel('../../CNV_screen/tier0_postthesis_40-45_v2.xlsx',
                                 index_col='Unnamed: 0'),
                   pd.read_excel('../../CNV_screen/tier0_ambiguousdeletions_v2_thru45.xlsx',
                                 index_col='Unnamed: 0'),
                                 pd.read_excel('../../CNV_screen/tier1.5 CNV screen/tier1.5_postthesis_40-45.xlsx',
                      sheetname = 'tier1.5_sindura.csv',index_col='Unnamed: 0')])

#hits = pd.concat([hits1,hits2])
hits = hits[hits['Validation result'].astype(str).str.contains('\+')].sort_values(
        ['Chromosome','from'])

##side note--pull the mutsigcv q values
mutsig=pd.read_csv('../Somatic_mutations/MutSigCV2015.sig_genes.txt',sep='\t',
                   index_col='gene')
muts = mutsig[mutsig.index.isin(hits.index)]
print(muts[['q']].reindex(hits.index).reset_index(drop=False).drop_duplicates())
print(muts[['q']].reindex(hits.index).reset_index(drop=False).drop_duplicates().q.to_csv(index=False))

genes = hits.index.unique().tolist()

amps = hits[hits['CNV type']=='amp'].index.unique().tolist()
dels = hits[hits['CNV type']=='del'].index.unique().tolist()

pan = pd.read_excel('1-s2.0-S0092867418302290-mmc1.xlsx',sheetname='TCGA-CDR')
brca = pan[pan.type=='BRCA'].set_index('bcr_patient_barcode')

#cnv = pd.read_csv('../../running scripts/CNVruns/new_TN_CNV_run/BRCA_CNVs_foldchange_TN_filtered2.csv',
#                  index_col='TCGA_ID')
cnv = pd.read_csv('../../running scripts/CNVruns/new_TN_CNV_run/BRCA_CNVs_foldchange_all_filtered2.csv',
                  index_col='TCGA_ID')
cnv = cnv[genes]

brca = brca.merge(cnv,how='right',left_index=True,right_index=True)
#brca.to_csv('BRCA_CNVgenes_TN_survivalendpoints_UPDATEASNEEDED.csv')
brca.to_csv('BRCA_CNVgenes_survivalendpoints_UPDATEASNEEDED.csv')

brca = pd.read_csv('BRCA_CNVgenes_survivalendpoints_UPDATEASNEEDED.csv',index_col=0,header=0)
#brca = pd.read_csv('BRCA_CNVgenes_TN_survivalendpoints_UPDATEASNEEDED.csv',index_col=0,header=0)

amptest = pd.DataFrame(index=amps,data=[np.nan]*len(amps))
deltest = pd.DataFrame(index=dels,data=[np.nan]*len(dels))
ampmeds = amptest.copy()
delmeds = deltest.copy()

points = ['OS','PFI','DSS','DFI'] #,'DSS','DFI'

##RANDOMLY SHUFFLE FOR COMPARISON - COMMENT THIS OUT WHEN GETTING REAL RESULTS
#for time in points:
#    np.random.shuffle(brca[time+'.time'])

def lr(row,cutoff,endpoint):
    
    gene = row.name
    
    #set endpoint to use
    endtime = endpoint+'.time'

    test = brca.dropna(subset=[endpoint,endtime,gene])

    lrtest = logrank_test(test[test[gene]>cutoff][endtime],
                                               test[test[gene]<=cutoff][endtime],
                                               test[test[gene]>cutoff][endpoint],
                                               test[test[gene]<=cutoff][endpoint],
                                               alpha=0.99)

    return lrtest.p_value

def med(row,cutoff,endpoint):
    gene = row.name
    
    #set endpoint to use
    endtime = endpoint+'.time'

    test = brca.dropna(subset=[endpoint,endtime,gene])

    kmf = KaplanMeierFitter()

    kmf.fit(test[test[gene]<=cutoff][endtime],
            event_observed=test[test[gene]<=cutoff][endpoint])
    less = kmf.median_survival_time_
    kmf.fit(test[test[gene]>cutoff][endtime],
            event_observed=test[test[gene]>cutoff][endpoint])
    greater = kmf.median_survival_time_
    return greater-less

def cox(row,cutoff,endpoint):
    gene = row.name
    
    #set endpoint to use
    endtime = endpoint+'.time'

    test = brca.dropna(subset=[endpoint,endtime,gene])[[endpoint,endtime,gene]]
    
    if cutoff>1:
        test[gene]=test[gene].apply(lambda x: x>cutoff).replace({True:1,False:0})
    else:
        test[gene]=test[gene].apply(lambda x: x<cutoff).replace({True:1,False:0})

    cph = CoxPHFitter()
    cph.fit(test,duration_col=endtime,event_col=endpoint)
    
    return cph.hazard_ratios_.ix[gene]#'exp(coef)',

#test each gene by itself
a = 2**0.3
d = 2**-0.3
for end in points:
    print(end)
    amptest[end+'_logrank_p'] = amptest.apply(lr,cutoff=a,endpoint=end,axis=1)
    deltest[end+'_logrank_p'] = deltest.apply(lr,cutoff=d,endpoint=end,axis=1)
    amptest[end+'_median_diff'] = amptest.apply(med,cutoff=a,endpoint=end,axis=1)
    deltest[end+'_median_diff'] = deltest.apply(med,cutoff=d,endpoint=end,axis=1)
    amptest[end+'_HR'] = amptest.apply(cox,cutoff=a,endpoint=end,axis=1)
    deltest[end+'_HR'] = -deltest.apply(cox,cutoff=d,endpoint=end,axis=1)

##    #produce an overall cox model:
##    print(end)
##    cph = CoxPHFitter()
##    cph.fit(test[[endpoint]+[endtime]+amps+dels].dropna(),
##            duration_col=endtime,event_col=endpoint,show_progress=True,step_size=0.1)
##    cph.print_summary()

rslt = pd.concat([amptest,deltest])[[point+'_logrank_p' for point in points]+
                                    [point+'_median_diff' for point in points]+
                                    [point+'_HR' for point in points]]
#fdr each test individually - reasonable
for end in points:
    rslt[end+'_logrank_fdr'] = mult(rslt[end+'_logrank_p'],method='fdr_bh')[1]
#rslt.to_csv('survival_medians_logranktests_cnvgenes_TN.csv')
rslt.reindex(genes).to_csv('survival_medians_logranktests_cnvgenes.csv')

#this is to fdr all 4 kinds of tests together - more conservative
##unstacked = pd.DataFrame()
##for end in points:
##    unstacked = pd.concat([unstacked,pvals[[end]].set_index(pvals.index+'-'+end).rename(
##        columns={end:'p'})])
##unstacked['fdr'] = mult(unstacked.p,method='fdr_bh')[1]

def kmplot(gene,cutoff,endpoint,colors,tag):

    endtime = endpoint+'.time'
    test = brca.dropna(subset=[endpoint,endtime,gene])
    
    kmf = KaplanMeierFitter()
    plt.figure()
    ax = plt.subplot(111)

    kmf.fit(test[test[gene]<=cutoff][endtime],
            event_observed=test[test[gene]<=cutoff][endpoint],label='≤'+str(cutoff)[:4])
    kmf.plot(ax=ax,color = colors[0])

    kmf.fit(test[test[gene]>cutoff][endtime],
            event_observed=test[test[gene]>cutoff][endpoint],label='>'+str(cutoff)[:4])
    kmf.plot(ax=ax, color = colors[1])

    plt.ylim([0,1])
    plt.title(endpoint+' for '+gene)
    plt.xlabel('Time (days)')
    if endpoint=='OS':
        plt.ylabel('Survival')
    else:
        plt.ylabel('Fraction event-free')
    #plt.show()
    plt.tight_layout()
    plt.savefig(gene+'_'+endpoint+'_'+tag+'.pdf')



cols = [col for col in rslt.columns if 'logrank_p' in col]
rslt[(rslt[cols]<0.1).any(axis=1)]


#TN:
#        OS_logrank_p  PFI_logrank_p  DSS_logrank_p  DFI_logrank_p  \
#E2F3        0.183591       0.086201       0.258809       0.219509   
#SOX4        0.175186       0.065615       0.208570       0.138703   
#STXBP1      0.012682       0.127284       0.122741       0.086020   
#PTEN        0.036200       0.231327       0.116317       0.160458
#   
##I don't really believe these plots. It looks like an effect of small numbers. so still with the whole dataset analysis
#kmplot('E2F3',a,'PFI',('gray','r'),'TN')
#kmplot('SOX4',a,'PFI',('gray','r'),'TN')
#kmplot('STXBP1',d,'OS',('b','gray'),'TN')
#kmplot('PTEN',d,'OS',('b','gray'),'TN')

#        PFI_logrank_fdr  DSS_logrank_fdr  DFI_logrank_fdr  
#E2F3           0.864501          0.98276         0.973847  
#SOX4           0.864501          0.98276         0.973847  
#STXBP1         0.864501          0.98276         0.973847  
#PTEN           0.864501          0.98276         0.973847  
# 
#           OS_HR    PFI_HR    DSS_HR    DFI_HR  
#E2F3    0.470532  0.351289  0.419432  0.394332  
#SOX4    0.463583  0.327940  0.382721  0.331641  
#STXBP1 -5.381682 -2.351322 -3.174719 -3.566318  
#PTEN   -4.371333 -2.003874 -4.667613 -2.926956 

#         OS_logrank_p  PFI_logrank_p  DSS_logrank_p  DFI_logrank_p  \
#AKIRIN1      0.176707       0.005876       0.062265       0.008160   
#HEYL         0.352724       0.023095       0.228750       0.017141   
#PPCS         0.773287       0.011658       0.184049       0.157269   
#GRHL1        0.091250       0.107012       0.080021       0.192757   
#KLF11        0.102229       0.111394       0.086857       0.192499   
#SOX4         0.664678       0.245921       0.736185       0.076216   
#MYB          0.043401       0.131926       0.024732       0.882551   
#ANKRD46      0.054446       0.484441       0.044372       0.147595   
#GRHL2        0.064398       0.366083       0.048345       0.077426   
#TRPS1        0.031825       0.237945       0.023735       0.084968   
#MYC          0.039543       0.360481       0.019710       0.089503   
#MSRA         0.017381       0.385540       0.008645       0.049973  

kmplot('AKIRIN1',a,'PFI',('gray','r'),'all')
kmplot('HEYL',a,'PFI',('gray','r'),'all')
kmplot('PPCS',a,'PFI',('gray','r'),'all')
kmplot('RBM34',a,'OS',('gray','r'),'all')
kmplot('GRHL1',a,'OS',('gray','r'),'all')
#skip KLF - DSS and DFI are unreliable
kmplot('C6orf203',a,'OS',('gray','r'),'all')
kmplot('MYB',a,'OS',('gray','r'),'all')
kmplot('ANKRD46',a,'OS',('gray','r'),'all')
kmplot('GRHL2',a,'OS',('gray','r'),'all')
kmplot('TRPS1',a,'OS',('gray','r'),'all')
kmplot('MYC',a,'OS',('gray','r'),'all')
kmplot('MSRA',d,'OS',('b','gray'),'all')
kmplot('SERPINB8',d,'OS',('b','gray'),'all')
#kmplot('MSRA',a,'OS',('gray','r'),'all_amp')
#kmplot('MSRA',a,'PFI',('gray','r'),'all_amp')

#         PFI_logrank_fdr  DSS_logrank_fdr  DFI_logrank_fdr  
#AKIRIN1         0.198187         0.302432         0.277445  
#HEYL            0.261748         0.555536         0.291398  
#PPCS            0.198187         0.521473         0.546143  
#GRHL1           0.541057         0.328127         0.546143  
#KLF11           0.541057         0.328127         0.546143  
#SOX4            0.648719         0.834343         0.434728  
#MYB             0.560683         0.210221         0.974558  
#ANKRD46         0.648719         0.273955         0.546143  
#GRHL2           0.648719         0.273955         0.434728  
#TRPS1           0.648719         0.210221         0.434728  
#MYC             0.648719         0.210221         0.434728  
#MSRA            0.648719         0.210221         0.434728   
#
#            OS_HR    PFI_HR    DSS_HR    DFI_HR  
#AKIRIN1  1.444811  2.019048  1.865499  2.393684  
#HEYL     1.310729  1.847518  1.563176  2.280106  
#PPCS     1.091283  1.906732  1.596665  1.643713  
#GRHL1    1.582633  1.572157  1.847055  1.619827  
#KLF11    1.559040  1.563259  1.822085  1.620272  
#SOX4     1.103948  0.740498  1.108164  0.502084  
#MYB      1.660638  1.491670  2.006785  0.943099  
#ANKRD46  1.401133  1.130754  1.623986  1.400625  
#GRHL2    1.387621  1.174565  1.623142  1.519926  
#TRPS1    1.469788  1.232945  1.746612  1.495852  
#MYC      1.448980  1.176948  1.796522  1.496191  
#MSRA    -0.662768 -0.859999 -0.534009 -0.638937  

#15 genes out of 34 were significant on some measure. but no significant fdrs
#12 significnat on some measure for whole dataset. 10 if you exclude DSS and DFI. 
###this is probably the metric I care most about.

#after randomizing patient endpoints:
#15 from whole analysis - 12 if you exclude DSS and DFI.
#8 from TN analysis
# skip the TN analysis, it looks lik garbage anyway
 