#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2018/05/15     **
#**  last modified: 2018/05/15     **
#************************************

import sys
sys.path.append('../../code/python/') 
import classModel as cmo
import classReaction as cre
import classMetabolite as cme
import classConsortium as cco
import classPlotter as plotter
import classConvertionFactors as ccf
import cobra
import pandas
from cobra.flux_analysis import parsimonious
import numpy as np
import matplotlib.pyplot as plt
import argparse
import plstyles
import random
import cPickle
import json
import copy
from matplotlib import gridspec

def plotVarma1994_twoEC(args, model, expcond, expcond2, tit):

    n1, n2 = model.dmodelsKeys
    
    bmn1 = 'biomass_%s' % n1
    bmn2 = 'biomass_%s' % n2

    # LOAD DATA from Varma1994
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1)
    # Top left
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
    ax1.set_ylabel('Q [gDW]')
    ax2.set_ylabel('C [mM]')
    ax2.set_xlabel('Time [hours]')
    if args.phitransition and args.psitransition:
        tit=tit+', w/ transition'
    else:
        tit=tit+', w/o transition'
    
    if args.runvarmafedbatch:
        dataAC=pandas.read_csv('../../ecoli/varma1994_data/varma10_ac.csv', sep=',',header=0)
        dataBM=pandas.read_csv('../../ecoli/varma1994_data/varma10_bm.csv', sep=',',header=0)
        dataGL=pandas.read_csv('../../ecoli/varma1994_data/varma10_gl.csv', sep=',',header=0)
        fignum = '_fig10'
    else:
        dataAC=pandas.read_csv('../../ecoli/varma1994_data/varma7_ac.csv', sep=',',header=0)
        dataBM=pandas.read_csv('../../ecoli/varma1994_data/varma7_bm.csv', sep=',',header=0)
        dataGL=pandas.read_csv('../../ecoli/varma1994_data/varma7_gl.csv', sep=',',header=0)
        fignum = '_fig7'
    ax1.set_title(tit)
    conctomM = 1000. if model.dmodels[n1].volumes.externalUnits == 'mL' else 1.

    x = np.array(model.dmodels[n1].T)
    x2=np.array(model.dmodels[n2].T)
    ybm1 = np.array(model.dmodels[n1].dmetabolites[bmn1].quantity)
    ygl1 = np.array(model.dmodels[n1].dmetabolites['ex_glucose'].concentration)*conctomM
    yac1 = np.array(model.dmodels[n1].dmetabolites['ex_acetate'].concentration)*conctomM
    #yam1 = np.array(model.dmodels[n1].dmetabolites['ex_ammonium'].quantity)
    ybm2 = np.array(model.dmodels[n2].dmetabolites[bmn2].quantity)
    yac2 = np.array(model.dmodels[n2].dmetabolites['ex_acetate'].concentration)*conctomM
    ygl2 = np.array(model.dmodels[n2].dmetabolites['ex_glucose'].concentration)*conctomM
    #yam2 = np.array(model.dmodels[n2].dmetabolites['ex_ammonium'].quantity)

    ax1.plot(x, ybm1, '#2980b9', label=bmn1.replace('_', ' '))
    ax1.plot(x, ybm2, '#16a085', label=bmn2.replace('_', ' '))
    ax1.plot(x, ybm1+ybm2, 'k', label='biomass Tot')
    #ax1.set_xlim(0., 1.)
    ll = ax1.legend(loc='center right', prop={'size':10})
    ax1.set_title(tit)
    ax2.plot(x, ygl1, '#c0392b', label='Glucose')
    ax2.plot(x, yac1, '#f39c12', label='Acetate')

    ax1.plot(dataBM['x'], dataBM['Curve1']*volExt, 'bs', label='Biomass (exp)')
    ax2.plot(dataGL['x'], dataGL['Curve1'], 'rs', label='Glucose (exp)')
    ax2.plot(dataAC['x'], dataAC['Curve1'], '#f1c40f',  linestyle='None', marker='s', label='Acetate (exp)')

    ll = ax1.legend(loc='best', prop={'size':10})
    ll2 = ax2.legend(loc='best', prop={'size':10})
    fig.savefig(args.pathout+'/varma1994_'+expcond2+fignum+'.png')
    
    return

def plotEnjalbert2015_growth2EC(transitions, model, expcond, expcond2, tit):

    n1, n2 = model.dmodelsKeys
    
    bmn1 = 'biomass_%s' % n1
    bmn2 = 'biomass_%s' % n2

    conctomM = 1000. if model.dmodels[n1].volumes.externalUnits == 'mL' else 1.

    if transitions:
        tit=tit+', w/ transition'
    else:
        tit=tit+', w/o transition'
    
    x = np.array(model.dmodels[n1].T)
    x2=np.array(model.dmodels[n2].T)
    ybm1 = np.array(model.dmodels[n1].dmetabolites[bmn1].quantity)
    ygl1 = np.array(model.dmodels[n1].dmetabolites['ex_glucose'].concentration)*conctomM
    yac1 = np.array(model.dmodels[n1].dmetabolites['ex_acetate'].concentration)*conctomM
    #yam1 = np.array(model.dmodels[n1].dmetabolites['ex_ammonium'].quantity)
    ybm2 = np.array(model.dmodels[n2].dmetabolites[bmn2].quantity)
    yac2 = np.array(model.dmodels[n2].dmetabolites['ex_acetate'].concentration)*conctomM
    ygl2 = np.array(model.dmodels[n2].dmetabolites['ex_glucose'].concentration)*conctomM
    #yam2 = np.array(model.dmodels[n2].dmetabolites['ex_ammonium'].quantity)

    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1)
    # Top left
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
    ax1.set_ylabel('Q [gDW]')
    ax2.set_ylabel('C [mM]')
    ax2.set_xlabel('Time [hours]')

    ## COLORS http://htmlcolorcodes.com/color-chart/
    ## blues: #2980b9 #3498db  #1abc9c #16a085
    ## yellows: #f4d03f f5b041 eb984e  #dc7633 
    ax1.plot(x, ybm1, '#2980b9', label=bmn1.replace('_', ' '))
    ax1.plot(x, ybm2, '#16a085', label=bmn2.replace('_', ' '))
    #ax1.set_xlim(0., 1.)
    ll = ax1.legend(loc='center right', prop={'size':10})
    ax1.set_title(tit)
    ax2.plot(x, ygl1, '#c0392b', label='Glucose')
    ax2.plot(x, yac1, '#f39c12', label='Acetate')
    # ax2.plot(x, ygl1, '#3498db', label='Glucose 1')
    # ax2.plot(x, yac1, '#1abc9c', label='Acetate 1')
    # ax2.plot(x, ygl2, '#f5b041', label='Glucose 2')
    # ax2.plot(x, yac2, '#eb984e', label='Acetate 2')
    ax1.plot(x, ybm1+ybm2, 'k', label='biomass Tot')

    fignum = ''
    if expcond == "batch_low_Glc":
        # dataAC=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2A_ac.csv', sep=',',header=None, names=['x','y'])
        # dataBM=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2A_bm.csv', sep=',',header=None, names=['x','y'])
        # dataGL=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2A_gl.csv', sep=',',header=None, names=['x','y'])
        # fignum = '_fig2A'
        # ax1.plot(dataBM['x'], dataBM['y']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        # ax2.plot(dataGL['x'], dataGL['y'], 'rs', label='Glucose (exp)')
        # ax2.plot(dataAC['x'], dataAC['y'], '#f1c40f',  linestyle='None', marker='s', label='Acetate (exp)')
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2_fromEnjalbert2015.csv', sep=',')
        ax1.errorbar(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, yerr=dataFIG['OD SD']*ODtoGDW*volExt, fmt='bs', label='Biomass (exp)')
        ax2.errorbar(dataFIG['Time'], dataFIG['Glucose mM'], yerr=dataFIG['Glucose SD'], fmt='rs', label='Glucose (exp)')
        ax2.errorbar(dataFIG['Time'], dataFIG['Acetate mM'], yerr=dataFIG['Acetate SD'], label='Acetate (exp)', color='#f1c40f',  linestyle='None', marker='s')
        fignum = '_fig2A'
    elif expcond == "fedbatch_low_Ac":
        # dataBM=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6a_bm.csv', sep=',',header=None, names=['x','y'])
        # ax1.plot(dataBM['x'], dataBM['y']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        #ax1.errorbar(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, yerr=dataFIG['OD SD']*ODtoGDW*volExt, fmt='bs', label='Biomass (exp)')
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6a_fromEnjalbert2015.csv', sep=',')
        ax1.plot(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        fignum = '_fig6A'
    elif expcond == "fedbatch_high_Ac":
        # dataBM=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6c_bm.csv', sep=',',header=None, names=['x','y'])
        # ax1.plot(dataBM['x'], dataBM['y']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6b_fromEnjalbert2015_4h.csv', sep=',')
        ax1.plot(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        fignum = '_fig6C'
    #ax.plot(x, ybm1, ':', label='BM1')
    #ax.plot(x, ybm2, ':', label='BM2')
    #ll = ax1.legend(loc='upper left', prop={'size':10})
    #ll2 = ax2.legend(loc='upper left', prop={'size':10})
    ll = ax1.legend(loc='best', prop={'size':10})
    ll2 = ax2.legend(loc='best', prop={'size':10})
    plt.setp(ax1.get_xticklabels(), visible=False)
    fig.savefig('enjalbert2015-%s-totBM%s.png' % (expcond2, fignum))

    #f1.set_title('E. Coli core Aerobic '+expcond.capitalize())
    #ll = f1.legend(loc='center left', prop={'size':10})
    #ll = f1b.legend(loc='center right', prop={'size':10})
    #fig1.savefig(args.pathout+'/enjalbert2015_'+expcond2+fignum+'.png')   
    return


ODtoGDW=0.33
volExt = 0.03
volUn = 'L'

mpath = 'endOfSimulation-ecoli_core-fedbatch_low_Ac-diauxic_shift_ECgl_ECac_0p95.p'

thrPc = 0.99999
pcECgl = 0.95
strratioecgl = ("%.2f" % pcECgl).replace('.','p')
if pcECgl > thrPc:
    pcECgl = thrPc
if pcECgl < (1 - thrPc):
    pcECgl = 1 - thrPc

expcond = 'batch_low_Glc'
glucose0 = 15.0 #mM
acetate0 = 0.0 #mM
fb=0.
afb=0.
ac_thr = 0.0
t_glc = 0.0

runhighacetate = runlowacetate = runbatchglc = runmixedacetate = runfedlowacetate = runfedhighacetate = runvarmabatch = runvarmafedbatch =False
runfedlowacetate = True

if runhighacetate:
    print('Simulating Batch Growth on 45mM Acetate as in Enjalbert2015')
    expcond = 'batch_high_Ac'
    glucose0 = 0.0
    acetate0 = 45.0
elif runlowacetate:
    print('Simulating Batch Growth on 4mM Acetate as in Enjalbert2015')
    expcond = 'batch_low_Ac'
    glucose0 = 0.0
    acetate0 = 4.0
elif runbatchglc:
    print('Simulating Batch Growth on 15mM Glucose as in Enjalbert2015')
    expcond = 'batch_low_Glc'
    glucose0 = 15.0 #mM
    acetate0 = 0.0 #mM
    fb=0.
    afb=0.
elif runmixedacetate:
    print('Simulating Batch Growth on 15mM Glucose 32mM Acetate as in Enjalbert2015')
    expcond = 'batch_mixed_Ac'
    glucose0 = 15.0
    acetate0 = 32.0
elif runfedlowacetate:
    print('Simulating Batch Growth on 15mM Glucose and constant 4mM Acetate as in Enjalbert2015')
    expcond = 'fedbatch_low_Ac'
    glucose0 = 15.0
    acetate0 = 0.0
    afb=1.0
    ac_thr = 4.0
    t_glc = 4.0
elif runfedhighacetate:
    print('Simulating Batch Growth on 15mM Glucose and constant 32mM Acetate as in Enjalbert2015')
    expcond = 'fedbatch_high_Ac'
    glucose0 = 15.0
    acetate0 = 32.0
    afb=1.0
    ac_thr = 32.0
    t_glc = 4.0
elif runvarmabatch:
    print('Simulating Batch Growth on Glucose as in Varma1994 Fig.7')
    expcond = 'varma_batch'
    glucose0 = 10.8
    acetate0 = 0.4
elif runvarmafedbatch:
    print('Simulating Fed Batch Growth on Glucose as in Varma1994 Fig.10')
    expcond = 'varma_fedbatch'
    glucose0 = 0.82
    acetate0 = 0.1
    #glucose is continuously provided
    #0.2 g/(L*hr)
    #1 g = 5.5 mmol
    #dGLC/dt = ... + fb ==> units = mmol/hr
    fb = 5.5*0.2*volExt #mmol/hr
else:
    print('Simulating custom conditions')

label='diauxic_shift'
ename = 'ecoli core'
exitname = '%s-%s-%s' % (ename, expcond, label) 

print 'Loading simulation file'
model = cPickle.load(open(mpath, 'r'))
exitname = mpath.split('/')[-1].split('.')[0].split('-')

titlab = 'Initial biomass %.0f%% ECgl and %.0f%% ECac' % (pcECgl*100, (1-pcECgl)*100)
if runvarmabatch or runvarmafedbatch:
    plotVarma1994_twoEC(args, model, exitname[2], '-'.join(exitname[1:]), titlab)
else:
    plotEnjalbert2015_growth2EC(True, model, exitname[2], '-'.join(exitname[1:]), titlab)
  

