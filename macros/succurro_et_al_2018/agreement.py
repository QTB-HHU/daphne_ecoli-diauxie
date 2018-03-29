#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2017/11/11     **
#**  last modified: 2018/01/15     **
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
from sklearn.metrics import auc
from matplotlib import gridspec
from scipy.stats import ks_2samp, linregress, chisquare

# http://bigg.ucsd.edu/models/e_coli_core
# http://bionumbers.hms.harvard.edu//bionumber.aspx?&id=107924&ver=8
### 0.47 gDW of E coli dissolved in 1L give a read of 1 OD600 unit
ODtoGDW=0.33
OD0= 0.219029
## VOLUME -- Enjalbert 2015 -- 50 mL and 200 mL flasks
volExt = 0.2
volUn = 'L'

## Fig6A : OD0 = 0.33
## Fig6C : OD0 = 0.5
expcondLabels = {'batch_high_Ac': 'grown on 45 mM Ac',
                 'batch_low_Ac': 'grown on 4 mM Ac',
                 'batch_mixed_Ac': 'grown on 15 mM Glc, 32 mM Ac',
                 'batch_low_Glc': 'grown on 15 mM Glc',
                 'fedbatch_low_Ac': 'grown on 15 mM Glc, fed 4 mM Ac',
                 'fedbatch_high_Ac': 'grown on 15 mM Glc, fed 32 mM Ac'}

def main():
    '''
    Reproduce results from
    ** Enjalbert et al 2015
    With Single EC model
    '''
    args = options()
    verbose = args.verbose

    biomass0 = float(args.biomassi)
    pcECgl = float(args.ratioecgl)
    thrPc = 0.99999
    if pcECgl > thrPc:
        pcECgl = thrPc
    if pcECgl < (1 - thrPc):
        pcECgl = 1 - thrPc
    strratioecgl = args.ratioecgl.replace('.','p')
    maxtime = float(args.tmax)
    expcond = 'batch_low_Glc'
    glucose0 = 15.0 #mM
    acetate0 = 0.0 #mM
    fb=0.
    afb=0.
    mpath = args.model
    print 'Loading simulation file'
    model = cPickle.load(open(mpath, 'r'))
    exitname = mpath.split('/')[-1].split('.')[0].split('-')
    expcond = exitname[2]

    agrBM = False
    agrGL = False
    agrAC = False
    stdd = False
    
    if args.runconsortium:
        n1, n2 = model.dmodelsKeys
        bmn1 = 'biomass_%s' % n1
        bmn2 = 'biomass_%s' % n2
        conctomM = 1000. if model.dmodels[n1].volumes.externalUnits == 'mL' else 1.
        ybm1 = np.array(model.dmodels[n1].dmetabolites[bmn1].quantity)
        ygl = np.array(model.dmodels[n1].dmetabolites['ex_glucose'].concentration)*conctomM
        yac = np.array(model.dmodels[n1].dmetabolites['ex_acetate'].concentration)*conctomM
        #yam1 = np.array(model.dmodels[n1].dmetabolites['ex_ammonium'].quantity)
        ybm2 = np.array(model.dmodels[n2].dmetabolites[bmn2].quantity)
        #yac2 = np.array(model.dmodels[n2].dmetabolites['ex_acetate'].concentration)*conctomM
        #ygl2 = np.array(model.dmodels[n2].dmetabolites['ex_glucose'].concentration)*conctomM
        ybm = ybm1+ybm2
        ty = np.array(model.dmodels[n1].T)
    else:
        conctomM = 1000. if model.volumes.externalUnits == 'mL' else 1.
        ybm = np.array(model.dmetabolites['biomass_ECgl'].quantity)
        x = np.array(model.T)
        ygl = np.array(model.dmetabolites['ex_glucose'].concentration)*conctomM
        yac = np.array(model.dmetabolites['ex_acetate'].concentration)*conctomM
        ty = np.array(model.T)
        
    if expcond == "varma_fedbatch":
        agrBM, agrGL, agrAC = True, True, True
        dataAC=pandas.read_csv('../../ecoli/varma1994_data/varma10_ac.csv', sep=',',header=0)
        dataBM=pandas.read_csv('../../ecoli/varma1994_data/varma10_bm.csv', sep=',',header=0)
        dataGL=pandas.read_csv('../../ecoli/varma1994_data/varma10_gl.csv', sep=',',header=0)
        ebm = np.array(dataBM['Curve1'])*volExt
        egl = np.array(dataGL['Curve1'])
        eac = np.array(dataAC['Curve1'])
        teb = np.array(dataBM['x'])
        teg = np.array(dataGL['x'])
        tea = np.array(dataAC['x'])
    elif expcond == "varma_batch":
        agrBM, agrGL, agrAC = True, True, True
        dataAC=pandas.read_csv('../../ecoli/varma1994_data/varma7_ac.csv', sep=',',header=0)
        dataBM=pandas.read_csv('../../ecoli/varma1994_data/varma7_bm.csv', sep=',',header=0)
        dataGL=pandas.read_csv('../../ecoli/varma1994_data/varma7_gl.csv', sep=',',header=0)
        ebm = np.array(dataBM['Curve1'])*volExt
        egl = np.array(dataGL['Curve1'])
        eac = np.array(dataAC['Curve1'])
        teb = np.array(dataBM['x'])
        teg = np.array(dataGL['x'])
        tea = np.array(dataAC['x'])
    elif expcond == "batch_low_Glc":
        agrBM, agrGL, agrAC, stdd = True, True, True, True
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2_fromEnjalbert2015.csv', sep=',')
        ebm = np.array(dataFIG['OD 600nm'])*ODtoGDW*volExt
        ebm_err = np.array(dataFIG['OD SD'])*ODtoGDW*volExt
        egl = np.array(dataFIG['Glucose mM'])
        egl_err = np.array(dataFIG['Glucose SD'])
        eac = np.array(dataFIG['Acetate mM'])
        eac_err = np.array(dataFIG['Acetate SD'])
        te = np.array(dataFIG['Time'])
    elif expcond == "fedbatch_low_Ac":
        agrBM = True
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6a_fromEnjalbert2015.csv', sep=',')
        ebm = np.array(dataFIG['OD 600nm'])*ODtoGDW*volExt
        te = np.array(dataFIG['Time'])
    elif expcond == "fedbatch_high_Ac":
        agrBM = True
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6b_fromEnjalbert2015_4h.csv', sep=',')
        ebm = np.array(dataFIG['OD 600nm'])*ODtoGDW*volExt
        te = np.array(dataFIG['Time'])

    # find time points in ty corresponding to te
    if "varma" in expcond:
        te = teb
    tidx = [np.argmin(np.abs(ty - i)) for i in te]
    #print ty[tidx]-te
    if agrBM:
        if "varma" in expcond:
            te = teb
            tidx = [np.argmin(np.abs(ty - i)) for i in te]
        sbm = ybm[tidx]
        nsim = sbm/auc(ty[tidx], sbm)
        nexp = ebm/auc(te, ebm)
        #print 'KS 2 sample test on Biomass: ', ks_2samp(ebm, sbm)
        #print 'KS 2 sample test on Biomass normalized to auc: ', ks_2samp(nexp, nsim)
        #print 'Chi sq for Biomass normalized to auc: ', chisquare(nexp, nsim)
        slope, intercept, r_value, p_value, std_err = linregress(nexp, nsim)
        print 'R sq for Biomass normalized to auc: ', r_value**2
        slope, intercept, r_value, p_value, std_err = linregress(ebm, sbm)
        print 'R sq for Biomass: ', r_value**2
        fig = plt.figure()
        plt.plot(ebm, sbm, "o", label="Simulated vs Measured data")
        plt.plot(ebm, intercept + slope*ebm, 'r', label='fitted line')
        plt.xlabel("Measured")
        plt.ylabel("Simulated")
        plt.legend()
        plt.title(r'Correlation for Biomass data $R^2$ = %.3f' % (r_value**2))
        fig.savefig('%s/rSq_%s-%s_bm.png'%(args.pathout, expcond, exitname[3]))
    if agrGL:
        if "varma" in expcond:
            te = teg
            tidx = [np.argmin(np.abs(ty - i)) for i in te]
        sgl = ygl[tidx]
        nsim = sgl/auc(ty[tidx], sgl)
        nexp = egl/auc(te, egl)
        #print 'KS 2 sample test on Glucose: ', ks_2samp(egl, sgl)
        #print 'KS 2 sample test on Glucose normalized to auc: ', ks_2samp(nexp, nsim)
        slope, intercept, r_value, p_value, std_err = linregress(nexp, nsim)
        print 'R sq for Glucose normalized to auc: ', r_value**2
        #print 'Chi sq for Glucose normalized to auc: ', chisquare(nexp, nsim)
        slope, intercept, r_value, p_value, std_err = linregress(egl, sgl)
        print 'R sq for Glucose: ', r_value**2
        fig = plt.figure()
        plt.plot(egl, sgl, "o", label="Simulated vs Measured data")
        plt.plot(egl, intercept + slope*egl, 'r', label='fitted line')
        plt.xlabel("Measured")
        plt.ylabel("Simulated")
        plt.legend()
        plt.title(r'Correlation for Glucose data $R^2$ = %.3f' % (r_value**2))
        fig.savefig('%s/rSq_%s-%s_gl.png'%(args.pathout, expcond, exitname[3]))
    if agrAC:
        if "varma" in expcond:
            te = tea
            tidx = [np.argmin(np.abs(ty - i)) for i in te]
        sac = yac[tidx]
        nsim = sac/auc(ty[tidx], sac)
        nexp = eac/auc(te, eac)
        #print 'KS 2 sample test on Acetate: ', ks_2samp(eac, sac)
        #print 'KS 2 sample test on Acetate normalized to auc: ', ks_2samp(nexp, nsim)
        #print 'Chi sq for Acetate normalized to auc: ', chisquare(nexp, nsim)
        slope, intercept, r_value, p_value, std_err = linregress(nexp, nsim)
        print 'R sq for Acetate normalized to auc: ', r_value**2
        slope, intercept, r_value, p_value, std_err = linregress(eac, sac)
        print 'R sq for Acetate: ', r_value**2
        fig = plt.figure()
        plt.plot(eac, sac, "o", label="Simulated vs Measured data")
        plt.plot(eac, intercept + slope*eac, 'r', label='fitted line')
        plt.xlabel("Measured")
        plt.ylabel("Simulated")
        plt.legend()
        plt.title(r'Correlation for Acetate data $R^2$ = %.3f' % (r_value**2))
        fig.savefig('%s/rSq_%s-%s_ac.png'%(args.pathout, expcond, exitname[3]))

    return

def chisq(obs, exp):
    
    return

def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-A', '--runhighacetate', help='batch 45mM Acetate conditions', action='store_true')
    parser.add_argument('-B', '--runlowacetate', help='batch 4mM Acetate conditions', action='store_true')
    parser.add_argument('-D', '--runglucose', help='batch 15mM Glucose conditions (default)', action='store_true')
    parser.add_argument('-G', '--runmixedacetate', help='batch 15mM Glucose 32mM Acetate conditions', action='store_true')
    parser.add_argument('-F', '--runfedlowacetate', help='batch 15mM Glucose fedbatch 4mM Acetate conditions', action='store_true')
    parser.add_argument('-H', '--runfedhighacetate', help='batch 15mM Glucose fedbatch 32mM Acetate conditions', action='store_true')
    parser.add_argument('-I', '--runvarmabatch', help='batch Glucose conditions as in Varma Fig 7', action='store_true')
    parser.add_argument('-K', '--runvarmafedbatch', help='fed-batch Glucose conditions as in Varma Fig 10', action='store_true')
    parser.add_argument('-E', '--extraplots', help='produce extra plots', action='store_true')
    parser.add_argument('-C', '--coremodel', help='use core model', action='store_true')
    parser.add_argument('-X', '--runecgl', help='run ECgl', action='store_true')
    parser.add_argument('-Y', '--runecac', help='run ECac', action='store_true')
    parser.add_argument('-Z', '--runconsortium', help='run consortium', action='store_true')
    parser.add_argument('-S', '--runsingle', help='run single', action='store_true')
    parser.add_argument('-T', '--psitransition', help='activate psi transition from Glc state to Ac state', action='store_true')
    parser.add_argument('-P', '--phitransition', help='activate phi transition from Ac state to Glc state', action='store_true')
    parser.add_argument('-R', '--run', help='run the dFBA', action='store_true')
    parser.add_argument('-J', '--death', help='add death rate', action='store_true')
    parser.add_argument('-M', '--moma', help='run dFBA with MOMA', action='store_true')
    parser.add_argument('-b', '--biomassi', help='initial biomass', default='0.001')
    parser.add_argument('-d', '--deathrate', help='death rate', default='0.01')
    parser.add_argument('-e', '--efftrans', help='transition efficiency', default='1.')
    parser.add_argument('-j', '--kmtranspsi', help='Km psi transition', default='0.')
    parser.add_argument('-k', '--kmtransphi', help='Km phi transition', default='0.')
    parser.add_argument('-l', '--label', help='output file label', default='diauxic_shift')
    parser.add_argument('-m', '--model', help='model path to load', default='../../ecoli/bigg/e_coli_core.xml')
    parser.add_argument('-n', '--nsteps', help='max number of steps', default='1000')
    parser.add_argument('-o', '--ofile', help='output file name', default='toymodel-test.xml')
    parser.add_argument('-p', '--pathout', help='path for outputs', default='/tmp')
    parser.add_argument('-r', '--ratioecgl', help='ratio of ECgl/ECac: BMECgl0 = X*biomass0', default='1.')
    parser.add_argument('-s', '--minstep', help='min step size', default='0.0')
    parser.add_argument('-t', '--tmax', help='max time', default='10.5')
    parser.add_argument('-x', '--vminexoxygen', help='V min EX O2', default='10.')
    parser.add_argument('-y', '--vmaxpsi', help='V max psi transition', default='0.0')
    parser.add_argument('-z', '--vmaxphi', help='V max phi transition', default='0.0')
    parser.add_argument('-w', '--psioffset', help='offset psi transition', default='0.05')
    parser.add_argument('-u', '--phioffset', help='offset phi transition', default='0.05')
    args = parser.parse_args()
    if args.verbose:
        print "verbosity turned on"
        print args
    return args

if __name__=="__main__":
    main()
