#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2017/10/04     **
#**  last modified: 2018/01/15     **
#************************************

import sys
sys.path.append('../../code/python/') 
import classModel as cmo
import classReaction as cre
import classMetabolite as cme
import classPlotter as plotter
import classConvertionFactors as ccf
import cobra
import pandas
from cobra.flux_analysis import parsimonious, flux_variability_analysis
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.pyplot import cm 
import argparse
import plstyles
import random
import cPickle
import json
import seaborn as sb
from scipy.optimize import curve_fit


def main():
    '''
    Test growth matrix data on different media
    '''
    args = options()
    verbose = args.verbose
    iflist = args.infile.split(' ')

    if args.verbose:
        return loadModelAndDF(iflist[0])

    if len(iflist) < 2:
        return 

    t = multiSimulationAnalyses(args, iflist)
    
    return t

def multiSimulationAnalyses(args, infiles):
    '''
    Run simulation comparisons
    '''
    maxt = []
    dfs = []
    fnames = []
    pcgs = []
    flabels = args.labels.split(' ')
    ptitle = args.title
    for f in infiles:
        a, b, t = loadModelAndDF(f)
        c, d = findModelNames(f)
        maxt.append(t)
        dfs.append(b)
        fnames.append(c)
        pcgs.append(c[-1].split('_')[-1])
        cond = c[2]
    tM = max(np.array(maxt))
    tax = np.linspace(0, tM, 100)
    mk = a.dmodelsKeys
    #plotVaryRatio(args, df, ):
    print(dfs[0].columns)
    plotVaryRatio(args, tax, dfs, pcgs, cond, mk)
    return dfs[0]

def plotVaryRatio(args, tax, df, pcgs, expcond, mk):

    fig = plt.figure()
    gs = gridspec.GridSpec(2,4)
    # Top left
    ax1 = fig.add_subplot(gs[0,0:3])
    ax2 = fig.add_subplot(gs[1,0:2])
    ax3 = fig.add_subplot(gs[1,2:4])
    ax1.set_xlabel('Time (hr)')
    ax2.set_xlabel('Time (hr)')
    ax3.set_xlabel('Time (hr)')
    ax1.set_ylabel('Q [gDW]')
    ax2.set_ylabel('Q [gDW]')
    ax3.set_ylabel('Q [gDW]')
    ax1.set_title('Total biomass')
    ax2.set_title('%s biomass' % mk[0])
    ax3.set_title('%s biomass' % mk[1])

    colors = cm.rainbow(np.linspace(0, 1, len(pcgs)))
    #symbols = ['*', 'd', 's', 'o', 'x']
    symbols = ['-', ':', '--', '-.', 'x']
    
    #ax1.plot(tax, len(tax)*[np.nan])
    #ax2.plot(tax, len(tax)*[np.nan])
    #ax3.plot(tax, len(tax)*[np.nan])
    for i in range(len(pcgs)):
        bm0 = np.array(df[i]['B_%s' % mk[0]])
        bm1 = np.array(df[i]['B_%s' % mk[1]])
        t0 = np.array(df[i]['T_%s' % mk[0]])
        t1 = np.array(df[i]['T_%s' % mk[1]])
        ax1.plot(t0, bm0+bm1, symbols[i], markersize=3, color=colors[i], label='BM0 = %d%% ECgl' % int(100*float(pcgs[i].replace('p','.'))))
        ax2.plot(df[i]['T_%s' % mk[0]], bm0, symbols[i], markersize=3, color=colors[i])
        ax3.plot(df[i]['T_%s' % mk[1]], bm1, symbols[i], markersize=3, color=colors[i])

    ODtoGDW=0.33
    volExt = 0.03
    if expcond == "batch_low_Glc":
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2_fromEnjalbert2015.csv', sep=',')
        ax1.errorbar(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, yerr=dataFIG['OD SD']*ODtoGDW*volExt, fmt='bo', markerfacecolor='None', label='Biomass (exp)')
    elif expcond == "fedbatch_low_Ac":
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6a_fromEnjalbert2015.csv', sep=',')
        ax1.plot(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, 'bo', markerfacecolor='None', label='Biomass (exp)')
    elif expcond == "fedbatch_high_Ac":
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6b_fromEnjalbert2015_4h.csv', sep=',')
        ax1.plot(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, 'bo', markerfacecolor='None', label='Biomass (exp)')
    ax1.legend(loc='upper left', bbox_to_anchor=(1., 0.6), prop={'size': 10}) #, transform=ax3.transAxes)
    gs.tight_layout(fig)
    fig.savefig(args.pathout+'/vary_initial_bm_ratio_'+expcond+'.png')

    return

def loadModelAndDF(fname):

    model = cPickle.load(open(fname, 'r'))
    n1, n2 = model.dmodelsKeys
    bmn1 = 'biomass_%s' % n1
    bmn2 = 'biomass_%s' % n2
    conctomM = 1000. if model.dmodels[n1].volumes.externalUnits == 'mL' else 1.
    x = np.array(model.dmodels[n1].T)
    x2=np.array(model.dmodels[n2].T)
    ybm1 = np.array(model.dmodels[n1].dmetabolites[bmn1].quantity)
    ygl1 = np.array(model.dmodels[n1].dmetabolites['ex_glucose'].concentration)*conctomM
    yac1 = np.array(model.dmodels[n1].dmetabolites['ex_acetate'].concentration)*conctomM
    ybm2 = np.array(model.dmodels[n2].dmetabolites[bmn2].quantity)
    yac2 = np.array(model.dmodels[n2].dmetabolites['ex_acetate'].concentration)*conctomM
    ygl2 = np.array(model.dmodels[n2].dmetabolites['ex_glucose'].concentration)*conctomM
    df = pandas.DataFrame({'T_%s' % n1: x, 'B_%s' % n1: ybm1, 'G_%s' % n1: ygl1, 'A_%s' % n1: yac1,
                           'T_%s' % n2: x2, 'B_%s' % n2: ybm2, 'G_%s' % n2: ygl2, 'A_%s' % n2: yac2})
    return model, df, x2[-1]

def findModelNames(fname):
    '''
    gsmname = ifn[1]
    expcon = ifn[2]
    bmfname = ifn[3]
    '''
    ifn = fname.split('/')[-1].split('.')[0].split('-')
    aname = 'fig'
    aname = '-'.join(ifn[1:4])
    return ifn, aname


def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-E', '--exchange', help='check all exchange reactions', action='store_true')
    parser.add_argument('-B', '--biomass', help='check all exchange reactions', action='store_true')
    parser.add_argument('-F', '--fva', help='run fva', action='store_true')
    parser.add_argument('-T', '--tstepdist', help='plot time step distribution', action='store_true')
    parser.add_argument('-M', '--maxv', help='plot max V reactions', action='store_true')
    parser.add_argument('-D', '--deltas', help='plot delta V', action='store_true')
    parser.add_argument('-C', '--compare', help='compare fluxes', action='store_true')
    parser.add_argument('-G', '--compareGE', help='compare fluxes before after GE', action='store_true')
    parser.add_argument('-P', '--plots', help='plot V', action='store_true')
    parser.add_argument('-i', '--infile', help='file path to load', default='../../outputs/spentmedia/endOfSimulation_e_coli_core_expSpentMedia_glucose_media_BIOMASS_Ecoli_core_w_GAM.p')
    parser.add_argument('-p', '--pathout', help='path for outputs', default='../../outputs/spentmedia/')
    parser.add_argument('-r', '--rxns', help='list of reactions to check, e.g. "EX_glc__D_e EX_nh4_e"', default='')
    parser.add_argument('-m', '--mets', help='list of metabolites to check', default='glc__D_e nh4_e')
    parser.add_argument('-l', '--labels', help='labels for the simulations to compare', default='EC')
    parser.add_argument('-t', '--title', help='title for plot', default='')
    parser.add_argument('-k', '--tskip', help='time steps to skip to measure expr after GE', default='0')
    parser.add_argument('-y', '--ylims', help='tple of y ranges', default='(None, None)')
    args = parser.parse_args()
    if args.verbose:
        print "verbosity turned on"
        print args
    return args

if __name__=="__main__":
    main()
