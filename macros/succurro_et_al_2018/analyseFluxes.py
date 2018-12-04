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
import argparse
import plstyles
import random
import cPickle
import json
import seaborn as sb
from scipy.optimize import curve_fit

def logisticGrowth(t, n0, r, k):
    n = k/(1 + (k-n0)*np.exp(-r*t)/n0)
    return n

def findGlExhTime(m0):
    ''' 
    Consider Glc exhausted when no more is imported
    '''
    #i0 = (np.abs(np.array(m0.dmetabolites['ex_glucose'].quantity) - 0.1)).argmin()
    #i1 = m0.dreactions['growth'].flux.index(0)

    ##Consider Glc exhausted when acetate starts getting consumed
    #acex = np.array(m0.dreactions['acetate_exchange'].flux)
    #i0 = np.argmax(acex < 0)

    i0 = np.argmin( np.abs(np.array(m0.dreactions['glucose_exchange'].flux)) )
    
    #if i0 != i1:
    #    print('ACHTUNG!! ',i0,i1)
    strpr = 'At simulation step %d (time %.3f) Glucose is %.3f mmol; glucose exchange is %.3f; acetate exchange is %.3f ' % (i0, m0.T[i0], m0.dmetabolites['ex_glucose'].quantity[i0], m0.dreactions['glucose_exchange'].flux[i0], m0.dreactions['acetate_exchange'].flux[i0])
    print(strpr)
    return i0
# http://bigg.ucsd.edu/models/e_coli_core
# ecoli = '../../ecoli/bigg/e_coli_core.xml'
# iJO1366 Escherichia coli
# http://bigg.ucsd.edu/models/iJO1366
# ecoli = '../../ecoli/bigg/iJO1366.xml'
# iJR904 Escherichia coli
# http://bigg.ucsd.edu/models/iJR904
# ecoli = '../../ecoli/bigg/iJR904.xml'

def main():
    '''
    Test growth matrix data on different media
    '''
    args = options()
    verbose = args.verbose
    iflist = args.infile.split(' ')

    if args.verbose:
        return loadModelAndDF(iflist[0])

    if len(iflist) > 1:
        return multiSimulationAnalyses(args, iflist)
        
        
    for ifl in iflist:
        singleSimulationAnalyses(args, ifl)
        
    return

def multiSimulationAnalyses(args, infiles):
    '''
    Run simulation comparisons
    '''
    maxt = []
    dfs = []
    fnames = []
    flabels = args.labels.split(' ')
    ptitle = args.title
    tGE = []
    #tiGE = []
    for f in infiles:
        a, b = loadModelAndDF(f)
        c, d = findModelNames(f)
        maxt.append(a.T[-1])
        dfs.append(b)
        fnames.append(c)
        if args.compareGE:
            iGE = findGlExhTime(a)
            tGE.append(a.T[iGE])
            if args.compareGE:
                tArr = np.array(a.T)
                iPreGE = (np.abs(tArr-(tGE[-1] - 1.9))).argmin()
                iPostGE =  (np.abs(tArr-(tGE[-1] + float(args.tskip)))).argmin()
                with open('%s/sol_%s_GE.json' % (args.pathout, d), 'w') as f:
                    json.dump(a.FBAsolutions[iGE].to_dict(), f)
                with open('%s/sol_%s_%dmin_before_GE.json' % (args.pathout, d, 1.9*60), 'w') as f:
                    json.dump(a.FBAsolutions[iPreGE].to_dict(), f)
                with open('%s/sol_%s_%dmin_after_GE.json' % (args.pathout, d, float(args.tskip)*60), 'w') as f:
                    json.dump(a.FBAsolutions[iPostGE].to_dict(), f)

            #print('GE time for ', f, ':', tGE[-1])
            #tiGE.append((np.abs(np.array(a.T)-tGE[-1])).argmin())
    tM = max(np.array(maxt))
    tax = np.linspace(0, tM, 100)

    rxnsid = args.rxns.split(' ')
    if len(rxnsid[0]) < 1:
        rxnsid = []

    # tind = (np.abs(np.array(a.T)-tGE)).argmin()
    # tPre = tind-5
    # tPost = tind+5
    # print('GE @ ',tGE,'; check times ', a.T[tPre], ' and ', a.T[tPost])
    
    ndfs = []
    for d in dfs:
        ndfs.append(alignDf(d, tax))

    if len(flabels) != len(ndfs):
        flabels = ['-'.join(x[1:4]) for x in fnames]
        
    rxndfs = []
    for r in rxnsid:
        rxndfs.append( mergeRxnVs(flabels, ndfs, r, tax) )
        #rxndfs.append( alignRxnVs(fnames, dfs, r, tax) )
        if args.plots:
            cnames = list(rxndfs[-1].columns)[1:]
            plotFromDF(rxndfs[-1], 'time', cnames, r, 'Time [hr]', 'V [mmol/(gDW*hr)]', '%s/v-%s_%s.png' % (args.pathout, r, '-'.join(cnames)))
    if args.plots and args.compare:
        ptit = ptitle if len(ptitle) > 0 else 'Fluxes'
        plotFromTwoDF(ndfs, flabels, 'time', rxnsid, ptit, 'Time [hr]', 'V [mmol/(gDW*hr)]', '%s/%s-%s_%s.png' % (args.pathout, ptit.replace(' ','-'), '-'.join(rxnsid), '-'.join(cnames)))
    if args.compare:
        #return rxndfs, ndfs
        rx = printVComp(rxnsid, rxndfs, flabels, args.pathout, args.compareGE, tGE, float(args.tskip))
        if args.compareGE:
            geneExprGE(args.pathout)
        else:
            geneExprFullGrowth(args.pathout)
    #growth_arr = np.array([ vgrowth[(np.abs(tmod-t)).argmin()] if t < tmod[-1] else 0. for t in tarr ])
    #(np.abs(df[xname]-t)).argmin()
    return ndfs

def printVComp(rs, dfs, ls, pt, doGE, tGE, tskip=0.0):
    '''
    Get a condition x reactions matrix
    '''
    #for d in dfs:
    #    d['Ratio'] = d[ls[0]]/d[ls[1]]
    tit = 'Normalized fluxes during exponential growth'
    if doGE:
        tit = 'Flux changes upon switch from glucose to acetate'
    mx = np.zeros((len(ls), len(rs)))
    for i in range(len(ls)):
        for j in range(len(rs)):
            if doGE:
                tiGE = (np.abs(np.array(dfs[j]['time'])-tGE[i])).argmin()
                ## Expr levels relative to GE - 115 min
                tiPreGE = (np.abs(np.array(dfs[j]['time'])-(tGE[i] - 1.9))).argmin()
                #tiPostGE = tiGE + tskip
                tiPostGE =  (np.abs(np.array(dfs[j]['time'])-(tGE[i] + tskip))).argmin()
                if j==0:
                    print tskip, np.array(dfs[j]['time'])[tiPreGE], np.array(dfs[j]['time'])[tiGE], np.array(dfs[j]['time'])[tiPostGE]
                    #print('For ', ls[i], ' GE time: ', tGE[i], ' after alignment: ', dfs[j].loc[tiGE, 'time'],'; Post GE time delta: ', dfs[j].loc[tiPostGE, 'time'] - tGE[i])
                VPostGE = dfs[j].loc[tiPostGE, ls[i]]
                VPreGE  = dfs[j].loc[tiPreGE, ls[i]]
                #mx[i,j] = dfs[j].loc[tiPostGE, ls[i]] - VPreGE
                #mx[i,j] = VPostGE/VPreGE
                dV = VPostGE - VPreGE
                absdV = np.abs(VPostGE) - np.abs(VPreGE)
                maxabsV = 0. #max(np.abs(VPostGE), np.abs(VPreGE))
                sV = VPostGE + VPreGE
                #print(rs[j], ls[i], VPreGE, VPostGE, dV, sV, tiPreGE, tiPostGE)
                mx[i,j] = absdV/maxabsV if maxabsV > 0 else absdV
            else:
                mx[i,j] = np.abs(dfs[j].loc[5,ls[i]])
    #plt.subplots(figsize=(20,5))
    # if doGE:
    #     mx = np.append([expupdo], mx, axis=0)
    #     #print(mx)
    #     ls = ['GeneExpData']+ls
    rx = pandas.DataFrame(mx, index=ls, columns=rs)
    if doGE:
        ax = sb.heatmap(rx, cmap="RdBu_r", center=0, linewidths=.5, cbar=True, cbar_kws={"shrink": .45}, square=True, vmin=-7, vmax=5)
        ax.set_title('Flux change from GE - 1.9 hr to GE + %.2f hr'%tskip)
        plt.yticks(rotation=0)
        plt.savefig('%s/geneexp_v_comp_%s_%s_tskip%.2f.png' % (pt, '-'.join(ls), '-'.join(rs), tskip))
        makePaperFigureGluExh(rx, tit, pt)
    else:
        normmx = mx/mx.sum(axis=0)[np.newaxis,:]
        #normmx = mx/(mx).max()#/mx[1,:]
        rxn = pandas.DataFrame(normmx, index=ls, columns=rs)
        ax = sb.heatmap(rxn, cmap="Greens", linewidths=.5, cbar=False, cbar_kws={"shrink": .45}, square=True)
        #ax = sb.heatmap(rxn, cmap="Greens", linewidths=.5, cbar=False, square=True)
        ax.set_title(tit)
        plt.yticks(rotation=0)
        plt.savefig('%s/geneexp_norm_v_comp_%s_%s_tskip%.2f.png' % (pt, '-'.join(ls), '-'.join(rs), tskip))
        makePaperFigureExpGro(rxn, tit, pt)
    plt.clf()
    return rx


def makePaperFigureGluExh(rxn, tit, pt):
    #rxns="ACKr PPCK PPS FBP ICL MALS PFK PYK PPC ICDHyr".split(' ')
    #rxns="ACKr PPCK PPS FBP ICL MALS PFK PYK PPC ICDHyr".split(' ')


    fig1 = plt.figure()
    gs = gridspec.GridSpec(2, 4)

    ax1 = fig1.add_subplot(gs[0,0:3])
    mx = np.array([[1, 1, 0, 0, 0, 0, -1, -1, -1, -1]])
    rxns="acs pck pps fbp icl mls pfkA pykF ppc icd".split(' ')
    ls=['']
    rxn0 = pandas.DataFrame(mx, index=ls, columns=rxns)
    ax1 = sb.heatmap(rxn0, cmap="RdBu_r", center=0, linewidths=.5, cbar=False, square=True)#, cbar=True, cbar_kws={"shrink": .25}, square=True)
    ax1.set_title('Qualitative gene expression after glucose exhaustion')
    ax1.set_yticklabels(['Experim'], rotation=0)
    ttl1 = ax1.title
    ttl1.set_position([.5, 1.3])


    ax3 = fig1.add_axes([.8, .3, .03, .4])

    ax2 = fig1.add_subplot(gs[1,0:3])
    sb.heatmap(rxn, ax=ax2, cmap="RdBu_r", center=0, linewidths=.5, cbar=True, square=True, cbar_ax=ax3)#, vmin=-5, vmax=5)
    ax2.set_title(tit)
    ax2.set_yticklabels(rxn.index, rotation=0)
    plt.xticks(rotation=20)
    ttl2 = ax2.title
    ttl2.set_position([.5, 1.3])
    
    plt.savefig('%s/gene_exp_glucose_exhaustion.png' % (pt))
    plt.savefig('../../outputs/pqfigures/fig2b.png', format='png', dpi=1500)

    return

def makePaperFigureExpGro(rxn, tit, pt):



    fig1 = plt.figure()
    gs = gridspec.GridSpec(2, 5)

    ax1 = fig1.add_subplot(gs[0,1:4])

    rxns="acs pck pps fbp icl mls pfkA pykF".split(' ')
    glucose = np.array(len(rxns)*[0.5])
    acetate = np.array([ 1, 1, 1, 1, 1, 1, 0, 0])
    ls=['Acetate (Exp)','Glucose (Exp)']
    mx = np.array([acetate, glucose])
    normmx = mx/mx.sum(axis=0)[np.newaxis,:]
    #normmx[normmx==0.] = 1./3
    #normmx[normmx==1.] = 2./3
    normmx[normmx==1./3] = 0.
    normmx[normmx==2./3] = 1.
    rxn0 = pandas.DataFrame(normmx, index=ls, columns=rxns)
    ax1 = sb.heatmap(rxn0, cmap="Greens", linewidths=.5, cbar=False, square=True)
    ax1.set_title('Genes expressed during full growth on glucose or acetate')
    ax1.set_yticklabels(ls, rotation=0)
    ttl1 = ax1.title
    ttl1.set_position([.5, 1.3])

    ax3 = fig1.add_axes([.85, .3, .03, .4])
    ax2 = fig1.add_subplot(gs[1,1:4])
    ax2 = sb.heatmap(rxn, ax=ax2, cmap="Greens", linewidths=.5, cbar=True, square=True, cbar_ax=ax3, vmin=0, vmax=1)
    ax2.set_title(tit)
    plt.xticks(rotation=20)
    ttl2 = ax2.title
    ttl2.set_position([.5, 1.3])

    ls2 = [r+' (pFBA)' for r in rxn.index]
    
    ax2.set_yticklabels(ls2, rotation=0)

    plt.savefig('%s/gene_exp_exponential_growth.png' % (pt))
    plt.savefig('../../outputs/pqfigures/fig2a.png', format='png', dpi=1500)

    return

def geneExprGE(pt):
    #rxns="ACKr PPCK PPS FBP ICL MALS PFK PYK PPC ICDHyr".split(' ')
    #rxns="ACKr PPCK PPS FBP ICL MALS PFK PYK PPC ICDHyr".split(' ')
    mx = np.array([[1, 1, 0, 0, 0, 0, -1, -1, -1, -1]])
    rxns="acs pck pps fbp icl mls pfkA pykF ppc icd".split(' ')
    ls=['']
    rxn = pandas.DataFrame(mx, index=ls, columns=rxns)
    ax = sb.heatmap(rxn, cmap="RdBu_r", center=0, linewidths=.5, cbar=False, square=True)#, cbar=True, cbar_kws={"shrink": .25}, square=True)
    ax.set_title('Gene expression after glucose exhaustion')
    plt.yticks(rotation=0)
    plt.savefig('%s/geneexp_comp_%s_%s.png' % (pt, '-'.join(ls), '-'.join(rxns)))
    plt.clf()
    return

def geneExprFullGrowth(pt):
    #rxns="ACKr PPCK PPS FBP ICL MALS PFK PYK PPC ICDHyr".split(' ')
    rxns="acs pck pps fbp icl mls pfkA pykF".split(' ')
    glucose = np.array(len(rxns)*[0.5])
    acetate = np.array([ 1, 1, 1, 1, 1, 1, 0, 0])
    ls=['Acetate','Glucose']
    mx = np.array([acetate, glucose])
    normmx = mx/mx.sum(axis=0)[np.newaxis,:]
    #normmx[normmx==0.] = 1./3
    #normmx[normmx==1.] = 2./3
    normmx[normmx==1./3] = 0.
    normmx[normmx==2./3] = 1.
    rxn = pandas.DataFrame(normmx, index=ls, columns=rxns)
    ax = sb.heatmap(rxn, cmap="Greens", linewidths=.5, cbar=False, square=True)
    ax.set_title('Genes expressed during full growth on glucose or acetate')
    plt.yticks(rotation=0)
    plt.savefig('%s/geneexp_comp_%s_%s.png' % (pt, '-'.join(ls), '-'.join(rxns)))
    return

def mergeRxnVs(fnames, ndfs, r, t):
    '''
    '''
    rdf = pandas.DataFrame(t, columns=['time'])
    for n, d in zip(fnames, ndfs):
        rdf[n] = d[r]
    return rdf

def alignRxnVs(fnames, dfs, r, t):
    '''
    '''
    rdf = pandas.DataFrame(t, columns=['time'])
    for n, d in zip(fnames, dfs):
        print('Aligning ', r, ' of simulation ', ' '.join(n))
        v = d[r]
        tsim = np.array(d['time'])
        nv = np.array( [v[(np.abs(tsim - tarr)).argmin()] if tarr < tsim[-1] else np.nan for tarr in t ] )
        cname = '-'.join(n[1:4])
        rdf[cname] = nv
    #return fnames, dfs, r, t
    return rdf

def alignDf(df, tx):
    '''
    '''
    ndf = pandas.DataFrame(tx, columns=['time'])
    tsim = np.array(df['time'])
    idx = [(np.abs(tsim - tarr)).argmin() if tarr < tsim[-1] else np.nan for tarr in tx ]
    ndf = df.loc[idx].rename(columns={'time':'OTime'}).reset_index()
    ndf['time'] = tx
    return ndf

def loadModelAndDF(fname):
    dmodel = cPickle.load(open(fname, 'r'))
    if dmodel.FBAsolutions[-1] is None:
        del dmodel.FBAsolutions[-1]
        del dmodel.T[-1]
    #dd = tempFix(dmodel.FBAsolutions)
    df = pandas.DataFrame(dmodel.FBAsolutions)
    #df = pandas.DataFrame(dd)
    df['time'] = dmodel.T
    df = df.set_index('time')
    df = df.reset_index()

    return dmodel, df

def tempFix(fsol):
    nfsol = []
    for i in range(len(fsol)):
        if type(fsol[i]) == dict:
            nfsol.append(pandas.Series(fsol[i]))
        else:
            nfsol.append(fsol[i])
    return nfsol

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

def singleSimulationAnalyses(args, infile):
    '''
    One simulation results
    '''

    dmodel, df = loadModelAndDF(infile)
    cbm = dmodel.cbModel

    ifn, aname = findModelNames(infile)
    
    ylims = eval(args.ylims)
    ymin = ylims[0]
    ymax = ylims[1]
    
    taxlab = 'Time [hr]'
    vaxlab = 'V [mmol/(gDW*hr)]'

    rxnsid = args.rxns.split(' ')
    if len(rxnsid[0]) < 1:
        rxnsid = []
    #print(rxnsid)

    if args.plots:
        if len(rxnsid) > 0:
            plotFromDF(df, 'time', rxnsid, 'Fluxes over time', taxlab, vaxlab, '%s/%s-v-%s.png' % (args.pathout, aname, '-'.join(rxnsid)))
            for r in rxnsid:
                plotFromDF(df, 'time', r, 'Fluxes over time', taxlab, vaxlab, '%s/%s-v-%s.png' % (args.pathout, aname, r))
        

    if args.tstepdist:
        tdist = pandas.DataFrame( df['time'] )
        tdist['Dt'] = tdist.diff(1)
        ax = tdist['Dt'].plot.hist(bins=100)
        ax.set_title('Time step distribution')
        ax.set_xlabel('Delta T [hr]')
        ax.set_yscale('log')
        plt.savefig('%s/%s-delta_t_hist.png' % (args.pathout, aname))
        ax = tdist.plot(subplots=True, title='Time evolution at simulation steps')
        #plt.suptitle()
        ax[1].set_xlabel('Simulation step')
        ax[0].set_ylabel(taxlab)
        ax[1].set_ylabel(taxlab)
        plt.savefig('%s/%s-t_dist.png' % (args.pathout, aname))
    
    maxrxnsid = []
    if args.maxv:
        del df['time']
        maxvdf = df.idxmax(axis=1)
        maxrxnsid = list(maxvdf.value_counts().index)
        print(maxrxnsid)
        mdf = df[maxrxnsid]
        mdf['time'] = dmodel.T
        plotFromDF(mdf, 'time', None, 'Maximal fluxes over time', taxlab, vaxlab, '%s/%s-max_v.png' % (args.pathout, aname))
        mdf['max'] = maxvdf
        dfs = []
        for i in maxrxnsid:
            dfs.append(mdf[mdf['max']==i][['time',i]])
        fig1 = plt.figure()
        f1 = fig1.add_subplot(111)
        f1.set_title('Maximal fluxes over time')
        plt.xlabel(taxlab)
        plt.ylabel(vaxlab)
        #f1.plot(dmodel.T)
        for i,d in zip(maxrxnsid, dfs):
            f1.plot(d['time'], d[i], label=i)
        l = f1.legend(loc='best', prop={'size':10})
        fig1.savefig('%s/%s-max_v_cut.png' % (args.pathout, aname))
        df['time'] = dmodel.T

    if args.deltas:
        # Compute the euclidean distance between the flux vectors at T and T+1
        # dV_t(j) = sqrt( Sum_i (v_i(t(j+t)) - v_i(t(t)))^2 )
        ddf = prepareDistDF(df, dmodel.T)#, cols=rxnsid, eucliddist=True)
        plotFromDF(ddf, 'step', 'deltaV', 'Euclidean distance of flux vector solutions', 'Simulation step', r'$\Delta \vec{V}$ [mmol/(gDW*hr)]',
                   '%s/%s-delta_v.png' % (args.pathout, aname), ymin=ymin, ymax=ymax, logy=True)

        plotFromDF(ddf, 'time-step', 'deltaV', 'Euclidean distance of flux vector solutions', taxlab, r'$\Delta \vec{V}$ [mmol/(gDW*hr)]',
                   '%s/%s-delta_v_2x.png' % (args.pathout, aname), ymin=ymin, ymax=ymax, logy=True)

        ddf = prepareDistDF(df, dmodel.T, eucliddist=False, sumcols=False)
        del ddf['time']
        del ddf['step']
        maxdvdf = ddf.idxmax(axis=1)
        maxdvid = list(maxdvdf.value_counts().index)
        print('Run next:   -P -r "%s"' % ' '.join(maxdvid))
        mddf = ddf[maxdvid]
        mddf['time'] = dmodel.T
        mddf['step'] = mddf.index
        plotFromDF(mddf, 'time-step', None, 'Maximal change in fluxes over time', taxlab, r'$|v_{t_{i+1}} - v_{t_{i}}|$ [mmol/(gDW*hr)]', '%s/%s-max_deltav.png' % (args.pathout, aname), ymin=ymin, ymax=ymax, logy=True)

        if len(rxnsid) > 0:
            ddf = prepareDistDF(df, dmodel.T, cols=rxnsid, eucliddist=False, sumcols=False)
            del ddf['time']
            plotFromDF(ddf, 'step', None, 'Euclidean distance of flux vector solutions', 'Simulation step', r'$\Delta \vec{V}$ [mmol/(gDW*hr)]',
                       '%s/%s-delta_v-%s.png' % (args.pathout, aname, '-'.join(rxnsid)), ymin=ymin, ymax=ymax, logy=True)

    if args.exchange:
        # Plot EX_rxn flux vs Time
        # Distinguish among EX_rxns == 0; >= 0; <= 0; changing
        exrxs_nul = []
        exrxs_mix = []
        exrxs_neg = []
        exrxs_pos = []
        exrxns = []
        metSeedID = []
        metDict = {}
        for r in cbm.reactions:
            if 'EX_' in r.id:
                metSeedID.append(cbm.metabolites.get_by_id(r.id[3:]).annotation.get('seed.compound'))
                metDict[r.id] = metSeedID[-1]
                exrxns.append(r)
                x = np.array(df[r.id])
                if min(x) >= 0:
                    if max(x) > 0:
                        exrxs_pos.append(r.id)
                    else:
                        exrxs_nul.append(r.id)
                else:
                    if max(x) > 0:
                        exrxs_mix.append(r.id)
                    else:
                        exrxs_neg.append(r.id)
                df.plot('time', r.id)
                plt.savefig('%s/%s-v_%s.png' % (args.pathout, aname, r.id))
        if len(exrxs_nul) > 0:
            plotFromDF(df, 'time', exrxs_nul, 'Off fluxes', taxlab, vaxlab, '%s/%s-v_allEXnul.png' % (args.pathout, aname))
        if len(exrxs_pos) > 0:
            plotFromDF(df, 'time', exrxs_pos, 'Out fluxes', taxlab, vaxlab, '%s/%s-v_allEXpos.png' % (args.pathout, aname))
        if len(exrxs_neg) > 0:
            plotFromDF(df, 'time', exrxs_neg, 'In fluxes', taxlab, vaxlab, '%s/%s-v_allEXneg.png' % (args.pathout, aname))
        if len(exrxs_mix) > 0:
            plotFromDF(df, 'time', exrxs_mix, 'Exchange fluxes', taxlab, vaxlab, '%s/%s-v_allEXmix.png' % (args.pathout, aname))

    if args.biomass:
        growthdf = pandas.read_csv('../../ecoli/richard_data/ecoli_spentmedia_glc.csv')
        plotBiomass(dmodel, (growthdf['Time'], growthdf['glucose']), '%s/%s-biomass_%s.png' % (args.pathout, aname, '_'.join(ifn[1:])), 'E. coli %s growth on glucose' % ifn[1])

    if args.fva:
        dmodel.setCbModelBounds(tindex = len(dmodel.T)/2)
        #solution = parsimonious.pfba(cbm, solver=m.solver)
        if args.exchange:
            md = pandas.DataFrame.from_dict(metDict, orient='index')
            md.columns = ['seedID']
            dt = flux_variability_analysis(cbm, exrxns)
            dt['secretable'] = (dt['maximum'])>0
            dt = pandas.concat([dt, md], axis=1)
            #dt['secretable'] = np.abs(dt['maximum'])>0
            #print(dt)
            #with open('%s/%s-fva.json' % (args.pathout, aname), 'w') as fp:
            #    json.dump(dt, fp)
            dt.index.name = 'rxn'
            dt[['seedID','secretable']].to_csv('%s/%s-v_fva.csv' % (args.pathout, aname))
            #glmap = json.load(open('map-%s-%s.json' % (self.name, '-'.join(rxnsid)), 'r'))

    return dmodel, df


def prepareDistDF(df, tcol, cols=[], eucliddist=True, sumcols=True):
    # Compute the euclidean distance between the flux vectors at T and T+1
    # dV_t(j) = sqrt( Sum_i (v_i(t(j+t)) - v_i(t(t)))^2 )
    if 'time' in df.columns:
        del df['time']
    if len(cols) > 0:
        df = df[cols]
    if eucliddist:
        ddf = np.square(df.diff(1))
        if sumcols:
            ddf = pandas.DataFrame( np.sqrt( ddf.sum(axis=1) ) )
    else:
        ddf = np.abs(df.diff(1))
        if sumcols:
            ddf = pandas.DataFrame( ddf.sum(axis=1) )
    if sumcols:
        ddf.columns = ['deltaV']
    ddf['time'] = tcol
    ddf = ddf.set_index('time')
    ddf = ddf.reset_index()
    ddf['step'] = ddf.index
    return ddf
    

def plotFromMultiDF(dfs, dfnames, xname, ynames, tlab, xlab, ylab, figname, pstyle='line', ymin=None, ymax=None):
    # 'time', rxnsid, 'Enjalbert 2013 GE', 'Time [hr]', 'V [mmol/(gDW*hr)]', '%s/v-%s_%s.png' % (args.pathout, '-'.join(rxnsid), '-'.join(cnames)))
    '''
    '''
    if len(dfs) < 2:
        print('Only one DF!')
        return
    ax = dfs[0].plot(xname, ynames, title=tlab, kind=pstyle)
    for i in range(1,len(dfs)):
        dfs[i].plot(xname, ynames, title=tlab, kind=pstyle, ax=ax)
    plt.savefig(figname)
    return


def f(s):
    if abs(s).max() < 0.00001:
        return s
    return s/abs(s).max()

def plotFromTwoDF(dfs, dfnames, xname, ynames, tlab, xlab, ylab, figname, pstyle='line', ymin=None, ymax=None, normdfs=True):
    # 'time', rxnsid, 'Enjalbert 2013 GE', 'Time [hr]', 'V [mmol/(gDW*hr)]', '%s/v-%s_%s.png' % (args.pathout, '-'.join(rxnsid), '-'.join(cnames)))
    '''
    '''
    if len(dfs) != 2:
        print('This works only of 2 DFs!')
        return
    if normdfs:
        dfs[0] = dfs[0].apply(f, axis=0)
        dfs[1] = dfs[1].apply(f, axis=0)

    ax = dfs[0].plot(xname, ynames, title=tlab, kind=pstyle)
    ax2 = ax.twinx()
    ax.plot(0, 0, 'k-', label=dfnames[0])
    ax.plot(0, 0, 'k--', label=dfnames[1])
    dfs[1].plot(xname, ynames, title=tlab, kind=pstyle, ax=ax2, ls='dashed', legend=False)
    ax.legend(loc='upper left')
    y1 = ax.get_ylim()
    y2 = ax2.get_ylim()
    ylims = (min(y1[0],y2[0]) if ymin is None else ymin, max(y1[1], y2[1])  if ymax is None else ymax)
    ax.set_ylim(ylims)
    ax2.set_ylim(ylims)
    ax2.get_yaxis().set_visible(False)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    #l = ax.get_legend_handles_labels()
    plt.savefig(figname)
    return

def plotFromDF(df, xname, ynames, tlab, xlab, ylab, figname, pstyle='line', ymin=None, ymax=None, logy=False):
    doublex = False
    if xname == 'time-step':
        doublex = True
        xname = 'time'
        x2name = 'step'
        x2vals = np.array(df[x2name])
        del df[x2name]
        x2lab = 'Simulation step'
        
    df.plot(xname, ynames, title=tlab, kind=pstyle) #'time', exrxs_pos, title='Out fluxes')
    plt.xlabel(xlab) #'Time [hr]')
    plt.ylabel(ylab)#'V [mmol/(gDW*hr)]')
    if ymax is not None:
        plt.ylim(ymax=ymax)
    if ymin is not None:
        plt.ylim(ymin=ymin)

    if doublex:
        ax = plt.gca()
        ax2 = ax.twiny()
        axt = ax.get_xticks()
        ax2t = [(np.abs(df[xname]-t)).argmin() for t in axt]
        ax2.set_xticks(axt)
        ax2.set_xticklabels(['{:g}'.format(t) for t in ax2t], y=0.92)
        ax2.set_xlim(ax.get_xlim())
        ax2.xaxis.set_label_coords(0.5,0.9)
        ax.legend(loc='center left')
        ax2.set_xlabel(x2lab, color='r')
        ax2.tick_params(direction='in', colors='r')
        if logy:
            ax.set_yscale('log')
    else:
        if logy:
           plt. yscale('log')

    plt.savefig(figname) #'%s/v_allEXpos.png' % (args.pathout))
    return
    


def plotSingleReaction():
    
    return

def plotBiomass(model, expdata, figname, title):
    
    ybm1 = np.array(model.dmetabolites['biomass'].quantity)
    ybmf = np.array(model.dreactions['growth'].ub)
    x = np.array(model.T)
    fig1 = plt.figure()
    f1 = fig1.add_subplot(111)
    f1.set_title(title)
    plt.xlabel('Time [hr]')
    plt.ylabel('OD600 / gDW')
    f1.plot(x, ybm1, 'b-', label='Biomass (sim)')
    #f1.plot(x, ybmf, 'y.', label='Growth (sim)')
    f1.plot(expdata[0], expdata[1], 'k.', label='Biomass (exp)')
    l = f1.legend(loc='best', prop={'size':10})
    fig1.savefig(figname)

    return

def plotFit(params, expdata, figname):
    
    t = np.linspace(0, 30, 100)
    n = logisticGrowth(t, *params)
    fig1 = plt.figure()
    f1 = fig1.add_subplot(111)
    f1.plot(expdata[0], expdata[1], 'k.', label='Biomass')
    f1.plot(t, n, 'r-', label='Logistic fit')
    #f1.plot(x, ybmf, 'y.', label='Growth (sim)')
    fig1.savefig(figname)

    return

def computeLagPhase(model, t1):
    '''
    Enjalbert2013
    (t1 - tm) = tlag = (t1 - t0) - ln(X1/X0)/mumax
    tm = t0 + ln(X1/X0)/mumax
    '''
    mumax = max(model.dreactions['growth'].flux)
    #tlag = t1 - 
    

def cobratest():
    model = cobra.io.read_sbml_model('../../ecoli/bigg/e_coli_core.xml')
    print(len(model.reactions))
    print(len(model.metabolites))
    print(len(model.genes))
    # for r in model.reactions:
    #     print(r.id)
    #     print('bounds: ', r.lower_bound, r.upper_bound)
    #     print('formula: ', r.reaction)
    #     print(r.check_mass_balance())
    return model

def solveStatPhaseFBA(x, y):
    '''
    The organism will enter a stationary phase (parsimonious FBA, switching off multiple objectives)
    After the run a random algorithm will decide if the organism will exit the stationary phase
    The probability of exiting should increase with the lenght of the stationary phase
    '''
    print('\n %d %d' % (x,y))
    s = x+y
    if s < 5:
        tryAgain = random.randint(0,4)
        if tryAgain:
            print 'Retry '
            x = solveStatPhaseFBA(x, y+1)
            print('try ', tryAgain, s)
            return x
        else:
            print 'Exit '
            return 0
    return 10


def fbaSol(m):
    if type(m) == cmo.DynamicModel:
        cm = m.cbModel
    else:
        cm = m
    solution = parsimonious.pfba(cm)
    sdf = pandas.DataFrame(solution.x_dict)
    sdf.columns = ['flux']
    sdf['rxn'] = sdf.index
    print('FBA solution')
    print(cm.summary())
    
    print('FVA=1. solution')
    print(cm.summary(fva=1.))
    
    return

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
