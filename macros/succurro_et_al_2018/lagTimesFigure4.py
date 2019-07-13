#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2018/01/08     **
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
import os

from scipy.optimize import curve_fit
from scipy.integrate import odeint

## WE USE GROWTH RATES MEASURED! THEY "INCLUDE" death already
deathrate=0.0

def logisticGrowth(t, n0, r, k):
    n = k/(1 + (k-n0)*np.exp(-r*t)/n0)
    return n

# http://bigg.ucsd.edu/models/e_coli_core
# ecoli = '../../ecoli/bigg/e_coli_core.xml'
# iJO1366 Escherichia coli
# http://bigg.ucsd.edu/models/iJO1366
# ecoli = '../../ecoli/bigg/iJO1366.xml'
# iJR904 Escherichia coli
# http://bigg.ucsd.edu/models/iJR904
# ecoli = '../../ecoli/bigg/iJR904.xml'

batchcond = {'batch_high_Ac' : 'M9A',
             'fedbatch_high_Ac': 'M9GA',
             'batch_mixed_Ac': 'M9GA',
             'batch_low_Glc' : 'M9G'}

def main():
    '''
    Test growth matrix data on different media
    '''
    args = options()
    verbose = args.verbose

    if args.memory:
        memo = 0.2
    else:
        memo = 0.0
    if args.storage:
        sto = 0.2
    else:
        sto = 0.0
    if args.debug:
        if args.scan:
            if args.usepapermu:
                mumaxGl = 0.57
                mumaxAc = 0.23
            elif args.usesamemu:
                mumaxGl = 0.7
                mumaxAc = 0.3
            else:
                mumaxGl = 0.549
                mumaxAc = 0.124
            scanLag(args.pathout, mumaxGl, mumaxAc, memo, sto, args.extreme)
        else:
            runDebug(args)
        return
    
    iflist = args.infile.split(' ')

    fn, an = findModelNames(iflist[0])
    ibc = batchcond[fn[2]]

    tarr = getInitialConditions(args, iflist[0], ibc)
    if args.getime:
        return
    if len(iflist) < 2 and not args.odesys:
        return

    flab = '_muSimu'
    if args.usepapermu:
        mumaxGl = 0.57
        mumaxAc = 0.23
        flab = '_muExp'
    elif args.usesamemu:
        mumaxGl = 0.7
        mumaxAc = 0.5
        flab = '_muHyp'
    else:
        mumaxGl = 0.549
        mumaxAc = 0.124

    print('*********************')
    print('Mu Max Gl: ', mumaxGl)
    print('Mu Max Ac: ', mumaxAc)
    print('*********************')
    
    tlag = {'M9A' : np.zeros(len(tarr)), 
            'M9G' : np.zeros(len(tarr))}
    tlag_up = {'M9A' : np.zeros(len(tarr)), 
               'M9G' : np.zeros(len(tarr))}
    tlag_do = {'M9A' : np.zeros(len(tarr)), 
               'M9G' : np.zeros(len(tarr))}
    tlag_err = {'M9A' : np.zeros([len(tarr), 2]),
                'M9G' : np.zeros([len(tarr), 2])}

    npyfname = {}

    bm0 = float(args.bminit)
    
    for k in tlag.keys():
        npyfname[k] = '%s/tlag_%s-%s%s.npy' % (args.pathout, ibc, k, flab)
        if not args.tlagcomp and os.path.isfile(npyfname[k]):
            tlag[k] = np.load(npyfname[k])

    if args.plots:
        ifl = iflist[0]
        fn, an = findModelNames(ifl)
        print(fn[-1])
        bc = batchcond[fn[2]]
        print('%s-%s' % (ibc, bc))
        dmodel = cPickle.load(open(ifl, 'r'))
        mumax = mumaxAc if bc == 'M9A' else mumaxGl
        plotBiomass(args.pathout, mumax, dmodel.dmodels['ECgl'], dmodel.dmodels['ECac'], 10., fn[-1], ibc, logy=True)
        
    #for i, t in enumerate(tarr):
    for ifl in iflist[1:]:
        fn, an = findModelNames(ifl)
        print(fn[-1])
        tidx = int( fn[-1].split('_')[0][7:] )
        bc = batchcond[fn[2]]
        print('%s-%s' % (ibc, bc))
        dmodel = cPickle.load(open(ifl, 'r'))
        mumax = mumaxAc if bc == 'M9A' else mumaxGl
        if args.plots:
            plotBiomass(args.pathout+'/'+ibc, mumax, dmodel.dmodels['ECgl'], dmodel.dmodels['ECac'], float(args.tval), fn[-1], '%s-%s' % (ibc, bc))
        print(fn, mumax)
        if args.tlagcomp:
            print('Compute Lag Time')
            tlag[bc][tidx] = getLagTime(mumax, dmodel.dmodels['ECgl'], dmodel.dmodels['ECac'], t1val=float(args.tval))

    if args.odesys:
        tmax = 10.
        vpsi, kpsi, vphi, kphi, eps = float(args.vmaxpsi), float(args.kmtranspsi), float(args.vmaxphi), float(args.kmtransphi), 0.9
        #vpsi, kpsi, vphi, kphi, eps = 1., 4.0, 1., 12.5, 0.9
        if args.scan:
            xGl = np.linspace(0.01, 0.99, 99)
            xAc = 1-xGl
        else:
            xGl = np.load('%s/xGl_%s.npy' % (args.pathout, ibc))
            xAc = np.load('%s/xAc_%s.npy' % (args.pathout, ibc))
        for bc in ['M9G', 'M9A']:
            print 'Daughter: ', bc
            X0 = [0., 0., 0., 0.] #[ECac0, ECgl0, gl0, ac0]
            mumax = 0
            if bc == 'M9G':
                muG = mumaxGl
                muA = 0.
                X0[2] = 15. #mM
                mumax = mumaxGl
                if ibc == 'M9GA' and args.memory:
                    print 'Adding a memory factor!'
                    muA = memo*mumaxAc
            else:
                muG = 0.
                muA = mumaxAc
                X0[3] = 40. #mM
                mumax = mumaxAc
                if ibc == 'M9GA' and args.memory:
                    print 'Adding a memory factor!'
                    muG = memo*mumaxGl
                if args.storage:
                    if muG < 0.001:
                        print 'Adding a storage effect'
                        muG = sto*mumaxGl
            for tidx in range(len(xGl)):
                X0[0] = bm0*xAc[tidx]
                X0[1] = bm0*xGl[tidx]
                #print(X0, tidx, bm0)
                ts, Xs = odesys(tmax, muG, muA, X0, vpsi, kpsi, vphi, kphi, eps, phi0=float(args.phioffset), psi0=float(args.psioffset))
                tlag[bc][tidx] = getLagTimeODE(mumax, ts, Xs, t1val=float(args.tval))#, inMinutes=True):
                print 'Computed Lag Time from ODE sys solution: ', tlag[bc][tidx]
                if args.tlagerrs:
                    X0[1] = bm0*min(0.999, xGl[tidx]*1.15)
                    X0[0] = bm0*(1.0 - X0[1])
                    ts, Xs = odesys(tmax, muG, muA, X0, vpsi, kpsi, vphi, kphi, eps, phi0=float(args.phioffset), psi0=float(args.psioffset))
                    tlag_up[bc][tidx] = getLagTimeODE(mumax, ts, Xs, t1val=float(args.tval))
                    print 'T starting with 15% more ECgl (', min(0.999, xGl[tidx]*1.15),'): ', tlag_up[bc][tidx]
                    X0[1] = bm0*min(0.999, xGl[tidx]*0.85)
                    X0[0] = bm0*(1.0- X0[1])
                    ts, Xs = odesys(tmax, muG, muA, X0, vpsi, kpsi, vphi, kphi, eps, phi0=float(args.phioffset), psi0=float(args.psioffset))
                    tlag_do[bc][tidx] = getLagTimeODE(mumax, ts, Xs, t1val=float(args.tval))
                    print 'T starting with 15% less ECgl (', min(0.999, xGl[tidx]*0.85),'): ', tlag_do[bc][tidx]
                    print [np.abs(tlag[bc][tidx] - tlag_do[bc][tidx]), np.abs(tlag[bc][tidx] - tlag_up[bc][tidx])]
                    tlag_err[bc][tidx] = np.array([np.abs(tlag[bc][tidx] - tlag_do[bc][tidx]), np.abs(tlag[bc][tidx] - tlag_up[bc][tidx])])
                #getLagTime(mumax, dmodel.dmodels['ECgl'], dmodel.dmodels['ECac'], t1val=float(args.tval))
    
    if args.tlagcomp:
        for k in tlag.keys():
            np.save(npyfname[k], tlag[k])
            if args.tlagerrs:
                np.save(npyfname[k].replace('tlag_', 'tlag_up_'), tlag_up[k])
                np.save(npyfname[k].replace('tlag_', 'tlag_do_'), tlag_do[k])

    if args.scan:
        plotLagTimesScan(args, tarr, tlag, ibc, flab)
        return

    #plotLagTimes(args, tarr, tlag, ibc, flab, tlag_err)
    plotLagTimes(args, tarr, tlag, ibc, flab, [tlag_do, tlag_up])
    
    return

def runDebug(args):
    #xGl = np.array([0.999, 0.5, 0.001])
    xGl = np.linspace(0.5, 0.999, 10)
    xAc = 1-xGl
    tlag = {'M9A' : np.zeros(len(xGl)), 
            'M9G' : np.zeros(len(xGl))}
    mumaxGl = 0.57
    mumaxAc = 0.23
    flab = '_debug'
    tmax = 10.
    vpsi, kpsi, vphi, kphi, eps = float(args.vmaxpsi), float(args.kmtranspsi), float(args.vmaxphi), float(args.kmtransphi), 0.9
    phi0=float(args.phioffset)
    psi0=float(args.psioffset)
    #vpsi, kpsi, vphi, kphi, eps = 1., 4.0, 1., 12.5, 0.9
    for bc in ['M9G', 'M9A']:
        X0 = [0., 0., 0., 0.] #[ECac0, ECgl0, gl0, ac0]
        mumax = 0
        if bc == 'M9G':
            muG = mumaxGl
            muA = 0.
            X0[2] = 15. #mM
            mumax = mumaxGl
        else:
            muG = 0.
            muA = mumaxAc
            X0[3] = 40. #mM
            mumax = mumaxAc
        for tidx in range(len(xGl)):
            X0[0] = xAc[tidx]
            X0[1] = xGl[tidx]
            print(X0, tidx)
            ts, Xs = odesys(tmax, muG, muA, X0, vpsi, kpsi, vphi, kphi, eps, phi0=phi0, psi0=psi0)
            print('Compute Lag Time from ODE sys solution')
            tlag[bc][tidx] = getLagTimeODE(mumax, ts, Xs, t1val=float(args.tval))#, inMinutes=True):
    plotLagTimes(args, xGl, tlag, 'Debug', flab, None, plotdata=False, convertToMin=False, titlestr='vpsi: %.2f; kpsi: %.2f; vphi: %.2f; kphi: %.2f; phi0: %.2f; psi0: %.2f'%(vpsi, kpsi, vphi, kphi, phi0, psi0))
    

def plotLagTimes(args, tarr, tlag, ibc, flab, tlag_err, plotdata=True, convertToMin=True, secondaxis=False, titlestr='Title'):

    cvf = 1.
    xlab = 'Sampling Time [hours]'
    if args.debug:
        xlab = 'initial BM ECgl / BM tot'
    ylab = 'Lag times [hours]'
    fig1 = plt.figure()
    if args.subplots:
        f1 = fig1.add_subplot(211)
        f2 = fig1.add_subplot(212)
    else:
        f1 = fig1.add_subplot(111)
    f1.set_title('Lag time after re-inoculation')#; GE @ %.3f hr' % tarr[4])
    if titlestr!='Title':
        f1.set_title(titlestr)
    if convertToMin:
        print('GE = ', tarr[4])
        tarr = 60*(tarr-tarr[4])
        cvf = 60.
        xlab = 'Time relative to Glc exhaustion [minutes]'
        ylab = 'Lag time [minutes]'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if args.tlagerrs:
        #yerrA = cvf*tlag_err['M9A'].transpose()
        #yerrG = cvf*tlag_err['M9G'].transpose()
        #f1.errorbar(tarr, cvf*tlag['M9A'], yerr=yerrA, marker='s', color='#5d6d7e', label=ibc+'-'+'M9A, sim')
        #f1.errorbar(tarr, cvf*tlag['M9G'], yerr=yerrG, marker='d', color='#808b96', label=ibc+'-'+'M9G, sim')
        f1.plot(tarr, cvf*tlag['M9A'], 's', color='#5d6d7e', label=ibc+'-'+'M9A, sim')
        f1.plot(tarr, cvf*tlag['M9G'], 'd', color='#808b96', label=ibc+'-'+'M9G, sim')
        tldo_m9a = tlag_err[0]['M9A']
        tlup_m9a = tlag_err[1]['M9A']
        tldo_m9g = tlag_err[0]['M9G']
        tlup_m9g = tlag_err[1]['M9G']
        for i,t in enumerate(tarr):
            f1.plot([t, t], [cvf*tldo_m9a[i], cvf*tlup_m9a[i]], '-', color='#5d6d7e')
            f1.plot([t, t], [cvf*tldo_m9g[i], cvf*tlup_m9g[i]], '-', color='#808b96')
    else:
        f1.plot(tarr, cvf*tlag['M9A'], 's', color='#5d6d7e', label=ibc+'-'+'M9A, sim')
        f1.plot(tarr, cvf*tlag['M9G'], 'd', color='#808b96', label=ibc+'-'+'M9G, sim')
    fign = 'fig4'+ibc
    if not args.debug:
        f1.axvline(0, linestyle=':')
    f1.axhline(0, linestyle=':')
    if plotdata:
        if ibc == 'M9G':
            #data_GA=pandas.read_csv('../../ecoli/enjalbert2015_data/fig4a_ga.csv', sep=',',header=None)
            #data_GG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig4a_gg.csv', sep=',',header=None)
            data_GA=pandas.read_csv('../../ecoli/enjalbert2015_data/fig4a_fromEnjalbert2015_ga.csv', sep=',')
            data_GG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig4a_fromEnjalbert2015_gg.csv', sep=',')
            fign='fig4a'
        else:
            #data_GA=pandas.read_csv('../../ecoli/enjalbert2015_data/fig4b_ga.csv', sep=',',header=None)
            #data_GG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig4b_gg.csv', sep=',',header=None)
            data_GA=pandas.read_csv('../../ecoli/enjalbert2015_data/fig4b_fromEnjalbert2015_gaa.csv', sep=',')
            data_GG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig4b_fromEnjalbert2015_gag.csv', sep=',')
            fign='fig4b'
        if secondaxis:
            ax = f1.twinx()
            #ax.plot(data_GA[0], data_GA[1], 'P', color='#e74c3c', label=ibc+'-'+'M9A, exp')
            #ax.plot(data_GG[0], data_GG[1], 'X', color='#f39c12', label=ibc+'-'+'M9G, exp')
            ax.errorbar(data_GA['Time'], data_GA['lag'], yerr=data_GA['SD'], fmt='P', color='#e74c3c', label=ibc+'-'+'M9A, exp')
            ax.errorbar(data_GG['Time'], data_GG['lag'], yerr=data_GG['SD'], fmt='X', color='#f39c12', label=ibc+'-'+'M9G, exp')
            handles0, labels0 = f1.get_legend_handles_labels()
            handles1, labels1 = ax.get_legend_handles_labels()
            ax.tick_params(axis='y', colors='#a93226')
            ax.spines['right'].set_color('#a93226')
            #f1.plot(0, -100, 'rP'
            #f1.plot(0, -100, 'rX'
            #print(handles0, labels0)
            l = f1.legend(handles0 + handles1, labels0+labels1)#, loc='best', prop={'size':10})
        else:
            f1.errorbar(data_GA['Time'], data_GA['lag'], yerr=data_GA['SD'], fmt='P', color='#e74c3c', label=ibc+'-'+'M9A, exp')
            f1.errorbar(data_GG['Time'], data_GG['lag'], yerr=data_GG['SD'], fmt='X', color='#f39c12', label=ibc+'-'+'M9G, exp')
            f1.set_ylim(-2., 100)
            #f1.plot(data_GA[0], data_GA[1], 'P', color='#e74c3c', label=ibc+'-'+'M9A, exp')
            #f1.plot(data_GG[0], data_GG[1], 'X', color='#f39c12', label=ibc+'-'+'M9G, exp')
            l = f1.legend(loc='best', prop={'size':10})
    else:
        l = f1.legend(loc='best', prop={'size':10})
    if args.subplots:
        xGl = np.load('%s/xGl_%s.npy' % (args.pathout, ibc))
        xAc = np.load('%s/xAc_%s.npy' % (args.pathout, ibc))
        f2.plot(tarr, xGl, label=ibc+' ECgl ratio')
        f2.plot(tarr, xGl+xAc, label=ibc+' ECgl+ECac')
        l2 = f2.legend(loc='best', prop={'size':10})

    
    if args.tlagerrs:
        flab = flab+'_werr'
    if args.memory:
        flab = flab+'_memo'
    if args.storage:
        flab = flab+'_store'

    fig1.savefig('%s/%s%s.png' % (args.pathout, fign, flab))
    #fig1.savefig('../../outputs/pqfigures/%s%s.png' % (fign, flab), format='png', dpi=1500)

    return

def plotBiomass(pathout, mumax, modGl, modAc, t1val, fname, expcon, logy=False):
    
    t1ind = (np.abs(np.array(modGl.T)-t1val)).argmin()
    t1 = modGl.T[t1ind]
    print('*******',t1val,t1)
    tarr = np.array(modGl.T[:t1ind+1])
    bmGl = np.array(modGl.dmetabolites['biomass_ECgl'].quantity[:t1ind+1])
    bmAc = np.array(modAc.dmetabolites['biomass_ECac'].quantity[:t1ind+1])
    bm0 = bmGl[0]+bmAc[0]

    fig1 = plt.figure()
    gs = gridspec.GridSpec(5, 1)
    # Top left
    ax1 = fig1.add_subplot(gs[0:2,0])
    ax2 = fig1.add_subplot(gs[2:4,0], sharex=ax1)
    ax3 = fig1.add_subplot(gs[4,0], sharex=ax1)
    tit = 'Growth for %s; #mu = %.3f' % (expcon, mumax)
    ax1.set_title(tit)
    ax1.plot(tarr, bm0*np.exp(mumax*tarr), '#e74c3c', label='Exp growth')
    ax1.plot(tarr, bmGl,'#16a085', label='BM ECgl')
    ax1.plot(tarr, bmAc,'#2980b9', label='BM ECac')
    ax1.plot(tarr, bmGl+bmAc, 'k', label='BM tot')
    if logy:
        ax1.set_yscale('log')
    ax1.legend(loc='best', prop={'size':10})
    ax1.set_ylabel('Q [gDW]')
    ax2.set_ylabel('C [mM]')
    ax3.set_xlabel('Time [hours]')
    ax2.plot(tarr, np.array(modGl.dmetabolites['ex_glucose'].concentration[:t1ind+1]), '#c0392b', label='Glucose')
    ax2.plot(tarr, np.array(modGl.dmetabolites['ex_acetate'].concentration[:t1ind+1]), '#f39c12', label='Acetate')
    ax2.legend(loc='best', prop={'size':10})
    ax3.set_ylabel('BM ECgl / BM tot')
    ax3.plot(tarr, bmGl/(bmAc + bmGl), 'k')
    ax3.set_ylim(0.0, 1.0)
    fig1.savefig('%s/biomass_%s_%s.png' % (pathout, expcon, fname))

    return

def getLagTime(mumax, modGl, modAc, t1val, inMinutes=False):
    t1ind = (np.abs(np.array(modGl.T)-t1val)).argmin()
    t1 = modGl.T[t1ind]
    x1 = modGl.dmetabolites['biomass_ECgl'].quantity[t1ind] + modAc.dmetabolites['biomass_ECac'].quantity[t1ind]
    x0 = modGl.dmetabolites['biomass_ECgl'].quantity[0]  + modAc.dmetabolites['biomass_ECac'].quantity[0]
    tm = np.log(x1/x0)/mumax
    tlag = t1 - tm
    if inMinutes:
        tlag = tlag*60.
    print(t1, x1, x0, tm, tlag)
    return tlag

def findGlExhTime(m0):
    i0 = np.array(m0.dmetabolites['ex_glucose'].quantity).argmin()
    i1 = m0.dreactions['growth'].flux.index(0)
    if i0 != i1:
        print('ACHTUNG!! ',i0,i1)
    return m0.T[i1]

def getInitialConditions(args, infile, ibc):
    dmodel = cPickle.load(open(infile, 'r'))
    mGl = dmodel.dmodels['ECgl']
    mAc = dmodel.dmodels['ECac']

    if args.scan:
        return np.linspace(0.01, 0.99, 99)

    bmGl = np.array(mGl.dmetabolites['biomass_ECgl'].quantity)
    bmAc = np.array(mAc.dmetabolites['biomass_ECac'].quantity)
    # Time in hours, from 4h -60 minutes to +90 minutes
    #t = numpy.linspace(-1, 1.5, 11)
    t = np.linspace(-1, 1.5, 11)
    tGE = findGlExhTime(mGl)
    t = t+tGE
    tidx = [(np.abs(np.array(mGl.T)-ti)).argmin() for ti in t]
    
    xGl = bmGl[tidx]/(bmGl[tidx]+bmAc[tidx])
    xAc = bmAc[tidx]/(bmGl[tidx]+bmAc[tidx])

    np.save(('%s/xGl_%s.npy' % (args.pathout, ibc)), xGl)
    np.save(('%s/xAc_%s.npy' % (args.pathout, ibc)), xAc)
    #print xGl+xAc
    #print xGl

    vpsi, kpsi, vphi, kphi = float(args.vmaxpsi), float(args.kmtranspsi), float(args.vmaxphi), float(args.kmtransphi)
    
    if args.getime:
        print(' %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f ' % (vpsi, kpsi, vphi, kphi, tGE, mGl.T[tidx[4]], xGl[4], xAc[4] ))
    
    for i, x in enumerate(xGl):
        #trans="--phitransition --psitransition --kmtransphi 12.5 --kmtranspsi 4 --vmaxpsi 0.2 --vmaxphi 0.2"
        #trans="--phitransition --psitransition --kmtransphi 4 --kmtranspsi 4 --vmaxpsi 0.2 --vmaxphi 0.2"
        #trans="--phitransition --psitransition --kmtransphi 1 --kmtranspsi 4 --vmaxpsi 0.2 --vmaxphi 0.2"
        trans="--phitransition --psitransition --kmtransphi %.3f --kmtranspsi %.3f --vmaxpsi %.3f --vmaxphi %.3f --phioffset %.3f --psioffset %.3f"%(kphi, kpsi,vpsi,vphi,float(args.phioffset), float(args.psioffset))
        baselinecomm = "python daughterCultures.py -p %s --run --runconsortium -b 0.0025 --ratioecgl '%.3f' -t 5 -x '-11.5' %s -e '0.9' -P -l sampleT%d" % (args.pathout, x, trans, i)
        if args.shellscript:
            # 15mM Glc
            print('%s_M9G --runglucose' % baselinecomm)
            # 45mM Ac
            print('%s_M9A --runhighacetate' % baselinecomm)
    
    return t

def loadModelAndDF(fname):
    dmodel = cPickle.load(open(fname, 'r'))
    if dmodel.FBAsolutions[-1] is None:
        del dmodel.FBAsolutions[-1]
        del dmodel.T[-1]
    df = pandas.DataFrame(dmodel.FBAsolutions)
    df['time'] = dmodel.T
    df = df.set_index('time')
    df = df.reset_index()
    return dmodel, df

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

def plotFromDF(df, xname, ynames, tlab, xlab, ylab, figname, pstyle='line', ymin=None, ymax=None):
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
    plt.savefig(figname) #'%s/v_allEXpos.png' % (args.pathout))
    return

def computeLagPhase(model, t1):
    '''
    Enjalbert2013
    (t1 - tm) = tlag = (t1 - t0) - ln(X1/X0)/mumax
    tm = t0 + ln(X1/X0)/mumax
    '''
    mumax = max(model.dreactions['growth'].flux)
    #tlag = t1 - 
    return

def plotLagTimesScan(args, tarr, tlag, ibc, flab, plotdata=True, convertToMin=True, secondaxis=False):

    cvf = 1.
    xlab = 'Initial ECgl ratio'
    ylab = 'Lag times [hours]'
    fig1 = plt.figure()
    f1 = fig1.add_subplot(111)
    f1.set_title('Lag time after re-inoculation')#; GE @ %.3f hr' % tarr[4])
    if convertToMin:
        print('GE = ', tarr[4])
        cvf = 60.
        ylab = 'Lag time [minutes]'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    f1.plot(tarr, cvf*tlag['M9A'], label='M9A, sim')
    f1.plot(tarr, cvf*tlag['M9G'], label='M9G, sim')
    l = f1.legend(loc='best', prop={'size':10})
    fig1.savefig('%s/scanIC%s.png' % (args.pathout, flab))
    
def odesys(tmax, muG, muA, X0, vpsi, kpsi, vphi, kphi, eps, phi0=0.0, psi0=0.0, h=5):
    '''
    float tmax = end time 
    float muG = growth rate on glc
    float muA = growth rate on ac
    list X0 = [bm0ECac, bm0ECgl, gl0, ac0]
    '''

    def dX_dt(X, t):
        return [X[0]*(muA - (phi0 + vphi*(X[2]**h)/(X[2]**h + kphi**h)) - deathrate) + X[1]*eps*(psi0 + vpsi*(X[3]**h)/(X[3]**h + kpsi**h)),
                X[1]*(muG - (psi0 + vpsi*(X[3]**h)/(X[3]**h + kpsi**h)) - deathrate) + X[0]*eps*(phi0 + vphi*(X[2]**h)/(X[2]**h + kphi**h)),
                0,
                0]

    ts = np.linspace(0, tmax, 100)
    Xs = odeint(dX_dt, X0, ts)
    
    ecac = Xs[:,0]
    ecgl = Xs[:,1]
    
    fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    plt.plot(ts, ecac, "+", label="EC ac")
    plt.plot(ts, ecgl, "x", label="EC gl")
    plt.xlabel("Time")
    plt.ylabel("Population")
    plt.legend()
    fig.savefig('/tmp/odesys_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f.png' % (muG, muA, X0[0], X0[1], X0[2], X0[3], vpsi, kpsi, vphi, kphi, eps))
    
    return ts, Xs
    
def tableLag(muG, muA, X0, psi, phi):
    '''
    float tmax = end time 
    float muG = growth rate on glc
    float muA = growth rate on ac
    list X0 = [bm0ECac, bm0ECgl]
    '''

    def dX_dt(X, t):
        a = muA - deathrate - phi
        b = 0.9*psi
        c = muG - deathrate - psi
        d = 0.9*phi
        return [X[0]*a + X[1]*b,
                X[1]*c + X[0]*d]

    tmax = 2.
    ts = np.linspace(0, tmax, 100)
    Xs = odeint(dX_dt, X0, ts)
    
    ecac = Xs[:,0]
    ecgl = Xs[:,1]

    t0 = 0
    t1i = (np.abs(ts-1.5)).argmin()
    t1 = ts[t1i]
    x0 = ecac[t0] + ecgl[t0]
    x1 = ecac[t1i] + ecgl[t1i]
    muM = max(muG, muA)
    tlag = t1 - (np.log(x1/x0))/muM
    return tlag*60


def scanLag(pathout, muG=0.683, muA=0.194, memory=0.0, stored=0.0, extreme=True):

    if extreme:
        X0G = [0.1, 0.9]
        xg = 90
        X0A = [0.6, 0.4]
        xa = 40
    else:
        X0G = [0.2, 0.8]
        xg = 80
        X0A = [0.5, 0.5]
        xa = 50

    xpsi = np.linspace(0., 1., 50)
    yphi = np.linspace(0., 1., 50)
    Xpsi = np.tile(np.array([xpsi]).transpose(), (1,len(xpsi)))
    Yphi = np.tile(np.array(yphi), (len(yphi),1))

    ZGG = np.zeros([len(xpsi), len(yphi)])
    ZAG = np.zeros([len(xpsi), len(yphi)])
    ZGA = np.zeros([len(xpsi), len(yphi)])
    ZAA = np.zeros([len(xpsi), len(yphi)])

    memorystr = ''
    if memory > 0:
        memorystr += 'withMemo_'
    if stored > 0:
        memorystr += 'withStored_'

    # M9G
    for i,psi in enumerate(xpsi):
        for j,phi in enumerate(yphi):
            ZGG[i,j] = tableLag(muG, 0.0, X0G, psi, phi)
            ZGA[i,j] = tableLag(muG*stored, muA, X0G, psi, phi)
            ZAG[i,j] = tableLag(muG, muA*memory, X0A, psi, phi)
            ZAA[i,j] = tableLag(max(muG*memory, muG*stored), muA, X0A, psi, phi)

    fig1 = plt.figure()
    h1 = plt.contourf(Xpsi, Yphi, ZGG, cmap=plt.cm.Greens)
    plt.colorbar()
    plt.title(('Lag time in M9G approximation, %d%% ECgl' % xg))
    plt.xlabel('Psi')
    plt.ylabel('Phi')
    fig1.savefig('%s/lagMap_%sM9G_%dECgl.png' % (pathout, memorystr, xg))

    fig2 = plt.figure()
    h2 = plt.contourf(Xpsi, Yphi, ZAG, cmap=plt.cm.Blues)
    plt.colorbar()
    plt.title(('Lag time in M9G approximation, %d%% ECgl' % xa))
    plt.xlabel('Psi')
    plt.ylabel('Phi')
    fig2.savefig('%s/lagMap_%sM9G_%dECgl.png' % (pathout, memorystr, xa))

    fig3 = plt.figure()
    h3 = plt.contourf(Xpsi, Yphi, ZGA, cmap=plt.cm.Reds)
    plt.colorbar()
    plt.title('Lag time in M9A approximation, %d%% ECgl' % xg)
    plt.xlabel('Psi')
    plt.ylabel('Phi')
    fig3.savefig('%s/lagMap_%sM9A_%dECgl.png' % (pathout, memorystr, xg))

    fig4 = plt.figure()
    h4 = plt.contourf(Xpsi, Yphi, ZAA, cmap=plt.cm.Purples)
    plt.colorbar()
    plt.title('Lag time in M9A approximation, %d%% ECgl' % xa)
    plt.xlabel('Psi')
    plt.ylabel('Phi')
    fig4.savefig('%s/lagMap_%sM9A_%dECgl.png' % (pathout, memorystr, xa))

    # fig1 = plt.figure()
    # f1 = fig1.add_subplot(221)
    # h1 = f1.contourf(Xpsi, Yphi, ZGG, cmap=plt.cm.Greens)
    # plt.colorbar(h1,cax=f1)
    # f1.set_title('Lag time in M9G approximation, 90% ECgl')
    # f1.set_xlabel('Psi')
    # f1.set_ylabel('Phi')
    # f2 = fig1.add_subplot(222)
    # h2 = plt.contourf(Xpsi, Yphi, ZAG, cmap=plt.cm.Blues)
    # plt.colorbar(h2,cax=f2)
    # f2.set_title('Lag time in M9G approximation, 40% ECgl')
    # f2.set_xlabel('Psi')
    # f2.set_ylabel('Phi')
    # f3 = fig1.add_subplot(223)
    # h3 = f3.contourf(Xpsi, Yphi, ZGA, cmap=plt.cm.Greens)
    # plt.colorbar(h3,cax=f3)
    # f3.set_title('Lag time in M9A approximation, 90% ECgl')
    # f3.set_xlabel('Psi')
    # f3.set_ylabel('Phi')
    # f4 = fig1.add_subplot(224)
    # h4 = plt.contourf(Xpsi, Yphi, ZAA, cmap=plt.cm.Blues)
    # plt.colorbar(h4,cax=f4)
    # f4.set_title('Lag time in M9A approximation, 40% ECgl')
    # f4.set_xlabel('Psi')
    # f4.set_ylabel('Phi')
    # fig1.savefig('%s/lagMap.png' % (pathout))



def scanLagB(pathout):
    muG = 0.683
    muA = 0.194
    X0G = [0.1, 0.9]
    X0A = [0.6, 0.4]

    xpsi = np.linspace(0., 1., 10)
    yphi = np.linspace(0., 1., 10)
    Xpsi = np.tile(np.array([xpsi]).transpose(), (1,len(xpsi)))
    Yphi = np.tile(np.array(yphi), (len(yphi),1))

    ZG = np.zeros([len(xpsi), len(yphi)])
    ZA = np.zeros([len(xpsi), len(yphi)])

    # M9G
    for i,psi in enumerate(xpsi):
        for j,phi in enumerate(yphi):
            ZG[i,j] = tableLag(muG, 0., X0G, psi, phi)
            ZA[i,j] = tableLag(0., muA, X0A, psi, phi)

    fig1 = plt.figure()
    h = plt.contourf(Xpsi, Yphi, ZG, cmap=plt.cm.Greens)
    plt.colorbar()
    plt.title('Lag time in M9G approximation')
    plt.xlabel('Psi')
    plt.ylabel('Phi')
    fig1.savefig('%s/lagMap_M9G.png' % (pathout))

    fig2 = plt.figure()
    h = plt.contourf(Xpsi, Yphi, ZA, cmap=plt.cm.Blues)
    plt.colorbar()
    plt.title('Lag time in M9A approximation')
    plt.xlabel('Psi')
    plt.ylabel('Phi')
    fig2.savefig('%s/lagMap_M9A.png' % (pathout))

def getLagTimeODE(mumax, ts, Xs, t1val, inMinutes=False):
    t1ind = (np.abs(ts-t1val)).argmin()
    t1 = ts[t1ind]
    x1 = Xs[t1ind, 0] + Xs[t1ind, 1]
    x0 = Xs[0, 0] +Xs[0, 1]
    tm = np.log(x1/x0)/mumax
    #print(t1, x1, x0, tm)
    tlag = t1 - tm
    if inMinutes:
        tlag = tlag*60.
    return tlag
    
def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-D', '--debug', help='debug', action='store_true')
    parser.add_argument('-E', '--tlagerrs', help='add errors', action='store_true')
    parser.add_argument('-G', '--getime', help='check ratios at glucose exhaustion point', action='store_true')
    parser.add_argument('-K', '--scan', help='scan initial conditions', action='store_true')
    parser.add_argument('-M', '--memory', help='add a memory effect', action='store_true')
    parser.add_argument('-N', '--storage', help='add a storage effect', action='store_true')
    parser.add_argument('-O', '--odesys', help='run simple ode system', action='store_true')
    parser.add_argument('-P', '--plots', help='plot V', action='store_true')
    parser.add_argument('-R', '--subplots', help='plot bm ratios', action='store_true')
    parser.add_argument('-S', '--shellscript', help='printout shell script', action='store_true')
    parser.add_argument('-T', '--tlagcomp', help='compute t lag', action='store_true')
    parser.add_argument('-U', '--usepapermu', help='use mu max from paper', action='store_true')
    parser.add_argument('-X', '--extreme', help='use 90% and 40% for lag scan', action='store_true')
    parser.add_argument('-Z', '--usesamemu', help='use same mu max', action='store_true')
    parser.add_argument('-i', '--infile', help='file path to load', default='../../outputs/twoEC_pfba_enjalbert2015_fig6C/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p')
    parser.add_argument('-p', '--pathout', help='path for outputs', default='/tmp/')
    parser.add_argument('-t', '--tval', help='value of t1', default='1.5')
    parser.add_argument('-y', '--ylims', help='tple of y ranges', default='(None, None)')
    parser.add_argument('-b', '--bminit', help='initial biomass', default='1.')
    parser.add_argument('-j', '--kmtranspsi', help='Km psi transition', default='0.')
    parser.add_argument('-k', '--kmtransphi', help='Km phi transition', default='0.')
    parser.add_argument('-l', '--vmaxpsi', help='V max psi transition', default='0.0')
    parser.add_argument('-m', '--vmaxphi', help='V max phi transition', default='0.0')
    parser.add_argument('-w', '--psioffset', help='offset psi transition', default='0.0')
    parser.add_argument('-u', '--phioffset', help='offset phi transition', default='0.0')
    args = parser.parse_args()
    if args.verbose:
        print "verbosity turned on"
        print args
    return args

if __name__=="__main__":
    main()
