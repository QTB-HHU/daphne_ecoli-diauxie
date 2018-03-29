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
from matplotlib.pyplot import cm 
import argparse
import plstyles
import random
import cPickle
import json
import os

from scipy.optimize import curve_fit
from scipy.integrate import odeint

deathrate=0.03

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

    if args.debug:
        xG  = 100*np.load('/home/succurro/repositories/gitlab/dfba-ode-framework/outputs/fig6_enjalbert2015_fig4/xGl_M9G.npy')
        xGA = 100*np.load('/home/succurro/repositories/gitlab/dfba-ode-framework/outputs/fig6_enjalbert2015_fig4/xGl_M9GA.npy')
        sG = 'M9G '
        sA = 'M9GA '
        for i in range(len(xG)):
            sG = sG+('& %.1f' % xG[i])
            sA = sA+('& %.1f' % xGA[i])
        print sG
        print sA
        return
    
    mumaxGl = 0.57
    mumaxAc = 0.23
    
    fig1 = plt.figure()

    gs = gridspec.GridSpec(2,3)
    # Top left
    ax1 = fig1.add_subplot(gs[0,0:2])
    ax2 = fig1.add_subplot(gs[1,0:2], sharex=ax1)
    #ax3 = fig1.add_subplot(gs[0:2,2])
    ax1.set_ylabel('% ECgl')
    ax2.set_ylabel('Lag time (min)')
    ax1.set_title('% ECgl value (mother cultures) and lag time (daughter cultures)')
    ax2.set_xlabel('T relative to GE (1/hr)')

    ibc = args.mother
    ax1.plot(np.nan, np.nan, 'kd', label=ibc+'-M9G')
    ax1.plot(np.nan, np.nan, 'k-', label=ibc+'-M9A')

    wgt=''
    if args.gthr:
        wgt='_glcthr'

    t = np.linspace(-1, 1.5, 11)

    muGA = 0.0
    muAA = mumaxAc

    muGG = mumaxGl
    muAG = 0.0

    psi0=0.04
    phi0=0.04
    eps=0.9
    kphi=5.0
    ks = [float(args.kmtranspsi)]
    vs = [0.2, 0.5, 0.8]
    colors = cm.rainbow(np.linspace(0, 1, len(vs)*len(vs)*len(ks)))

    j=0
    for vpsi in vs:
        for vphi in vs:
            for kpsi in ks:
                tlag = {'M9A' : np.zeros(len(t)), 
                        'M9G' : np.zeros(len(t))}
                parstr = 'vpsi%.1f_kpsi%.1f_vphi%.1f_kphi%.1f' % (vpsi, kpsi, vphi, kphi)
                parlab = 'psi (%.1f, %.1f) phi (%.1f, %.1f)' % (vpsi, kpsi, vphi, kphi)
                xGl = np.load('%s/%s%s/xGl_%s.npy' % (args.pathout, parstr, wgt, ibc))
                xAc = np.load('%s/%s%s/xAc_%s.npy' % (args.pathout, parstr, wgt, ibc))
                for i in range(len(xGl)):
                    tlag['M9A'][i] = getlag(muGA, muAA, [xAc[i], xGl[i], 0., 40.], vpsi, kpsi, vphi, kphi, eps, phi0, psi0)
                    tlag['M9G'][i] = getlag(muGG, muAG, [xAc[i], xGl[i], 15., 0.], vpsi, kpsi, vphi, kphi, eps, phi0, psi0)
                ax1.plot(t, 100*xGl, 'x', color=colors[j])
                ax1.plot(np.nan, np.nan, 'o', color=colors[j], label=parlab)
                ax2.plot(t, tlag['M9A'], '-', color=colors[j])
                ax2.plot(t, tlag['M9G'], 'd', color=colors[j])
                j += 1

    ax1.legend(loc='upper left', bbox_to_anchor=(1., 0.6), prop={'size': 8}) #, transform=ax3.transAxes)
    gs.tight_layout(fig1)
    fig1.savefig('%s/kpsi%d/ratios_%s%s.png' % (args.pathout, int(ks[0]), ibc, wgt))

    return

    
def getlag(muG, muA, X0, vpsi, kpsi, vphi, kphi, eps, phi0=0.0, psi0=0.0, h=5):
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

    
def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-D', '--debug', help='debug', action='store_true')
    parser.add_argument('-E', '--tlagerrs', help='add errors', action='store_true')
    parser.add_argument('-G', '--gthr', help='simulation with glc threshold', action='store_true')
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
    parser.add_argument('-m', '--mother', help='mother condition', default='M9G')
    parser.add_argument('-p', '--pathout', help='path for outputs', default='/home/succurro/repositories/gitlab/dfba-ode-framework/outputs/lagscans/debug/')
    parser.add_argument('-a', '--kmtranspsi', help='Km psi transition', default='0.')
    parser.add_argument('-b', '--kmtransphi', help='Km phi transition', default='0.')
    parser.add_argument('-c', '--vmaxpsi', help='V max psi transition', default='0.0')
    parser.add_argument('-d', '--vmaxphi', help='V max phi transition', default='0.0')
    parser.add_argument('-e', '--psioffset', help='offset psi transition', default='0.0')
    parser.add_argument('-f', '--phioffset', help='offset phi transition', default='0.0')
    args = parser.parse_args()
    if args.verbose:
        print "verbosity turned on"
        print args
    return args

if __name__=="__main__":
    main()
