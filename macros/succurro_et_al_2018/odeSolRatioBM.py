import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.pyplot import cm 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from plotly import __version__
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import Scatter, Figure, Layout, Surface
from matplotlib import gridspec
import csv

def computeLambda(muG, muA, psi, phi, eps):
    '''
    
    '''
    a = 1.
    b = psi + phi - muG - muA
    c = muG*muA - muA*psi - muG*phi + (1 - eps**2)*psi*phi
    delta = b**2 - 4*a*c
    if delta < 0:
        return nan, nan
    lambda1 = (np.sqrt(delta) - b)/(2*a)
    lambda2 = (-1*np.sqrt(delta) - b)/(2*a)
    return lambda1, lambda2

def computeGammaPos(muG, muA, psi, phi, eps):
    '''
    
    '''
    a = -1*eps*psi
    b = muG - muA - psi + phi
    c = eps*phi
    delta = b**2 - 4*a*c
    #if delta < 0:
    #    return nan, nan
    gamma = (np.sqrt(delta) - b)/(2*a)
    #gamma[gamma < 0] = np.nan
    return gamma

def computeGammaNeg(muG, muA, psi, phi, eps):
    '''
    
    '''
    a = -1*eps*psi
    b = muG - muA - psi + phi
    c = eps*phi
    delta = b**2 - 4*a*c
    #if delta < 0:
    #    return nan, nan
    gamma = (-1*np.sqrt(delta) - b)/(2*a)
    #gamma[gamma < 0] = np.nan
    return gamma

def findOkGammas(pathout, muG, muA):
    eps = 0.9
    x = np.linspace(0., 1., 50)
    y = np.linspace(0., 1., 50)
    X = np.tile(np.array([x]).transpose(), (1,len(x)))
    Y = np.tile(np.array(y), (len(y),1))
    Za = computeGammaNeg(muG, 0.0, X, Y, eps)
    Zb = computeGammaNeg(0.0, muA, X, Y, eps)
    Zc = computeGammaNeg(muG, muA, X, Y, eps)
    Ga = Za > 1
    Gb = Zb < 1
    Gc = Zc > 1
    G = Ga*Gb*Gc

    fig1 = plt.figure()
    h = plt.contourf(X, Y, G, cmap=plt.cm.Greens)
    plt.colorbar()
    plt.title('Gamma conditions')
    plt.xlabel('Psi')
    plt.ylabel('Phi')
    fig1.savefig('%s/booleanMap.png' % (pathout))
    

    fig2 = plt.figure()
    h = plt.contourf(X, Y, np.log(Za), cmap=plt.cm.Reds)
    plt.colorbar()
    plt.title('log(Gamma_glc)')
    plt.xlabel('Psi')
    plt.ylabel('Phi')
    fig2.savefig('%s/logGammaGlc.png' % (pathout))

    fig3 = plt.figure()
    h = plt.contourf(X, Y, np.log(Zb), cmap=plt.cm.Blues)
    plt.colorbar()
    plt.title('log(Gamma_ac)')
    plt.xlabel('Psi')
    plt.ylabel('Phi')
    fig3.savefig('%s/logGammaAc.png' % (pathout))

    fig4 = plt.figure()
    h = plt.contourf(X, Y, np.log(Zc), cmap=plt.cm.Purples)
    plt.colorbar()
    plt.title('log(Gamma_mix)')
    plt.xlabel('Psi')
    plt.ylabel('Phi')
    fig4.savefig('%s/logGammaMix.png' % (pathout))

    
    fig5 = plt.figure()
    h = plt.contourf(X, Y, np.log(Za)*np.log(Zb), cmap=plt.cm.Greens)
    plt.colorbar()
    plt.title('log(Gamma_glc) * log(Gamma_ac)')
    plt.xlabel('Psi')
    plt.ylabel('Phi')
    fig5.savefig('%s/logGammaProductAcGlc.png' % (pathout))

    return G
    

def odesys(tmax, muG, muA, psi, phi, eps, bm0, r0):
    '''

    '''

    def dX_dt(X, t):
        return [X[0]*(muA - phi) + X[1]*eps*psi,
                X[1]*(muG - psi) + X[0]*eps*phi ]

    ts = np.linspace(0, tmax, 100)
    X0 = [bm0*(1-r0), bm0*r0]
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
    fig.savefig('/tmp/odesys_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f.png' % (muG, muA, psi, phi, eps, bm0, r0))
    
    return Xs

def main():

    args = options()
    verbose = args.verbose


    if args.scan:
        plotLagTimesScan(args)
        return
    
    muG = 0
    muA = 0

    muGref = 0.57
    muAref = 0.23

    if args.okgammas:
        G =  findOkGammas(args.pathout, muGref, muAref)
        x = np.linspace(0., 1., 50)
        phi0 = 0.2
        psis = x[G[:,np.argmin(np.abs(x-phi0))]]
        print 'Psi & (a) Glucose only & (b) Acetate only & (c) Mixed \\\\'
        for i in psis:
            s = '%.3f & %.3f & %.3f & %.3f \\\\' % (i, np.log(computeGammaNeg(muGref, 0., i, 0.2, 0.9)), np.log(computeGammaNeg(0., muAref, i, 0.2, 0.9)), np.log(computeGammaNeg(muGref, muAref, i, 0.2, 0.9)) )
            print s
        return

    if args.runhighacetate:
        muA = muAref
        title = 'acetate'
    elif args.runglucose:
        muG = muGref
        title = 'glucose'
    elif args.runmixedacetate:
        muA = muAref
        muG = muGref
        title = 'mixed'        

    biomass0 = float(args.biomassi)
    pcECgl = float(args.ratioecgl)
    thrPc = 0.99999
    if pcECgl > thrPc:
        pcECgl = thrPc
    if pcECgl < (1 - thrPc):
        pcECgl = 1 - thrPc
    strratioecgl = args.ratioecgl.replace('.','p')

    tmax = 12
    if args.steady:
        tmax = 100
    phi = 0.2
    psi = 0.2
    eps = float(args.efftrans)
    gammas = {'Xphi': [],
              'Xpsi': [],
              'Xeps': [],
              'Xbmr': []}
    if args.odesys:
        odesys(tmax, muG, muA, psi, phi, eps, biomass0, pcECgl)
        fig = plt.figure()
        gs = gridspec.GridSpec(4, 5)
        # Top left
        #fig.tight_layout()
        ax1 = fig.add_subplot(gs[0:2,0:2])
        ax2 = fig.add_subplot(gs[2:4,0:2])#, sharex=ax1)
        ax3 = fig.add_subplot(gs[0:2,2:4])#, sharex=ax1)
        ax4 = fig.add_subplot(gs[2:4,2:4])#, sharex=ax1)
        #ax5 = fig.add_subplot(gs[0,4])#, sharex=ax1)
        #eca="EC ac"
        #ecg="EC gl"
        eca="X_A"
        ecg="X_G"
        ax1.set_ylabel(eca+' cells')
        ax2.set_ylabel(eca+' cells')
        ax3.set_ylabel(eca+' cells')
        ax2.set_xlabel(ecg+' cells')
        ax4.set_xlabel(ecg+' cells')
        ax4.set_xlabel(ecg+' cells')
        ax1.set_title('Vary phi')
        ax2.set_title('Vary psi')
        ax3.set_title('Vary eps')
        ax4.set_title('Vary bm0 ratio')
        ax3.text(1.2, 0.8, 'Condition:\n%s;\nValues:' % title, transform=ax3.transAxes, fontsize=12, bbox=dict(facecolor='yellow', alpha=0.5))
        vpar = np.linspace(0.0, 1.0, 11)
        for p in vpar:
            Xphi = odesys(tmax, muG, muA, psi, p, eps, biomass0, pcECgl)
            Xpsi = odesys(tmax, muG, muA, p, phi, eps, biomass0, pcECgl)
            Xeps = odesys(tmax, muG, muA, psi, phi, p, biomass0, pcECgl)
            pcP = p
            if p < 0.001:
                pcP = 0.001
            if p > 0.999:
                pcP = 0.999
            Xbmr = odesys(tmax, muG, muA, psi, phi, eps, biomass0, pcP)
            ax1.plot(Xphi[:,1], Xphi[:,0], '-')
            ax2.plot(Xpsi[:,1], Xpsi[:,0], '-')
            ax3.plot(Xeps[:,1], Xeps[:,0], '-', label='%.2f' % p)
            ax4.plot(Xbmr[:,1], Xbmr[:,0], '-')
            ax4.plot(Xbmr[0,1], Xbmr[0,0], 'k*')
            ax4.set_xscale('log')
            ax4.set_yscale('log')
            ax3.legend(loc='upper left', bbox_to_anchor=(1.15, 0.6)) #, transform=ax3.transAxes)
            gs.tight_layout(fig)
            #fig.suptitle(title)
            gammas['Xphi'].append(Xphi[-1,1]/Xphi[-1,0])
            gammas['Xpsi'].append(Xpsi[-1,1]/Xpsi[-1,0])
            gammas['Xeps'].append(Xeps[-1,1]/Xeps[-1,0])
            gammas['Xbmr'].append(Xbmr[-1,1]/Xbmr[-1,0])
            fig.savefig('/tmp/vary_odesys_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f.png' % (muG, muA, psi, phi, eps, biomass0, pcECgl))
        ax1.plot(Xphi[0,1], Xphi[0,0], 'k*')
        ax3.plot(Xeps[0,1], Xeps[0,0], 'k*')
        ax2.plot(Xpsi[0,1], Xpsi[0,0], 'k*')
        fig.savefig('%s/vary_odesys_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f.png' % (args.pathout, muG, muA, psi, phi, eps, biomass0, pcECgl))
        with open('%s/final_gamma_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f.csv' % (args.pathout, muG, muA, psi, phi, eps, biomass0, pcECgl), 'w') as f:
            #w = csv.DictWriter(f, gammas.keys())
            writer = csv.writer(f)
            writer.writerow(gammas.keys())
            writer.writerows(zip(*gammas.values()))
            
        return
        
    #plotGamma(muG, muA, float(args.efftrans), args.pathout, title, args.logscale)
    l1, l2 = computeLambda(muG, muA, 0.3, 0.14, float(args.efftrans))
    print(l1, l2)
    
    return

def plotGamma(muG, muA, eps, pathout, title, doLog=False):

    x = np.linspace(0., 1., 50)
    y = np.linspace(0., 1., 50)

    X = np.tile(np.array([x]).transpose(), (1,len(x)))
    Y = np.tile(np.array(y), (len(y),1))

    Zpos =  computeGammaPos(muG, muA, X, Y, eps)
    Zneg =  computeGammaNeg(muG, muA, X, Y, eps)
    gammaTit = '$\\Gamma$'

    cdist = np.zeros(Zpos.shape)

    if doLog:
        Zpos = np.log(Zpos)
        Zneg = np.log(Zneg)
        gammaTit = '$\\log\\Gamma$'
    
    #trace1 = Surface()
    data = [dict(x=X, y=Y, z=Zpos, opacity=0.9, showscale=False, surfacecolor=cdist,
                 colorscale=[[0.0, 'rgb(0, 0, 180)'], [1.0, 'rgb(0, 0, 180)']],
                 type='surface', name='Positive root'),
            dict(x=X, y=Y, z=Zneg, opacity=0.9, showscale=False, surfacecolor=cdist,
                 colorscale=[[0.0, 'rgb(0, 180, 180)'], [1.0, 'rgb(0, 180, 180)']],
                 type='surface', name='Negative root')]

    data0 = [dict(x=X, y=Y, z=Zpos, opacity=0.9, showscale=False, surfacecolor=cdist,
                  colorscale=[[0.0, 'rgb(0, 0, 180)'], [1.0, 'rgb(0, 0, 180)']],
                  type='surface', name='Positive root')]
    data1 = [dict(x=X, y=Y, z=Zneg, opacity=0.9, showscale=False, surfacecolor=cdist,
                  colorscale=[[0.0, 'rgb(0, 180, 0)'], [1.0, 'rgb(0, 180, 0)']],
                  type='surface', name='Negative root')]
    
    fn = '%s/gamma_%s_eff%d' % (pathout, title, 100*eps)

    layout = Layout(scene = dict(xaxis = dict(title='$\\psi$'),
                                 yaxis = dict(title='$\\phi$'),
                                 zaxis = dict(title=gammaTit),))
    fig = Figure(data=data, layout=layout)
    plot(fig, filename=fn)
    return

def plotGamma2(muG, muA, eps, pathout, title):

    x = np.linspace(0., 1., 50)
    y = np.linspace(0., 1., 50)

    X = np.tile(np.array([x]).transpose(), (1,len(x)))
    Y = np.tile(np.array(y), (len(y),1))

    Zpos =  computeGammaPos(muG, muA, X, Y, eps)
    Zneg =  computeGammaNeg(muG, muA, X, Y, eps)

    fn = '%s/gamma_%s_eff%d' % (pathout, title, 100*eps)

    layout = Layout(scene = dict(xaxis = dict(title='$\psi$'),
                                 yaxis = dict(title='$\phi$'),
                                 zaxis = dict(title='$\Gamma$'),),
                    width=700,
                    margin=dict(r=20, b=10,l=10, t=10))
    fig = Figure(layout=layout)
    plot(fig, [dict(z=Zpos, type='surface'),
               dict(z=Zneg, showscale=False, opacity=0.9, type='surface')],
         filename=fn )
    return

def previous():
    
    x = np.linspace(0., 1., 50)
    y = np.linspace(0., 1., 50)

    X = np.tile(np.array([x]).transpose(), (1,len(x)))
    Y = np.tile(np.array(y), (len(y),1))

    Zpos =  computeGammaPos(muG, muA, X, Y, eps)
    Zneg =  computeGammaNeg(muG, muA, X, Y, eps)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(X, Y, Zpos, rstride=1, cstride=1, cmap=cm.coolwarm)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    
    ax.set_title('Media = %s' % (title))
    ax.set_xlabel('#psi')
    ax.set_ylabel('#phi')
    ax.set_zlabel('#Gamma')
    fig.savefig('%s/gamma_%s_eff%d' % (pathout, title, 100*eps))

    #plt.show()
    return

def doMM1(X, Y):

    Zpos =  1 - (mmkin(X, K_vit, 4))*(mmkin(Y, K_fe, 4))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(X, Y, Zpos, rstride=1, cstride=1, cmap=cm.coolwarm)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    
    ax.set_title('')
    ax.set_xlabel('[Vit]')
    ax.set_ylabel('[Fe]')
    ax.set_zlabel('  f')
    fig.savefig('../../outputs/fiona/test_hills_doca1.png')

    #plt.show()
    return

def doMM2(X, Y):

    Zpos =  (mmkin1(X, K_vit, 4))*(mmkin1(Y, K_fe, 4))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(X, Y, Zpos, rstride=1, cstride=1, cmap=cm.coolwarm)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    
    ax.set_title('')
    ax.set_xlabel('[Vit]')
    ax.set_ylabel('[Fe]')
    ax.set_zlabel('  f')
    fig.savefig('../../outputs/fiona/test_hills_doca2.png')

    #plt.show()
    return

def varH(k):

    hills = range(1,10)
    colorK = cm.rainbow(np.linspace(0,1,len(hills)))
    fig1, f1 = plt.subplots(2, sharex=True)
    f1[1].set_xlabel('Nutrients')
    f1[1].set_ylabel('MM')
    f1[0].set_ylabel('1 - MM')
    x = np.linspace(0., 10., 100)
    for i, h in enumerate(hills):
        f1[0].plot(x, 1-mmkin(x, k, h), '-', color=colorK[i], label=('h = %d'%h))
        f1[1].plot(x, mmkin(x, k, h), '-', color=colorK[i], label=('h = %d'%h))
    l0 = f1[0].legend(loc='best', prop={'size':10})
    l1 = f1[1].legend(loc='best', prop={'size':10})
    fig1.savefig('../../outputs/fiona/test-mm-h-k%d.png' % k)
    return

def varK(h):

    kcoeff = np.linspace(0.1,1,10)
    colorK = cm.rainbow(np.linspace(0,1,len(kcoeff)))
    fig1, f1 = plt.subplots(2, sharex=True)
    f1[1].set_xlabel('Nutrients')
    f1[1].set_ylabel('MM')
    f1[0].set_ylabel('1 - MM')
    x = np.linspace(0., 10., 100)
    for i,k in enumerate(kcoeff):
        f1[0].plot(x, 1-mmkin(x, k, h), '-', color=colorK[i], label=('K = %.2f'%k))
        f1[1].plot(x, mmkin(x, k, h), '-', color=colorK[i], label=('K = %.2f'%k))
    l0 = f1[0].legend(loc='best', prop={'size':10})
    l1 = f1[1].legend(loc='best', prop={'size':10})
    fig1.savefig('../../outputs/fiona/test-mm-k-h%d.png' % h)
    return


def createDict(len1, len2):
    d = {'M9G': {},
         'M9A': {}}
    for x in range(len1):
        lab='P%d' % x
        d['M9G'][lab] = np.zeros(len2)
        d['M9A'][lab] = np.zeros(len2)
    return d

def plotLagTimesScan(args):#, tarr, tlag, ibc, flab, plotdata=True, convertToMin=True, secondaxis=False):

    xGl = np.linspace(0.01, 0.99, 11)
    xAc = 1-xGl

    vpar = np.linspace(0.0, 1.0, 11)
    kpar = np.linspace(0.0, 40.0, 11)

    tmax = 10.
    vpsi, kpsi, vphi, kphi, eps = float(args.vmaxpsi), float(args.kmtranspsi), float(args.vmaxphi), float(args.kmtransphi), 0.9

    Xvpsi, Xkpsi, Xvphi, Xkphi = scanLagTimes(float(args.tval), tmax, vpsi, kpsi, vphi, kphi, eps, h=int(args.hillcoeff), bias=float(args.bias), inMinutes=True)
    
    ylab = 'Lag time [minutes]'
    xlab = 'Initial ECgl ratio'
    fig = plt.figure()
    gs = gridspec.GridSpec(4, 5)
    ax1 = fig.add_subplot(gs[0:2,0:2])
    ax2 = fig.add_subplot(gs[2:4,0:2], sharex=ax1)
    ax3 = fig.add_subplot(gs[0:2,2:4])#, sharey=ax1)
    ax4 = fig.add_subplot(gs[2:4,2:4], sharex=ax2)#, sharey=ax2)
    ax1.set_ylabel(ylab)
    ax2.set_ylabel(ylab)
    ax3.set_ylabel(ylab)
    ax4.set_ylabel(ylab)
    ax2.set_xlabel(xlab)
    ax4.set_xlabel(xlab)
    ax1.set_title('Vary Vphi')
    ax2.set_title('Vary Vpsi')
    ax3.set_title('Vary Kphi')
    ax4.set_title('Vary Kpsi')
    #ax3.text(1.2, 0.8, 'Condition:\n%s;\nValues:' % title, transform=ax3.transAxes, fontsize=12, bbox=dict(facecolor='yellow', alpha=0.5))
    colors = cm.rainbow(np.linspace(0, 1, len(vpar)))

    ax3.plot(np.nan, np.nan, 'kd', label='M9G')
    ax3.plot(np.nan, np.nan, 'k+', label='M9A')
    for i in range(len(vpar)):
        parid='P%d' % i
        ax1.plot(xGl, Xvphi['M9G'][parid], 'd', color=colors[i], markersize=3)
        ax1.plot(xGl, Xvphi['M9A'][parid], '+', color=colors[i], markersize=3)
        ax2.plot(xGl, Xvpsi['M9G'][parid], 'd', color=colors[i], markersize=3)
        ax2.plot(xGl, Xvpsi['M9A'][parid], '+', color=colors[i], markersize=3)
        ax3.plot(xGl, Xkphi['M9G'][parid], 'd', color=colors[i], markersize=3)
        ax3.plot(xGl, Xkphi['M9A'][parid], '+', color=colors[i], markersize=3)
        ax4.plot(xGl, Xkpsi['M9G'][parid], 'd', color=colors[i], markersize=3)
        ax4.plot(xGl, Xkpsi['M9A'][parid], '+', color=colors[i], markersize=3)
        ax3.plot(np.nan, np.nan, 'o', color=colors[i], label='V=%.2f/K=%.2f' % (vpar[i], kpar[i]))
    ax3.legend(loc='upper left', bbox_to_anchor=(1.15, 0.6), prop={'size': 6}) #, transform=ax3.transAxes)
    gs.tight_layout(fig)
    #fig.suptitle('Lag time of daughter cultures')
    #ax1.legend(loc='best')
    fig.savefig('%s/scan_odesyslag_%.3f_%.3f_%.3f_%.3f_%.3f_h%d_b%.2f.png' % (args.pathout, vpsi, kpsi, vphi, kphi, eps, int(args.hillcoeff), float(args.bias)))
    

def scanLagTimes(t1, tmax, vpsi, kpsi, vphi, kphi, eps, h=5, bias=0.0, inMinutes=False):

    xGl = np.linspace(0.01, 0.99, 11)
    xAc = 1-xGl

    vpar = np.linspace(0.0, 1.0, 11)
    kpar = np.linspace(0.0, 40.0, 11)

    tlagXvphi = createDict(len(vpar), len(xGl))
    tlagXvpsi = createDict(len(vpar), len(xGl))
    tlagXkphi = createDict(len(kpar), len(xGl))
    tlagXkpsi = createDict(len(kpar), len(xGl))

    for bc in ['M9G', 'M9A']:
        X0 = [0., 0., 0., 0.] #[ECac0, ECgl0, gl0, ac0]
        mumax = 0
        if bc == 'M9G':
            muG = 0.57
            muA = 0.
            X0[2] = 15. #mM
            mumax = muG
        else:
            muG = 0.
            muA = 0.23
            X0[3] = 40. #mM
            mumax = muA 

        for bm0idx in range(len(xGl)):
            X0[0] = xAc[bm0idx]
            X0[1] = xGl[bm0idx]
            print(X0, bm0idx)
            for i,p in enumerate(vpar):
                parid='P%d' % i
                #V phi
                ts, Xs = odesyslag(tmax, muG, muA, X0, vpsi, kpsi, p, kphi, eps, h, bias)
                tlagXvphi[bc][parid][bm0idx] = getLagTimeODE(mumax, ts, Xs, t1, inMinutes)
                #V psi
                ts, Xs = odesyslag(tmax, muG, muA, X0, p, kpsi, vphi, kphi, eps, h, bias)
                tlagXvpsi[bc][parid][bm0idx] = getLagTimeODE(mumax, ts, Xs, t1, inMinutes)
            for i,p in enumerate(kpar):
                parid='P%d' % i
                ts, Xs = odesyslag(tmax, muG, muA, X0, vpsi, kpsi, vphi, p, eps, h, bias)
                tlagXkphi[bc][parid][bm0idx] = getLagTimeODE(mumax, ts, Xs, t1, inMinutes)
                ts, Xs = odesyslag(tmax, muG, muA, X0, vpsi, p, vphi, kphi, eps, h, bias)
                tlagXkpsi[bc][parid][bm0idx] = getLagTimeODE(mumax, ts, Xs, t1, inMinutes)


    #return {'vphi': tlagXvphi, 'kpsi': tlagXkpsi, 'vpsi': tlagXvpsi, 'kphi': tlagXkphi}
    return tlagXvpsi, tlagXkpsi, tlagXvphi, tlagXkphi

def odesyslag(tmax, muG, muA, X0, vpsi, kpsi, vphi, kphi, eps, h=5, b=0.0):
    '''
    float tmax = end time 
    float muG = growth rate on glc
    float muA = growth rate on ac
    list X0 = [bm0ECac, bm0ECgl, gl0, ac0]
    '''

    def dX_dt(X, t):
        deathrate=0.03
        return [X[0]*(muA - (b + vphi*(X[2]**h)/(X[2]**h + kphi**h)) - deathrate) + X[1]*eps*(b + vpsi*(X[3]**h)/(X[3]**h + kpsi**h)),
                X[1]*(muG - (b + vpsi*(X[3]**h)/(X[3]**h + kpsi**h)) - deathrate) + X[0]*eps*(b + vphi*(X[2]**h)/(X[2]**h + kphi**h)),
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
    parser.add_argument('-S', '--steady', help='up to steady state', action='store_true')
    parser.add_argument('-L', '--logscale', help='plot log gamma', action='store_true')
    parser.add_argument('-O', '--odesys', help='solve ode sys', action='store_true')
    parser.add_argument('-A', '--runhighacetate', help='batch 45mM Acetate conditions', action='store_true')
    parser.add_argument('-D', '--runglucose', help='batch 15mM Glucose conditions (default)', action='store_true')
    parser.add_argument('-G', '--runmixedacetate', help='batch 15mM Glucose 32mM Acetate conditions', action='store_true')
    parser.add_argument('-K', '--okgammas', help='make gamma plots', action='store_true')
    parser.add_argument('-N', '--scan', help='scan transition parameters', action='store_true')
    parser.add_argument('-o', '--ofile', help='output file name', default='toymodel-test.xml')
    parser.add_argument('-l', '--label', help='output file label', default='diauxic_shift')
    parser.add_argument('-p', '--pathout', help='path for outputs', default='/tmp')
    parser.add_argument('-b', '--biomassi', help='initial biomass', default='0.001')
    parser.add_argument('-r', '--ratioecgl', help='ratio of ECgl/ECac: BMECgl0 = X*biomass0', default='1.')
    parser.add_argument('-d', '--deathrate', help='death rate', default='0.')
    parser.add_argument('-e', '--efftrans', help='transition efficiency', default='0.9')
    parser.add_argument('-j', '--kmtranspsi', help='Km psi transition', default='4.')
    parser.add_argument('-k', '--kmtransphi', help='Km phi transition', default='12.5')
    parser.add_argument('-n', '--vmaxpsi', help='V max psi transition', default='0.2')
    parser.add_argument('-m', '--vmaxphi', help='V max phi transition', default='0.2')
    parser.add_argument('-q', '--hillcoeff', help='hill coefficient', default='5')
    parser.add_argument('-z', '--bias', help='bias', default='0.0')
    parser.add_argument('-t', '--tval', help='value of t1', default='1.5')
    args = parser.parse_args()
    if args.verbose:
        print "verbosity turned on"
        print args
    return args

if __name__=="__main__":
    main()
