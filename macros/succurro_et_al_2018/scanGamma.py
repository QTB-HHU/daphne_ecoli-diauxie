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

def findOkGammas(pathout, muG, muA, noaconly=False, noglonly=False):
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
    if noaconly:
        G = Ga*Gc
    elif noglonly:
        G = Gb*Gc
    else:
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

    muG = 0
    muA = 0

    muGref = 0.57
    muAref = 0.23

    if args.okgammas:
        x = np.linspace(0., 1., 50)
        
        G =  findOkGammas(args.pathout, muGref, muAref, noglonly=True)

        
        phi0 = float(args.minphi)
        print 'Minimal phi = ',phi0,' (no Glucose)'
        psis = x[G[:,np.argmin(np.abs(x-phi0))]]
        print 'Psi & (b) Acetate only & (c) Mixed \\\\'
        for i in psis:
            s = '%.3f & %.3f & %.3f \\\\' % (i, np.log(computeGammaNeg(0., muAref, i, phi0, 0.9)), np.log(computeGammaNeg(muGref, muAref, i, phi0, 0.9)) )
            print s
            
        psi0 = float(args.maxpsi)
        print 'Maximal psi = ',psi0,' (high Acetate)'
        phis = x[G[np.argmin(np.abs(x-psi0)),:]]
        print 'Phi & (b) Acetate only & (c) Mixed \\\\'
        for i in phis:
            s = '%.3f & %.3f & %.3f \\\\' % (i, np.log(computeGammaNeg(0., muAref, psi0, i, 0.9)), np.log(computeGammaNeg(muGref, muAref, psi0, i, 0.9)) )
            print s

        G =  findOkGammas(args.pathout, muGref, muAref, noaconly=True)

        psi0 = float(args.minpsi)
        print 'Minimal psi = ',psi0,' (no Acetate)'
        phis = x[G[np.argmin(np.abs(x-psi0)),:]]
        print 'Phi & (a) Glucose only & (c) Mixed \\\\'
        for i in phis:
            s = '%.3f & %.3f & %.3f \\\\' % (i, np.log(computeGammaNeg(muGref, 0., psi0, i, 0.9)), np.log(computeGammaNeg(muGref, muAref, psi0, i, 0.9)) )
            print s

        phi0 = float(args.maxphi)
        print 'Maximal phi = ',phi0,' (high Glucose)'
        phi0 = 0.2
        psis = x[G[:,np.argmin(np.abs(x-phi0))]]
        print 'Psi & (a) Glucose only & (c) Mixed \\\\'
        for i in psis:
            s = '%.3f & %.3f & %.3f \\\\' % (i, np.log(computeGammaNeg(muGref, 0., i, phi0, 0.9)), np.log(computeGammaNeg(muGref, muAref, i, phi0, 0.9)) )
            print s

    return

def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-K', '--okgammas', help='make gamma plots', action='store_true')
    parser.add_argument('-p', '--pathout', help='path for outputs', default='/tmp')
    parser.add_argument('-a', '--minphi', help='minimal phi', default='0.04')
    parser.add_argument('-b', '--maxphi', help='maximal phi', default='0.2')
    parser.add_argument('-c', '--minpsi', help='minimal psi', default='0.04')
    parser.add_argument('-d', '--maxpsi', help='maximal psi', default='0.2')
    args = parser.parse_args()
    if args.verbose:
        print "verbosity turned on"
        print args
    return args

if __name__=="__main__":
    main()
