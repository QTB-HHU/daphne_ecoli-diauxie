#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  last modified: 2016/11/04     **
#************************************

import sys
sys.path.append('../../code/python/')
import time
import classModel as cmo
import classReaction as cre
import classMetabolite as cme
import classPlotter as plotter
import classConvertionFactors as ccf
import cobra
import numpy as np
import matplotlib.pyplot as plt
import optparse
import plstyles

def main():
    '''
    Test dynamicModel on a pure ODE system
    '''
    opts, args = options()
    verbose = opts.verbose
    timeunits = opts.timeunits
    outpath = opts.pathout

    volExt = 1
    volUn = 'mL'
    exitname = "lotkavolterra"

    #parameters = {'a': 1., 'b': 0.1, 'c': 1.5, 'd': 0.75}
    parameters = {'a': float(opts.apar),
                  'b': float(opts.bpar),
                  'c': float(opts.cpar),
                  'd': float(opts.dpar)}
                  
    y0 = eval(opts.yicond)
    frabbit = 'self.p["a"]*RABBITS - self.p["b"]*RABBITS*FOXES'
    ffoxes = '-self.p["c"]*FOXES + self.p["d"]*self.p["b"]*RABBITS*FOXES'

    
    #ODE system
    r_biomass = cme.Biomass([y0[0]], {'r_growth': [(1, 'r_biomass')], 'r_death': [(-1, 'r_biomass', 'f_biomass')]})
    f_biomass = cme.Biomass([y0[1]], {'f_death': [(-1, 'f_biomass')], 'f_growth': [(1, 'r_biomass', 'f_biomass')]})
    
    r_growth = cre.DReaction('r_growth', cre.FixedBound(0., parameters['a']), isODE=True)
    r_death = cre.DReaction('r_death', cre.FixedBound(0., parameters['b']), isODE=True)
    f_growth = cre.DReaction('f_growth', cre.FixedBound(0., (parameters['d']*parameters['b'])), isODE=True)
    f_death = cre.DReaction('f_death', cre.FixedBound(0., parameters['c']), isODE=True)

    dyRxn = {'r_growth': r_growth,
             'r_death':  r_death,
             'f_growth': f_growth,
             'f_death':  f_death}

    dyMet = {'r_biomass': r_biomass,
             'f_biomass': f_biomass}

    model = cmo.DynamicModel(dyRxn, dyMet, None, volExt, volUn, 'optlang-glpk', exitname, timeU=timeunits, savePath=outpath)

    model.initializeConcentrations()
    model.runDynamicFBA(18, 'vode', (0, 0., 1., 1000), False)

    mname1 = exitname.replace(' ','_')
    mname = mname1.split('_')[0]

    x = np.array(model.T)
    bmf = np.array(model.dmetabolites['f_biomass'].quantity)
    bmr = np.array(model.dmetabolites['r_biomass'].quantity)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(x, bmf, 'b-', label='Foxes')
    ax1.plot(x, bmr, 'r-', label='Rabbits')
    ax1.set_ylabel('Population size')
    ax1.set_xlim(0, 18)
    ax1.grid()
    ll = ax1.legend(loc='best', prop={'size':8})
    fig1.savefig('%s/%s.png' % (outpath, mname1))

    return


def options():
    '''define here in-line arguments'''
    parser = optparse.OptionParser(description='Parsing options')
    parser.add_option('-v', '--verbose', dest='verbose', help='increase output verbosity', action='store_true')
    parser.add_option('-a', '--apar', dest='apar', help='parameter a', default='1')
    parser.add_option('-b', '--bpar', dest='bpar', help='parameter b', default='0.1')
    parser.add_option('-c', '--cpar', dest='cpar', help='parameter c', default='1.5')
    parser.add_option('-d', '--dpar', dest='dpar', help='parameter d', default='0.75')
    parser.add_option('-y', '--yicond', dest='yicond', help='initial conditions', default='[10., 5.]')
    parser.add_option('-t', '--timeunits', dest='timeunits', help='time units', default='hr')
    parser.add_option('-p', '--pathout', dest='pathout', help='path for outputs', default='./')
    opts, args = parser.parse_args()
    if opts.verbose:
        print "verbosity turned on"
        print opts
        print args
    return opts, args

if __name__=="__main__":
    inittime = time.time()
    main()
    elapsedtime = time.time() - inittime
    h = divmod(elapsedtime,3600)  # hours
    m = divmod(h[1],60)  # minutes
    s = m[1]  # seconds
    print "Elapsed time: %d hours, %d minutes, %f seconds" % (h[0],m[0],s)
