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
import matplotlib.pyplot as plt
from matplotlib import gridspec
import argparse
import plstyles
import cPickle
import json
import seaborn as sb
import numpy as np

def main():
    '''
    Test growth matrix data on different media
    '''
    args = options()
    verbose = args.verbose

    dmodel = cPickle.load(open(args.fname, 'r'))
    if dmodel.FBAsolutions[-1] is None:
        del dmodel.FBAsolutions[-1]
        del dmodel.T[-1]

    mplot = plotter.Plotter({args.name: dmodel}, args.pathout)

    if args.timex:
        xvars = (np.array(dmodel.T), 'Time', dmodel.Tunits)
    else:
        xvars = eval(args.xvars)
    mplot.plotReactions(args.name, eval(args.yvars), xvars, args.oname, eval(args.lims))
    return


def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-T', '--timex', help='x axis is time', action='store_true')
    parser.add_argument('-R', '--relflux', help='plot relative fluxes', action='store_true')
    parser.add_argument('-p', '--pathout', help='path for outputs', default='../../outputs/mytests/')
    parser.add_argument('-n', '--name', help='plot name', default='varma1994_batch')
    parser.add_argument('-o', '--oname', help='out file name', default='varma1994_batch')
    parser.add_argument('-f', '--fname', help='file name', default='../../outputs/mytests/endOfSimulation-ecoli_core-varma_batch-varma_o2_11p5_ac_3_ECgl.p')
    parser.add_argument('-y', '--yvars', help='y fluxes in plotter format', default="[['glucose_exchange', ('flux','lb')],['acetate_exchange', ('flux','lb','ub')], ['oxygen_exchange', ('flux',)]]")
    parser.add_argument('-x', '--xvars', help='x variable in plotter format', default="(np.array(dmodel.T), 'Time', dmodel.Tunits)")
    parser.add_argument('-l', '--lims', help='x limit range', default="(7,8.5)")
    args = parser.parse_args()
    if args.verbose:
        print "verbosity turned on"
        print args
    return args

if __name__=="__main__":
    main()
