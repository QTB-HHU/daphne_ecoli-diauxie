import numpy as np
import argparse

def main():
    '''
    
    '''
    args = options()
    verbose = args.verbose
    od = float(args.odval)
    vol = float(args.vol)
    bm = float(args.bm)
    ODtoGDW = float(args.conv)

    gdw = od*vol*ODtoGDW
    print('%.3f OD corresponds to %.4f gDW in %.3f L using as conversion factor 1 OD = %.3f gDW/L' % (od, gdw, vol, ODtoGDW))
    

def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-o', '--odval', help='od600 value', default='0.213')
    parser.add_argument('-v', '--vol', help='volume value in L', default='0.03')
    parser.add_argument('-b', '--bm', help='biomass value', default='0.0')
    parser.add_argument('-c', '--conv', help='conversion factor value in gDW/(L*OD)', default='0.47')
    args = parser.parse_args()
    if args.verbose:
        print "verbosity turned on"
        print args
    return args

if __name__=="__main__":
    main()
