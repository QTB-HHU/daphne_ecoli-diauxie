import numpy as np
import cobra
import pandas
from cobra.flux_analysis import parsimonious
import argparse
from matplotlib import gridspec

def main():
    '''
    Reproduce results from
    ** Enjalbert et al 2015
    With Single EC model
    '''
    args = options()
    verbose = args.verbose

    bmname = 'BIOMASS_Ecoli_core_w_GAM'
    m = cobra.io.read_sbml_model(args.model)
    exrxnname = ['EX_glc__D_e', 'EX_ac_e', 'EX_o2_e']
    exrxn = []

    otname = args.pathout+'tab_v.csv'
    if args.constraintAcSec:
        m.reactions.get_by_id('EX_ac_e').upper_bound = 3
        otname = args.pathout+'maxAcOut_tab_v.csv'
        
    for n in exrxnname:
        exrxn.append(m.reactions.get_by_id(n))

    vsp = np.linspace(-10, 0, 11)
    for r in exrxn:
        r.lower_bound = vsp[0]

    msol = parsimonious.pfba(m)
    df_tmp = msol.x_dict[exrxnname+[bmname]]

    idxlab = '_'.join(3*['%.1f' % vsp[0]])
    df = pandas.DataFrame(data=df_tmp).T
    df.index = [idxlab]

    # with glucose
    for r in exrxn:
        for v in vsp:
            r.lower_bound = v
            msol = parsimonious.pfba(m)
            df_tmp = msol.x_dict[exrxnname+[bmname]]
            lab = []
            for rl in exrxn:
                lab.append('%.1f' % rl.lower_bound)
            df_tmp.name = '_'.join(lab)
            df = df.append(df_tmp)
        r.lower_bound = vsp[0]

    # without glucose
    exrxn[0].lower_bound = 0
    for r in exrxn[1:]:
        for v in vsp:
            r.lower_bound = v
            try:
                msol = parsimonious.pfba(m)
                df_tmp = msol.x_dict[exrxnname+[bmname]]
            except:
                msol = None
                df_tmp = pandas.Series(dict(zip(exrxnname+[bmname], 4*[np.nan])))
            lab = []
            for rl in exrxn:
                lab.append('%.1f' % rl.lower_bound)
            df_tmp.name = '_'.join(lab)
            df = df.append(df_tmp)
        r.lower_bound = vsp[0]

    df.to_csv(otname)
    return

def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-C', '--constraintAcSec', help='max acetate secretion', action='store_true')
    parser.add_argument('-p', '--pathout', help='path for outputs', default='../../outputs/mytests/')
    parser.add_argument('-r', '--rxnlist', help='list of reactions to constraint', default='')
    parser.add_argument('-m', '--model', help='model path to load', default='../../ecoli/bigg/e_coli_core.xml')
    args = parser.parse_args()
    if args.verbose:
        print "verbosity turned on"
        print args
    return args

if __name__=="__main__":
    main()
