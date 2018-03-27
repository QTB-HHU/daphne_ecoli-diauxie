#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  last modified: 2015/08/05     **
#************************************

import cobra as cb
import numpy as np
import optparse

def buildToyModel(toyv=0, cofbm=0, costq3=10., verbose=False):
    '''
    function to build a toymodel
    '''

    toymodel = cb.Model('toymodel')

    rxnNames = ['EX_subs(e)', 'EX_cofa(e)', 'STRANSP', 'CTRANSP', 'Q1SYN', 'Q2SYN', 'Biomass', 'IN_cofa(i)']
    rxnBounds = {'EX_subs(e)': (-10, 1000), 'EX_cofa(e)': (-10, 1000), 'IN_cofa(i)':(0, 0.1),
                 'Q1SYN':(0, 10), 'Q2SYN':(0, 10), 'Biomass':(0, 1000)}
    metNames = ['subs[e]', 'cofa[e]', 'subs[i]', 'cofa[i]', 'q1[i]', 'q2[i]']

    rxns = {}
    mets = {}

    compartments = {'e': 'external', 'i': 'internal'}
    
    addRxnNames = []
    addRxnBounds = []
    addMetNames = []

    if toyv == 2:
        addRxnNames = ['IESYN', 'DESYN']
        addRxnBounds = [(0, 10.), (0, 10.)]
    elif toyv == 3:
        addRxnNames = ['ST_cofa(i)', 'US_cofa(s)', 'IN_cofa(s)']
        addRxnBounds = [(0, 0.0), (0, 0.0), (0, 0.0)]
        addMetNames = ['cofa[s]']
    elif toyv == 4:
        addRxnNames = ['ST_cofa(i)', 'US_cofa(s)', 'IN_cofa(s)', 'Q3SYN']
        addRxnBounds = [(0, 0.1), (0, 0.1), (0, 0.1), (0, 10)]
        addMetNames = ['cofa[s]', 'q3[i]']

    for r in rxnNames:
        rxns[r] = cb.Reaction(r)
        rxns[r].name = r
        bounds = rxnBounds.get(r, (-1000., 1000.))
        rxns[r].lower_bound = bounds[0]
        rxns[r].upper_bound = bounds[1]

    for r,b in zip(addRxnNames, addRxnBounds):
        rxns[r] = cb.Reaction(r)
        rxns[r].name = r
        rxns[r].lower_bound = b[0]
        rxns[r].upper_bound = b[1]

    if toyv > 2:
        compartments['s'] = 'storage'

    for m in metNames+addMetNames:
        mets[m] = cb.Metabolite(m, name = m, compartment = m[m.find('[')+1])
        if verbose:
            print mets[m].compartment

    rxns['EX_subs(e)'].add_metabolites({mets['subs[e]'] : -1.0})
    rxns['EX_cofa(e)'].add_metabolites({mets['cofa[e]'] : -1.0})
    rxns['STRANSP'].add_metabolites({mets['subs[e]'] : -1.0,
                                     mets['subs[i]'] : 1.0})
    rxns['CTRANSP'].add_metabolites({mets['cofa[e]'] : -1.0,
                                     mets['cofa[i]'] : 1.0})
    rxns['Q1SYN'].add_metabolites({mets['subs[i]'] : -1.0,
                                   mets['q1[i]'] : 1.0,})
    rxns['Q2SYN'].add_metabolites({mets['subs[i]'] : -1.0,
                                   mets['q2[i]'] : 1.0,})
    if cofbm < 0:
        rxns['Biomass'].add_metabolites({mets['q2[i]'] : -1.0,
                                        mets['q1[i]'] : -1.0,
                                        mets['cofa[i]'] : cofbm})
    else:
        rxns['Biomass'].add_metabolites({mets['q2[i]'] : -1.0,
                                        mets['q1[i]'] : -1.0,})
    rxns['IN_cofa(i)'].add_metabolites({mets['cofa[i]'] : -1.0})


    if toyv == 1:
        rxns['IN_cofa(i)'].add_metabolites({mets['q2[i]'] : -1.0})
    elif toyv == 2:
        rxns['IESYN'].add_metabolites({mets['q1[i]']: -0.5,
                               mets['q2[i]'] : -0.5})
        rxns['DESYN'].add_metabolites({mets['q1[i]']: -0.1,
                               mets['q2[i]'] : -0.1,
                               mets['cofa[i]'] : -0.05})
    elif toyv == 3:
        rxns['ST_cofa(i)'].add_metabolites({mets['cofa[i]'] : -1.0,
                                            mets['cofa[s]'] : 1.0})
        rxns['US_cofa(s)'].add_metabolites({mets['cofa[s]'] : -1.0,
                                            mets['cofa[i]'] : 1.0})
        rxns['IN_cofa(s)'].add_metabolites({mets['cofa[s]'] : -1.0})
    elif toyv == 4:
        rxns['Q3SYN'].add_metabolites({mets['subs[i]'] : -1.0,
                                        mets['q3[i]'] : 1.0,})
        rxns['ST_cofa(i)'].add_metabolites({mets['cofa[i]'] : -1.0,
                                            mets['q3[i]'] : -1*costq3,#SUBS
                                            mets['cofa[s]'] : 1.0})
        rxns['US_cofa(s)'].add_metabolites({mets['cofa[s]'] : -1.0,
                                            mets['cofa[i]'] : 1.0})
        rxns['IN_cofa(s)'].add_metabolites({mets['cofa[s]'] : -1.0})

    for r in rxnNames+addRxnNames:
        toymodel.add_reaction(rxns[r])
        if verbose:
            print rxns[r].reaction
            print rxns[r].objective_coefficient

    #toymodel.objective = 'Biomass + IN_cofa(i)'
    rxns['Biomass'].objective_coefficient = 1
    rxns['IN_cofa(i)'].objective_coefficient = 1

    toymodel.compartments = compartments
    return toymodel

def options():
    '''define here in-line arguments'''
    parser = optparse.OptionParser(description='Parsing options')
    parser.add_option('-v', '--verbose', dest='verbose', help='increase output verbosity', action='store_true')
    parser.add_option('-m', '--modelversion', dest='modelversion', help='integer for model version', default='0')
    parser.add_option('-o', '--ofile', dest='ofile', help='output file name', default='toymodel-test.xml')
    opts, args = parser.parse_args()
    if opts.verbose:
        print "verbosity turned on"
        print opts
        print args
    return opts, args

def main():
    opts, args = options()
    tm = buildToyModel(int(opts.modelversion), verbose=opts.verbose)
    #cb.io.write_sbml_model(tm, opts.ofile, use_fbc_package=False)
    cb.io.write_sbml_model(tm, opts.ofile, use_fbc_package=False)

if __name__=="__main__":
    main()
