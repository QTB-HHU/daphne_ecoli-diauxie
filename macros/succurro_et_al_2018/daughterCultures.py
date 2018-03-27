#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2018/03/22     **
#**  last modified: 2018/03/22     **
#************************************

import sys
sys.path.append('../../code/python/') 
import classModel as cmo
import classReaction as cre
import classMetabolite as cme
import classConsortium as cco
import classPlotter as plotter
import classConvertionFactors as ccf
import cobra
import pandas
from cobra.flux_analysis import parsimonious
import numpy as np
import matplotlib.pyplot as plt
import argparse
import plstyles
import random
import cPickle
import json
import copy
from matplotlib import gridspec

# http://bigg.ucsd.edu/models/e_coli_core
# http://bionumbers.hms.harvard.edu//bionumber.aspx?&id=107924&ver=8
### 0.47 gDW of E coli dissolved in 1L give a read of 1 OD600 unit
ODtoGDW=0.33
#ODtoGDW=0.47
OD0= 0.219029
## VOLUME -- Enjalbert 2015 -- 50 mL and 200 mL flasks filled at 15% --> 30mL
volExt = 0.03
volUn = 'L'

## Fig2  : OD0 = 0.213
## Fig6A : OD0 = 0.3
## Fig6C : OD0 = 0.48
expcondLabels = {'batch_high_Ac': 'grown on 45 mM Ac',
                 'batch_low_Ac': 'grown on 4 mM Ac',
                 'batch_mixed_Ac': 'grown on 15 mM Glc, 32 mM Ac',
                 'batch_low_Glc': 'grown on 15 mM Glc',
                 'fedbatch_low_Ac': 'grown on 15 mM Glc, fed 4 mM Ac',
                 'fedbatch_high_Ac': 'grown on 15 mM Glc, fed 32 mM Ac'}

def main():
    '''
    Reproduce results from
    ** Enjalbert et al 2015
    With Single EC model
    '''
    args = options()
    verbose = args.verbose

    biomass0 = float(args.biomassi)
    pcECgl = float(args.ratioecgl)
    thrPc = 0.99999
    if pcECgl > thrPc:
        pcECgl = thrPc
    if pcECgl < (1 - thrPc):
        pcECgl = 1 - thrPc
    strratioecgl = args.ratioecgl.replace('.','p')
    maxtime = float(args.tmax)
    expcond = 'batch_low_Glc'
    glucose0 = 15.0 #mM
    acetate0 = 0.0 #mM
    fb=0.
    afb=0.
    mpath = args.model
    if 'xml' in mpath:
        print 'Running dFBA simulation'
    else:
        print 'Loading simulation file'
        model = cPickle.load(open(mpath, 'r'))
        return model
    
    ac_thr = 0.0
    t_glc = 0.0
    if args.runhighacetate:
        print('Simulating Batch Growth on 45mM Acetate as in Enjalbert2015')
        expcond = 'batch_high_Ac'
        glucose0 = 0.0
        acetate0 = 45.0
    elif args.runlowacetate:
        print('Simulating Batch Growth on 4mM Acetate as in Enjalbert2015')
        expcond = 'batch_low_Ac'
        glucose0 = 0.0
        acetate0 = 4.0
    elif args.runmixedacetate:
        print('Simulating Batch Growth on 15mM Glucose 32mM Acetate as in Enjalbert2015')
        expcond = 'batch_mixed_Ac'
        glucose0 = 15.0
        acetate0 = 32.0
    elif args.runfedlowacetate:
        print('Simulating Batch Growth on 15mM Glucose and constant 4mM Acetate as in Enjalbert2015')
        expcond = 'fedbatch_low_Ac'
        glucose0 = 15.0
        acetate0 = 0.0
        afb=1.0
        ac_thr = 4.0
        t_glc = 4.0
    elif args.runfedhighacetate:
        print('Simulating Batch Growth on 15mM Glucose and constant 32mM Acetate as in Enjalbert2015')
        expcond = 'fedbatch_high_Ac'
        glucose0 = 15.0
        acetate0 = 32.0
        afb=1.0
        ac_thr = 32.0
        t_glc = 4.0
    elif args.runvarmabatch:
        print('Simulating Batch Growth on Glucose as in Varma1994 Fig.7')
        expcond = 'varma_batch'
        glucose0 = 10.8
        acetate0 = 0.4
    elif args.runvarmafedbatch:
        print('Simulating Fed Batch Growth on Glucose as in Varma1994 Fig.10')
        expcond = 'varma_fedbatch'
        glucose0 = 0.82
        acetate0 = 0.1
        #glucose is continuously provided
        #0.2 g/(L*hr)
        #1 g = 5.5 mmol
        #dGLC/dt = ... + fb ==> units = mmol/hr
        fb = 5.5*0.2*volExt #mmol/hr
    else:
        print('Simulating Batch Growth on 15mM Glucose as in Enjalbert2015')

    ename = 'ecoli core'
    exitname = '%s-%s-%s' % (ename, expcond, args.label) 
        
    bmf = 'BIOMASS_Ecoli_core_w_GAM'
    rxnnames = {'EX_glc_e':'EX_glc__D_e', 
                'EX_o2_e': 'EX_o2_e',
                'EX_ac_e':'EX_ac_e'}

    ## biomass0 is already gDW
    ## glucose0 and acetate0 are concentrations (mM)
    ## dMetabolites are initialized with quantities (mmol)
    acetate0 = volExt*acetate0 # mM*L = mmol
    glucose0 = volExt*glucose0 # mM*L = mmol

    #biomass dilution rate in chemostat, 1/hr
    ch = 0. 

    ## Uptake parameters
    ## Uptake Vmax 10 mmol/g/hr (Gosset, 2005)
    ## Uptake Km 10 muM = 0.01 mM (Gosset, 2005)
    vmaxexglc = 10. #mmol/g/hr
    vmaxexace = 10. #mmol/g/hr
    kmuptake = 0.01 #mM
    vminoxygen = float(args.vminexoxygen)
    ubexace = 0.3*vmaxexace
    ubexace = 3.0

    # cm.reactions.get_by_id('ATPM').lower_bound --> 8.39
    #NGAM = 8.39
    death_rate = -1*float(args.deathrate)

    psi0 = 0.
    psi_transition_rate = 0.
    psi_transition_KM = float(args.kmtranspsi) #mM; Ac concentration @ transition_rate/2
    ### ACHTUNG with units!! concentrations in the model!!
    if args.psitransition:
        psi_transition_rate = float(args.vmaxpsi)
        psi0 = float(args.psioffset)
    phi0 = 0.
    phi_transition_rate = 0.
    phi_transition_KM = float(args.kmtransphi) #mM; Glc concentration @ transition_rate/2
    ### ACHTUNG with units!! concentrations in the model!!
    if args.phitransition:
        phi_transition_rate = float(args.vmaxphi)
        phi0 = float(args.phioffset)
    transition_efficiency = float(args.efftrans)
        
    ### ECgl model
    
    biomass_ECgl = cme.Biomass([biomass0*pcECgl], {'growth': [(1, 'biomass_ECgl')], 'dilution': [(ch, None)], 'death': [(death_rate, 'biomass_ECgl')],
                                                   'psi_transition': [(-1., 'biomass_ECgl')], 'phi_transition': [(transition_efficiency, 'biomass_ECac')]})
    ex_glucose_ECgl = cme.DMetabolite('glc_D_e', [glucose0], False, {'glucose_exchange': [(1, 'biomass_ECgl')], 'glucose_fed': [(fb, None)] })
    ex_acetate_ECgl = cme.DMetabolite('ac_e', [acetate0], False, {'acetate_exchange': [(1, 'biomass_ECgl')]})#, 'acetate_fed': [(afb, None)]})

    
    oxygen_exchange_ECgl = cre.DReaction(rxnnames['EX_o2_e'], cre.FixedBound(vminoxygen, 0), False)
    ### Only Ac secretion
    acetate_exchange_ECgl = cre.DReaction(rxnnames['EX_ac_e'], cre.FixedBound(0., ubexace), True)
    if args.runsingle:
        biomass_ECgl = cme.Biomass([biomass0*pcECgl], {'growth': [(1, 'biomass_ECgl')], 'dilution': [(ch, None)], 'death': [(death_rate, 'biomass_ECgl')]})
        acetate_exchange_ECgl = cre.DReaction(rxnnames['EX_ac_e'], cre.MichaelisMenten1(ex_acetate_ECgl, 1, vmaxexace, kmuptake, 1, upperBound=ubexace), False)
    glucose_exchange_ECgl = cre.DReaction(rxnnames['EX_glc_e'], cre.MichaelisMenten1(ex_glucose_ECgl, 1, vmaxexglc, kmuptake, 1, upperBound=15.), False)
    growth_ECgl = cre.DReaction(bmf, cre.FixedBound(0., 1000.))


    ### ECac model
    biomass_ECac = cme.Biomass([biomass0*(1-pcECgl)], {'growth': [(1, 'biomass_ECac')], 'dilution': [(ch, None)], 'death': [(death_rate, 'biomass_ECac')],
                                                       'psi_transition': [(transition_efficiency, 'biomass_ECgl')], 'phi_transition': [(-1., 'biomass_ECac')]})
    ex_glucose_ECac = cme.DMetabolite('glc_D_e', [glucose0], False, {'glucose_exchange': [(1, 'biomass_ECac')] })
    ex_acetate_ECac = cme.DMetabolite('ac_e', [acetate0], False, {'acetate_exchange': [(1, 'biomass_ECac')], 'acetate_fed': [(afb, 'biomass_ECac')]})
    
    oxygen_exchange_ECac = cre.DReaction(rxnnames['EX_o2_e'], cre.FixedBound(vminoxygen, 0), False)
    acetate_exchange_ECac = cre.DReaction(rxnnames['EX_ac_e'], cre.MichaelisMenten1(ex_acetate_ECac, 1, vmaxexace, kmuptake, 1, upperBound=ubexace), False)
    ### EX_glc off
    glucose_exchange_ECac = cre.DReaction(rxnnames['EX_glc_e'], cre.FixedBound(0.,0.), False)
    growth_ECac = cre.DReaction(bmf, cre.FixedBound(0., 1000.))

    #acetate_fed_ECac = cre.DReaction(None, cre.ConcentrationMaintenanceFunction(ex_acetate_ECac, ac_thr, t_glc, linkedReactions=[acetate_exchange_ECac]), True, isODE=True)
    acetate_fed_ECac = cre.DReaction(None, cre.SquareWave(9.1, 20, 1, t_glc), True, isODE=True)


    #########################
    ## Transitions
    #########################

    #lrxn = (growth_ECgl, 1.)
    lrxn = None
    #lrxn2 = (growth_ECac, 1.)
    lrxn2 = None

    hctrans = 5
    #hctrans = 1

    
    if args.glthreshold:
        glccthr = 13.0
        ## phi = Glc -> Ac ##
        ## ECglc senses Ac and transitions
        biomass_psi_transition_ECgl = cre.DReaction(None, cre.MichaelisMentenLinkedSubstrate(ex_acetate_ECgl, lowerBound=0.,
                                                                                             maxVelocity=psi_transition_rate, hillCoeff=hctrans,
                                                                                             mmConstant=psi_transition_KM, linkedSubstrate=ex_glucose_ECgl,
                                                                                             onThrC=0., offThrC=glccthr, offset=psi0), True, isODE=True)
        # Need a copy to avoid duplicate entries; TODO: add option for class model to not update object dMetabolite/dReaction
        #biomass_psi_transition_ECac = copy.copy(biomass_psi_transition_ECgl)
        biomass_psi_transition_ECac = cre.DReaction(None, cre.MichaelisMentenLinkedSubstrate(ex_acetate_ECac, lowerBound=0.,
                                                                                             maxVelocity=psi_transition_rate, hillCoeff=hctrans,
                                                                                             mmConstant=psi_transition_KM, linkedSubstrate=ex_glucose_ECac,
                                                                                             onThrC=0., offThrC=glccthr, offset=psi0), True, isODE=True)
    else:
        ## phi = Glc -> Ac ##
        ## ECglc senses Ac and transitions
        biomass_psi_transition_ECgl = cre.DReaction(None, cre.MichaelisMentenLinked(ex_acetate_ECgl, lowerBound=0.,
                                                                                    maxVelocity=psi_transition_rate, hillCoeff=hctrans,
                                                                                    mmConstant=psi_transition_KM, linkedReaction=lrxn, onThr=0., offset=psi0), True, isODE=True)
        # Need a copy to avoid duplicate entries; TODO: add option for class model to not update object dMetabolite/dReaction
        #biomass_psi_transition_ECac = copy.copy(biomass_psi_transition_ECgl)
        biomass_psi_transition_ECac = cre.DReaction(None, cre.MichaelisMentenLinked(ex_acetate_ECgl, lowerBound=0.,
                                                                                    maxVelocity=psi_transition_rate, hillCoeff=hctrans,
                                                                                    mmConstant=psi_transition_KM, linkedReaction=lrxn, onThr=0., offset=psi0), True, isODE=True)

    ## phi = Ac -> Glc ##
    ## ECac senses Glc and transitions
    biomass_phi_transition_ECac = cre.DReaction(None, cre.MichaelisMentenLinked(ex_glucose_ECac, lowerBound=0.,
                                                                                maxVelocity=phi_transition_rate, hillCoeff=hctrans,
                                                                                mmConstant=phi_transition_KM, linkedReaction=lrxn2, onThr=0., offset=phi0), True, isODE=True)
    # Need a copy to avoid duplicate entries; TODO: add option for class model to not update object dMetabolite/dReaction
    #biomass_phi_transition_ECgl = copy.copy(biomass_phi_transition_ECac)
    biomass_phi_transition_ECgl = cre.DReaction(None, cre.MichaelisMentenLinked(ex_glucose_ECac, lowerBound=0.,
                                                                                maxVelocity=phi_transition_rate, hillCoeff=hctrans,
                                                                                mmConstant=phi_transition_KM, linkedReaction=lrxn2, onThr=0., offset=phi0), True, isODE=True)
    #acetate_fed_ECgl = cre.DReaction(None, cre.ConcentrationMaintenanceFunction(ex_acetate_ECgl, ac_thr, t_glc, linkedReactions=[acetate_exchange_ECgl]), True, isODE=True)
    #acetate_fed_ECgl = cre.DReaction(None, cre.SquareWave(4, 20, 1, t_glc), True, isODE=True)
    
    dyRxn_ECgl = {'glucose_exchange': glucose_exchange_ECgl, 
                  'oxygen_exchange': oxygen_exchange_ECgl,
                  #'acetate_fed': acetate_fed_ECgl,
                  'acetate_exchange': acetate_exchange_ECgl,
                  'psi_transition': biomass_psi_transition_ECgl,
                  'phi_transition': biomass_phi_transition_ECgl,
                  'growth': growth_ECgl}
    dyMet_ECgl = {'biomass_ECgl': biomass_ECgl,
                  'ex_glucose': ex_glucose_ECgl,
                  'ex_acetate': ex_acetate_ECgl}

    model_ECgl = cmo.DynamicModel(dyRxn_ECgl, dyMet_ECgl, mpath, volExt, volUn, 'optlang-glpk', exitname+'_ECgl', savePath=args.pathout)

    
    dyRxn_ECac = {'glucose_exchange': glucose_exchange_ECac, 
                  'oxygen_exchange': oxygen_exchange_ECac,
                  'acetate_fed': acetate_fed_ECac,
                  'acetate_exchange': acetate_exchange_ECac,
                  'psi_transition': biomass_psi_transition_ECac,
                  'phi_transition': biomass_phi_transition_ECac,
                  'growth': growth_ECac}
    dyMet_ECac = {'biomass_ECac': biomass_ECac,
                  'ex_glucose': ex_glucose_ECac,
                  'ex_acetate': ex_acetate_ECac}

    model_ECac = cmo.DynamicModel(dyRxn_ECac, dyMet_ECac, mpath, volExt, volUn, 'optlang-glpk', exitname+'_ECac', savePath=args.pathout)

    model_ECgl.loadSBMLModel()
    model_ECgl.resetObjective()
    model_ECgl.setObjective(growth_ECgl.modName, 1.)
    model_ECgl.constrainSources(element='C', count=0, uptakeOff=True)
    model_ECgl.initializeConcentrations()
    model_ECgl.setParsimoniousFBA(True)
    model_ECgl.setMinimalMaintenance(growth_ECgl.modName, -0.15)
    model_ECgl.setQuitOnFBAError(False)
    model_ECgl.setBoundsThreshold(0.8*vmaxexglc)
    model_ECgl.cbModel.reactions.get_by_id('ATPM').upper_bound = 8.39
    model_ECgl.cbModel.reactions.get_by_id('ATPM').lower_bound = 8.39
    
    model_ECac.loadSBMLModel()
    model_ECac.resetObjective()
    model_ECac.setObjective(growth_ECac.modName, 1.)
    model_ECac.constrainSources(element='C', count=0, uptakeOff=True)
    model_ECac.initializeConcentrations()
    model_ECac.setParsimoniousFBA(True)
    model_ECac.setMinimalMaintenance(growth_ECac.modName, -0.15)
    model_ECac.setQuitOnFBAError(False)
    model_ECac.setBoundsThreshold(0.8*vmaxexace)
    model_ECac.cbModel.reactions.get_by_id('ATPM').upper_bound = 8.39
    model_ECac.cbModel.reactions.get_by_id('ATPM').lower_bound = 8.39
    
    if args.moma:
        model_ECgl.setMOMA(True)
        model_ECac.setMOMA(True)

    model_comm = cco.Consortium({'ECgl': model_ECgl, 'ECac': model_ECac},
                                {'biomass_ECgl': ['ECgl'], 
                                 'biomass_ECac': ['ECac'],
                                 'ex_glucose': ['ECgl', 'ECac'],
                                 'ex_acetate': ['ECgl', 'ECac']
                                }, 'optlang-glpk', name=exitname+'_ECgl_ECac_'+strratioecgl, savePath=args.pathout)

    print('<<<>>> Simulation parameters: ')
    print('<<<>>> Run type: consortium | single | ecgl | ecac ')
    print('<<<>>> \t %d | %d | %d | %d ' % (args.runconsortium, args.runsingle, args.runecgl, args.runecac))
    print('<<<>>> Exp condition: high acetate | low acetate | mixed acetate | fed low acetate | fed high acetate')
    print('<<<>>> \t %d | %d | %d | %d | %d ' % (args.runhighacetate, args.runlowacetate, args.runmixedacetate, args.runfedlowacetate, args.runfedhighacetate))
    print('<<<>>> BM0  \t ECgl \t ECac \t Glc0 \t Ac0 \t death rate')
    print('<<<>>> %.3f \t %.3f \t %.3f \t %.3f \t %.3f \t %.3f ' % (biomass0, pcECgl, (1-pcECgl), glucose0, acetate0, death_rate))
    print('<<<>>> VMGlc \t KMGlc \t VMAc \t KMAc \t XI \t ZETA')
    print('<<<>>> %.3f \t %.3f \t %.3f \t %.3f \t %.3f \t %.3f ' % (vmaxexglc, kmuptake, vmaxexace, kmuptake, fb, afb))
    print('<<<>>> PSI transition %d ' % args.psitransition)
    print('<<<>>> Psi0 \t VMpsi \t KMpsi \t epsilonpsi')
    print('<<<>>> %.3f \t %.3f \t %.3f \t %.3f ' % (psi0, psi_transition_rate, psi_transition_KM, transition_efficiency))
    print('<<<>>> PHI transition %d ' % args.phitransition)
    print('<<<>>> Phi0 \t VMphi \t KMphi \t epsilonphi')
    print('<<<>>> %.3f \t %.3f \t %.3f \t %.3f ' % (phi0, phi_transition_rate, phi_transition_KM, transition_efficiency))

    if verbose:
        for m in model_comm.dmodelsKeys:
            for r in model_comm.dmodels[m].dreactionsKeys:
                print(m, r)
                attrs = vars(model_comm.dmodels[m].dreactions[r].kinF)
                print ', '.join("%s: %s" % item for item in attrs.items())
    
    if not args.run:
        return model_comm
        
    if args.runconsortium:
        model_comm.runConsortiumDynamicFBA(1.5, 'dopri5', (1, 0., 1.5, 1081), verbose=False)
        #model_comm.runConsortiumDynamicFBA(maxtime, 'vode', (0, float(args.minstep), 1., int(args.nsteps)), verbose=False)
    elif args.runecgl:
        model_ECgl.runDynamicFBA(1.5, 'dopri5', (1, 0., 1.5, 1081), verbose=False, quitOnFBAError=False)

    return

def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-A', '--runhighacetate', help='batch 45mM Acetate conditions', action='store_true')
    parser.add_argument('-B', '--runlowacetate', help='batch 4mM Acetate conditions', action='store_true')
    parser.add_argument('-C', '--coremodel', help='use core model', action='store_true')
    parser.add_argument('-D', '--runglucose', help='batch 15mM Glucose conditions (default)', action='store_true')
    parser.add_argument('-E', '--extraplots', help='produce extra plots', action='store_true')
    parser.add_argument('-F', '--runfedlowacetate', help='batch 15mM Glucose fedbatch 4mM Acetate conditions', action='store_true')
    parser.add_argument('-G', '--runmixedacetate', help='batch 15mM Glucose 32mM Acetate conditions', action='store_true')
    parser.add_argument('-H', '--runfedhighacetate', help='batch 15mM Glucose fedbatch 32mM Acetate conditions', action='store_true')
    parser.add_argument('-I', '--runvarmabatch', help='batch Glucose conditions as in Varma Fig 7', action='store_true')
    parser.add_argument('-J', '--death', help='add death rate', action='store_true')
    parser.add_argument('-K', '--runvarmafedbatch', help='fed-batch Glucose conditions as in Varma Fig 10', action='store_true')
    parser.add_argument('-M', '--moma', help='run dFBA with MOMA', action='store_true')
    parser.add_argument('-P', '--phitransition', help='activate phi transition from Ac state to Glc state', action='store_true')
    parser.add_argument('-Q', '--glthreshold', help='Glc concentration threshold for psi transition', action='store_true')
    parser.add_argument('-R', '--run', help='run the dFBA', action='store_true')
    parser.add_argument('-S', '--runsingle', help='run single', action='store_true')
    parser.add_argument('-T', '--psitransition', help='activate psi transition from Glc state to Ac state', action='store_true')
    parser.add_argument('-X', '--runecgl', help='run ECgl', action='store_true')
    parser.add_argument('-Y', '--runecac', help='run ECac', action='store_true')
    parser.add_argument('-Z', '--runconsortium', help='run consortium', action='store_true')
    parser.add_argument('-b', '--biomassi', help='initial biomass', default='0.001')
    parser.add_argument('-d', '--deathrate', help='death rate', default='0.03')
    parser.add_argument('-e', '--efftrans', help='transition efficiency', default='1.')
    parser.add_argument('-j', '--kmtranspsi', help='Km psi transition', default='0.')
    parser.add_argument('-k', '--kmtransphi', help='Km phi transition', default='0.')
    parser.add_argument('-l', '--label', help='output file label', default='diauxic_shift')
    parser.add_argument('-m', '--model', help='model path to load', default='../../ecoli/bigg/e_coli_core.xml')
    parser.add_argument('-n', '--nsteps', help='max number of steps', default='1000')
    parser.add_argument('-o', '--ofile', help='output file name', default='toymodel-test.xml')
    parser.add_argument('-p', '--pathout', help='path for outputs', default='/tmp')
    parser.add_argument('-r', '--ratioecgl', help='ratio of ECgl/ECac: BMECgl0 = X*biomass0', default='1.')
    parser.add_argument('-s', '--minstep', help='min step size', default='0.0')
    parser.add_argument('-t', '--tmax', help='max time', default='10.5')
    parser.add_argument('-x', '--vminexoxygen', help='V min EX O2', default='10.')
    parser.add_argument('-y', '--vmaxpsi', help='V max psi transition', default='0.0')
    parser.add_argument('-z', '--vmaxphi', help='V max phi transition', default='0.0')
    parser.add_argument('-w', '--psioffset', help='offset psi transition', default='0.0')
    parser.add_argument('-u', '--phioffset', help='offset phi transition', default='0.0')
    args = parser.parse_args()
    if args.verbose:
        print "verbosity turned on"
        print args
    return args

if __name__=="__main__":
    main()
