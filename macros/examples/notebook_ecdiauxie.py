#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2018/05/15     **
#**  last modified: 2018/05/15     **
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


expcondLabels = {'batch_high_Ac': 'grown on 45 mM Ac',
                 'batch_low_Ac': 'grown on 4 mM Ac',
                 'batch_mixed_Ac': 'grown on 15 mM Glc, 32 mM Ac',
                 'batch_low_Glc': 'grown on 15 mM Glc',
                 'fedbatch_low_Ac': 'grown on 15 mM Glc, fed 4 mM Ac',
                 'fedbatch_high_Ac': 'grown on 15 mM Glc, fed 32 mM Ac'}

# http://bigg.ucsd.edu/models/e_coli_core
ODtoGDW=0.33
volExt = 0.03
volUn = 'L'

maxtime = 10.5
biomass0 = 0.0038
mpath = '../../ecoli/bigg/e_coli_core.xml'

thrPc = 0.99999
pcECgl = 0.95
strratioecgl = ("%.2f" % pcECgl).replace('.','p')
if pcECgl > thrPc:
    pcECgl = thrPc
if pcECgl < (1 - thrPc):
    pcECgl = 1 - thrPc

expcond = 'batch_low_Glc'
glucose0 = 15.0 #mM
acetate0 = 0.0 #mM
fb=0.
afb=0.
ac_thr = 0.0
t_glc = 0.0

runhighacetate = runlowacetate = runbatchglc = runmixedacetate = runfedlowacetate = runfedhighacetate = runvarmabatch = runvarmafedbatch =False
runfedlowacetate = True

if runhighacetate:
    print('Simulating Batch Growth on 45mM Acetate as in Enjalbert2015')
    expcond = 'batch_high_Ac'
    glucose0 = 0.0
    acetate0 = 45.0
elif runlowacetate:
    print('Simulating Batch Growth on 4mM Acetate as in Enjalbert2015')
    expcond = 'batch_low_Ac'
    glucose0 = 0.0
    acetate0 = 4.0
elif runbatchglc:
    print('Simulating Batch Growth on 15mM Glucose as in Enjalbert2015')
    expcond = 'batch_low_Glc'
    glucose0 = 15.0 #mM
    acetate0 = 0.0 #mM
    fb=0.
    afb=0.
elif runmixedacetate:
    print('Simulating Batch Growth on 15mM Glucose 32mM Acetate as in Enjalbert2015')
    expcond = 'batch_mixed_Ac'
    glucose0 = 15.0
    acetate0 = 32.0
elif runfedlowacetate:
    print('Simulating Batch Growth on 15mM Glucose and constant 4mM Acetate as in Enjalbert2015')
    expcond = 'fedbatch_low_Ac'
    glucose0 = 15.0
    acetate0 = 0.0
    afb=1.0
    ac_thr = 4.0
    t_glc = 4.0
elif runfedhighacetate:
    print('Simulating Batch Growth on 15mM Glucose and constant 32mM Acetate as in Enjalbert2015')
    expcond = 'fedbatch_high_Ac'
    glucose0 = 15.0
    acetate0 = 32.0
    afb=1.0
    ac_thr = 32.0
    t_glc = 4.0
elif runvarmabatch:
    print('Simulating Batch Growth on Glucose as in Varma1994 Fig.7')
    expcond = 'varma_batch'
    glucose0 = 10.8
    acetate0 = 0.4
elif runvarmafedbatch:
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
    print('Simulating custom conditions')

label='diauxic_shift'
ename = 'ecoli core'
exitname = '%s-%s-%s' % (ename, expcond, label) 
        
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

death_rate = -0.03

## Uptake parameters
## Uptake Vmax 10 mmol/g/hr (Gosset, 2005)
## Uptake Km 10 muM = 0.01 mM (Gosset, 2005)
vmaxexglc = 10. #mmol/g/hr
vmaxexace = 10. #mmol/g/hr
kmuptake = 0.01 #mM
## Parametrized with Varma 1994
vminoxygen = 11.5
ubexace = 3.0


transition_efficiency = 0.9

# No transition
psi0 = 0.
psi_transition_rate = 0.
psi_transition_KM = 0.
phi0 = 0.
phi_transition_rate = 0.
phi_transition_KM = 0.

# Noise transition
psi0 = 0.04
phi0 = 0.04
psi_transition_rate = 0.
psi_transition_KM = 0.
phi_transition_rate = 0.
phi_transition_KM = 0.

# Concentration dependent transition
psi0 = 0.04
phi0 = 0.04
psi_transition_rate = 0.2
phi_transition_rate = 0.2
psi_transition_KM = 30.0
phi_transition_KM = 5.0

        
### ECgl model
biomass_ECgl = cme.Biomass([biomass0*pcECgl], {'growth': [(1, 'biomass_ECgl')], 'dilution': [(ch, None)], 'death': [(death_rate, 'biomass_ECgl')],
                                               'psi_transition': [(-1., 'biomass_ECgl')], 'phi_transition': [(transition_efficiency, 'biomass_ECac')]})
ex_glucose_ECgl = cme.DMetabolite('glc_D_e', [glucose0], False, {'glucose_exchange': [(1, 'biomass_ECgl')], 'glucose_fed': [(fb, None)] })
ex_acetate_ECgl = cme.DMetabolite('ac_e', [acetate0], False, {'acetate_exchange': [(1, 'biomass_ECgl')]})#, 'acetate_fed': [(afb, None)]})
    
oxygen_exchange_ECgl = cre.DReaction(rxnnames['EX_o2_e'], cre.FixedBound(vminoxygen, 0), False)
### Only Ac secretion
acetate_exchange_ECgl = cre.DReaction(rxnnames['EX_ac_e'], cre.FixedBound(0., ubexace), True)

# if args.runsingle:
#     biomass_ECgl = cme.Biomass([biomass0*pcECgl], {'growth': [(1, 'biomass_ECgl')], 'dilution': [(ch, None)], 'death': [(death_rate, 'biomass_ECgl')]})
#     acetate_exchange_ECgl = cre.DReaction(rxnnames['EX_ac_e'], cre.MichaelisMenten1(ex_acetate_ECgl, 1, vmaxexace, kmuptake, 1, upperBound=ubexace), False)

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

lrxn = None
lrxn2 = None

hctrans = 5

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

model_ECgl = cmo.DynamicModel(dyRxn_ECgl, dyMet_ECgl, mpath, volExt, volUn, 'optlang-glpk', exitname+'_ECgl', savePath='./')

    
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

model_ECac = cmo.DynamicModel(dyRxn_ECac, dyMet_ECac, mpath, volExt, volUn, 'optlang-glpk', exitname+'_ECac', savePath='./')

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
    
#if args.moma:
#    model_ECgl.setMOMA(True)
#    model_ECac.setMOMA(True)

model_comm = cco.Consortium({'ECgl': model_ECgl, 'ECac': model_ECac},
                            {'biomass_ECgl': ['ECgl'], 
                             'biomass_ECac': ['ECac'],
                             'ex_glucose': ['ECgl', 'ECac'],
                             'ex_acetate': ['ECgl', 'ECac']
                            }, 'optlang-glpk', name=exitname+'_ECgl_ECac_'+strratioecgl, savePath='./')

print('<<<>>> Simulation parameters: ')
print('<<<>>> Exp condition: high acetate | low acetate | mixed acetate | fed low acetate | fed high acetate')
print('<<<>>> \t %d | %d | %d | %d | %d ' % (runhighacetate, runlowacetate, runmixedacetate, runfedlowacetate, runfedhighacetate))
print('<<<>>> BM0  \t ECgl \t ECac \t Glc0 \t Ac0 \t death rate')
print('<<<>>> %.3f \t %.3f \t %.3f \t %.3f \t %.3f \t %.3f ' % (biomass0, pcECgl, (1-pcECgl), glucose0, acetate0, death_rate))
print('<<<>>> VMGlc \t KMGlc \t VMAc \t KMAc \t XI \t ZETA')
print('<<<>>> %.3f \t %.3f \t %.3f \t %.3f \t %.3f \t %.3f ' % (vmaxexglc, kmuptake, vmaxexace, kmuptake, fb, afb))
print('<<<>>> PSI transition ')
print('<<<>>> Psi0 \t VMpsi \t KMpsi \t epsilonpsi')
print('<<<>>> %.3f \t %.3f \t %.3f \t %.3f ' % (psi0, psi_transition_rate, psi_transition_KM, transition_efficiency))
print('<<<>>> PHI transition')
print('<<<>>> Phi0 \t VMphi \t KMphi \t epsilonphi')
print('<<<>>> %.3f \t %.3f \t %.3f \t %.3f ' % (phi0, phi_transition_rate, phi_transition_KM, transition_efficiency))


model_comm.runConsortiumDynamicFBA(maxtime, 'vode', (0, 0.0, 1., 10000), verbose=False)

#if not args.run:
#    return model_comm
        
#if args.runconsortium:
#    model_comm.runConsortiumDynamicFBA(maxtime, 'vode', (0, float(args.minstep), 1., int(args.nsteps)), verbose=False)
#elif args.runecgl:
#    model_ECgl.runDynamicFBA(maxtime, 'vode', (0, float(args.minstep), 1., int(args.nsteps)), False, quitOnFBAError=False)
#elif args.runecac:
#    model_ECac.runDynamicFBA(maxtime, 'vode', (0, float(args.minstep), 1., int(args.nsteps)), False, quitOnFBAError=False)


def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-M', '--moma', help='run dFBA with MOMA', action='store_true')
    parser.add_argument('-P', '--phitransition', help='activate phi transition from Ac state to Glc state', action='store_true')
    parser.add_argument('-Q', '--glthreshold', help='Glc concentration threshold for psi transition', action='store_true')
    parser.add_argument('-R', '--run', help='run the dFBA', action='store_true')
    parser.add_argument('-S', '--runsingle', help='run single', action='store_true')
    parser.add_argument('-T', '--psitransition', help='activate psi transition from Glc state to Ac state', action='store_true')
    parser.add_argument('-X', '--runecgl', help='run ECgl', action='store_true')
    parser.add_argument('-Y', '--runecac', help='run ECac', action='store_true')
    parser.add_argument('-Z', '--runconsortium', help='run consortium', action='store_true')
    parser.add_argument('-n', '--nsteps', help='max number of steps', default='1000')
    parser.add_argument('-o', '--ofile', help='output file name', default='toymodel-test.xml')
    parser.add_argument('-p', '--pathout', help='path for outputs', default='/tmp')
    parser.add_argument('-r', '--ratioecgl', help='ratio of ECgl/ECac: BMECgl0 = X*biomass0', default='1.')
    parser.add_argument('-s', '--minstep', help='min step size', default='0.0')
    args = parser.parse_args()
    if args.verbose:
        print "verbosity turned on"
        print args
    return args
