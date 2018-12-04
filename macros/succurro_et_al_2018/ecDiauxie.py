#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2017/11/11     **
#**  last modified: 2018/01/15     **
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
#ODtoGDW=0.35
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
        exitname = mpath.split('/')[-1].split('.')[0].split('-')
        if args.runconsortium:
            #return model
            titlab = 'Initial biomass %.0f%% ECgl and %.0f%% ECac' % (pcECgl*100, (1-pcECgl)*100)
            if args.runvarmabatch or args.runvarmafedbatch:
                plotVarma1994_twoEC(args, model, exitname[2], '-'.join(exitname[1:]), titlab)
                return
            plotEnjalbert2015_growth2EC(args, model, exitname[2], '-'.join(exitname[1:]), titlab)
            if args.runecgl:
                plotEnjalbert2015_growth1EC(args, model.dmodels['ECgl'], exitname[2], '-'.join(exitname[1:]))
            if args.runecac:
                plotEnjalbert2015_growth1EC(args, model.dmodels['Ecac'], exitname[2], '-'.join(exitname[1:]))
        if args.runsingle:
            if args.runvarmabatch or args.runvarmafedbatch:
                plotVarma1994(args, model, exitname[2], '-'.join(exitname[1:]))
            if args.runecgl:
                plotEnjalbert2015_growth1EC(args, model, exitname[2], '-'.join(exitname[1:]))
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
        model_comm.runConsortiumDynamicFBA(maxtime, 'vode', (0, float(args.minstep), 1., int(args.nsteps)), verbose=False)
    elif args.runecgl:
        model_ECgl.runDynamicFBA(maxtime, 'vode', (0, float(args.minstep), 1., int(args.nsteps)), False, quitOnFBAError=False)
    elif args.runecac:
        model_ECac.runDynamicFBA(maxtime, 'vode', (0, float(args.minstep), 1., int(args.nsteps)), False, quitOnFBAError=False)

    return


def plotVarma1994(args, model, expcond, expcond2):
    # LOAD DATA from Varma1994
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1)
    # Top left
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
    ax1.set_ylabel('Q [gDW]')
    ax2.set_ylabel('C [mM]')
    ax2.set_xlabel('Time [hours]')
    tit = 'E. Coli core aerobic'
    if args.runvarmafedbatch:
        dataAC=pandas.read_csv('../../ecoli/varma1994_data/varma10_ac.csv', sep=',',header=0)
        dataBM=pandas.read_csv('../../ecoli/varma1994_data/varma10_bm.csv', sep=',',header=0)
        dataGL=pandas.read_csv('../../ecoli/varma1994_data/varma10_gl.csv', sep=',',header=0)
        fignum = '_fig10'
        tit = tit+' fed-batch'
    else:
        dataAC=pandas.read_csv('../../ecoli/varma1994_data/varma7_ac.csv', sep=',',header=0)
        dataBM=pandas.read_csv('../../ecoli/varma1994_data/varma7_bm.csv', sep=',',header=0)
        dataGL=pandas.read_csv('../../ecoli/varma1994_data/varma7_gl.csv', sep=',',header=0)
        fignum = '_fig7'
        tit = tit+' batch'
    if args.moma:
        tit = tit+', MOMA'
    else:
        tit = tit+', pFBA'
    ax1.set_title(tit)
    conctomM = 1000. if model.volumes.externalUnits == 'mL' else 1.
    ybm1 = np.array(model.dmetabolites['biomass_ECgl'].quantity)
    x = np.array(model.T)
    ygl1 = np.array(model.dmetabolites['ex_glucose'].concentration)*conctomM
    yac1 = np.array(model.dmetabolites['ex_acetate'].concentration)*conctomM

    ax1.plot(x, ybm1, '#2980b9', label='Biomass EC')
    #ax1.set_xlim(0., 1.)
    ax2.plot(x, ygl1, '#c0392b', label='Glucose')
    ax2.plot(x, yac1,  '#f39c12', label='Acetate')

    ax1.plot(dataBM['x'], dataBM['Curve1']*volExt, 'bs', label='Biomass (exp)')
    ax2.plot(dataGL['x'], dataGL['Curve1'], 'rs', label='Glucose (exp)')
    ax2.plot(dataAC['x'], dataAC['Curve1'], '#f1c40f',  linestyle='None', marker='s', label='Acetate (exp)')

    ll = ax1.legend(loc='best', prop={'size':10})
    ll2 = ax2.legend(loc='best', prop={'size':10})
    fig.savefig(args.pathout+'/varma1994_'+expcond2+fignum+'.png')
    
    return


def plotVarma1994_twoEC(args, model, expcond, expcond2, tit):

    n1, n2 = model.dmodelsKeys
    
    bmn1 = 'biomass_%s' % n1
    bmn2 = 'biomass_%s' % n2

    # LOAD DATA from Varma1994
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1)
    # Top left
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
    ax1.set_ylabel('Q [gDW]')
    ax2.set_ylabel('C [mM]')
    ax2.set_xlabel('Time [hours]')
    if args.phitransition and args.psitransition:
        tit=tit+', w/ transition'
    else:
        tit=tit+', w/o transition'
    
    if args.runvarmafedbatch:
        dataAC=pandas.read_csv('../../ecoli/varma1994_data/varma10_ac.csv', sep=',',header=0)
        dataBM=pandas.read_csv('../../ecoli/varma1994_data/varma10_bm.csv', sep=',',header=0)
        dataGL=pandas.read_csv('../../ecoli/varma1994_data/varma10_gl.csv', sep=',',header=0)
        fignum = '_fig10'
    else:
        dataAC=pandas.read_csv('../../ecoli/varma1994_data/varma7_ac.csv', sep=',',header=0)
        dataBM=pandas.read_csv('../../ecoli/varma1994_data/varma7_bm.csv', sep=',',header=0)
        dataGL=pandas.read_csv('../../ecoli/varma1994_data/varma7_gl.csv', sep=',',header=0)
        fignum = '_fig7'
    ax1.set_title(tit)
    conctomM = 1000. if model.dmodels[n1].volumes.externalUnits == 'mL' else 1.

    x = np.array(model.dmodels[n1].T)
    x2=np.array(model.dmodels[n2].T)
    ybm1 = np.array(model.dmodels[n1].dmetabolites[bmn1].quantity)
    ygl1 = np.array(model.dmodels[n1].dmetabolites['ex_glucose'].concentration)*conctomM
    yac1 = np.array(model.dmodels[n1].dmetabolites['ex_acetate'].concentration)*conctomM
    #yam1 = np.array(model.dmodels[n1].dmetabolites['ex_ammonium'].quantity)
    ybm2 = np.array(model.dmodels[n2].dmetabolites[bmn2].quantity)
    yac2 = np.array(model.dmodels[n2].dmetabolites['ex_acetate'].concentration)*conctomM
    ygl2 = np.array(model.dmodels[n2].dmetabolites['ex_glucose'].concentration)*conctomM
    #yam2 = np.array(model.dmodels[n2].dmetabolites['ex_ammonium'].quantity)

    ax1.plot(x, ybm1, '#2980b9', label=bmn1.replace('_', ' '))
    ax1.plot(x, ybm2, '#16a085', label=bmn2.replace('_', ' '))
    ax1.plot(x, ybm1+ybm2, 'k', label='biomass Tot')
    #ax1.set_xlim(0., 1.)
    ll = ax1.legend(loc='center right', prop={'size':10})
    ax1.set_title(tit)
    ax2.plot(x, ygl1, '#c0392b', label='Glucose')
    ax2.plot(x, yac1, '#f39c12', label='Acetate')

    ax1.plot(dataBM['x'], dataBM['Curve1']*volExt, 'bs', label='Biomass (exp)')
    ax2.plot(dataGL['x'], dataGL['Curve1'], 'rs', label='Glucose (exp)')
    ax2.plot(dataAC['x'], dataAC['Curve1'], '#f1c40f',  linestyle='None', marker='s', label='Acetate (exp)')

    ll = ax1.legend(loc='best', prop={'size':10})
    ll2 = ax2.legend(loc='best', prop={'size':10})
    fig.savefig(args.pathout+'/varma1994_'+expcond2+fignum+'.png')
    
    return

def plotEnjalbert2015_growth1EC(args, model, expcond, expcond2):
    #fig to compare with consortium fig
    n1 = 'ECgl' if args.runecgl else 'Ecac'
    n2 = 'E. coli' if args.runsingle else n1
    bmn1 = 'biomass_%s' % n1
    conctomM = 1000. if model.volumes.externalUnits == 'mL' else 1.
    ybm1 = np.array(model.dmetabolites[bmn1].quantity)
    x = np.array(model.T)
    ygl1 = np.array(model.dmetabolites['ex_glucose'].concentration)*conctomM
    yac1 = np.array(model.dmetabolites['ex_acetate'].concentration)*conctomM

    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1)
    # Top left
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
    ax1.set_ylabel('Q [gDW]')
    ax2.set_ylabel('C [mM]')
    ax2.set_xlabel('Time [hours]')
    tit = '%s %s' % (n2, expcondLabels[expcond])
    ## COLORS http://htmlcolorcodes.com/color-chart/
    ## blues: #2980b9 #3498db  #1abc9c #16a085
    ## yellows: #f4d03f f5b041 eb984e  #dc7633
    if args.runsingle:
        bmn1 = 'biomass_EC'
        #n1 = 'EC'
        if args.moma:
            tit = tit+', MOMA'
        else:
            tit = tit+', pFBA'

    ax1.plot(x, ybm1, '#2980b9', label=bmn1.replace('_', ' '))
    #ax1.set_xlim(0., 1.)
    ll = ax1.legend(loc='center right', prop={'size':10})
    ax1.set_title(tit)
    ax2.plot(x, ygl1, '#c0392b', label='Glucose')
    ax2.plot(x, yac1,  '#f39c12', label='Acetate')

    fignum = ''
    if expcond == "batch_low_Glc":
        # dataAC=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2A_ac.csv', sep=',',header=None, names=['x','y'])
        # dataBM=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2A_bm.csv', sep=',',header=None, names=['x','y'])
        # dataGL=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2A_gl.csv', sep=',',header=None, names=['x','y'])
        # ax1.plot(dataBM['x'], dataBM['y']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        # ax2.plot(dataGL['x'], dataGL['y'], 'rs', label='Glucose (exp)')
        # ax2.plot(dataAC['x'], dataAC['y'], '#f1c40f',  linestyle='None', marker='s', label='Acetate (exp)')
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2_fromEnjalbert2015.csv', sep=',')
        ax1.errorbar(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, yerr=dataFIG['OD SD']*ODtoGDW*volExt, fmt='bs', label='Biomass (exp)')
        ax2.errorbar(dataFIG['Time'], dataFIG['Glucose mM'], yerr=dataFIG['Glucose SD'], fmt='rs', label='Glucose (exp)')
        ax2.errorbar(dataFIG['Time'], dataFIG['Acetate mM'], yerr=dataFIG['Acetate SD'], label='Acetate (exp)', color='#f1c40f',  linestyle='None', marker='s')
        fignum = '_fig2A'
    elif expcond == "fedbatch_low_Ac":
        # dataBM=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6a_bm.csv', sep=',',header=None, names=['x','y'])
        # ax1.plot(dataBM['x'], dataBM['y']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        #ax1.errorbar(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, yerr=dataFIG['OD SD']*ODtoGDW*volExt, fmt='bs', label='Biomass (exp)')
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6a_fromEnjalbert2015.csv', sep=',')
        ax1.plot(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        fignum = '_fig6A'
    elif expcond == "fedbatch_high_Ac":
        # dataBM=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6c_bm.csv', sep=',',header=None, names=['x','y'])
        # ax1.plot(dataBM['x'], dataBM['y']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6b_fromEnjalbert2015_4h.csv', sep=',')
        ax1.plot(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        fignum = '_fig6C'

    #ax.plot(x, ybm1, ':', label='BM1')
    #ax.plot(x, ybm2, ':', label='BM2')
    #l = ax1.legend(loc='best', prop={'size':10})
    ll = ax1.legend(loc='best', prop={'size':10})
    ll2 = ax2.legend(loc='best', prop={'size':10})
    plt.setp(ax1.get_xticklabels(), visible=False)
    fig.savefig('%s/enjalbert2015-%s-%s%s.png' % (args.pathout, expcond2, n1, fignum))
    fig.savefig('../../outputs/pqfigures/enjalbert2015-%s-%s%s.png' % (expcond2, n1, fignum), format='png', dpi=1500)

    return

def plotEnjalbert2015_growth2EC(args, model, expcond, expcond2, tit):

    n1, n2 = model.dmodelsKeys
    
    bmn1 = 'biomass_%s' % n1
    bmn2 = 'biomass_%s' % n2

    conctomM = 1000. if model.dmodels[n1].volumes.externalUnits == 'mL' else 1.

    if args.phitransition and args.psitransition:
        tit=tit+', w/ transition'
    else:
        tit=tit+', w/o transition'
    
    x = np.array(model.dmodels[n1].T)
    x2=np.array(model.dmodels[n2].T)
    ybm1 = np.array(model.dmodels[n1].dmetabolites[bmn1].quantity)
    ygl1 = np.array(model.dmodels[n1].dmetabolites['ex_glucose'].concentration)*conctomM
    yac1 = np.array(model.dmodels[n1].dmetabolites['ex_acetate'].concentration)*conctomM
    #yam1 = np.array(model.dmodels[n1].dmetabolites['ex_ammonium'].quantity)
    ybm2 = np.array(model.dmodels[n2].dmetabolites[bmn2].quantity)
    yac2 = np.array(model.dmodels[n2].dmetabolites['ex_acetate'].concentration)*conctomM
    ygl2 = np.array(model.dmodels[n2].dmetabolites['ex_glucose'].concentration)*conctomM
    #yam2 = np.array(model.dmodels[n2].dmetabolites['ex_ammonium'].quantity)

    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1)
    # Top left
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
    ax1.set_ylabel('Q [gDW]')
    ax2.set_ylabel('C [mM]')
    ax2.set_xlabel('Time [hours]')

    ## COLORS http://htmlcolorcodes.com/color-chart/
    ## blues: #2980b9 #3498db  #1abc9c #16a085
    ## yellows: #f4d03f f5b041 eb984e  #dc7633 
    ax1.plot(x, ybm1, '#2980b9', label=bmn1.replace('_', ' '))
    ax1.plot(x, ybm2, '#16a085', label=bmn2.replace('_', ' '))
    #ax1.set_xlim(0., 1.)
    ll = ax1.legend(loc='center right', prop={'size':10})
    ax1.set_title(tit)
    ax2.plot(x, ygl1, '#c0392b', label='Glucose')
    ax2.plot(x, yac1, '#f39c12', label='Acetate')
    # ax2.plot(x, ygl1, '#3498db', label='Glucose 1')
    # ax2.plot(x, yac1, '#1abc9c', label='Acetate 1')
    # ax2.plot(x, ygl2, '#f5b041', label='Glucose 2')
    # ax2.plot(x, yac2, '#eb984e', label='Acetate 2')
    ax1.plot(x, ybm1+ybm2, 'k', label='biomass Tot')

    fignum = ''
    if expcond == "batch_low_Glc":
        # dataAC=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2A_ac.csv', sep=',',header=None, names=['x','y'])
        # dataBM=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2A_bm.csv', sep=',',header=None, names=['x','y'])
        # dataGL=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2A_gl.csv', sep=',',header=None, names=['x','y'])
        # fignum = '_fig2A'
        # ax1.plot(dataBM['x'], dataBM['y']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        # ax2.plot(dataGL['x'], dataGL['y'], 'rs', label='Glucose (exp)')
        # ax2.plot(dataAC['x'], dataAC['y'], '#f1c40f',  linestyle='None', marker='s', label='Acetate (exp)')
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig2_fromEnjalbert2015.csv', sep=',')
        ax1.errorbar(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, yerr=dataFIG['OD SD']*ODtoGDW*volExt, fmt='bs', label='Biomass (exp)')
        ax2.errorbar(dataFIG['Time'], dataFIG['Glucose mM'], yerr=dataFIG['Glucose SD'], fmt='rs', label='Glucose (exp)')
        ax2.errorbar(dataFIG['Time'], dataFIG['Acetate mM'], yerr=dataFIG['Acetate SD'], label='Acetate (exp)', color='#f1c40f',  linestyle='None', marker='s')
        fignum = '_fig2A'
    elif expcond == "fedbatch_low_Ac":
        # dataBM=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6a_bm.csv', sep=',',header=None, names=['x','y'])
        # ax1.plot(dataBM['x'], dataBM['y']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        #ax1.errorbar(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, yerr=dataFIG['OD SD']*ODtoGDW*volExt, fmt='bs', label='Biomass (exp)')
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6a_fromEnjalbert2015.csv', sep=',')
        ax1.plot(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        fignum = '_fig6A'
    elif expcond == "fedbatch_high_Ac":
        # dataBM=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6c_bm.csv', sep=',',header=None, names=['x','y'])
        # ax1.plot(dataBM['x'], dataBM['y']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        dataFIG=pandas.read_csv('../../ecoli/enjalbert2015_data/fig6b_fromEnjalbert2015_4h.csv', sep=',')
        if args.shiftgetime:
            ax1.plot(dataFIG['Time']+0.8, dataFIG['OD 600nm']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        else:
            ax1.plot(dataFIG['Time'], dataFIG['OD 600nm']*ODtoGDW*volExt, 'bs', label='Biomass (exp)')
        fignum = '_fig6C'
    #ax.plot(x, ybm1, ':', label='BM1')
    #ax.plot(x, ybm2, ':', label='BM2')
    #ll = ax1.legend(loc='upper left', prop={'size':10})
    #ll2 = ax2.legend(loc='upper left', prop={'size':10})
    ll = ax1.legend(loc='best', prop={'size':10})
    ll2 = ax2.legend(loc='best', prop={'size':10})
    plt.setp(ax1.get_xticklabels(), visible=False)
    fig.savefig('%s/enjalbert2015-%s-totBM%s.png' % (args.pathout, expcond2, fignum))
    fig.savefig('../../outputs/pqfigures/enjalbert2015-%s-totBM%s.png' % (expcond2, fignum), format='png', dpi=1500)

    #f1.set_title('E. Coli core Aerobic '+expcond.capitalize())
    #ll = f1.legend(loc='center left', prop={'size':10})
    #ll = f1b.legend(loc='center right', prop={'size':10})
    #fig1.savefig(args.pathout+'/enjalbert2015_'+expcond2+fignum+'.png')   
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
    parser.add_argument('-L', '--shiftgetime', help='shift GE time', action='store_true')
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
