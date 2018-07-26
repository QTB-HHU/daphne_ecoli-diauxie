#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  last modified: 2015/11/12     **
#************************************

from __future__ import print_function
from six import iteritems
import numpy as np
from scipy.integrate import ode, odeint
import warnings
import json, cPickle
import random
import pandas
import cobra
from cobra.flux_analysis import parsimonious, moma
from buildToyModel import buildToyModel as btm
import classConvertionFactors

DEBUGVS = False
DEBUGQD = False
DEBUGCS = False
# From http://www.nature.com/nprot/journal/v2/n3/full/nprot.2007.99.html
#  Variables -- The values that are computed by the toolbox are fluxes,
#   which can be best understood as reaction rates. The units for fluxes
#   used throughout this protocol are mmol gDW- 1 h- 1, where gDW is the
#   dry weight of the cell in grams. The biomass functions included in each
#   SBML file are weighted combinations of molecules that are required for
#   cellular growth and reproduction and are scaled such that the units are h-1.
#   Concentrations are expressed in units of mmol gDW-1 ml-1, where the unit of volume is arbitrary.

'''
Begin helper code
'''
class NoSBMLNameException(Exception):
    '''
    NoSBMLNameException is raised when a dynamicModel initialized without a sbmlName attribute
    tries to load a sbml model - the basic toy model is built instead
    '''
class NoBiomassException(Exception):
    '''
    NoBiomassException is raised when a dynamicModel.reactions attribute does not
    contain 'biomass' as key
    Biomass is needed to follow growth!
    '''
class NotOptimalFBASolutionException(Exception):
    '''
    NotOptimalFBASolutionException is raised when the FBA Solution from
    the cobra.optimize method is not optimal
    '''
'''
End helper code
'''

class Volume(object):
    '''
    Very simple volume object, internal and external attributes
    '''
    def __init__(self, bm, external=None, externalUnits=None):
        '''
        Internal volume is biomass volume
        Init with biomass dmetabolite object and the external volume value and unit
        Since fluxes in FBA are in mmol/(gDW*hr), the flux boundaries must also be in the same units
        and therefore all concentrations must be expressed in mmol/gDW
        So default external volume units are gDW
        TODO: converter if units are different
        '''
        self.internal = bm.quantity[0]
        self.bioavailable = bm.quantity[0]
        self.internalUnits = bm.quantityUnits
        self.external = external if external is not None else 1
        self.externalUnits = externalUnits if externalUnits is not None else 'mL'
        self.bm = bm
        self.convertionFactors = classConvertionFactors.ConvertionFactors()
        self.convertToStandardUnits()

    def update(self):
        '''
        Update the internal volume with growth
        '''
        self.internal = max(self.bm.quantity)
        self.bioavailable = self.bm.quantity[-1]
        return

    def convertToStandardUnits(self):
        '''
        Update the volumes to the standard units
        '''
        if self.internalUnits != 'gDW':
            print('WARNING! Something went wrong, biomass is not in gDW')
        # if self.externalUnits != 'gDW':
        #     print 'Converting ', self.external, self.externalUnits, ' in gDW:'
        #     self.external = self.external*self.convertionFactors.AtoB(self.externalUnits, 'gDW')
        #     self.externalUnits = 'gDW'
        #     print '\t', self.external,self.externalUnits
        return


class DynamicModel(object):
    '''
    Model for the dynamic FBA
    Uses COBRApy to solve FBA problems
    '''
    def __init__(self, reactions, metabolites, sbmlName=None, externalV=None, externalVU=None, solver=None, name=None, timeU=None, savePath=None, startT=0.):
        '''
        Init with dictionaries of DReactions and DMetabolites, SBML file (if none, toymodel is built),
        eternal volume with units, the solver to use (cglpk by default) and a name for the dFBA model
        '''
        self.dreactions = reactions
        self.dmetabolites = metabolites
        self.dmetabolitesKeys = self.dmetabolites.keys()
        self.biomassKey = self.getBiomassKey()
        self.dmetabolitesKeys.insert(0, self.dmetabolitesKeys.pop(self.dmetabolitesKeys.index(self.biomassKey)))
        self.dreactionsKeys = self.dreactions.keys()
        self.metCb2DMod = {self.dmetabolites[k].modName: k for k in self.dmetabolitesKeys if self.dmetabolites[k].modName is not None}
        self.rxnCb2DMod = {self.dreactions[k].modName: k for k in self.dreactionsKeys if self.dreactions[k].modName is not None}
        self.smblName = sbmlName if sbmlName is not None else None
        self.volumes = Volume(self.dmetabolites[self.biomassKey], externalV, externalVU) if self.biomassKey is not 'na' else None
        self.T = [startT]
        self.inStatPhase = []
        self.Tunits = timeU if timeU is not None else 'hr'
        self.FBAsolutions = []
        self.FBAsol = None
        self.solver = solver if solver is not None else 'optlang-glpk'
        self.name = name.replace(' ','_') if name is not None else 'a_dynamic_model'
        self.savePath = savePath if savePath is not None else '/tmp/'
        self.retryFlag = 0
        self.statPhaseID = 0
        self.statPhaseTimeI = -1
        self.statPhaseObj = []
        self.objective_sense = ''
        self.thresholdForStatPhase = None
        self.growthRxnId = None
        self.growthMaintenance = None
        self.offDReactions = []
        self.cbModel = None
        self.pFBA = False
        self.moma = False
        self.quitOnFBAError = True
        self.boundsThreshold = 10.
        #self.solver = solver if solver is not None else 'gurobi'
        #.initializeValueStorage()
        
        # Reactions DataFrame supposed to contain one DF per time point
        # Format:
        #      | Upper | Lower | Sol
        # RXN1 |       |       |
        # RXN2 |       |       |
        self.rxnDF = []

    def getBiomassKey(self):
        '''
        Simple method to find the biomass key amongst the DMetabolites dictionary
        TODO: What if more than one biomass key? Should throw an error?
        '''
        bmkeys = filter(lambda x: 'biomass' in x, self.dmetabolitesKeys)
        if len(bmkeys) > 0:
            return bmkeys[0]
        else:
            return 'na'

    def setParsimoniousFBA(self, usePFBA=True):
        if usePFBA:
            self.pFBA = True
        else:
            self.pFBA = False
        return

    def setMOMA(self, useMOMA=True):
        if useMOMA:
            self.moma = True
        else:
            self.moma = False
        return

    def setCbModelMedium(self, mdf):
        '''
        mdf = dataframe with index 'metabolite' and cols ['concentration', 'mM', 'limiting']
        sets the lower bounds of cbmodel exchange reactions according to media composition
        NB everything else is set to 0
        '''
        cm2dm = dict( map( lambda i: (i[1].modName, i[0]), iteritems(self.dreactions) ) )
        for r in self.cbModel.reactions:
            if 'EX_' not in r.id:
                continue
            m = r.id[3:]
            if m in mdf.index:
                ks = '%s in medium composition, setting lower bound as' % m
                if not mdf.loc[m, 'limiting']:
                    ks = ks + ' unlimited: '
                    r.lower_bound = -1000.
                else:
                    # medium concentration [mmol/L] * 0.001 * external volume [L] / biomass [gDW] ==> [mmol/gDW]
                    c = mdf.loc[m, 'concentration']*mdf.loc[m, 'mM']*self.volumes.external*0.001/self.volumes.bioavailable
                    dr = cm2dm.get(r.id, None)
                    if dr is not None:
                        # NB the medium is external ==> volume is bioavailable biomass
                        lb = -1 * self.dreactions[dr].kinF.computeRateFromC(c)
                        ks = '%s constrained by the kinF of DReaction %s: ' % (ks, dr)
                    else:
                        lb = -1 * c/24.
                        ks = ks + ' constrained considering 24h for exhaustion: '
                    r.lower_bound = lb
                print(ks, r.lower_bound)
            else:
                r.lower_bound = 0.
        return
            
        # for m in mdf.metabolite:
        #     # check it exists in cbmodel
        #     if m not in self.cbModel.metabolites:
        #         print(m, ' not in cobra model metabolites, ignoring')
        #         continue
        #     r = 'EX_'+m
        #     if r not in self.cbModel.reactions:
        #         print('ACHTUNG something wrong!')
        #         continue
            
        return
        
    
    def removeFluxes(self, itemi=-1):
        '''
        Update the fluxes of the DReactions
        If DReaction.isODE just append the upperBound
        else use a dictionary of flux solution (vsol) from FBA COBRA
        The flux attribute of a DReaction is a list
        TODO: to keep dimension consistent between time and fluxes, if at the last step of integration no FBA solution is found (vsol is None), just copy the last flux - but this is not correct!
        '''
        for r in self.dreactionsKeys:
            del self.dreactions[r].flux[itemi]
            del self.dreactions[r].ub[itemi]
            del self.dreactions[r].lb[itemi]
        del self.rxnDF[itemi]
        return
    
    def updateFluxes(self, vsol):
        '''
        Update the fluxes of the DReactions
        If DReaction.isODE just append the upperBound
        else use a dictionary of flux solution (vsol) from FBA COBRA
        The flux attribute of a DReaction is a list
        TODO: to keep dimension consistent between time and fluxes, if at the last step of integration no FBA solution is found (vsol is None), just copy the last flux - but this is not correct!
        '''
        for r in self.dreactionsKeys:
            if self.dreactions[r].isODE:
                self.dreactions[r].flux.append(self.dreactions[r].upperBound)
            else:
                if vsol is None:
                    #self.dreactions[r].flux.append(self.dreactions[r].flux[-1])
                    self.dreactions[r].flux.append(0.)
                else:
                    self.dreactions[r].flux.append(vsol[self.dreactions[r].modName])
            self.dreactions[r].appendBounds()
            self.rxnDF[-1].loc[r, 'rate'] = self.dreactions[r].flux[-1]
        return


    def updateQuantitiesDQ(self, dq):
        '''
        Update the DMetabolites quantities with an increment dq (list of quantity variation)
        '''
        for m in self.dmetabolitesKeys:
            self.dmetabolites[m].updateQuantityDQ(dq, self.volumes)
        return

    def updateQuantitiesQ(self, q):
        '''
        Update the DMetabolites quantities with q (list of quantity values)
        '''
        for m in self.dmetabolitesKeys:
            self.dmetabolites[m].updateQuantityQ(dq, self.volumes)
        return

    def initCBModel(self):
        '''
        Setups for cobrapy model
        '''
        # Hardcoding max number of calls to glpk simplex
        if self.solver == 'optlang-glpk':
            self.cbModel.solver.configuration._smcp.it_lim = 50000
        else:
            print('not optlang-glpk')
        return


    def loadSBMLModel(self):
        '''
        Use the COBRApy method to load the SBML model named according the attribute smblName
        If no smblName is defined build a toy model of type 0 (basic)
        TODO: differences in the sbml loading methods?
        TODO: if Tunits is not hr, convert all pre-existing flux boundaries!
        '''
        if self.smblName is not None:
            try:
                #self.cbModel = cobra.io.read_legacy_sbml(self.smblName)
                self.cbModel = cobra.io.read_sbml_model(self.smblName)
                #self.cbModel = cobra.io.sbml.create_cobra_model_from_sbml_file(self.smblName)
            except NoSBMLNameException:
                raise
        else:
            self.cbModel = btm(0)
        self.initCBModel()
        return

    def loadToyModel(self, tt=0, cc=0, qq=0, vv=False):
        '''
        Build a toy model of type tt (see buildToyModel function)
        '''
        self.cbModel = btm(tt, cc, qq, vv)
        self.initCBModel()
        return

    def updateVolume(self):
        '''
        Call the update method of the volume class
        '''
        self.volumes.update()
        return

    def initializeConcentrations(self):
        '''
        Initialize the DMetabolites concentrations by updating the volume and setting concentrations and concentration units
        '''
        self.updateVolume()
        for m in self.dmetabolitesKeys:
            self.dmetabolites[m].updateConcentration(self.volumes)
            self.dmetabolites[m].setConcentrationUnits(self.volumes)
        return

    def updateConcentrations(self):
        '''
        Update the DMetabolites concentrations after updating the volume
        '''
        self.updateVolume()
        for m in self.dmetabolitesKeys:
            self.dmetabolites[m].updateConcentration(self.volumes)
        return

    def constrainSources(self, element='C', count=0, uptakeOff=False, exportOff=False):
        '''
        Set the upper (exportOff==True) and/or lower (uptakeOff==True) bounds to 0
        for exchange reactions involving at least *count*+1 atoms of *element*
        If both uptakeOff and exportOff are False, this function will return the list
        of all the sources of a certain element and its contraints
        '''
        exrxns = [r for r in self.cbModel.reactions if 'EX_' in r.id and r.reactants[0].elements.get(element, 0) > count]
        if uptakeOff or exportOff:
            for r in exrxns:
                if uptakeOff:
                    r.lower_bound = 0.
                if exportOff:
                    r.upper_bound = 0.
            
        return exrxns

    
    def resetObjective(self, hardReset=True):
        '''
        Set all coefficients of the cobra model objective vector to 0
        '''
        self.objective_dict = {}
        if hardReset:
            for r in self.cbModel.reactions:
                r.objective_coefficient = 0.
            self.objective_ids = []
        else:
            for rid in self.objective_ids:
                print('Setting reaction %s to 0 from %.2f coeff' % (rid))
                self.cbModel.reactions.get_by_id(rid[0]).objective_coefficient = 0.
        return

    def changeTemporaryObjective(self, objectiveList):
        '''
        change the cbModel objective but not the object original attributes
        '''
        self.changeObjective(objectiveList, True)

    def changeObjective(self, objectiveList, temporary=False):
        '''
        change the cbModel objective but not the object original attributes
        NB need to have always the same reactions in list!!
        '''
        if self.pFBA and len(objectiveList) > 1:
            print('NEW COBRApy function pfba supports models with multiple objective function, not switching off pFBA!')
            #self.setParsimoniousFBA(False)
        previousRxns = self.objective_dict.keys()
        self.resetObjective(hardReset=(not temporary))
        for objtpl in objectiveList:
            reactionName = objtpl[0]
            objectiveCoefficient = objtpl[1]
            objective_sense='maximize'
            if len(objtpl) > 2:
                objective_sense = objtpl[2]
            self.setObjective(reactionName, objectiveCoefficient, objective_sense, temporary)
        if sorted(self.objective_dict.keys()) != sorted(previousRxns):
            print('WARNING: the objective did not change all the reactions coefficients!')
        return


    def setStatPhaseObjective(self, reactionName, objectiveCoefficient, objective_sense='maximize'):
        '''
        Set the objective in the cobra model
        More reactions can be added as objectives
        '''
        if self.objective_sense and self.objective_sense != objective_sense:
            print('Adding a stat phase objective function with different sense, inverting objective coefficient')
            objectiveCoefficient = -1*objectiveCoefficient
        self.statPhaseObj.append((reactionName, objectiveCoefficient))
        return

    def setObjective(self, reactionName, objectiveCoefficient, objective_sense='maximize', temporary=False):
        '''
        Set the objective in the cobra model
        More reactions can be added as objectives
        '''
        if self.objective_sense and self.objective_sense != objective_sense:
            print('Adding an objective function with different sense, inverting objective coefficient')
            objectiveCoefficient = -1*objectiveCoefficient
        else:
            self.objective_sense = objective_sense
        self.objective_dict[self.cbModel.reactions.get_by_id(reactionName)] = objectiveCoefficient
        #self.cbModel.reactions.get_by_id(self.dreactions[reactionName].modName).objective_coefficient = objectiveCoefficient
        self.cbModel.reactions.get_by_id(reactionName).objective_coefficient = objectiveCoefficient
        if temporary:
            return
        self.objective_ids.append((reactionName, objectiveCoefficient))
        if self.thresholdForStatPhase is None:
            self.thresholdForStatPhase = (reactionName, 0.001)
        return

    def growthLimitationMap(self, search=True, save=True, plot=True, nsteps=20):
        '''
        Fill rxns list with EX reactions importing C, N, P, S and with -1000 < lb < 0 
        rxnsid dict 'RxnID : (#C, #N, #P, #S)'
        '''
        rxns = [r for r in self.cbModel.reactions if 'EX_' in r.id and
                r.lower_bound > -1000. and r.lower_bound < 0. and
                (r.reactants[0].elements.get('C', 0) > 0 or r.reactants[0].elements.get('N', 0) > 0 or r.reactants[0].elements.get('P', 0) > 0 or r.reactants[0].elements.get('S', 0) > 0)]
        rxnCNPS = {}
        dmetsCScaled0 = {}
        rids = []
        rxnKinF = {}
        for r in rxns:
            dmetname = self.metCb2DMod.get( r.reactants[0].id, '' )
            drxnname = self.rxnCb2DMod.get( r.id, '' )
            if len(dmetname) > 0 and len(drxnname) > 0:
                rids.append(r.id)
                rxnCNPS[r.id] = (r.reactants[0].elements.get('C', 0), r.reactants[0].elements.get('N', 0), r.reactants[0].elements.get('P', 0), r.reactants[0].elements.get('S', 0) )
                dmetsCScaled0[r.id] = self.dmetabolites.get( dmetname ).concentrationScaled[0]
                rxnKinF[r.id] = self.dreactions.get( drxnname ).kinF
        if len(rids) > 2:
            print('ACHTUNG: 3d maps not yet implemented!')
            rids = rids[:2]
        try:
            glmap = json.load(open('map-%s-%s.json' % (self.name, '-'.join(rxnsid)), 'r'))
        except:
            glmap = self.buildGrowthLimMap(rids, rxnCNPS, dmetsCScaled0, rxnKinF, plot, nsteps)
        return glmap

    def buildGrowthLimMap(self, rids, CNPS, c0, kin, draw=True, xysize=20):
        '''
        rids = list of cbModel reactions
        CNPS = dict of rids: (#C, #N, #P, #S)
        c0   = dict of rids: initial conc scaled of substrate
        kin  = dict of rids: kinetic function
        Compute lower bounds for reactions in interval (0, c0, 100)
        '''
        if self.growthRxnId is None:
            print('ERROR no growthRxnId attribute')
            return
        rX = rids[0]
        rY = rids[1]
        #cX = np.linspace(0, min(c0[rX], kin[rX].mmConstant+abs(c0[rX] - kin[rX].mmConstant)/2), xysize)
        #cY = np.linspace(0, min(c0[rY], kin[rY].mmConstant+abs(c0[rY] - kin[rY].mmConstant)/2), xysize)
        cX = np.linspace(0, min(c0[rX], kin[rX].mmConstant), xysize)
        cY = np.linspace(0, min(c0[rY], kin[rY].mmConstant), xysize)
        vX = kin[rX].computeRateFromC(cX)
        vY = kin[rY].computeRateFromC(cY)
        # This is the correct ordering of the 2D arrays to then fill Z(i,j)
        X = np.tile(np.array([vX]).transpose(), (1, xysize))
        Y = np.tile(vY, (xysize,1))
        ##X = np.tile(np.array(Xax), (xsize,1))
        ##Y = np.tile(np.array([np.array(Yax)]).transpose(), (1,ysize))
        Z = np.zeros((xysize, xysize))
        # set cbModel objective to max biomass
        for cbr, oc in iteritems(self.objective_dict):
            if cbr.id == self.growthRxnId:
                cbr.objective_coefficient = 1.
            else:
                cbr.objective_coefficient = 0
        # compute max achievable growth setting to c0s
        self.cbModel.reactions.get_by_id(rX).lower_bound = -1*kin[rX].computeRateFromC(c0[rX])
        self.cbModel.reactions.get_by_id(rY).lower_bound = -1*kin[rY].computeRateFromC(c0[rY])
        msol = parsimonious.pfba(self.cbModel) #, objective=self.objective_dict, solver=self.solver)
        maxvbm = msol.x_dict[self.growthRxnId]
        
        for i, x in enumerate(-1*vX):
            #self.cbModel.reactions.get_by_id(rX).upper_bound = x
            self.cbModel.reactions.get_by_id(rX).lower_bound = x
            for j, y in enumerate(-1*vY):
                #self.cbModel.reactions.get_by_id(rY).upper_bound = y
                self.cbModel.reactions.get_by_id(rY).lower_bound = y
                # solve pFBA
                try:
                    sol = parsimonious.pfba(self.cbModel) #, objective=self.objective_dict, solver=self.solver)
                    s = sol.x_dict[self.growthRxnId]
                except:
                    print('In buildGrowthLimMap, ',i,'-',j,'th entry: LP Solver called by pfba() raised an exception')
                    s = 0.
                Z[i,j] = s/maxvbm

        # reset cbModel objective 
        for cbr, oc in iteritems(self.objective_dict):
            cbr.objective_coefficient = oc

        if draw:
            import matplotlib.pyplot as plt
            from matplotlib import cm, gridspec
            from mpl_toolkits.mplot3d import Axes3D
            from matplotlib.ticker import LinearLocator, FormatStrFormatter
            fig0 = plt.figure()
            ax = fig0.add_subplot(111, projection='3d')
            surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm)
            ax.set_xlabel(rX)
            ax.set_ylabel(rY)
            fig0.savefig('./3d-map-%s-%s.png' % (self.name, '-'.join(rids)))

            fig = plt.figure()
            gs = gridspec.GridSpec(2, 3, width_ratios=[1,3,0.3], height_ratios=[1,3])
            ax1 = fig.add_subplot(gs[1,1])
            ax2 = fig.add_subplot(gs[0,1], sharex=ax1)
            ax3 = fig.add_subplot(gs[1,0], sharey=ax1)
            ax1b = fig.add_subplot(gs[1,2])
            im = ax1.pcolor(X, Y, Z, cmap=cm.jet)
            ax1.set_title('Relative max growth rate')
            ax1.set_xlabel('-1*L.B. %s (mmol/gDw/hr)' % (rX))
            ax3.set_ylabel('-1*L.B. %s (mmol/gDw/hr)' % (rY))
            ax2.plot(vX, cX)
            ax3.plot(cY, vY)
            ax2.set_title(rX)
            ax3.set_title(rY)
            ax2.set_ylabel('Conc %s' % (self.cbModel.reactions.get_by_id(rX).reactants[0]).id)
            ax3.set_xlabel('Conc %s' % (self.cbModel.reactions.get_by_id(rY).reactants[0]).id)
            plt.colorbar(im, cax=ax1b)
            plt.setp(ax1.get_yticklabels(), visible=False)
            plt.setp(ax2.get_xticklabels(), visible=False)
            fig.savefig('./growth-limiting-map-%s-%s.png' % (self.name, '-'.join(rids)))

            fig = plt.figure()
            gs = gridspec.GridSpec(2, 3, width_ratios=[1,3,0.3], height_ratios=[1,3])
            ax1 = fig.add_subplot(gs[1,1])
            ax2 = fig.add_subplot(gs[0,1], sharex=ax1)
            ax3 = fig.add_subplot(gs[1,0], sharey=ax1)
            ax1b = fig.add_subplot(gs[1,2])
            im = ax1.pcolor(cX, cY, Z, cmap=cm.jet)
            ax1.set_title('Relative max growth rate')
            ax1.set_xlabel('Conc %s (mmol/gDw)' % (self.cbModel.reactions.get_by_id(rX).reactants[0]).id)
            ax3.set_ylabel('Conc %s (mmol/gDw)' % (self.cbModel.reactions.get_by_id(rY).reactants[0]).id)
            ax2.plot(cX, vX)
            ax3.plot(vY, cY)
            ax2.set_title(rX)
            ax3.set_title(rY)
            ax2.set_ylabel('-1*L.B.')#'L.B. (mmol/gDw/hr)')
            ax3.set_xlabel('-1*L.B.')#'L.B. (mmol/gDw/hr)')
            plt.colorbar(im, cax=ax1b)
            plt.setp(ax1.get_yticklabels(), visible=False)
            plt.setp(ax2.get_xticklabels(), visible=False)
            fig.savefig('./growth-limiting-map-%s-%s-%s.png' % (self.name, (self.cbModel.reactions.get_by_id(rX).reactants[0]).id, (self.cbModel.reactions.get_by_id(rY).reactants[0]).id))
        
        return Z

    def fillRxnDF(self):
        '''
        Append the current time reaction Data Frame
        '''
        tmpdf = pandas.DataFrame(index=self.dreactionsKeys)
        for r in self.dreactionsKeys:
            self.dreactions[r].updateBounds(time=self.T)
            tmpdf.loc[r, 'upperBound'] = self.dreactions[r].upperBound
            tmpdf.loc[r, 'lowerBound'] = self.dreactions[r].lowerBound
            if not self.dreactions[r].isODE:
                self.cbModel.reactions.get_by_id(self.dreactions[r].modName).upper_bound = tmpdf.loc[r, 'upperBound']
                self.cbModel.reactions.get_by_id(self.dreactions[r].modName).lower_bound = tmpdf.loc[r, 'lowerBound']
        
        self.rxnDF.append(tmpdf)
        return

    
    def setBounds(self):
        '''
        Loop over the DReactions and updates the bounds of the dFBA and the cobra models
        The updateBounds() method of DReaction objects calls the computeBounds() method of KineticLaw objects (user-defined) and copies the lowerBound, upperBound attributes
            The computeBounds() method of KineticLaw object updates the lowerBound, upperBound attributes
        '''
        for r in self.dreactionsKeys:
            if self.dreactions[r].isODE:
                self.dreactions[r].updateBounds(time=self.T)
            else:
                self.dreactions[r].updateBounds(time=self.T)
                self.cbModel.reactions.get_by_id(self.dreactions[r].modName).upper_bound = self.dreactions[r].upperBound
                self.cbModel.reactions.get_by_id(self.dreactions[r].modName).lower_bound = self.dreactions[r].lowerBound
        return
    
    def setCbModelBounds(self, tindex):
        '''
        Loop over the DReactions and updates the bounds of the dFBA and the cobra models
        The updateBounds() method of DReaction objects calls the computeBounds() method of KineticLaw objects (user-defined) and copies the lowerBound, upperBound attributes
            The computeBounds() method of KineticLaw object updates the lowerBound, upperBound attributes
        '''
        for r in self.dreactionsKeys:
            if self.dreactions[r].isODE:
                continue
            else:
                self.cbModel.reactions.get_by_id(self.dreactions[r].modName).upper_bound = self.dreactions[r].ub[tindex]
                self.cbModel.reactions.get_by_id(self.dreactions[r].modName).lower_bound = self.dreactions[r].lb[tindex]
        return

    def noMetabolitesExhausted(self, qThr=1e-8):
        '''
        Checks if there are DMetabolites with quantity below a threshold
        '''
        return len(filter(lambda x: self.dmetabolites[x].quantity[-1] < qThr, self.dmetabolitesKeys)) < 1

    def checkForNegativeSolutions(self, qsol, qThr=0):
        '''
        Checks if some DMetabolite quantity went below zero
        Return boolean and the negative DMetabolite (None in case all positive)
        '''
        minQ = min(qsol)
        isNeg = minQ < qThr
        if isNeg:
            metNeg = np.where(qsol==minQ)[0][0]
            #metNeg = self.dmetabolitesKeys[np.where(qsol==minQ)[0][0]]
            return isNeg, metNeg
        else:
            return isNeg, None

    def exitSimulation(self, allGood):
        print('End of dFBA simulation\n\t stopDFBA message:\t%s\n\t stopFBA message:\t%s\n' % (self.stopDFBA[1], self.stopFBA[1]))
        ofile = open(self.savePath+'/endOfSimulation-'+self.name+'.p', 'w')
        cPickle.dump(self, ofile)
        ofile.close()
        if not allGood:
            print('Something not good!')
            exit(1)
        return

    def setQuitOnFBAError(self, quitOnFBAError):
        self.quitOnFBAError = quitOnFBAError
        return

    def setBoundsThreshold(self, thr):
        self.boundsThreshold = thr

    def returnFBASolOnError(self, message):
        '''
        Manage what is returned if FBA is not solvable
        '''
        if self.quitOnFBAError:
            self.stopDFBA = (True, ('%s at T[%d] = %.2f' % (message, len(self.T)-1, self.T[-1]) ))
            print(self.stopDFBA[1])
            return None
        else:
            self.stopFBA = (True, ('Continue without FBA: %s at T[%d] = %.2f' % (message, len(self.T)-1, self.T[-1]) ), len(self.T)-1 )
            print(self.stopFBA[1])
            # COBRApy returns now pandas Series, not anymore a dict!
            #return dict( zip( self.cobraReactions,  len(self.cobraReactions)*[0.0] ) )
            return pandas.Series( len(self.cobraReactions)*[0.0], index=self.cobraReactions )
        return None
        
    def testDeltaBounds(self):
        dbs = np.array( map(lambda r: ( (r.lb[self.stopFBA[2]] - r.lb[-1])**2 , (r.ub[self.stopFBA[2]] - r.ub[-1])**2 ) if not r.isODE else (0.,0.), self.dreactions.values() ) )
        #print('Delta B @ ', len(self.T)-1, ': ', np.sqrt(dbs.sum()) )
        return  np.sqrt(dbs.sum()) > self.boundsThreshold

    def setGrowthMaintenance(self, rxnid, gthr):
        self.growthRxnId = rxnid
        self.growthMaintenance = gthr
        return

    def setMinimalMaintenance(self, rxnid, vthr):
        '''
        Find what is the minimal flux for FBA feasible solutions in stationary phase
        '''
        self.thresholdForStatPhase = (rxnid, vthr)
        return

    def relaxConstraints(self):
        '''
        Look for reactions with ub == lb and change it to either (0, ub) or (lb, 0), depending on sign
        '''
        retry = False
        rxnsToChange = filter(lambda x: (not self.dreactions[x].isODE) and abs(self.dreactions[x].lowerBound - self.dreactions[x].upperBound) < 0.00001, self.dreactionsKeys)
        # print rxnsToChange
        if len(rxnsToChange) < 1:
            return retry
        for rid in rxnsToChange:
            if abs(self.dreactions[rid].lowerBound) < 1e-10:
                continue
            # print rid, self.dreactions[rid].lowerBound, self.dreactions[rid].upperBound
            if self.dreactions[rid].lowerBound > 0:
                self.dreactions[rid].lowerBound = 0.
                self.cbModel.reactions.get_by_id(self.dreactions[rid].modName).lower_bound = 0.
                self.rxnDF[-1].loc[rid, 'lowerBound'] = 0.
            else:
                self.dreactions[rid].upperBound = 0.
                self.cbModel.reactions.get_by_id(self.dreactions[rid].modName).upper_bound = 0.
                self.rxnDF[-1].loc[rid, 'upperBound'] = 0.
            # print rid, self.dreactions[rid].lowerBound, self.dreactions[rid].upperBound
            retry = True
        return retry

    def initFBA(self):
        '''
        First FBA run to initialize
           Prevent T=0 discontinuities! 
        First call to initialize the reaction rates:
           run FBA, fill constraints, re-run, remove initial entries and update constraints
        '''
        if self.cbModel is None:
            self.FBAsolutions.append(self.solveFBA(init=True))
            self.updateFluxes(self.FBAsolutions[-1])
            return
        self.cobraReactions = map(lambda x: x.id, self.cbModel.reactions )
        self.FBAsolutions.append(self.solveFBA(init=True))
        #XXXX +self.odefluxes[-1]
        self.updateFluxes(self.FBAsolutions[-1])
        if len(self.FBAsolutions) > 1:
            print('WARNING! Initializing solveFBA from non empty solution list!')
        self.FBAsolutions[0] = self.solveFBA(init=True)
        self.removeFluxes()
        self.updateFluxes(self.FBAsolutions[-1])
        if len(self.FBAsolutions) > 1:
            print('WARNING! Initializing solveFBA returns too many entries!')


    def solveFBA(self, firstcall=True, init=False):
        '''
        Solve the FBA problem:
        1) update boundary conditions using the setBounds method
        2) solve FBA using the COBRApy optimize method
        3) return the solution fluxes and the status
        TODO: lp parameters, setBounds parameters etc
        http://cobrapy.readthedocs.io/en/latest/_modules/cobra/flux_analysis/parsimonious.html
        '''
        if self.stopFBA[0]:
            # the model is "dead"
            self.setBounds()
            self.fillRxnDF()
            resumeFBA = self.testDeltaBounds()
            if not resumeFBA:
                return self.FBAsolutions[-1]
            else:
                self.stopFBA = (False, 'Resume FBA at T[%d] = %.2f' % (len(self.T)-1, self.T[-1]), -1)
                print(self.stopFBA[1])
        elif firstcall:
            self.setBounds()
            self.fillRxnDF()
        if self.cbModel is None:
            return
        self.FBAsolStatus = 'error'
        self.inStatPhase.append(self.statPhaseID > 0)
        if self.statPhaseID > 0:
            print('Model %s in StatPhase %d; TO CHANGE: evolve objective with increased stay in stat phase' % (self.name, self.statPhaseID))
            if self.statPhaseID < 2:
                self.changeTemporaryObjective(self.statPhaseObj)
            #print('TO CHANGE: no different solver, but different objective => YOU SHOULD NOT BE HERE')
            #return self.solveStatPhaseFBA()
        if self.moma and firstcall:
            if not init:
                # Switch off pFBA, not needed
                #self.pFBA = False
                # Obtain a solution that is minimally different from previous one (contained in self.FBAsol)
                # N.B. self.cbModel has to stay with usual constraints: copy, modify, copy, delete
                self.cbmCopy = self.cbModel.copy()
                self.cbModel = moma.add_moma(self.cbmCopy, solution=self.FBAsol, linear=True, runcopy=True)
        if self.pFBA:
            try:
                self.FBAsol = parsimonious.pfba(self.cbModel) #, objective=self.objective_dict, solver=self.solver)#, fraction_of_optimum = 0.95)
            except:
                if firstcall:
                    print('LP Solver called by pfba() raised an exception at T[%d] = %.2f; not quitting yet' % (len(self.T)-1, self.T[-1]))
                    self.FBAsol = None
                else:
                    return self.returnFBASolOnError('LP Solver called by pfba() raised an exception')
        else:
            try:
                self.FBAsol = self.cbModel.optimize(objective_sense = self.objective_sense, solver=self.solver)
            except:
                return self.returnFBASolOnError('LP Solver called by optimize() raised an exception')
        #optimize_minimal_flux(self.cbModel, objective_sense = self.objective_sense, solver=self.solver)
        if self.FBAsol is not None:
            self.FBAsolStatus = self.FBAsol.status
        if self.FBAsolStatus is not 'optimal':
            print('No FBA Solution: status ', self.FBAsolStatus)
            tryAgain = self.relaxConstraints()
            if tryAgain:
                print('Retry after relaxing constraints')
                tmpsol = self.solveFBA(False)
                return tmpsol
            else:
                return self.returnFBASolOnError('FBA %s' % (self.FBAsolStatus))
        if self.moma:
            if 'moma_old_objective' in self.cbModel.solver.variables:
                self.cbModel = self.cbmCopy.copy()
                del self.cbmCopy
        return self.FBAsol.x_dict

    def solveStatPhaseFBA(self, firstcall=True):
        '''
        The organism will enter a stationary phase (parsimonious FBA, switching off multiple objectives)
        After the run a random algorithm will decide if the organism will exit the stationary phase
        The probability of exiting should increase with the lenght of the stationary phase
        '''
        if self.statPhaseID == 1:
            if len(self.statPhaseObj) > 0:
                # TODO!!
                print('YOU SHOULD NOT BE HERE!')
                self.changeTemporaryObjective(self.statPhaseObj)
            # self.fop = 0.1
            # if self.growthRxnId:
            #     self.changeTemporaryObjective([(self.growthRxnId, 1)])
            #     if self.growthMaintenance and abs(self.FBAsolutions[-1][self.growthRxnId]) > 0.:
            #         self.fop = self.growthMaintenance/abs(self.FBAsolutions[-1][self.growthRxnId])
            # else:
            #     self.changeTemporaryObjective([self.objective_ids[0]])
        # if len(self.T) > 5650:
        #     ofile1 = open('/tmp/stopDFBAcobra.p', 'w')
        #     ofile2 = open('/tmp/stopDFBAdynamicmodel.p', 'w')
        #     cPickle.dump(self.cbModel, ofile1)
        #     cPickle.dump(self, ofile2)
        #     ofile1.close()
        #     ofile2.close()
        self.FBAsol = self.cbModel.optimize(objective_sense = self.objective_sense, solver=self.solver)
        if self.FBAsol.status is not 'optimal':
            print('In StatPhaseFBA, FBA not feasible (', self.FBAsol.status,'), relaxing constraints')
            tryAgain = self.relaxConstraints()
            if tryAgain:
                print('Retry after relaxing constraints, iteration: ',len(self.T))
                tmpsol = self.solveStatPhaseFBA(False)
                return tmpsol
            else:
                return self.returnFBASolOnError('FBA %s' % (self.FBAsol.status))
        self.FBAsol = parsimonious.pfba(self.cbModel) #, objective=self.objective_dict, solver=self.solver)#, fraction_of_optimum = self.fop)
        if self.FBAsol.status is not 'optimal':
            return self.returnFBASolOnError('No StatPhaseFBA Solution; pFBA %s' % (self.FBAsol.status))
        return self.FBAsol.x_dict

    def updateStatPhase(self):
        '''
        random.random() returns a random float 0.0 <= x < 1.0
        if x > 0.5/self.statPhaseID reset StatPhase
        '''
        if self.statPhaseID > 1:
            if random.random() > 1: # - 0.01*self.statPhaseID:
                print('Randomly exit StatPhase')
                self.changeTemporaryObjective(self.objective_ids)
                self.statPhaseID = 0
                return
        if self.statPhaseID > 5000:
            print(self.name, ' exit stationary phase')
            self.changeTemporaryObjective(self.objective_ids)
            self.statPhaseID = 0
            return

        rid = self.thresholdForStatPhase[0]
        rfl = self.thresholdForStatPhase[1]
        if abs(self.FBAsolutions[-1][rid]) < rfl:
            if self.statPhaseID < 1:
                self.statPhaseTimeI = len(self.T)-1
                print('at Time ', self.T[-1] ,' enter stationary phase, flux ',  rid, ' is ',self.FBAsolutions[-1][rid], ' < ', rfl)
            elif DEBUGVS:
                print('at Time ', self.T[-1] ,' stay in stationary phase, ID: ',self.statPhaseID)
            self.statPhaseID += 1
        else:
            if self.statPhaseID > 0:
                print(self.name, ' resetting stationary phase id, ', rid, ' was ',self.FBAsolutions[-1][rid])
            self.statPhaseID = 0
        return

    def getInitialConditions(self):
        '''
        Return last quantity and last time values to use in the next FBA step as initial conditions
        '''
        q0 = map(lambda x: self.dmetabolites[x].quantity[-1], self.dmetabolitesKeys)
        t0 = self.T[-1]
        return q0, t0

    def switchOffDReaction(self, metIdx):
        '''
        If a metabolite is exhausted, switch off the reaction
        '''
        m = self.dmetabolitesKeys[metIdx]
        for r in self.dmetabolites[m].metReactions.keys():
            if r in self.dreactionsKeys:
                if self.dreactions[r].flux[-1] < 0:
                    return
        self.offDReactions.append()
        return

    def differentialFBA(self, y, verbose=False):
        '''
        Define the differential problem with dQdt ODEs
        The rate can also depend on fixed components not in the dModel --> either take the flux or 1.
        1) call to solveFBA (set ub,lb and solve FBA)
             if no solution terminate loop and set the stopDFBA attribute ('No FBA Solution' message), save the cbModel (cPickle.dump) and return dQdt
             else
             TODO: change objective rendomly?!
        2) append the solution to the FBAsolutions attribute
             since differentialFBA is called multiple times before solution is used, solutions are store in self.FBAsolutions and fluxes update only once
        3) define ODEs for metabolites v = dqdt = sum_i f_i*c_i*q(m_i)
            loop over the DMetabolites and retrieve the flux solution f_i for the DReactions involved in dQ
            the corresponding metReactions attribute of DMetabolite gives c_i and m_i and the quantity of m_i is retrieved
        4) save dQdt as self.dQdt and return dQdt
        TODO: dirty way to exit in case no FBA solution... do smth smarter
        '''
        sizeQ = len(self.dmetabolitesKeys)
        dQdt = np.zeros((sizeQ,1))

        for i in range(sizeQ):
            m = self.dmetabolitesKeys[i]
            dQdt[i] = 0
            for r in self.dmetabolites[m].metReactions.keys():
                if DEBUGVS:
                    print('DEBUGVS r ', r, 'm ', m, ' at time step ', len(self.T))
                    if r in self.dreactionsKeys:
                        if not self.dreactions[r].isODE:
                            if self.FBAsolutions[-1] is None:
                                print('\t FBAsolution: None')
                            else:
                                print('\t FBAsolution: ', self.FBAsolutions[-1][self.dreactions[r].modName])
                        print('\t DReaction V: ', self.dreactions[r].flux[-1], ' UB: ', self.dreactions[r].upperBound, ' LB: ', self.dreactions[r].lowerBound)
                    else:
                        print('\t 1.')
                rate = self.dreactions[r].flux[-1] if r in self.dreactionsKeys else 1.
                metQ = 1
                first = True
                #AS 2015/12/02 example of possible need of sum
                #dGlc/dt = v_ex*(Biomass - smth)
                #{'glucose_exchange': [(1, 'biomass'), (-1, smth)]}
                for n in self.dmetabolites[m].metReactions[r]:
                    mq = 1.
                    for ni in range(1, len(n)):
                        if n[ni] is not None:
                            mq = mq*y[self.dmetabolitesKeys.index(n[ni])]
                    if first:
                        metQ = metQ*n[0]*mq
                    else:
                        metQ += n[0]*mq
                    first = False
                dQdt[i] += rate*metQ
            if verbose:
                print(m, ' ODE is: ')
                for r in self.dmetabolites[m].metReactions.keys():
                    print(self.dmetabolites[m].metReactions[r], '*v(',r,')')
        if verbose:
            print('dQdt:', dQdt)
        self.dQdt = dQdt
        return dQdt

    def runDynamicFBA(self, maxT=10, integrator='vode', stepChoiceLevel=(0, 0., 10., 1000), verbose=False, silent=True, restarting=False, quitOnFBAError=True):
        '''
        Run the dynamic FBA as an integration problem using ode imported from scipy
        choose a different automatization of integration method via the stepChoiceLevel flag
        stepChoiceLevel=(type, min T step, max T step, max N steps)
        TODO: for now only 0 choice, should look more into it
        stepChoiceLevel==0
            1) set the solver using ode (scipy) and the integrator (defaulted to 'vode')
            2) 'dopri5' integrator has nsteps setting
            3) start integration loop while self.stopDFBA is False, time is lower than maxT and step counter is lower than 10000. (anti-infinity condition)
            4) set initial conditions and integrate
            5) use the self.checkForNegativeSolutions method to check for negative DMetabolite quantities
                if so, go back to find the zero crossing point and restart from there
            6) update time, DMetabolites and DReaction fluxes using the last self.FBAsolutions[-1]
            7) stopDFBA, save the cbModel (cPickle.dump) and return
        Links of reference:
            http://stackoverflow.com/questions/12926393/using-adaptive-step-sizes-with-scipy-integrate-ode
            http://vlsd.blogspot.de/2012/11/adaptive-time-step-integrators-in-scipy.html
        self.stopDFBA = (bool, int) the integer corresponds to an exit flag
        TODO: better event handler? see http://scicomp.stackexchange.com/questions/16325/dynamically-ending-ode-integration-in-scipy
        '''
        def dqdt(t, y, mod, verbose=False):
            '''
            Return the dQdt system of ODEs defined from the flux solutions
            '''
            return mod.differentialFBA(y, verbose)

        smethod = 'bdf'
        smax_step = 1
        self.stopDFBA = (False, 'dFBA simulation OK')
        self.stopFBA = (False, 'FBA problem OK', -1)
        self.quitOnFBAError = quitOnFBAError

        if verbose:
            print('\n start runDynamicFBA \n')

        integratorSet = False
        
        if stepChoiceLevel[0] == 0:
            nMaxSteps = stepChoiceLevel[3]
            solver = ode(dqdt).set_integrator(integrator)
            if integrator == 'dopri5':
                solver.set_integrator(integrator, nsteps=1, max_step=stepChoiceLevel[2])
                integratorSet = True
            elif integrator in 'vode':
                solver.set_integrator(integrator, min_step=stepChoiceLevel[1], max_step=stepChoiceLevel[2], method = 'bdf', order = 5)
                integratorSet = True
            else:
                # 'lsoda'
                print('Option 0 is available only for vode and dopri5 SciPy solvers!')
                return
        elif stepChoiceLevel[0] == 1:
            nMaxSteps = stepChoiceLevel[3] - 1
            grid_t = np.linspace(stepChoiceLevel[1], stepChoiceLevel[2], stepChoiceLevel[3])
            grid_dt = grid_t[1] - grid_t[0]
            if integrator == 'dopri5':
                solver = ode(dqdt).set_integrator(integrator, nsteps=1, max_step=grid_dt)
                integratorSet =True
            else:
                # 'lsoda'
                print('Option 1 is available only for dopri5 SciPy solvers!')
                return

        if integratorSet:
            #set the parameters of the differential function dqdt: model and verbosity
            solver.set_f_params(self, verbose)
            # suppress Fortran-printed warning
            solver._integrator.iwork[2] = -1
            warnings.filterwarnings("ignore", category=UserWarning)

            q0, t0 = self.getInitialConditions()
            solver.set_initial_value(q0, t0)
            if verbose:
                print('Initial Conditions to start integration: ')
                print('\t t0: ', t0)
                print('\t q0: ', q0)

            #print len(self.FBAsolutions)
            #solveFBA() calls setBounds() which distinguishes ODE and COBRA reactions
            if not restarting:
                self.initFBA()

            step = 0
            substeps = 0
            eventIdx = 1

            while not self.stopDFBA[0] and self.T[-1]<maxT and step < nMaxSteps:
                step += 1
                if self.cbModel is not None:
                    self.updateStatPhase()
                if stepChoiceLevel[0] == 1:
                    solver.integrate(solver.t + grid_dt)
                else:
                    solver.integrate(maxT, step=True)
                if verbose:
                    print('solver time: ', solver.t, '; step : ',step)
                    print('solver solution: ')
                    for i in range(len(solver.y)):
                        print('\t', self.dmetabolitesKeys[i],': ',solver.y[i])
                elif not silent:
                    if step%10==0:
                        print('Step: ', step, ', Substep: ', substeps,', Time:', solver.t)
                handleZeroEvent, metIdx = self.checkForNegativeSolutions(solver.y)
                if handleZeroEvent:
                    if self.dmetabolitesKeys[metIdx] == self.biomassKey:
                        print('in handleZeroEvent for ', self.dmetabolitesKeys[metIdx], ', Exit ')
                        self.stopDFBA = (True, 'Negative Biomass at time %f' % self.T[-1])
                        continue
                    if verbose:
                        print('****************************\n!in handleZeroEvent!')
                    substeps += 1
                    #slope = (min(solver.y) - self.dmetabolites[metId].quantity[-1])/(solver.t - self.T[-1])
                    slope = float(self.dQdt[metIdx])
                    if slope == 0:
                        self.stopDFBA = (True, 'Zero slope')
                        continue
                    zeroDT = 0.-self.dmetabolites[self.dmetabolitesKeys[metIdx]].quantity[-1]/slope
                    if zeroDT < 0:
                        print('ACHTUNG!!! ',self.dmetabolitesKeys[metIdx], ' previous Q: ', self.dmetabolites[self.dmetabolitesKeys[metIdx]].quantity[-1] ,' slope ', slope, ' negative delta time!!!!', zeroDT)
                    for i in range(len(self.dmetabolitesKeys)):
                        dq = float(self.dQdt[i])*zeroDT
                        if i == metIdx:
                            self.dmetabolites[self.dmetabolitesKeys[metIdx]].updateQuantityQ(0., self.volumes)
                        else:
                            self.dmetabolites[self.dmetabolitesKeys[i]].updateQuantityDQ(dq, self.volumes)
                    self.T.append(self.T[-1]+zeroDT)
                    q0, t0 = self.getInitialConditions()
                    solver.set_initial_value(q0, t0)
                    #self.switchOffDReaction(metIdx)
                    if verbose:
                        print('new initial conditions: \n\t t0: ', t0, ', solver t: ',solver.t,'\n\t q0: ', q0)
                    else:
                        print('in handleZeroEvent for ', self.dmetabolitesKeys[metIdx], '\n\t ', self.dmetabolitesKeys,
                              '\n\t new q0 ',q0, '; zeroDQ', dq,'\n\t new t0 ', t0, '; zeroDT ',zeroDT)
                    self.statPhaseID = 0
                else:
                    for i in range(len(self.dmetabolitesKeys)):
                        self.dmetabolites[self.dmetabolitesKeys[i]].updateQuantityQ(solver.y[i], self.volumes)
                    if verbose:
                        print('    dt: ', solver.t-self.T[-1],'\n    dQdt: ', self.dQdt,'\n    dQ: ', self.dQdt*(solver.t-self.T[-1]), '\n', [self.dmetabolites[m].quantity[-2:] for m in self.dmetabolitesKeys])
                    self.T.append(solver.t)
                self.FBAsolutions.append(self.solveFBA())
                #XXXX +self.odefluxes[-1]
                self.updateFluxes(self.FBAsolutions[-1])
                #print len(self.FBAsolutions)

            warnings.resetwarnings()

        if verbose:
            print('\n end runDynamicFBA \n Total time steps: ',len(self.T))

        #print len(self.FBAsolutions)
        #time and q got one step more than Fluxes
        # del self.T[-1]
        # for m in self.dmetabolitesKeys:
        #     del self.dmetabolites[m].quantity[-1]
        #     del self.dmetabolites[m].concentration[-1]
        if self.T[-1] >= maxT:
            self.stopDFBA = (False, 'Reached max time %f' % self.T[-1])
        elif step >= nMaxSteps:
            self.stopDFBA = (False, 'Reached max steps %d' % step)
        if self.cbModel and self.stopDFBA[0]: #self.FBAsol.status is not 'optimal':
            print('Pickling cobra model and dynamic model in /tmp/')
            dftmp = pandas.DataFrame(self.FBAsolutions[-10:-1], index=self.T[-10:-1])
            dftmp.to_csv('/tmp/stopDFBAcobra.log')
            ofile1 = open('/tmp/stopDFBAcobra.p', 'w')
            ofile2 = open('/tmp/stopDFBAdynamicmodel.p', 'w')
            cPickle.dump(self.cbModel, ofile1)
            cPickle.dump(self, ofile2)
            ofile1.close()
            ofile2.close()

        self.exitSimulation(True)
        return
