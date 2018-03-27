#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  last modified: 2015/11/12     **
#************************************

#from __future__ import print_function
import cobra
#from cobra.flux_analysis.parsimonious import optimize_minimal_flux
import numpy as np
from scipy.integrate import ode, odeint
import random
from buildToyModel import buildToyModel as btm
import classConvertionFactors
import warnings
import json, cPickle


DEBUGVS = True
DEBUGVS = False
DEBUGQD = False
DEBUGCS = False


class Consortium(object):
    '''
    Model to simulate consortia of DynamicModels with dynamic FBA
    Implement different strategies for objectives
    '''
    def __init__(self, dmodels, dmetabolites, solver=None, name=None, timeU=None, savePath=None):
        '''
        Init with dictionary of DynamicModels,
        dictionary of metabolites where the key is the DMetabolite's key and the value is a list of dmodel keys,
        the solver to use (cglpk by default) and a name for the dFBA model
        '''
        self.dmodels = dmodels
        self.dmodelsKeys = self.dmodels.keys()
        self.dmetabolites = dmetabolites
        self.dmetabolitesKeys = self.dmetabolites.keys()
        self.T = [0.]
        self.Tunits = timeU if timeU is not None else 'hr'
        self.solver = solver if solver is not None else 'cglpk'
        self.name = name.replace(' ','_') if name is not None else 'a_dynamic_consortium'
        self.deaddmod = []
        self.savePath = savePath if savePath is not None else '/tmp/'
        #self.solver = solver if solver is not None else 'gurobi'
        #.initializeValueStorage()

    def addDynamicModel(self, dmod, dmodname):
        if dmodname in self.dmodelsKeys:
            print "Going to overwrite DynamicModel!"
        else:
            self.dmodelsKeys.append(dmodname)
        self.dmodels[dmodname] = dmod
        return


    def consortiumDifferentialFBA(self, y, verbose=False):
        '''
        Define the differential problem with dQdt ODEs
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
        #self.dmetabolitesKeys is initialized with the object consortium
        #it contains the keys of the shared (e.g. glucose, acetate...) and unique (e.g. biomass1, biomass2...) metabolites
        sizeQ = len(self.dmetabolitesKeys)
        dQdt = np.zeros((sizeQ,1))

        for i in range(sizeQ):
            m = self.dmetabolitesKeys[i]
            dQdt[i] = 0
            if DEBUGVS:
                print 'DEBUGVS met ', m, ' at time step ', len(self.T)
            #Loop over dmodels and use their dmetabolites
            for org in self.dmodelsKeys:
                #Try if the dmodel has the metabolite, else return none
                omet = self.dmodels[org].dmetabolites.get(m, None)
                if omet is None:
                    continue
                for r in omet.metReactions.keys():
                    if DEBUGVS:
                        print '\tDEBUGVS o ', org, ' r ', r
                        if r in self.dmodels[org].dreactionsKeys:
                            if not self.dmodels[org].dreactions[r].isODE:
                                if self.dmodels[org].FBAsolutions[-1] is None:
                                    print '\t\tFBAsolution: None'
                                else:
                                    print '\t\tFBAsolution: ', self.dmodels[org].FBAsolutions[-1][self.dmodels[org].dreactions[r].modName]
                            print '\t\tDReaction V: ', self.dmodels[org].dreactions[r].flux[-1]
                        else:
                            print '\t\t1.'
                    rate = self.dmodels[org].dreactions[r].flux[-1] if r in self.dmodels[org].dreactionsKeys else 1.
                    if DEBUGVS:
                        print '\t\tDEBUGVS rate ',  rate
                    metQ = 1
                    first = True
                    for n in omet.metReactions[r]:
                        mq = 1.
                        for ni in range(1, len(n)):
                            if n[ni] is not None:
                                mq = mq*y[self.dmetabolitesKeys.index(n[ni])]
                                if DEBUGVS:
                                    print '\t\tDEBUGVS mq ', mq, ' n ', n, ' ns index ', self.dmetabolitesKeys.index(n[1])
                        if first:
                            metQ = metQ*n[0]*mq
                        else:
                            metQ += n[0]*mq
                        first = False
                    dQdt[i] += rate*metQ
                    if DEBUGVS:
                        print '\t\t***DEBUGVS dQdt[i] += rate*metQ :: dQdt[i] += ', rate*metQ 
                if verbose:
                    print m, ' ODE for org ', org, ' is: '
                    for r in omet.metReactions.keys():
                        print omet.metReactions[r], '*v(',r,')'
                    print '\t dQdt is: ', dQdt[i]
        if verbose:
            print 'dQdt:', dQdt
        self.dQdt = dQdt
        return dQdt

    def getInitialConditions(self):
        '''
        Return last quantity and last time values to use in the next FBA step as initial conditions
        '''
        q0 = []
        for m in self.dmetabolitesKeys:
            dmodk = self.dmetabolites[m][0]
            omet = self.dmodels[dmodk].dmetabolites.get(m, None)
            if omet is None:
                print 'ACHTUNG! metabolite ',m, ' not found!'
                return
            q0.append(omet.quantity[-1])
        t0 = self.T[-1]
        return q0, t0

    def runConsortiumDynamicFBA(self, maxT=10, integrator='vode', stepChoiceLevel=(0, 0., 10., 1000), verbose=False, silent=True, restarting=False):
        '''
        Run the dynamic FBA as an integration problem using ode imported from scipy
        choose a different automatization of integration method via the stepChoiceLevel flag
        stepChoiceLevel=(type, max T step, max N steps)
        Randomize:
            * The first model to run
            * The objective function
        '''
        def dqdt(t, y, mod, verbose=False):
            '''
            Return the dQdt system of ODEs defined from the flux solutions
            '''
            return mod.consortiumDifferentialFBA(y, verbose)

        smethod = 'bdf'
        smax_step = 1
        self.stopDFBA = (False, 'Everything OK')

        if verbose:
            print '\n start runDynamicFBA \n'

        integratorSet = False
        
        if stepChoiceLevel[0] == 0:
            nMaxSteps = stepChoiceLevel[3]
            if verbose:
                nMaxSteps = 5
            solver = ode(dqdt).set_integrator(integrator)
            if integrator == 'dopri5':
                solver.set_integrator(integrator, nsteps=1, max_step=stepChoiceLevel[2])
                integratorSet =True
            elif integrator == 'vode':
                solver.set_integrator(integrator, min_step=stepChoiceLevel[1], max_step=stepChoiceLevel[2], method = 'bdf', order = 5)
                integratorSet =True
            else:
                print 'Option 0 is available only for vode and dopri5 SciPy solvers!'
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

            # XXX TODO:
            #choose randomly the organism to start with
            #dmk = self.dmodelsKeys[:]
            #random.shuffle(dmk)

            #use the first
            #set the parameters of the differential function dqdt: model and verbosity
            solver.set_f_params(self, verbose)
            # suppress Fortran-printed warning
            solver._integrator.iwork[2] = -1
            warnings.filterwarnings("ignore", category=UserWarning)
            q0, t0 = self.getInitialConditions()
            solver.set_initial_value(q0, t0)
            if verbose:
                print 'Initial Conditions to start integration: '
                print '\t t0: ', t0
                print '\t q0: ', q0

            #solve FBA with each dmodel's initial conditions
            #NB these are the fluxes the organisms would use if alone!
            #print len(self.FBAsolutions)
            #for dm in dmk:
            #    self.dmodels[dm].FBAsolutions.append(self.dmodels[dm].solveFBA())
            #    self.dmodels[dm].updateFluxes(self.dmodels[dm].FBAsolutions[-1])
            if not restarting:
                self.updateFBASolFluxesT(False)

            step = 0
            substeps = 0
            eventIdx = 1

            #now need to integrate the whole ODE system and update all organisms
            while not self.stopDFBA[0] and self.T[-1]<maxT and step < nMaxSteps:
                step += 1
                ## statPhase
                for dm in self.dmodelsKeys:
                    if self.dmodels[dm].cbModel is not None:
                        self.dmodels[dm].updateStatPhase()
                if stepChoiceLevel[0] == 1:
                    solver.integrate(solver.t + grid_dt)
                else:
                    solver.integrate(maxT, step=True)
                if verbose:
                    print 'solver time: ', solver.t
                    print 'solver solution: '
                    for i in range(len(solver.y)):
                        print '\t', self.dmetabolitesKeys[i],': ',solver.y[i]
                elif not silent:
                    if step%10==0:
                        print 'Step: ', step, ', Substep: ', substeps,', Time:', solver.t
                handleZeroEvent, metIdx = self.checkForNegativeSolutions(solver.y)
                if handleZeroEvent:
                    if verbose:
                        print '****************************'
                        print '!in handleZeroEvent!'
                    substeps += 1
                    #slope = (min(solver.y) - self.dmetabolites[metId].quantity[-1])/(solver.t - self.T[-1])
                    slope = float(self.dQdt[metIdx])
                    if slope == 0:
                        self.stopDFBA = (True, 'Zero slope')
                        continue
                    ##metIdx is the index in the consortium! use the key for the models
                    metKey = self.dmetabolitesKeys[metIdx]

                    zeroDT = 0.-self.dmodels[self.dmetabolites[metKey][0]].dmetabolites[metKey].quantity[-1]/slope
                    for i in range(len(self.dmetabolitesKeys)):
                        dq = float(self.dQdt[i])*zeroDT
                        if i == metIdx:
                            self.updateQuantity(metKey, 0., False)
                        else:
                            self.updateQuantity(self.dmetabolitesKeys[i], dq, True)
                    self.T.append(self.T[-1]+zeroDT)
                    q0, t0 = self.getInitialConditions()
                    solver.set_initial_value(q0, t0)
                    if verbose:
                        print 'new initial conditions: '
                        print '\t t0: ', t0, ', solver t: ',solver.t
                        print '\t q0: ', q0
                    else:
                        print 'in handleZeroEvent for ', self.dmetabolitesKeys[metIdx]
                        print '\t ', self.dmetabolitesKeys
                        print '\t new q0 ',q0, '; zeroDQ', dq,'\n\t new t0 ', t0, '; zeroDT ',zeroDT
                    #StatPhase NOT WORKING YET
                    #print 'Set lag-phase for model ',self.dmetabolites[metKey][0]
                    #self.dmodels[self.dmetabolites[metKey][0]].statPhaseID += 1
                else:
                    for i in range(len(self.dmetabolitesKeys)):
                        self.updateQuantity(self.dmetabolitesKeys[i], solver.y[i], False)
                    if verbose:
                        print '    dt: ', solver.t-self.T[-1]
                        print '    dQdt: ', self.dQdt
                        print '    dQ: ', self.dQdt*(solver.t-self.T[-1])
                        print [self.dmodels[self.dmetabolites[m][0]].dmetabolites[m].quantity[-2:] for m in self.dmetabolitesKeys]
                    self.T.append(solver.t)
                self.updateFBASolFluxesT()
                #print len(self.FBAsolutions)
                #self.stopDFBA = (True, 'Debugging')

            warnings.resetwarnings()

        if verbose:
            print '\n end runDynamicFBA \n Total time steps: ',len(self.T)

        #print len(self.FBAsolutions)
        #time and q got one step more than Fluxes
        # del self.T[-1]
        # for m in self.dmetabolitesKeys:
        #     del self.dmetabolites[m].quantity[-1]
        #     del self.dmetabolites[m].concentration[-1]
        self.exitSimulation(True)
        return

    def updateQuantity(self, dmetk, dq, isincrem=True):
        '''
        Retrieve the dmodel's dmetabolite and calls the updateq/updatedq method
        '''
        for dmod in self.dmetabolites[dmetk]:
            if isincrem:
                self.dmodels[dmod].dmetabolites[dmetk].updateQuantityDQ(dq, self.dmodels[dmod].volumes)
            else:
                self.dmodels[dmod].dmetabolites[dmetk].updateQuantityQ(dq, self.dmodels[dmod].volumes)
        return

    def updateFBASolFluxesT(self, notFirstCall=True):
        for dm in self.dmodelsKeys:
            if notFirstCall:
                self.dmodels[dm].T.append(self.T[-1])
            else:
                self.dmodels[dm].stopDFBA = self.stopDFBA
                self.dmodels[dm].stopFBA = (False, 'FBA problem OK')
                self.dmodels[dm].initFBA()
                continue
            if dm in self.deaddmod:
                continue
            self.dmodels[dm].FBAsolutions.append(self.dmodels[dm].solveFBA())
            self.dmodels[dm].updateFluxes(self.dmodels[dm].FBAsolutions[-1])
            #print 'Model ', dm, ' v EX_ac: ', self.dmodels[dm].FBAsolutions[-1]['EX_ac_e_']
            if self.dmodels[dm].stopDFBA[0]:
                # Set all fluxes to zero and remove from dFBA
                self.dmodels[dm].FBAsolutions[-1] = dict( zip( self.dmodels[dm].FBAsolutions[-2].keys(), len(self.dmodels[dm].FBAsolutions[-2].keys())*[0.0] ) )
                self.deaddmod.append(dm)
                #self.stopDFBA = (True, 'In model %s %s' % (dm, self.dmodels[dm].stopDFBA[1]))
        return

#Functions from DynamicModel class
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
        return

    def exitSimulation(self, allGood):
        print 'End of dFBA simulation, stopDFBA message:\n\t', self.stopDFBA[1]
        ofile = open(self.savePath+'endOfSimulation-'+self.name+'.p', 'w')
        cPickle.dump(self, ofile)
        ofile.close()
        if not allGood:
            exit(1)
        return
