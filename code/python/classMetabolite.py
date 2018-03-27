#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  last modified: 2015/11/12     **
#************************************

import numpy as np
import classConvertionFactors

class DMetabolite(object):
    '''
    General dynamic metabolite for dynamic FBA
    * modName: name of the metabolite in the SBML/cobra model
    * quantity: 1-dim-list with initial quantity
    * isInternal: True/False for internal/external metabolites
    * metReactions: dict of reaction names and coefficient for their rates in the differential equation
         e.g.: {'oxygen_exchange': [(1, 'biomass')], 'oxygen_mass_transfer1': [(kLa*oxygen_gas_phase_C, None)], 'oxygen_mass_transfer2': [(-1*kLa, 'ex_oxygen')]}
               means:          dQdt = 1*biomass.concentration[-1]*oxygen_exchange.fluxSOLUTION + kLa*oxygen_gas_phase_C - 1*kLa*ex_oxygen.concentration[-1]
               NB 'oxygen_mass_transfer1' and 'oxygen_mass_transfer2' do not enter the model reactions dictionary since they do not correspond to any FBA flux!
    * quantityUnits: units
    * concentration: empty list to be initialized and then updated together with the quantities, but needs the volume information from the model
    * concentrationUnits: the info is taken from the volume object provided
    TODO: implement a Warning for concentrations update - currently just a printout
    '''
    def __init__(self, modName, quantity, isInternal, metReactions, quantityUnits=None):
        '''
        Init with the name of the metabolite in the SBML model, the starting quantity as one element list,
        if it's internal, a dictionary of reaction names, the units of quantity (deafult mmol)
        '''
        self.modName = modName
        self.quantity = quantity
        self.isInternal = isInternal
        self.metReactions = metReactions
        self.quantityUnits = quantityUnits if quantityUnits is not None else 'mmol'
        self.concentration = []
        self.concentrationUnits = None
        self.concentrationScaled = []
        self.concentrationScaledUnits = None

    def updateConcentration(self, vol):
        '''
        Update the metabolite concentration using the last element in quantity and the volume object passed (vol)
        concentrationScaled is the concentration scaled to cell density (concentration of substrate available per biomass unit) and it's used to compute the reaction rates
        Ref. Varma and Palsson 1994
        '''
        #make sure to use the last value of biomass as internal volume
        vol.update()
        v = vol.internal if self.isInternal else vol.external
        # vscal = vol.internal if self.isInternal else vol.bioavailable
        # NOT SURE should internal metabolites' conc dependent reaction rates be computed according to bioavailable biomass?
        vscal = vol.internal if self.isInternal else vol.bioavailable 
        self.concentration.append( self.quantity[-1]/v )
        if vol.internal < 0.00000000001:
            self.concentrationScaled.append( float('nan') )
        else:
            self.concentrationScaled.append( self.quantity[-1]/vscal )
        #self.concentrationScaled.append( self.quantity[-1]/v )
        return

    def setConcentrationUnits(self, vol):
        '''
        Set the units for metabolite concentration checking the volume object passed (vol)
        '''
        vu = vol.internalUnits if self.isInternal else vol.externalUnits
        self.concentrationUnits = self.quantityUnits+'/('+vu+')'
        self.concentrationScaledUnits = self.quantityUnits+'/('+vol.internalUnits+')'
        return

    def updateQuantityDQ(self, dQ, vol=None, nullThr=1e-10):
        '''
        Update the quantity by a variation dQ
        If volume is provided concentrations are updated
        '''
        Q = self.quantity[-1] + dQ
        if Q < nullThr:
            Q = 0
        self.quantity.append(Q)
        if vol is None:
            print 'Warning: concentrations are not going to be updated!'
            return
        else:
            self.updateConcentration(vol)
            return

    def updateQuantityQ(self, Q, vol=None, nullThr=1e-10):
        '''
        Update the quantity to the passed value Q
        If volume is provided concentrations are updated
        '''
        if Q < nullThr:
            Q = 0
        self.quantity.append(Q)
        if vol is None:
            print 'Warning: concentrations are not going to be updated!'
            return
        else:
            self.updateConcentration(vol)
            return

class Biomass(DMetabolite):
    '''
    Biomass is a special metabolite
    '''
    def __init__(self, quantity, metReactions, quantityUnits=None):
        '''
        Init as external DMetabolite w/o a corresponding name in the SBML model (in SBML biomass is a reaction and sum of metabolites)
        '''
        self.convertionFactors = classConvertionFactors.ConvertionFactors()
        self.quantityUnits = quantityUnits if quantityUnits is not None else 'gDW'
        q0 = self.convertToGDW(quantity[0])
        DMetabolite.__init__(self, None, [q0], False, metReactions, 'gDW')
        self.concentrationScaledUnits = 'a.u.'

    def convertToGDW(self, qty):
        '''
        Biomass must be in gDW
        Convert accordingly
        '''
        cfct = 1.
        if self.quantityUnits != 'gDW':
            cfct = self.convertionFactors.AtoB(self.quantityUnits, 'gDW')
        return qty*cfct

    def updateConcentration(self, vol):
        '''
        Update the biomass concentration
        The first time the method is called, the external volume is converted into biomass units and stored
        Concentration scaled in the case of biomass is Q (gDW)/ Vext (gDW)
        '''
        self.concentration.append( self.quantity[-1]/vol.external )
        self.concentrationScaled.append( self.quantity[-1]/(vol.external*self.convertionFactors.AtoB(vol.externalUnits, 'gDW')))
        #self.concentrationScaled.append( self.quantity[-1]/(vol.external))
        return

    def spaceLimitReached(self, volThr=1e-8):
        '''
        Method to check that biomass volume is contained in the externally defined volume within an arbitrary threshold
        TODO: not tested?
        '''
        return self.concentrationScaled > 1 - volThr
