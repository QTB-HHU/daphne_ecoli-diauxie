#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  last modified: 2015/11/12     **
#************************************

import numpy as np
#from numpy import log,exp

'''
Begin helper code
'''
class NoKineticLawException(Exception):
    '''
    NoKineticLawException is raised when a DReaction tries to
    update its bounds without having a KineticLaw defined
    '''
'''
End helper code
'''


class KineticLaw(object):
    '''
    General Kinetic Law class
    Optional DMetabolite attribute: substrate
    Optional float attributes (defaulted to 1.0): catalyticEff; maxVelocity; mmConstant; hillCoeff
    If lowerBound/upperBound is specified, it will stay fixed and will never be overwrote
    '''
    def __init__(self, substrate=None, catalyticEff=None, maxVelocity=None, lowerBound=None, upperBound=None, useQs=False):
        '''
        Init with a DMetabolite as substrate, catalytic efficiency, max V, Michaelis-Menten constant, Hill coefficient
        '''
        self.catalyticEff = catalyticEff if catalyticEff is not None else 1.0
        self.maxVelocity = maxVelocity if maxVelocity is not None else 1.0
        self.substrate = substrate
        #self.lowerBound = lowerBound if lowerBound is not None else -1000.
        #self.upperBound = upperBound if upperBound is not None else 1000.
        self.lowerBound = lowerBound
        self.upperBound = upperBound
        self.useQs = useQs
        # TO DO: default to scaled or not??? right now is scaled!
        #, useCscaled=False
        #self.useCscaled = useCscaled

    def computeBounds(self, val=-1, time=None, rSign=(-1., 1.)):
        '''
        Call computeRate and return lowerBound and upperBound
        '''
        lowR, uppR = self.lowerBound, self.upperBound
        r = self.computeRate(val, time)
        if r < 1e-10:
            r = 0.
        if lowR is None:
            lowR = r
        if uppR is None:
            uppR = r
        lb = rSign[0]*abs(lowR)
        ub = rSign[1]*abs(uppR)
        return lb, ub

    def computeRate(self, val=-1, time=None):
        '''
        Function specific to kinetics, will be overwrote in the following classes, here needed for computeBounds
        '''
        return 1.

    def setZeroLowerBound(self, qThr=1e-4):
        '''
        Check if the lower bound should be set to 0
        Returns True if the substrate is defined and its concentration is lower than the threshold (defaulted to 1e-4)
        '''
        if self.substrate is None:
            return False
        else:
            #return self.substrate.concentrationScaled[-1] < qThr
            return self.substrate.concentration[-1] < qThr

    def heavisideStepFunction(self, val):
        '''
        Simple step function: 1 or 0 for positive or negative values respectively
        '''
        if val > 0:
            return 1
        else:
            return 0

    def mmrate(self, params, val=-1):
        '''
        Compute the reaction rate according to 1st order Henry-Michaelis-Menten equation
        params is an ntuple (conc, mmk, maxv=1, hic=1)
        '''
        #conc = params[0].quantity[val]
        conc = params[0]
        mmk = params[1]
        if len(params) < 3:
            maxv = 1.
        else:
            maxv = params[2]
        if len(params) < 4:
            hic = 1.
        else:
            hic = params[3]
        if mmk < 1e-5:
            rate = maxv
        else:
            rate = maxv*(conc**hic)/(conc**hic + mmk**hic)
        return rate

    def massActionRate(self, substrates, kcat, val=-1):
        concs = 1.
        for s in substrates:
            concs = concs*(s[0].concentration[val]**s[1])
        return kcat*concs

    def rateControlFunction(self, desiredConc, currentQty, deltaT, linkedRate=0., isGrowthLinked=False):
        '''
        XXXXX FIX XXXXX Compute the rate needed to keep the substrate concentration at the desiredConc level
        '''
        if isGrowthLinked:
            # the concentration gets diluted - scale this factor in order to keep up
            linkedRate = linkedRate*desiredConc
        rate = linkedRate + (desiredConc - currentQty)/deltaT
        return rate


class ConcentrationMaintenanceFunction(KineticLaw):
    '''
    Function to reproduce experimental setups where concentrations are artificially maintained
    '''
    def __init__(self, substrate, cthreshold, timedelay=None, linkedReactions=None):
        '''
        Init as KineticLaw with the DMetabolite to be maintained as substrate
        '''
        KineticLaw.__init__(self, substrate)
        self.conc = cthreshold
        self.deltat = timedelay if timedelay is not None else 0.
        self.linkedReactions=linkedReactions

    def computeRate(self, val=-1, time=None):
        rate = 0
        if time is not None and (time[val] - time[0] > self.deltat) and (self.conc - self.substrate.concentration[val]) > 0:
            rate = (self.conc - self.substrate.concentration[val])#/(time[val]-time[val-1])
            if self.linkedReactions is not None:
                for r in self.linkedReactions:
                    rate += -1*r.flux[val]
        #print('AAAAAA', rate)
        return rate

    def computeRateFromC(self, conc, val=-1, time=None):
        rate = 0
        if time is not None and time[val] - time[0] > self.deltat and self.conc > conc:
            rate = (self.conc - conc)/(time[val]-time[val-1])
        return rate
    
class MichaelisMenten1(KineticLaw):
    '''
    HenryMichaelisMenten function for 1st order reaction
    Inherits from KineticLaw, same attributes
    '''
    def __init__(self, substrate, catalyticEff=None, maxVelocity=None, mmConstant=None, hillCoeff=None, lowerBound=None, upperBound=None, minRate=0):
        '''
        Init as KineticLaw with a DMetabolite as substrate, catalytic efficiency (1/time), max V (concentration/time), Michaelis-Menten constant (concentration), Hill coefficient (dimensionless)
        '''
        KineticLaw.__init__(self, substrate, catalyticEff, maxVelocity, lowerBound, upperBound)
        self.mmConstant = mmConstant if mmConstant is not None else 1.0
        self.hillCoeff = hillCoeff if hillCoeff is not None else 1.0
        self.minRate = minRate

    def computeRate(self, val=-1, time=None):
        '''
        Compute the reaction rate according to 1st order Henry-Michaelis-Menten equation
        '''
        #conc = self.substrate.concentrationScaled[val]
        conc = self.substrate.concentration[val]
        rate = self.mmrate((conc, self.mmConstant, self.maxVelocity, self.hillCoeff))
        #rate = self.maxVelocity*(conc**self.hillCoeff)/(conc**self.hillCoeff + self.mmConstant**self.hillCoeff)
        return self.minRate + rate

    def computeRateFromC(self, conc):
        '''
        Compute the reaction rate according to 1st order Henry-Michaelis-Menten equation
        '''
        #rate = self.maxVelocity*(conc**self.hillCoeff)/(conc**self.hillCoeff + self.mmConstant**self.hillCoeff)
        rate = self.mmrate((conc, self.mmConstant, self.maxVelocity, self.hillCoeff))
        return self.minRate + rate

class MichaelisMentenMaintenance(KineticLaw):
    '''
    Rate function for substrate_storage --> substrate_internal
    Want to move substrate from the storage pool only if the concentration level is below the maintenance level:
        * A * Reaction off if [int] > maintenanceThr
        * B * Reaction on (HenryMichaelisMenten like) if [int] < maintenanceThr
    '''
    def __init__(self, substrate, control, maintenanceThr, linkedReaction=None, catalyticEff=None, maxVelocity=None, mmConstant=None, hillCoeff=None, lowerBound=None, upperBound=None, substratesWithSt=[]):
        '''
        Init as KineticLaw with a DMetabolite as substrate, catalytic efficiency (1/time), max V (concentration/time), Michaelis-Menten constant (concentration), Hill coefficient (dimensionless)
        '''
        KineticLaw.__init__(self, substrate, catalyticEff, maxVelocity, lowerBound, upperBound)
        self.control = control
        self.mmConstant = mmConstant if mmConstant is not None else 1.0
        self.hillCoeff = hillCoeff if hillCoeff is not None else 1.0
        self.maintenanceThr = maintenanceThr
        self.linkedReaction = linkedReaction[0] if linkedReaction is not None else None
        self.linkedReactionSF = linkedReaction[1] if linkedReaction is not None and len(linkedReaction) > 1 else 1.
        self.substratesWithSt = substratesWithSt

    def computeRate(self, val=-1, time=None):
        '''
        Reaction off if [int] > maintenanceThr
        Compute the reaction rate according to 1st order Henry-Michaelis-Menten equation
        '''
        #conc = self.substrate.concentrationScaled[val]
        conc = self.substrate.concentration[val]
        #concKM = self.control.concentrationScaled[val]
        concKM = self.control.concentration[val]
        ft = -1
        if len(self.substrate.quantity)>1 and self.linkedReaction is not None:
            ft = self.linkedReaction.flux[val]
        return self.computeRateFromC(conc, concKM, ft)

    def computeRateFromC(self, conc, concKM, linkedFlux=-1):
        '''
        Calls computeRate
        '''
        mmsub = self.mmrate((conc, self.mmConstant, self.maxVelocity, 1.))
        mmcon = self.mmrate((concKM, self.maintenanceThr, -1., 15))
        #print concKM, mmcon
        linkedrate = 999.
        if linkedFlux >= 0:
            linkedrate = self.maintenanceThr*self.linkedReactionSF*linkedFlux
        return min(linkedrate, mmsub*(1 + mmcon))

class MichaelisMentenStorage(MichaelisMentenMaintenance):
    '''
    Rate function for substrate_internal --> substrate_storage
    Want to move substrate into the storage pool only if the concentration_i level is above the maintenance level and the concentration_s is below saturation:
        * Reaction off if [int] < maintenanceThr
        * Reaction on (HenryMichaelisMenten like) if [sto] < saturationThr
        * Reaction on (following growth rate to maintain constant concentration)  if [sto] >= saturationThr
    '''
    def __init__(self, substrate, control, maintenanceThr, saturationThr, linkedReaction=None, catalyticEff=None, maxVelocity=None, mmConstant=None, hillCoeff=None, lowerBound=None, upperBound=None, substratesWithSt=[]):
        '''
        Init as MichaelisMentenSaturation, override computeRate, inherits computeRateFromC
        '''
        MichaelisMentenMaintenance.__init__(self, substrate, control, maintenanceThr, linkedReaction, catalyticEff, maxVelocity, mmConstant, hillCoeff, lowerBound, upperBound, substratesWithSt)
        self.saturationThr = saturationThr

    def computeRate(self, val=-1, time=None):
        '''
        Storage:
            1) switch off storing if [int] < maintenanceThr
            2) if saturated, control flux to keep up with growth
        '''
        #conc = self.substrate.concentrationScaled[val]
        conc = self.substrate.concentration[val]
        #concKM = self.control.concentrationScaled[val]
        concKM = self.control.concentration[val]
        ft = -1
        if len(self.substrate.quantity)>1 and self.linkedReaction is not None:
            ft = self.linkedReaction.flux[val]
        return self.computeRateFromC(conc, concKM, ft)

    def computeRateFromC(self, conc, concKM, linkedFlux=-1):
        '''
        Calls computeRate
        '''
        concShift = conc - self.maintenanceThr
        if concShift < 0:
            concShift = 0.
        #mmsub = self.mmrate((conc, self.mmConstant, self.maxVelocity, 5))
        mmsub = self.mmrate((concShift, self.mmConstant+self.maintenanceThr, self.maxVelocity, 1))
        mmcon = self.mmrate((concKM, self.saturationThr, -1., 15))
        #print concKM, mmcon
        linkedrate = 999.
        if linkedFlux >= 0:
            linkedrate = self.saturationThr*self.linkedReactionSF*linkedFlux
        return min(linkedrate, mmsub*(1 + mmcon))


class MichaelisMentenOnOff(KineticLaw):
    '''
    Rate function with MM Kinetics, ON only if onThr < linkedFlux < offThr
    E.g. active only if growth below maintenance threshold:  onThr=0., offThr=maintenance
    '''
    def __init__(self, substrate, catalyticEff=None, maxVelocity=None, mmConstant=None, hillCoeff=None, lowerBound=None, upperBound=None, linkedReaction=None, onThr=-1000., offThr=1000.):
        '''
        Init as KineticLaw with a DMetabolite as substrate, catalytic efficiency (1/time), max V (concentration/time), Michaelis-Menten constant (concentration), Hill coefficient (dimensionless)
        '''
        KineticLaw.__init__(self, substrate, catalyticEff, maxVelocity, lowerBound, upperBound)
        self.mmConstant = mmConstant if mmConstant is not None else 1.0
        self.hillCoeff = hillCoeff if hillCoeff is not None else 1.0
        self.linkedReaction = linkedReaction[0] if linkedReaction is not None else None
        self.linkedReactionSF = linkedReaction[1] if linkedReaction is not None and len(linkedReaction) > 1 else 1.
        self.onThr = onThr
        self.offThr = offThr

    def computeRate(self, val=-1, time=None):
        '''
        Reaction off if [int] > maintenanceThr
        Compute the reaction rate according to 1st order Henry-Michaelis-Menten equation
        '''
        #conc = self.substrate.concentrationScaled[val]
        conc = self.substrate.concentration[val]
        lf = 0.
        if len(self.substrate.quantity)>1 and self.linkedReaction is not None:
            lf = self.linkedReaction.flux[-1]*self.linkedReactionSF
        rate = self.computeRateFromC(conc, lf)
        #print(rate, lf)
        return rate


    def computeRateFromC(self, conc, linkedFlux=0):
        '''
        Calls computeRate
        '''
        if linkedFlux > self.onThr and linkedFlux < self.offThr:
            return self.mmrate((conc, self.mmConstant, self.maxVelocity, self.hillCoeff))
        return 0.


class MichaelisMentenLinked(KineticLaw):
    '''
    Rate function with MM Kinetics, ON only if onThr < linkedFlux < offThr
    E.g. active only if growth below maintenance threshold:  onThr=0., offThr=maintenance
    '''
    def __init__(self, substrate, catalyticEff=None, maxVelocity=None, mmConstant=None, hillCoeff=None, lowerBound=None, upperBound=None, linkedReaction=None, onThr=-1000., offThr=1000., offset=0.):
        '''
        Init as KineticLaw with a DMetabolite as substrate, catalytic efficiency (1/time), max V (concentration/time), Michaelis-Menten constant (concentration), Hill coefficient (dimensionless)
        '''
        KineticLaw.__init__(self, substrate, catalyticEff, maxVelocity, lowerBound, upperBound)
        self.mmConstant = mmConstant if mmConstant is not None else 1.0
        self.hillCoeff = hillCoeff if hillCoeff is not None else 1.0
        self.linkedReaction = linkedReaction[0] if linkedReaction is not None else None
        self.linkedReactionSF = linkedReaction[1] if linkedReaction is not None and len(linkedReaction) > 1 else 1.
        self.onThr = onThr
        self.offThr = offThr
        self.offset = offset

    def computeRate(self, val=-1, time=None):
        '''
        Reaction off if [int] > maintenanceThr
        Compute the reaction rate according to 1st order Henry-Michaelis-Menten equation
        '''
        #conc = self.substrate.concentrationScaled[val]
        conc = self.substrate.concentration[val]
        lf = 1.
        if len(self.substrate.quantity)>1 and self.linkedReaction is not None:
            lf = self.linkedReaction.flux[-1]*self.linkedReactionSF
        rate = self.computeRateFromC(conc, lf)
        #print(rate, lf)
        return rate
        #return rate*lf


    def computeRateFromC(self, conc, linkedFlux=0):
        '''
        Calls computeRate
        '''
        if linkedFlux > self.onThr and linkedFlux < self.offThr:
            return self.offset + self.mmrate((conc, self.mmConstant, self.maxVelocity, self.hillCoeff))
        return self.offset


class MichaelisMentenLinkedSubstrate(KineticLaw):
    '''
    Rate function with MM Kinetics, ON only if onThr < linkedFlux < offThr
    E.g. active only if growth below maintenance threshold:  onThr=0., offThr=maintenance
    '''
    def __init__(self, substrate, catalyticEff=None, maxVelocity=None, mmConstant=None, hillCoeff=None, lowerBound=None, upperBound=None, linkedSubstrate=None, onThrC=-1000., offThrC=1000., offset=0.):
        '''
        Init as KineticLaw with a DMetabolite as substrate, catalytic efficiency (1/time), max V (concentration/time), Michaelis-Menten constant (concentration), Hill coefficient (dimensionless)
        '''
        KineticLaw.__init__(self, substrate, catalyticEff, maxVelocity, lowerBound, upperBound)
        self.mmConstant = mmConstant if mmConstant is not None else 1.0
        self.hillCoeff = hillCoeff if hillCoeff is not None else 1.0
        self.linkedSubstrate = linkedSubstrate if linkedSubstrate is not None else None
        self.onThrC = onThrC
        self.offThrC = offThrC
        self.offset = offset

    def computeRate(self, val=-1, time=None):
        '''
        Reaction off if [int] > maintenanceThr
        Compute the reaction rate according to 1st order Henry-Michaelis-Menten equation
        '''
        #conc = self.substrate.concentrationScaled[val]
        conc = self.substrate.concentration[val]
        lf = -1
        if self.linkedSubstrate is not None:
            lf = self.linkedSubstrate.concentration[val]
        rate = self.computeRateFromC(conc, lf)
        #print(rate, lf)
        return rate
        #return rate*lf


    def computeRateFromC(self, conc, lconc=-1):
        '''
        Calls computeRate
        '''
        SF = 1.
        if lconc > 0:
            SF = 1 - max(1e-5, self.mmrate((lconc, self.offThrC, 1., 10)))
            #print(lconc, SF)
        return self.offset + self.mmrate((conc, self.mmConstant, self.maxVelocity, self.hillCoeff))*SF
        
    
class Enzymatic(KineticLaw):
    '''
    Simple Enzymatic rate law v = k_cat*[E]
    '''
    def __init__(self, enzyme, catalyticEff=None, maxVelocity=10., lowerBound=None, upperBound=None):
        '''
        Init with the enzyme as the substrate and the catalytic efficiency
        '''
        KineticLaw.__init__(self, enzyme, catalyticEff, maxVelocity, lowerBound, upperBound)

    def computeRate(self, val=-1, time=None):
        '''
        Compute the reaction rate according to v = k_cat*[E]
        The enzyme concentration is evaluated by default as the last updated value
        '''
        #rate = self.catalyticEff*self.substrate.concentrationScaled[val]
        rate = self.catalyticEff*self.substrate.concentration[val]
        return rate

    def computeRateFromC(self, conc):
        '''
        Compute the reaction rate according to v = k_cat*[E]
        The enzyme concentration is evaluated by default as the last updated value
        '''
        rate = self.catalyticEff*conc
        return rate

class StepCofactor(KineticLaw):
    '''
    Simple Step function of velocity based on cofactor concentration higher or lower than a threshold
    '''
    def __init__(self, cofactor, maxVelocity=10., cofactorThr=0.001, lowerBound=None, upperBound=None):
        '''
        Init with the cofactor as the substrate and the catalytic efficiency
        '''
        KineticLaw.__init__(self, cofactor, None, maxVelocity, lowerBound, upperBound)
        self.cofactorThr = cofactorThr

    def computeRate(self, val=-1, time=None):
        '''
        Compute the rate as maxV*(1+heaviside([C]-threshold))
        '''
        #rate = self.maxVelocity*(1 + self.heavisideStepFunction(self.substrate.concentrationScaled[val]-self.cofactorThr))
        rate = self.maxVelocity*(1 + self.heavisideStepFunction(self.substrate.concentration[val]-self.cofactorThr))
        return rate

    def computeRateFromC(self, conc):
        '''
        Compute the rate as maxV*(1+heaviside([C]-threshold))
        '''
        rate = self.maxVelocity*(1 + self.heavisideStepFunction(conc-self.cofactorThr))
        return rate

class StepSaturation(KineticLaw):
    '''
    Simple 1/0 saturation condition
    '''
    def __init__(self, substrate, maxVelocity=1000., concentrationThr=1., substrate2=None, lowerBound=None, upperBound=None):
        '''
        Init with the substrate, the max velocity and the concentration threshold
        concentrationThr can be an absolute number or a factor for the substrate2 concentration
        Units must be consistent!
        '''
        KineticLaw.__init__(self, substrate, None, maxVelocity, lowerBound, upperBound)
        self.concentrationThr = concentrationThr
        self.substrate2 = substrate2

    def computeRate(self, val=-1, time=None):
        '''
        Compute the rate as maxV*heaviside(thr - [C])
        '''
        threshold = self.concentrationThr
        if self.substrate2 is not None:
            #threshold = self.concentrationThr*self.substrate2.concentrationScaled[val]
            threshold = self.concentrationThr*self.substrate2.concentration[val]
        #rate = self.maxVelocity*self.heavisideStepFunction(threshold - self.substrate.concentrationScaled[val])
        rate = self.maxVelocity*self.heavisideStepFunction(threshold - self.substrate.concentration[val])
        return rate

    def computeRateFromC(self, conc, conc2):
        '''
        Compute the rate as maxV*heaviside(thr - [C])
        '''
        threshold = self.concentrationThr
        if self.substrate2 is not None:
            threshold = self.concentrationThr*conc2
        rate = self.maxVelocity*self.heavisideStepFunction(threshold - conc)
        return rate


class FixedBound(KineticLaw):
    '''
    No kinetic law
    '''
    def __init__(self, lowerBound=None, upperBound=None, substrate=None, linkedReaction=None):
        '''
        If lowerBound and/or upperBound are not None, they overwrite any other bound
        linkedReaction = (DReaction object, Optional Scale Factor, Optional Value to use as symmetric constraint)
        '''
        self.linkedReaction = linkedReaction[0] if linkedReaction is not None else None
        self.linkedReactionSF = linkedReaction[1] if linkedReaction is not None and len(linkedReaction) > 1 else 1.
        self.symBound = linkedReaction[2] if linkedReaction is not None and len(linkedReaction) > 2 else 1000.
        KineticLaw.__init__(self, substrate, None, None, lowerBound, upperBound)

    def computeRate(self, val=-1, time=None):
        '''
        Return rate zero if a substrate is defined and its concentration at val is lower than 1e-3
        Else Return the maximal rate
        '''
        conc = 1.
        if self.substrate is not None:
            #if self.substrate.concentrationScaled[val] < 1e-3:
            if self.substrate.concentration[val] < 1e-3:
                return 0.
            else:
                #conc = self.substrate.concentrationScaled[val]
                conc = self.substrate.concentration[val]
        if self.linkedReaction is not None:
            if len(self.linkedReaction.flux) > 0:
                return min(abs(self.symBound), abs(self.linkedReactionSF*self.linkedReaction.flux[-1]))
            else:
                return 0. #abs(self.symBound)
        return conc*self.upperBound

    def computeRateFromC(self, conc):
        '''
        Fake compute rate
        '''
        return self.computeRate()

class MinimalFixedBound(KineticLaw):
    '''
    No kinetic law
    '''
    def __init__(self, lowerBound=None, upperBound=None, substrate=None, linkedReaction=None):
        '''
        Ensure a minimal rate (self.symBound)
        Choose max(self.symBound, linked flux)
        If lowerBound and/or upperBound are not None, they overwrite any other bound
        linkedReaction = (DReaction object, Optional Scale Factor, Optional Value to use as symmetric constraint)
        '''
        self.linkedReaction = linkedReaction[0] if linkedReaction is not None else None
        self.linkedReactionSF = linkedReaction[1] if linkedReaction is not None and len(linkedReaction) > 1 else 1.
        self.symBound = linkedReaction[2] if linkedReaction is not None and len(linkedReaction) > 2 else 1000.
        KineticLaw.__init__(self, substrate, None, None, lowerBound, upperBound)

    def computeRate(self, val=-1, time=None):
        '''
        Return rate zero if a substrate is defined and its concentration at val is lower than 1e-3
        Else Return the maximal rate
        '''
        conc = 1.
        if self.substrate is not None:
            #if self.substrate.concentrationScaled[val] < 1e-3:
            if self.substrate.concentration[val] < 1e-3:
                return 0.
            else:
                #conc = self.substrate.concentrationScaled[val]
                conc = self.substrate.concentration[val]
        if self.linkedReaction is not None:
            if len(self.linkedReaction.flux) > 0:
                return max(abs(self.symBound), abs(self.linkedReactionSF*self.linkedReaction.flux[-1]))
            else:
                return 0. #abs(self.symBound)
        return conc*self.upperBound

    def computeRateFromC(self, conc):
        '''
        Fake compute rate
        '''
        return self.computeRate()

class SubstrateBound(KineticLaw):
    '''
    No kinetic law
    '''
    def __init__(self, lowerBound=None, upperBound=None, substrate=None, substratesWithFc=[], useQ=False):
        '''
        Init with a substrate and a maximal rate defaulted to 1000.
        computeBounds() is inherited
        '''
        KineticLaw.__init__(self, substrate, None, None, lowerBound, upperBound)
        self.substratesWithFc = substratesWithFc
        self.useQ = useQ

    def computeRate(self, val=-1, time=None):
        '''
        Return rate zero if a substrate is defined and its concentration at val is lower than 1e-3
        Else Return the maximal rate
        '''
        #if self.substrate is not None and self.substrate.concentrationScaled[val] < 1e-3:
        if self.substrate is not None and self.substrate.concentration[val] < 1e-3:
            return 0.
        k = 1.
        for sb in self.substratesWithFc:
            if self.useQ:
                k = k*sb[0].quantity[val]*sb[1]
            else:
                #k = k*sb[0].concentrationScaled[val]*sb[1]
                k = k*sb[0].concentration[val]*sb[1]
        return k

    def computeRateFromC(self, conc):
        '''
        Fake compute rate
        '''
        return self.computeRate()

class Logistic(KineticLaw):
    '''
    Logistic kinetics
    '''
    def __init__(self, substrate, params, useQs = True):
        '''
        Init: params[0] -> growth params[1] -> death
        '''
        KineticLaw.__init__(self, substrate, useQs=useQs)
        self.params = params

    def computeRate(self, val=-1, time=None):
        if not self.useQs:
            print 'ACHTUNG'
        rate = self.params[0] - self.substrate.quantity[val]*self.params[1]
        return rate

class GrowthCurve(KineticLaw):
    '''
    Growth kinetics as from:
    x = K/(1 + ((K-X0)/X0)*exp(-R*t) )
    dx/dt = Rx*(1 - x/K)
    rate = R*(1 - x[-1]/K)
    '''
    def __init__(self, substrate, params, useQs = True):
        '''
        params[0] -> R (growth rate)
        params[1] -> K (carrying capacity)
        '''
        KineticLaw.__init__(self, substrate, useQs=useQs)
        self.params = params

    def computeRate(self, val=-1, time=None):
        if not self.useQs:
            print 'ACHTUNG'
        rate = self.params[0]*(1 - self.substrate.quantity[val]/self.params[1])
        #print(self.params[0], self.substrate.quantity[val], self.params[1], rate)
        return rate

    def computePopulation(self, val=-1, time=None):
        if time is None:
            print 'ACHTUNG'
            return 0.
        R = self.params[0]
        X0 = self.substrate.quantity[0]
        K = self.params[1]
        Xt = K/(1 + ((K-X0)/X0)*np.exp(-R*t) )
        return Xt

    def computeBounds(self, val=-1, time=None, rSign=(1., 1.)):
        '''
        Call computeRate and return lowerBound and upperBound
        '''
        r = self.computeRate(val, time)
        return r, r

class ProductionConsumption(KineticLaw):
    '''
    Logistic kinetics
    '''
    def __init__(self, growth, ub=1., f1params=[], typef1=0, useQs = True):#, f2params, f3params
        '''
        Init: substrate is the organism's Biomass
        growth = growth kin function
        '''
        KineticLaw.__init__(self, upperBound=ub, useQs=useQs)
        self.growth = growth
        self.f1params = f1params
        self.typef1 = typef1

    def computeRate(self, val=-1, time=None):
        #if growth is positive production/consumption in possible, else not!
        rate = self.heavisideStepFunction(self.growth.upperBound)*self.growth.upperBound*self.upperBound
        if len(self.f1params) > 0:
            f1 = 1.
            for p in self.f1params:
                f1 = f1*self.mmrate(p)
            if self.typef1 == 1:
                f1 = (1 - f1)
            rate = rate*(1 + f1)
        return rate

    def computeBounds(self, val=-1, time=None, rSign=(-1., 1.)):
        '''
        Call computeRate and return lowerBound and upperBound
        '''
        r = self.computeRate(val)
        lb = -1*r
        ub = r
        return lb, ub

class DeathRateFunction(KineticLaw):
    '''
    Init:
    '''
    def __init__(self, vd, f1params, f2params, useQs = True, timedep=(False, 0.)):
        KineticLaw.__init__(self, upperBound=vd, useQs=useQs)
        self.f1params = f1params
        self.f2params = f2params
        self.vd = vd
        self.timedep = timedep

    def computeRate(self, val=-1, time=None):
        f1 = 1.
        f2 = 1.
        f3 = 1.
        if len(self.f1params) > 0:
            f1 = 1 + self.mmrate(self.f1params, val)
        if len(self.f2params) > 0:
            f2sum = 0.
            for i in self.f2params:
                if type(i) == int or type(i) == float:
                    f2sum += i
                else:
                    f2sum += i.upperBound
            f2 = 1/f2sum
        if self.timedep[0]:
            if time is not None:
                #f3 = np.log(10*time+1)
                f3 = self.timedep[1]*time[val]/(time[val]+1.)
                #f3 = np.log(time+1)*self.timedep[1]*time/(time+0.1)
                #print f3
        rate = self.vd*f1*f2*f3
        #print 'death:', rate
        return rate

    def computeBounds(self, val=-1, time=None, rSign=(-1., 1.)):
        '''
        Call computeRate and return lowerBound and upperBound
        '''
        r = self.computeRate(val, time)
        lb = -1*r
        ub = r
        return lb, ub

class LogisticComp(KineticLaw):
    '''
    Logistic kinetics
    '''
    def __init__(self, substrate, vg, cc, f1params, f2params=[], useQs = True, timedep=(False, 0.)):#, f2params, f3params
        '''
        Init: substrate is the organism's Biomass
        vg = growth rate
        cc = carrying capacity
        params are list of ntuples
        f1 = function multiplying the growth rate, will depend on positive substances, list of ntuples
        f2 = function entering the inhibition factor, will depend on other organisms and carrying capacity
        f3 = nothing
        '''
        KineticLaw.__init__(self, substrate, upperBound=vg, useQs=useQs)
        self.f1params = f1params
        self.f2params = f2params
        #self.f3params = f3params
        self.vg = vg
        self.cc = cc
        self.timedep = timedep

    def computeBounds(self, val=-1, time=None, rSign=(-1., 1.)):
        '''
        Call computeRate and return lowerBound and upperBound
        '''
        r = self.computeRate(val, time)
        lb = -1*r
        ub = r
        #print 'bounds growth: ', lb, ub
        return lb, ub

    def carryCapacityRate(self, val=-1, time=None):
        '''
        todo: make more general
        params is an ntuple with metabolite and carry capacity
        '''
        cc = self.cc
        tf = 1.
        if self.timedep[0] and time is not None:
            #tf = time/(time + self.timedep[1])
            tf = np.log(time[val]/self.timedep[1] + 1)
        rate = tf*self.substrate.quantity[val]/cc
        return rate

    def multiMMrate(self, plist, val=-1):
        '''
        Compute the reaction rate according to 1st order Henry-Michaelis-Menten equation
        params is a list of tples (conc, mmk, maxv=1, hic=1)
        '''
        conc = 0.
        concscal = 0.
        mmk = 0.
        hic = 1.
        for params in plist:
            maxv = 1.
            if len(params) > 2:
                maxv = params[2]
            concscal += maxv*params[0].quantity[-1]
            conc += params[0].quantity[val]
            mmk += params[1]
            DEBUGMMM = False
            if DEBUGMMM:
                print params[0].modName
                print len(params[0].quantity)
                print '\tC: ', conc
                print '\tC scaled: ', concscal
                print '\tVm: ', maxv
                print '\tKmm: ', mmk
            #if len(params) > 3:
            #    hic = hic*params[3]
        rate = (concscal**hic)/(conc**hic + mmk)
        if DEBUGMMM:
            print '\tRATE: ', rate
        return rate

    def computeRate(self, val=-1, time=None):
        f2 = 1.
        if not self.useQs:
            #f2 = self.substrate.concentrationScaled[val]/self.substrate.quantity[val]
            f2 = self.substrate.concentration[val]/self.substrate.quantity[val]
        f1 = 1.
        for i in range(len(self.f1params)):
            if len(self.f2params) == len(self.f1params):
                f1 = f1*self.multiMMrate([self.f1params[i], self.f2params[i]])
            else:
                f1 = f1*self.mmrate(self.f1params[i])
        #f2 = self.carryCapacityRate(val, time)
        f2 = f2*self.substrate.quantity[val]/self.cc
        #f3 = self.timeDependence()
        rate = self.vg*f1*(1 - f2)
        #print 'growth: ', rate
        return rate



class SquareWave(KineticLaw):
    '''
    Square wave impulse
    '''
    def __init__(self, amplitude, period, bias=1., delay=0.):
        '''
        Init:
        '''
        KineticLaw.__init__(self)
        self.amplitude = amplitude
        self.period = period
        self.bias = bias
        self.delay = delay

    def computeRate(self, val=-1, time=None):
        rate = 0.
        if time is not None:
            tt = (time[val] - self.delay)*2*np.pi/self.period
            rate = 0.5*self.amplitude*(self.bias + np.sign(np.sin(tt)))
        return rate

class cambridgeCC(KineticLaw):
    '''
    Kinetci Law like in Grant et al, 2014
    '''
    def __init__(self):
        '''
        Init with a substrate and a maximal rate defaulted to 1000.
        '''
        KineticLaw.__init__(self)
        #self.

    def computeRate(self, val=-1, time=None):
        '''
        Return rate zero if a substrate is defined and its concentration at val is lower than 1e-3
        Else Return the maximal rate
        '''
        rate = 1 - 1
        return self.upperBound


class DReaction(object):
    '''
    General reaction for dynamic FBA
    The only mandatory value to be initiated is the name of the reaction in the SBML model
    '''
    def __init__(self, modName, kinF=None, isIrrev=None, enzyme=None, cofactors=None, fluxUnits=None, reactionSign=(-1.,1.), isODE=False):
        '''
        Init with the name of the reaction in the SBML model
        Optionally pass a kinetic law, if the reaction is irreversible, the enzyme (a DMetabolite),
        the cofactor (a DMetabolite) and reaction velocity units (mmol/(gDW*hr) by default)
        TODO: Allow more cofactors?
        '''
        self.modName = modName
        self.kinF = kinF if kinF is not None else None
        self.isIrrev = isIrrev if isIrrev is not None else True
        self.enzyme = enzyme if enzyme is not None else 'ND'
        self.cofactors = cofactors if cofactors is not None else []
        self.fluxUnits = fluxUnits if fluxUnits is not None else 'mmol/(gDW*hr)'
        self.flux = []
        self.upperBound = 1000.0
        self.lowerBound = 0 if self.isIrrev else -1000.0
        self.reactionSign = reactionSign
        self.isODE = isODE
        #? add smth else?
        self.ub = []
        self.lb = []

    def updateUpperBound(self, val=None):
        '''
        For general reactions simply set the value directly
        '''
        if self.kinF is None:
            raise NoKineticLawException
        else:
            self.upperBound = self.kinF.computeRate(val)
        return

    def updateLowerBound(self, val=None):
        '''
        For general reactions simply set the value directly
        '''
        if self.kinF is None:
            raise NoKineticLawException
        elif self.isIrrev or self.kinF.setZeroLowerBound():
            self.lowerBound = 0.
        else:
            self.lowerBound = -1*self.kinF.computeRate(val)
        return

    def appendBounds(self):
        '''
        append ub and lb to self lists
        '''
        self.ub.append(self.upperBound)
        self.lb.append(self.lowerBound)
        return

    def updateBounds(self, valU=None, valL=None, time=None):
        '''
        Update upper and lower bounds
        Priority given to direct upper and lower bounds values given (valU and valL respectively)
        If not given, check existence of a kinetic law (if not raise NoKineticLawException)
        and compute rate according to the kinetic law
        TODO: not sure about symmetric boundaries! should lowerBound always be -upperBound?
        '''
        done = False
        # print self.modName, type(self)
        if valU is not None:
            self.upperBound = valU
            done = True
        if valL is not None:
            self.lowerBound = valL
            done = True
        if done:
            #self.appendBounds(self.upperBound, self.lowerBound)
            return

        if self.kinF is None:
            raise NoKineticLawException
            return

        self.lowerBound, self.upperBound = self.kinF.computeBounds(time=time, rSign=self.reactionSign)
        # if self.setFlux:
            # print 'XXX ', self.modName, self.lowerBound, self.upperBound
        return

    def updateBoundsV1(self, valU=None, valL=None):
        '''
        Previous version
        TODO: remove?
        '''
        self.updateUpperBound(valU)
        if self.isIrrev:
            return
        else:
            self.updateLowerBound(valL)
        return

    def appendBoundsV1(self, ub, lb):
        '''
        append ub and lb to self lists
        TODO: remove?
        '''
        self.ub.append(ub)
        self.lb.append(lb)
        return


class FixedRateReaction(DReaction):
    '''
    Simple fixed rate reaction, inherits from DReaction
    '''
    def __init__(self, modName, rate=1.):
        '''
        Init with the name in SBML model only
        '''
        DReaction.__init__(self, modName, FixedBound(lowerBound=rate, upperBound=rate))

    def updateUpperBound(self, val=None):
        '''
        Do nothing
        '''
        return

    def updateLowerBound(self, val=None):
        '''
        Do nothing
        '''
        return

    def updateBounds(self, val=None, time=None):
        '''
        Do nothing
        '''
        self.lowerBound, self.upperBound = self.kinF.computeBounds(time=time)
        return

    
class CBModelReaction(DReaction):
    '''
    Reaction of the COBRA model, inherits from DReaction
    No bound updates of any kind
    '''
    def __init__(self, modName):
        '''
        Init with the name in SBML model only
        '''
        DReaction.__init__(self, modName)

    def updateUpperBound(self, val=None):
        '''
        Do nothing
        '''
        return

    def updateLowerBound(self, val=None):
        '''
        Do nothing
        '''
        return

    def updateBounds(self, val=None, time=None):
        '''
        Do nothing
        '''
        return
