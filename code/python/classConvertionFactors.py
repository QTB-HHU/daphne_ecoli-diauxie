#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  last modified: 2015/08/13     **
#************************************

class ConvertionFactors(object):
    '''
    Class to help managing units
    '''
    def __init__(self):
        '''
        Init by creating the convertion factors dictionary
        '''
        self.cf = self.initConvertionFactors()

    def initConvertionFactors(self):
        '''
        Return the convertion factors dictionary
        '''
        cf = {'gDW-mL': 5., 'hr-min': 60., 'hr-sec': 3600., 'min-sec': 60., 'mL-L': 0.001}

        cf['gDW-L'] = 0.001*cf['gDW-mL']
        cf['gDW-muL'] = 1000.*cf['gDW-mL']
        cf['gDW-uL'] = 1000.*cf['gDW-mL']
        cf['gDW-nL'] = 1000.*cf['gDW-uL']
        cf['gDW-pL'] = 1000.*cf['gDW-nL']

        cf['mL-gDW'] = 1./cf['gDW-mL']
        cf['muL-gDW'] = 1./cf['gDW-muL']
        cf['uL-gDW'] = 1./cf['gDW-uL']
        cf['nL-gDW'] = 1./cf['gDW-nL']
        cf['pL-gDW'] = 1./cf['gDW-pL']
        cf['L-gDW'] = 1./cf['gDW-L']

        cf['min-hr'] = 1./cf['hr-min']
        cf['sec-min'] = 1./cf['min-sec']
        cf['sec-hr'] = 1./cf['hr-sec']

        return cf

    def AtoB(self, A, B):
        '''
        Return the factor to convert A to B
        X [B u.] = cf * X [A u.]
        '''
        return self.cf[A+'-'+B]
