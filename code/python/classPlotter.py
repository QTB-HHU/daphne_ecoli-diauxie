#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  last modified: 2015/08/13     **
#************************************

import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rcParams['legend.numpoints'] = 1

class Plotter(object):
    '''
    Object that draws plots from a DynamicModel
    '''
    def __init__(self, model, outpath='/tmp/'):
        '''
        model is a dictionary of models
        '''
        self.model = model
        self.modelKeys = self.model.keys()
        self.outpath = outpath
        self.colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k')
        self.linestyles = ('-', '--', ':', '-.')
        self.pointstyles = ('.', 's', 'd', '^', 'v')
        self.plstyle = [y+z for y in self.linestyles for z in self.colors]+['.'+z for z in self.colors]+['*'+z for z in self.colors]+['D'+z for z in self.colors]+['+'+z for z in self.colors]
        self.ppstyle = ['.'+z for z in self.colors]
        self.fluxstyles = {'flux': '.k', 'ub': '-b', 'lb': '-g'}

    def plotHistogramTSteps(self):
        '''
        Plot the distribution of time steps
        '''
        for k in self.modelKeys:
            tarr = np.array(self.model[k].T)
            tsteps = np.diff(tarr)
            ti = np.amin(tarr)
            tf = np.amax(tarr)
            dtmin = np.amin(tsteps)
            dtmax = np.amax(tsteps)
            afig, af = self.createPlot()
            plt.hist(tsteps, label=('Tot steps: %d'%len(tarr)))
            txtstr = 'T interval = (%.2f,%.2f)\nMin $\Delta$T = %.3f\nMax $\Delta$T int = %.3f' % (ti, tf, dtmin, dtmax)
            af.text(0.45, 0.95, txtstr, transform=af.transAxes, fontsize=12, verticalalignment='top')
            self.stylePlot(afig, af, k+'-time-steps-hist.png', 'Simulation Time Steps Distribution', ('Time Steps ('+self.model[k].Tunits+')', False), ('Count', False))
            self.closePlot(afig)
            if dtmax/dtmin > 100.:
                #newarr = tsteps/dtmin
                dtmean = np.mean(tsteps)
                delstr = '$\mu$=%.3f\nKeep $\Delta$T < 1.8*$\mu$\nDeleted elements: '%dtmean
                ntsteps = []
                for i in range(len(tsteps)):
                    if tsteps[i] > 1.8*dtmean:
                        delstr = delstr + ('\n$\Delta$T[%d] = %.2f @ T = %.2f -- %.2f' % (i, tsteps[i], tarr[i], tarr[i+1]))
                    else:
                        ntsteps.append(tsteps[i])
                afig, af = self.createPlot()
                plt.hist(np.array(ntsteps), label=('Tot steps: %d'%len(ntsteps)))
                af.text(0.05, 0.95, delstr, transform=af.transAxes, fontsize=10, verticalalignment='top')
                self.stylePlot(afig, af, k+'-time-steps-hist2.png', 'Simulation Time Steps Distribution w/o outliers', ('Time Steps ('+self.model[k].Tunits+')', False), ('Count', False))
                self.closePlot(afig)
        return

    def showStyles(self):
        for s in self.plstyle:
            print 'Style ', self.plstyle.index(s),': ', s
        return

    def plotAllMetabolitesQuantitiesVsTime(self, t='', logX=False, logY=False, doSinglePlots=False):
        for k in self.modelKeys:
            afig, af = self.createPlot()
            ytit = self.plotModelContent(afig, af, map(lambda x: (k, x, 1., 'Q'), self.model[k].dmetabolites.keys()), self.model[k].T, doSinglePlots)
            self.stylePlot(afig, af, k+'-'+t+'all_metabolites_quantities_vs_time.png', 'Metabolites Quantities Vs Time', ('Time ('+self.model[k].Tunits+')', logX), ('Quantity', logY))
            self.closePlot(afig)
        return

    def plotAllMetabolitesConcentrationsVsTime(self, t='', logX=False, logY=True, doSinglePlots=False):
        for k in self.modelKeys:
            afig, af = self.createPlot()
            ytit = self.plotModelContent(afig, af, map(lambda x: (k, x, 1., 'C'), self.model[k].dmetabolites.keys()), self.model[k].T, doSinglePlots)
            self.stylePlot(afig, af, k+'-'+t+'all_metabolites_concentrations_vs_time.png', 'Metabolites Concentrations Vs Time', ('Time ('+self.model[k].Tunits+')', logX), ('Concentration', logY))
            self.closePlot(afig)
        return

    def plotAllReactionsFluxesVsTime(self, t='', logX=False, logY=False, doSinglePlots=False):
        for k in self.modelKeys:
            afig, af = self.createPlot()
            ytit = self.plotModelContent(afig, af, map(lambda x: (k, x, 1., 'V'), self.model[k].dreactions.keys()), self.model[k].T, doSinglePlots, 28)
            self.stylePlot(afig, af, k+'-'+t+'all_reactions_fluxes_vs_time.png', 'Reactions Fluxes Vs Time', ('Time ('+self.model[k].Tunits+')', logX), ('Reaction Flux', logY))
            self.closePlot(afig)
        return

    def createPlot(self):
        '''
        create a plt.figure that is going to be updated
        '''
        fig0 = plt.figure()
        f0 = fig0.add_subplot(111)
        return fig0, f0

    def closePlot(self, pfig):
        pfig.clear()
        plt.close(pfig)
        return

    def stylePlot(self, pfig, pf, plotname='plot.png', plottitle='', xaxstyle=('', False), yaxstyle=('', True), xaxrange = (False, 0, 0), yaxrange = (False, 0, 0)):
        pf.set_title(plottitle)
        xtitle = xaxstyle[0]
        logX = xaxstyle[1]
        ytitle = yaxstyle[0]
        logY = yaxstyle[1]
        pf.set_xlabel(xtitle)
        pf.set_ylabel(ytitle)
        if xaxrange[0]:
            if xaxrange[2]-xaxrange[1] > 0:
                pf.set_xlim([xaxrange[1], xaxrange[2]])
        if yaxrange[0]:
            if yaxrange[2]-yaxrange[1] > 0:
                pf.set_ylim([yaxrange[1], yaxrange[2]])
            else:
                yl = pf.get_ylim()
                if logY:
                    if yl[0] > 0:
                        y0 = 0.01*yl[0]
                    else:
                        y0 = 0
                    y1 = 1.2*yl[1]
                else:
                    y0 = 0.9*yl[0]
                    y1 = 1.21*yl[1]
                pf.set_ylim(y0, y1)
        if logY:
            pf.set_yscale('log')
        if logX:
            pf.set_xscale('log')
        pf.legend(loc='best', prop={'size':10}, numpoints=1, fancybox=True, framealpha=0.5)
        #the following line sometimes throw the warning /usr/lib/python2.7/dist-packages/numpy/ma/core.py:3847 UserWarning: Warning: converting a masked element to nan.
        #probably related to axis issues - when getting ylim sometimes the values are strange...
        pfig.savefig(self.outpath+plotname)
        self.closePlot(pfig)
        return

    def updatePlot(self):
        '''
        gfig is an open plt.figure that is going to be updated
        '''
        return

    def plotArrays(self, pfig, pf, py, px, insty=0, plabel=''): #, plotname='plot.png', plottitle='', xaxstyle=('', False), yaxstyle=('', True), units=True, kinetics=True, xaxrange = (False, 0, 0), yaxrange = (False, 0, 0)):
        '''
        plot generic py vs px
        '''
        pf.plot(px, py, self.plstyle[insty], label=plabel)
        return

    def plotModelContent(self, pfig, pf, py, px=None, doSinglePlots=False, insty=0, kinetics=False, unitsOnLegend=True, interval=('full',0,0)): #, plotname='plot.png', plottitle='', xaxstyle=('', False), yaxstyle=('', True), units=True, kinetics=True, xaxrange = (False, 0, 0), yaxrange = (False, 0, 0)):
        '''
        py is a list of 4 element tuples:
            el[0] model
            el[1] dictionary name
            el[2] scale factor
            el[3] q or c or f ...
        '''
        if doSinglePlots:
            fig1 = plt.figure()
        plotBounds = False
        kmodel = list(set(map(lambda x: x[0], py)))
        style = insty
        for el in py:
            tmodel = self.model[el[0]]
            dname = el[1] #dictionary name
            sf = el[2] #scale factor
            variable = el[3] #q or c or f ...
            if px is None:
                xarray = np.array(tmodel.T)
            else:
                xarray = np.array(px)
            yarray = sf
            units = ''
            if variable.lower() == 'c':
                yarray = np.array(tmodel.dmetabolites[dname].concentration)
                units = ' ('+tmodel.dmetabolites[dname].concentrationUnits+')'
                ytitle = 'Concentration'+units
            elif variable.lower() == 'q':
                yarray = np.array(tmodel.dmetabolites[dname].quantity)
                units = ' ('+tmodel.dmetabolites[dname].quantityUnits+')'
                ytitle = 'Quantity'+units
            elif variable.lower() == 'ub':
                yarray = np.array(tmodel.dreactions[dname].ub)
                if max(yarray) > 999: #and logY
                    if sf == 1:
                        sf = 0.001
                units = ' ('+tmodel.dreactions[dname].fluxUnits+')'
                ytitle = 'U.B. on flux'+units
            elif variable.lower() == 'lb':
                yarray = np.array(tmodel.dreactions[dname].lb)
                if min(yarray) < -999: #and logY
                    if sf == 1:
                        sf = 0.001
                units = ' ('+tmodel.dreactions[dname].fluxUnits+')'
                ytitle = 'L.B. on flux'+units
            elif variable.lower() == 'f' or variable.lower() == 'v':
                yarray = np.array(tmodel.dreactions[dname].flux)
                if max(yarray) <= 0: #and logY
                    if sf > 0:
                        sf = -1*sf
                units = ' ('+tmodel.dreactions[dname].fluxUnits+')'
                ytitle = 'Flux'+units
                if doSinglePlots:
                    plotBounds = True
                if tmodel.dreactions[dname].kinF is None:
                    kinetics = False
                elif tmodel.dreactions[dname].kinF.substrate is None:
                    kinetics = False
                if kinetics:
                    ksty = style+7
                    yc = []
                    #for i in range(len(tmodel.dreactions[dname].kinF.substrate.concentration)):
                        #computeRate takes the index of the list!
                        #yc.append(tmodel.dreactions[dname].kinF.computeRate(i))
                    for i in tmodel.dreactions[dname].kinF.substrate.concentration:
                        yc.append(tmodel.dreactions[dname].kinF.computeRateFromC(i))
            plabel = ''
            #print len(yarray), len(xarray)
            if sf != 1:
                plabel = '%.3f * '%(sf)
            plabel = el[0]+': '+plabel+variable+' '+dname
            if unitsOnLegend:
                plabel = plabel+units
            if interval[0]=='first':
                xarray = xarray[:interval[1]]
                yarray = yarray[:interval[1]]
                if kinetics:
                    yc = yc[:interval[1]]
            elif interval[0]=='last':
                xarray = xarray[-interval[1]:]
                yarray = yarray[-interval[1]:]
                if kinetics:
                    yc = yc[-interval[1]:]
            pf.plot(xarray, yarray*sf, self.plstyle[style], label=plabel)
            if doSinglePlots:
                f1 = fig1.add_subplot(111)
                #f1.set_title(plottitle)
                if kinetics:
                    f1.plot(xarray, np.array(yc),  self.plstyle[ksty], label='boundary')
                    f1.plot(xarray, yarray, self.plstyle[style], label='solution')
                else:
                    f1.plot(xarray, yarray, self.plstyle[style], label=plabel)
                if plotBounds:
                    if max(tmodel.dreactions[dname].ub) < 1000.:
                        f1.plot(xarray, np.array(tmodel.dreactions[dname].ub), self.plstyle[15], label='U.B.')
                    if min(tmodel.dreactions[dname].lb) > -1000.:
                        f1.plot(xarray, np.array(tmodel.dreactions[dname].lb), self.plstyle[16], label='L.B.')
                #f1.set_xlabel(xtitle)
                f1.set_ylabel(ytitle)
                f1.legend(loc='best', prop={'size':10}, numpoints=1, fancybox=True, framealpha=0.5)
                fig1.savefig(self.outpath+variable+'_'+dname+'_'+el[0]+'.png')
                #clear fig1
                plt.clf()
            else:
                if kinetics:
                    pf.plot(xarray, np.array(yc),  self.plstyle[ksty], label='boundary')
            style += 1
        return ytitle

    def plotReactions(self, mod, rxns, xvar, plotname, timeint):
        '''
        rxns: list like [['rxnname1', ('flux')], ['rxnname2', ('flux', 'ub')], ...]
        xvar: list like [xarray, 'xlabel', 'xunits']
        '''
        #fluxes = dict(map(lambda x: (x, np.array(self.model[x].flux)), rxns))
        nax = len(rxns)
        xarr = xvar[0]
        fig, axarr = plt.subplots(nax, sharex=True)
        plt.suptitle(mod)
        axarr[-1].set_xlabel(xvar[1]+' ('+xvar[2]+')')
        axarr[nax-1].set_xlim(timeint)
        for i,r in enumerate(rxns):
            rid = r[0]
            toplot = r[1]
            drxn = self.model[mod].dreactions[rid]
            axarr[i].set_title(rid)
            for j,v in enumerate(toplot):
                yarr = np.array(getattr(drxn, v))
                axarr[i].plot(xarr, yarr, self.fluxstyles[v], label=v)
            axarr[i].set_ylabel('Rate ('+drxn.fluxUnits+')')
            if len(toplot) > 1:
                ll = axarr[i].legend(loc='best', prop={'size':10})
        fig.savefig(self.outpath+plotname)
        self.closePlot(fig)
        return

    def plotMetabolites(self, mod, mets, xvar, plotname, xaxint=None, yaxint=None, axorg=None, colorcodes=None):
        nax = len(mets)
        xarr = xvar[0]
        twinax = axorg[0]
        fig, axarr = plt.subplots(nax, sharex=True)
        twinaxarr = [None]*nax
        for a in range(nax):
            if twinax[a]:
                twinaxarr[a] = axarr[a].twinx()
        plt.suptitle(mod)
        axarr[-1].set_xlabel(xvar[1]+' ('+xvar[2]+')')
        if xaxint:
            axarr[nax-1].set_xlim(xaxint)
        titles = ['']*nax
        i = -1
        for axi,toplot in enumerate(mets):
            for j,t in enumerate(toplot):
                i += 1
                r = t.split('+')
                rid = r[0]
                if rid == 'HLINE':
                    self.plotHLine(twinaxarr[axi] if twinax[axi] else axarr[axi], r[1:])
                    continue
                v = r[1]
                vlab = v.replace('quantity', 'Q').replace('concentration', 'C')
                drxn = self.model[mod].dmetabolites[rid]
                jcolor = self.colors[i]
                if colorcodes:
                    jcolor = colorcodes[rid]
                if rid not in titles[axi]:
                    titles[axi] = titles[axi]+rid+', '
                yarr = np.array(getattr(drxn, v))
                jlab = rid
                if twinax[axi]:
                    jlab = rid+' '+vlab
                units = getattr(drxn, v+'Units')
                ylab = axarr[axi].get_ylabel()
                if ylab == '':
                    axarr[axi].set_ylabel(vlab+' ('+units+')')
                    ylab = axarr[axi].get_ylabel()
                if twinax[axi] and units not in ylab:
                    if twinaxarr[axi].get_ylabel() == '':
                        twinaxarr[axi].set_ylabel(vlab+' ('+units+')')
                    twinaxarr[axi].plot(xarr, yarr, jcolor, label=jlab)
                else:
                    axarr[axi].plot(xarr, yarr, jcolor, label=jlab)
        for axi in range(nax):
            if twinax[axi]:
                lines, labels = axarr[axi].get_legend_handles_labels()
                lines2, labels2 = twinaxarr[axi].get_legend_handles_labels()
                twinaxarr[axi].legend(lines + lines2, labels + labels2, loc='best', prop={'size':10})
                #ll = axarr[axi].legend(loc='center left', prop={'size':10})
                #ll = twinaxarr[axi].legend(loc='right', prop={'size':10})
            else:
                ll = axarr[axi].legend(loc='best', prop={'size':10})

        for i in range(nax):
            axarr[i].set_title(titles[i][:-2])
        fig.savefig(self.outpath+plotname)
        self.closePlot(fig)

        return

    def plotHLine(self, ax, argu):
        xy = eval(argu[0])
        ax.axhline(y=xy[1], color='k', linestyle=':')
        if len(argu) > 1:
            ax.text(xy[0], xy[1], argu[1], va='bottom')
        return

    def plotShareX(self, stit, mets, xvar, plotname, xaxint=None, yaxint=None, axorg=None, colorcodes=None):
        nax = len(mets)
        xarr = xvar[0]
        twinax = map(lambda x: len(x) > 1 and x[1], axorg)
        fig, axarr = plt.subplots(nax, sharex=True)
        twinaxarr = [None]*nax
        for a in range(nax):
            if twinax[a]:
                twinaxarr[a] = axarr[a].twinx()
        plt.suptitle(stit)
        axarr[-1].set_xlabel(xvar[1]+' ('+xvar[2]+')')
        if xaxint:
            axarr[nax-1].set_xlim(xaxint)
        titles = ['']*nax
        i = -1
        for axi,toplot in enumerate(mets):
            for j,t in enumerate(toplot):
                i += 1
                r = t.split('+')
                if len(r) < 4:
                    print 'Please give the correct input ',t
                    continue
                mod = r[0]
                att = r[1]
                rid = r[2]
                if rid == 'HLINE':
                    self.plotHLine(twinaxarr[axi] if twinax[axi] else axarr[axi], r[3:])
                    continue
                v = r[3]
                vlab = v.replace('quantity', 'Q').replace('concentration', 'C').replace('flux', 'V')
                vulab = v.replace('lb','flux').replace('ub','flux')+'Units'
                drxn = getattr(self.model[mod], att)[rid]
                yarr = np.array(getattr(drxn, v))
                if len(r) > 7:
                    print 'Array operations'
                    sign = r[4]
                    att2 = r[5]
                    rid2 = r[6]
                    v2 = r[7]
                    sf=1.
                    if len(r) > 8:
                        sf = float(r[8])
                    # print sf
                    if len(r) > 9:
                        rid = r[9]
                    drxn2 = getattr(self.model[mod], att2)[rid2]
                    yarr2 = sf*np.array(getattr(drxn2, v2))
                    if sign == 'minus':
                        # itime = (np.abs(xarr-0.8)).argmin()
                        # print yarr[itime:itime+10]
                        yarr = yarr - yarr2
                        # print yarr[itime:itime+10]
                    elif sign == 'plus':
                        yarr = yarr + yarr2
                    elif sign == 'times':
                        yarr = yarr * yarr2
                    elif sign == 'div':
                        yarr = yarr / yarr2
                    else:
                        print 'invalid operation'
                jcolor = self.colors[i]
                if colorcodes:
                    jcolor = colorcodes.get(rid+vlab, 'k.')
                if rid not in titles[axi]:
                    titles[axi] = titles[axi]+rid+', '
                jlab = rid
                if twinax[axi] or 'flux' in vulab:
                    jlab = rid+' '+vlab
                units = getattr(drxn, vulab)
                ontwin = False
                ylab = axarr[axi].get_ylabel()
                if ylab == '':
                    axarr[axi].set_ylabel(vlab+' ('+units+')')
                    ylab = axarr[axi].get_ylabel()
                if vlab not in ylab or (len(r) > 4 and r[4]=='TWIN'):
                    ontwin = True
                elif units not in ylab:
                    axarr[axi].set_ylabel(vlab+' (a.u.)')
                    ylab = axarr[axi].get_ylabel()
                if twinax[axi] and ontwin:
                    twinylab = twinaxarr[axi].get_ylabel()
                    if twinylab == '':
                        twinaxarr[axi].set_ylabel(vlab+' ('+units+')')
                    if units not in twinylab:
                        twinaxarr[axi].set_ylabel(vlab+' (a.u.)')
                        twinylab = twinaxarr[axi].get_ylabel()
                    twinaxarr[axi].plot(xarr, yarr, jcolor, label=jlab)
                else:
                    axarr[axi].plot(xarr, yarr, jcolor, label=jlab)
        for axi in range(nax):
            legendloc = 'best'
            if len(axorg[axi]) > 2:
                legendloc = axorg[axi][2]
            if axorg[axi][0]:
                axarr[axi].set_ylim(axorg[axi][0])
            if len(axorg[axi]) > 3:
                if 'ylog' in axorg[axi][3]:
                    axarr[axi].set_yscale('log')
            if twinax[axi]:
                if axorg[axi][1]:
                    twinaxarr[axi].set_ylim(axorg[axi][1])
                if len(axorg[axi]) > 4:
                    if 'ylog' in axorg[axi][4]:
                        twinaxarr[axi].set_yscale('log')
                lines, labels = axarr[axi].get_legend_handles_labels()
                lines2, labels2 = twinaxarr[axi].get_legend_handles_labels()
                twinaxarr[axi].legend(lines + lines2, labels + labels2, loc=legendloc, prop={'size':10})
                #ll = axarr[axi].legend(loc='center left', prop={'size':10})
                #ll = twinaxarr[axi].legend(loc='right', prop={'size':10})
            else:
                ll = axarr[axi].legend(loc=legendloc, prop={'size':10})

        for i in range(nax):
            axarr[i].set_title(titles[i][:-2])
        fig.savefig(self.outpath+plotname)
        self.closePlot(fig)

        return

def plotQuantities():
    return
