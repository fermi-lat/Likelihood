"""
Interface to SWIG-wrapped C++ classes.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
#
# $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/python/SrcAnalysis.py,v 1.15 2005/03/03 00:46:53 jchiang Exp $
#
import os
import numarray as num
import pyLike
from SrcModel import SourceModel

_funcFactory = pyLike.SourceFactory_funcFactory()

class SrcAnalysis(object):
    def __init__(self, srcModel, eventFile, scFile, expMap=None, irfs='TEST',
                 optimizer='Minuit'):
        self.optimizer = optimizer
        self.respFuncs = pyLike.ResponseFunctions()
        self.respFuncs.load(irfs)
        self.expMap = pyLike.ExposureMap()
        if expMap is not None:
            self.expMap.readExposureFile(expMap)
        self.scData = pyLike.ScData()
        self.expCube = pyLike.ExposureCube()
        observation = pyLike.Observation(self.respFuncs,
                                         self.scData,
                                         pyLike.RoiCuts_instance(),
                                         self.expCube,
                                         self.expMap)
        self.observation = observation
        self.logLike = pyLike.LogLike(self.observation)
        self._readData(scFile, eventFile)
        self.events = self.logLike.events();
        self.logLike.readXml(srcModel, _funcFactory)
        self.logLike.computeEventResponses()
        self.model = SourceModel(self.logLike)
        eMin, eMax = pyLike.RoiCuts_instance().getEnergyCuts()
        nee = 21
        estep = num.log(eMax/eMin)/(nee-1)
        self.energies = eMin*num.exp(estep*num.arange(nee, type=num.Float))
        self.disp = None
        self.resids = None
    def __call__(self):
        return -self.logLike.value()
    def _fileList(self, files):
        if isinstance(files, str):
            return [files]
        else:
            return files
    def _readData(self, scFile, eventFile):
        self._readScData(scFile)
        self._readEvents(eventFile)
    def _readScData(self, scFile):
        scFiles = self._fileList(scFile)
        self.scData.readData(scFiles[0], True)
        for file in scFiles[1:]:
            self.scData.readData(scFile)
    def _readEvents(self, eventFile):
        eventFiles = self._fileList(eventFile)
        pyLike.RoiCuts_instance().readCuts(eventFiles[0])
        for file in eventFiles:
            self.logLike.getEvents(file)
    def _Nobs(self, emin, emax):
        nobs = 0
        for event in self.events:
            if emin < event.getEnergy() < emax:
                nobs += 1
        return nobs
    def plotData(self, yrange=None):
        import hippoplotter as plot
        nt = plot.newNTuple(([], [], []),
                            ('energy', 'nobs', 'nobs_err'))
        self.data_nt = nt
        for emin, emax in zip(self.energies[:-1], self.energies[1:]):
            nobs = self._Nobs(emin, emax)
            nt.addRow((num.sqrt(emin*emax), nobs, num.sqrt(nobs)))
        self.disp = plot.XYPlot(nt, 'energy', 'nobs', yerr='nobs_err',
                                xlog=1, ylog=1)
        self.disp.getDataRep().setSymbol('filled_square', 2)
        if yrange is None:
            yrange = (0.5, max(nt.getColumn('nobs'))*1.5)
        self.disp.setRange('y', yrange[0], yrange[1])
    def _srcCnts(self, source):
        cnts = []
        for emin, emax in zip(self.energies[:-1], self.energies[1:]):
            cnts.append(source.Npred(emin, emax))
        return num.array(cnts)
    def _plot_model(self, source, yrange=None, lineStyle="Dot", oplot=False,
                    color='black'):
        import hippoplotter as plot
        if oplot and self.disp is not None:
            plot.canvas.selectDisplay(self.disp)
        else:
            self.plotData(yrange)
            
        self.col = lambda x: num.array(self.data_nt.getColumn(x))
        
        energies = self.col('energy')
        try:
            model = self._srcCnts(self.logLike.getSource(source))
        except:
            model = source
        plot.scatter(energies, model, oplot=True, pointRep='Line',
                     lineStyle=lineStyle, color=color)
        return model
    def _plot_residuals(self, model, oplot=None, color='black'):
        import hippoplotter as plot
        resid = (self.col('nobs') - model)/model
        resid_err = self.col('nobs_err')/model

        energies = self.col('energy')
        nt = plot.newNTuple((energies, resid, resid_err),
                            ('energy', 'residuals', 'resid_err'))
        if oplot and self.resids is not None:
            plot.canvas.selectDisplay(self.resids)
            rep = plot.XYPlot(nt, 'energy', 'residuals', yerr='resid_err',
                              xlog=1, oplot=1, color=color)
            rep.setSymbol('filled_square', 2)
        else:
            self.resids = plot.XYPlot(nt, 'energy', 'residuals',
                                      yerr='resid_err', xlog=1, color=color,
                                      yrange=(-1, 1))
            self.resids.getDataRep().setSymbol('filled_square', 2)
    def plot(self, srcs=None, oplot=False, yrange=None, color='black'):
        import hippoplotter as plot
        if oplot:
            color = 'red'
        if isinstance(srcs, str):
            total = self._plot_model(srcs, yrange=yrange, color=color, 
                                     oplot=oplot, lineStyle='Solid')
        else:
            if srcs is None:
                srcs = pyLike.StringVector()
                self.logLike.getSrcNames(srcs)
            total = self._plot_model(srcs[0], yrange=yrange, color=color,
                                     oplot=oplot)
            if len(srcs) > 1:
                for src in list(srcs[1:]):
                    total += self._plot_model(src, oplot=True, color=color)
            self._plot_model(total, color=color, oplot=True, lineStyle='Solid')
        self._plot_residuals(total, oplot=oplot, color=color)
    def fit(self, verbosity=3, tol=1e-5, optimizer=None):
        errors = self._errors(optimizer, verbosity, tol)
        return -self.logLike.value()
    def _errors(self, optimizer=None, verbosity=0, tol=1e-5):
        self.logLike.syncParams()
        if optimizer is None:
            optimizer = self.optimizer
        myOpt = eval("self.logLike.%s()" % optimizer)
        myOpt.find_min(verbosity, tol)
        errors = myOpt.getUncertainty()
        j = 0
        for i in range(len(self.model.params)):
            if self.model[i].isFree():
                self.model[i].setError(errors[j])
                j += 1
        return errors

if __name__ == '__main__':
    srcAnalysis = SrcAnalysis('galdiffuse_model.xml',
                              'galdiffuse_events_0000.fits',
                              'galdiffuse_scData_0000.fits',
                              'expMap_test.fits')
    srcAnalysis.plot('Galactic Diffuse')
