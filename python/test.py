#!/usr/bin/env python
import sys
from load_data import *

def makeFactories():
    srcFactory = Likelihood.SourceFactory()
    xmlFile = LikelihoodRoot + "/xml/A1_Sources.xml"
    funcFactory = Likelihood.SourceFactory_funcFactory()
    srcFactory.readXml(xmlFile, funcFactory )
    return srcFactory, funcFactory

def fitStatistic(binned=False):

#    srcFactory, funcFactory = makeFactories()
#    ourGalaxy = srcFactory.create("Galactic Diffuse Emission")
#    extragalactic = srcFactory.create("Extragalactic Diffuse Emission")
#    Crab = srcFactory.create("Bright Point Source")
#    Crab.setDir(ra0, dec0, 1, 0)
#    Crab.setName("Crab")

    funcFactory = Likelihood.SourceFactory_funcFactory()
    srcModel = "ptsrcModel.xml"

    if not binned:
        logLike = Likelihood.LogLike()
        eventFile = "test_events_0000.fits"
        logLike.getEvents(eventFile)
        logLike.readXml(srcModel, funcFactory)
        logLike.computeEventResponses()
    else:
        srcMapsFile = "sourceMaps.fits"
        dataMap = Likelihood.CountsMap(srcMapsFile)
        logLike = Likelihood.BinnedLikelihood(dataMap, srcMapsFile)
#        logLike = Likelihood.BinnedLikelihood_create(srcMapsFile)
        logLike.readXml(srcModel, funcFactory)

    print logLike.value()

#    myOpt = logLike.Drmngb()
#    myOpt = logLike.Lbfgs()
    myOpt = logLike.Minuit()

    verbose = 3
    tol = 1e-5
    myOpt.find_min(verbose, tol)

    return (logLike, myOpt)

if __name__ == "__main__":
    binned = "-b" in sys.argv
    (logLike, myOpt) = fitStatistic(binned)
    logLike.print_source_params()
