#!/usr/bin/env python

import os, sys
optimizersRoot = os.getenv("OPTIMIZERSROOT")
latResponseRoot = os.getenv("LATRESPONSEROOT")
LikelihoodRoot = os.getenv("LIKELIHOODROOT")
sys.path.append(optimizersRoot + "/python")
sys.path.append(latResponseRoot + "/python")
sys.path.append(LikelihoodRoot + "/python")

import optimizers, latResponse, Likelihood

irfsFactory = latResponse.irfsFactory()

ra0 = 83.57
dec0 = 22.01

Likelihood.RoiCuts_setCuts(LikelihoodRoot+"/data/RoiCuts.xml")

roiCuts = Likelihood.RoiCuts_instance()

sc_file = LikelihoodRoot + "/data/single_src_scData_0000.fits"
sc_hdu = 2
Likelihood.ScData_readData( sc_file, sc_hdu )

Likelihood.ResponseFunctions_addRespPtr(0, irfsFactory.create('DC1::Front'))
Likelihood.ResponseFunctions_addRespPtr(1, irfsFactory.create('DC1::Back'))

#expFile = LikelihoodRoot + "/src/test/Data/exp_" + obs_root + "_new.fits"
#Likelihood.ExposureMap_readExposureFile( expFile )

def make_SourceFactory():
    funcFactory = optimizers.FunctionFactory()

    funcFactory.addFunc("SkyDirFunction", Likelihood.SkyDirFunction())
    funcFactory.addFunc("SpatialMap", Likelihood.SpatialMap())

    srcFactory = Likelihood.SourceFactory()
    
    xmlFile = LikelihoodRoot + "/xml/A1_Sources.xml"

    srcFactory.readXml( xmlFile, funcFactory )

    return srcFactory

def fitStatistic():

    srcFactory = make_SourceFactory()

#    ourGalaxy = srcFactory.create("Galactic Diffuse Emission")
#    extragalactic = srcFactory.create("Extragalactic Diffuse Emission")
    Crab = srcFactory.create("Bright Point Source")
    Crab.setDir(ra0, dec0, 0)
    Crab.setName("Crab")

    logLike = Likelihood.LogLike()

#    logLike.addSource(ourGalaxy)
#    logLike.addSource(extragalactic)
    logLike.addSource(Crab)

    eventFile = LikelihoodRoot + "/data/single_src_events_0000.fits"

    logLike.getEvents(eventFile, 2)

    logLike.computeEventResponses()

#    myOpt = optimizers.Lbfgs(logLike)
    myOpt = optimizers.Minuit(logLike)

    verbose = 3
    myOpt.find_min(verbose)

    return (logLike, myOpt)

if __name__ == "__main__":
    (logLike, myOpt) = fitStatistic()
    logLike.print_source_params()
