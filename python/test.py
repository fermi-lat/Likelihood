#!/usr/bin/env python

import os, sys
optimizersRoot = os.getenv("OPTIMIZERSROOT")
latResponseRoot = os.getenv("LATRESPONSEROOT")
LikelihoodRoot = os.getenv("LIKELIHOODROOT")
sys.path.append(optimizersRoot + "/python")
sys.path.append(latResponseRoot + "/python")
sys.path.append(LikelihoodRoot + "/python")

import optimizers, latResponse, Likelihood

caldbPath = LikelihoodRoot + "/src/test/CALDB"

ra0 = 193.98
dec0 = -5.82
#Likelihood.RoiCuts_setCuts(ra0, dec0, 20.)

Likelihood.RoiCuts_setCuts(LikelihoodRoot+"/xml/RoiCuts.xml")
#ra0 = 0.; dec0 = 0.
#Likelihood.RoiCuts_getRaDec(ra0, dec0)
#print ra0, dec0

roiCuts = Likelihood.RoiCuts_instance()
#roiCuts.getRaDec( ra0, dec0 )
#print ra0, dec0
print roiCuts.extractionRegion().radius()

obs_root = "diffuse_test_5"

sc_file = LikelihoodRoot + "/src/test/Data/" + obs_root + "_sc_0000"
sc_hdu = 2
Likelihood.ScData_readData( sc_file, sc_hdu )

expFile = LikelihoodRoot + "/src/test/Data/exp_" + obs_root + "_new.fits"
Likelihood.ExposureMap_readExposureFile( expFile )

respFunctions = Likelihood.ResponseFunctions_instance()
respFunctions.addGlast25Resp(LikelihoodRoot + "/src/test/CALDB/", 4)

def make_SourceFactory():
    funcFactory = optimizers.FunctionFactory()

    funcFactory.addFunc("PowerLaw", optimizers.PowerLaw())
    funcFactory.addFunc("Gaussian", optimizers.Gaussian())
    funcFactory.addFunc("AbsEdge", optimizers.AbsEdge())

    funcFactory.addFunc("SkyDirFunction", Likelihood.SkyDirFunction())
    funcFactory.addFunc("ConstantValue", Likelihood.ConstantValue())
    funcFactory.addFunc("SpatialMap", Likelihood.SpatialMap())

    srcFactory = Likelihood.SourceFactory()
    
    xmlFile = LikelihoodRoot + "/xml/A1_Sources.xml"

    srcFactory.readXml( xmlFile, funcFactory )

    return srcFactory

def fitStatistic():

    srcFactory = make_SourceFactory()

    ourGalaxy = srcFactory.create("Galactic Diffuse Emission")
    extragalactic = srcFactory.create("Extragalactic Diffuse Emission")
    _3c279 = srcFactory.create("Bright Point Source")
    _3c279.setDir(ra0, dec0)
    _3c279.setName("3C 279")

    logLike = Likelihood.LogLike()

    logLike.addSource(ourGalaxy)
    logLike.addSource(extragalactic)
    logLike.addSource(_3c279)

    eventFile = LikelihoodRoot + "/src/test/Data/" + obs_root + "_0000"

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
