#!/usr/bin/env python

import os, sys
optimizersRoot = os.getenv("OPTIMIZERSROOT")
LikelihoodRoot = os.getenv("LIKELIHOODROOT")
sys.path.append(optimizersRoot + "/python")

import optimizers, Likelihood

caldbPath = LikelihoodRoot + "/src/test/CALDB"

string = optimizers.Function.string
vector = optimizers.DoubleVector

ra0 = 193.98
dec0 = -5.82
Likelihood.RoiCuts_setCuts(ra0, dec0, 20.)

obs_root = "diffuse_test_5"

sc_file = LikelihoodRoot + "/src/test/Data/" + obs_root + "_sc_0000"
sc_hdu = 2
Likelihood.ScData_readData( string(sc_file), sc_hdu )

expFile = LikelihoodRoot + "/src/test/Data/exp_" + obs_root + "_new.fits"
Likelihood.ExposureMap_readExposureFile( string(expFile) )

psf = Likelihood.Psf_instance()
psfFile = LikelihoodRoot + "/src/test/CALDB/psf_lat.fits"
psf.readPsfData( string(psfFile), Likelihood.Response.Combined )

aeff = Likelihood.Aeff_instance()
aeffFile = LikelihoodRoot + "/src/test/CALDB/aeff_lat.fits"
aeff.readAeffData( string(aeffFile), Likelihood.Response.Combined )

def make_SourceFactory():
    funcFactory = optimizers.FunctionFactory()

    funcFactory.addFunc(string("PowerLaw"), optimizers.PowerLaw())
    funcFactory.addFunc(string("Gaussian"), optimizers.Gaussian())
    funcFactory.addFunc(string("AbsEdge"), optimizers.AbsEdge())

    funcFactory.addFunc(string("SkyDirFunction"), Likelihood.SkyDirFunction())
    funcFactory.addFunc(string("ConstantValue"), Likelihood.ConstantValue())
    funcFactory.addFunc(string("SpatialMap"), Likelihood.SpatialMap())

    srcFactory = Likelihood.SourceFactory()
    
    xmlFile = LikelihoodRoot + "/xml/A1_Sources.xml"

    srcFactory.readXml( string(xmlFile), funcFactory )

    return srcFactory

def fitStatistic():

    srcFactory = make_SourceFactory()

    ourGalaxy = srcFactory.create("Galactic Diffuse Emission")
    extragalactic = srcFactory.create("Extragalactic Diffuse Emission")
    _3c279 = srcFactory.create("Bright Point Source")
    _3c279.setDir(ra0, dec0)
    _3c279.setName("3C 279")

    logLike = Likelihood.logLike_ptsrc()

    logLike.addSource(ourGalaxy)
    logLike.addSource(extragalactic)
    logLike.addSource(_3c279)

    eventFile = LikelihoodRoot + "/src/test/Data/" + obs_root + "_0000"

    logLike.getEvents(eventFile, 2)

    logLike.computeEventResponses()

    myOpt = optimizers.Minuit(logLike)

    verbose = 3
    myOpt.find_min(verbose)

    return (logLike, myOpt)

if __name__ == "__main__":
    (logLike, myOpt) = fitStatistic()
    logLike.print_source_params()
