from load_data import *

def make_SourceFactory():
#    funcFactory = optimizers.FunctionFactory()
#
#    funcFactory.addFunc("SkyDirFunction", Likelihood.SkyDirFunction())
#    funcFactory.addFunc("SpatialMap", Likelihood.SpatialMap())

    srcFactory = Likelihood.SourceFactory()
    
    xmlFile = LikelihoodRoot + "/xml/A1_Sources.xml"

    funcFactory = Likelihood.SourceFactory.funcFactory()
    srcFactory.readXml( xmlFile, funcFactory )

    return srcFactory

def fitStatistic():

#    srcFactory = make_SourceFactory()
#
#    ourGalaxy = srcFactory.create("Galactic Diffuse Emission")
#    extragalactic = srcFactory.create("Extragalactic Diffuse Emission")
#    Crab = srcFactory.create("Bright Point Source")
#    Crab.setDir(ra0, dec0, 1, 0)
#    Crab.setName("Crab")
#
#    funcFactory = optimizers.FunctionFactory()
#
#    funcFactory.addFunc("SkyDirFunction", Likelihood.SkyDirFunction())
#    funcFactory.addFunc("SpatialMap", Likelihood.SpatialMap())

    srcFactory = Likelihood.SourceFactory()
    funcFactory = Likelihood.SourceFactory.funcFactory(srcFactory)
    
#    logLike = Likelihood.LogLike()
    srcMapsFile = "sourceMaps.fits"
    logLike = Likelihood.BinnedLikelihood_create(srcMapsFile)
    srcModel = "anticenter_model.xml"
    logLike.readXml(srcModel, funcFactory)

    eventFile = "test_events_0000.fits"

    logLike.getEvents(eventFile)

    print logLike.value()

#    logLike.computeEventResponses()

#    myOpt = logLike.Drmngb()
#    myOpt = logLike.Lbfgs()
    myOpt = logLike.Minuit()

    verbose = 3
    tol = 1e-5
    myOpt.find_min(verbose, tol)
#
    return (logLike, myOpt)

if __name__ == "__main__":
    (logLike, myOpt) = fitStatistic()
    logLike.print_source_params()
