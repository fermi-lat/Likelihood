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
    
    logLike = Likelihood.LogLike()
    srcModel = LikelihoodRoot + "/data/anticenter_model.xml"
    logLike.readXml(srcModel, funcFactory)

    eventFile = LikelihoodRoot + "/data/single_src_events_0000.fits"

    logLike.getEvents(eventFile, 2)

#    print logLike.value()

#    logLike.computeEventResponses()

#    myOpt = optimizers.Lbfgs(logLike)
#    myOpt = optimizers.Minuit(logLike)

#    verbose = 3
#    myOpt.find_min(verbose)

#    return (logLike, myOpt)
    return (logLike, 1)

if __name__ == "__main__":
    (logLike, myOpt) = fitStatistic()
    logLike.print_source_params()
