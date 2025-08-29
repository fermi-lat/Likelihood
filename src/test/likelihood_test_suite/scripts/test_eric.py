import os
import pyLikelihood as pyLike
import BinnedAnalysis as ba

NOMINAL_MEMORY = 88
NOMINAL_TIME = 0.65

def test_eric(dirs):
    # test data is in dirs['LOCAL_DATA']

    _funcFactory = pyLike.SourceFactory_funcFactory()


    kw = dict(srcMaps=dirs['LOCAL_DATA'] + 'srcmap_00.fits',
            expCube=dirs['LOCAL_DATA'] + 'ltcube_239557414_428903014_z100_r180_gti.fits',
            binnedExpMap=dirs['LOCAL_DATA'] + 'bexpmap_00.fits',
            irfs='P8R3_SOURCE_V3')

    obs = ba.BinnedObs(**kw)

    srcModel = dirs['LOCAL_DATA'] + 'srcmdl_00.xml'
    optimizer = 'MINUIT'
    wmap = None
    edisp_bins = -1

    binned_config = ba.BinnedConfig(use_edisp=edisp_bins != 0,
                                    optimizer='MINUIT',
                                    applyPsfCorrections=True,
                                    performConvolutiion=True,
                                    resample=True,
                                    resamp_factor=2,
                                    minbinsz=0.05, 
                                    verbose=False,
    #                                save_all_srcmaps=True,
                                    edisp_bins=edisp_bins)
    print(binned_config.delete_local_fixed(),binned_config.save_all_srcmaps())
    print("1")
    logLike = pyLike.BinnedLikelihood(obs.countsMap, obs.observation, binned_config, obs.srcMaps, wmap)
    print("2")

    logLike.initOutputStreams()
    print("3")
    logLike.readXml(dirs['LOCAL_DATA'] + '/srcmdl_00.xml', _funcFactory, False, True, False)
    #logLike.readXml('/Users/echarles/software/fermi/packages/fermipy/fermipy_test_draco/srcmdl_00.xml', _funcFactory, False, True, False)
    print("4")
    nobs = logLike.countsSpectrum()
    print("5")
    logLike.buildFixedModelWts()
    print("6")
    print (logLike.npred())
    print ("Done.")
