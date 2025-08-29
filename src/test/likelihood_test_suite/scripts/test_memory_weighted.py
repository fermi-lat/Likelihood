from BinnedAnalysis import *
from SummedLikelihood import *

NOMINAL_MEMORY = 1265
NOMINAL_TIME = 28

def test_memory(dirs):
    roi = "_1001"
    suffix = "_E13"
    ltname = dirs['LOCAL_DATA'] + "ltcube_10years_zmax105.fits"
    srcmap = dirs['LOCAL_DATA'] + "srcMap" + roi + suffix + ".fits"
    bexpmap = dirs['LOCAL_DATA'] + "binned_exposure" + suffix + ".fits"
    srcmodel = dirs['LOCAL_DATA'] + "srcModel" + roi + "_memory.xml"
    wmap = dirs['LOCAL_DATA'] + "DataE2XWgts_3percent_P305_zmax105_PSF3.fits"

    like = SummedLikelihood()

    for k in range(20):
        obs = BinnedObs(irfs="CALDB", expCube=ltname, srcMaps=srcmap, binnedExpMap=bexpmap)
        analysis = BinnedAnalysis(obs, srcmodel, wmap = wmap)
        like.addComponent(analysis)

    like.optimizer = "Minuit"

    like.fit()