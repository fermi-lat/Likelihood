import pyLikelihood as pyLike
import BinnedAnalysis as ba

NOMINAL_MEMORY = 95
NOMINAL_TIME = 0.85

def test_lhAnalysis4(dirs):
    # test data is in dirs['LOCAL_DATA']

    #Test code goes here
    roi = "_4"
    bexpmap = dirs['LOCAL_DATA'] + "binned_exposure.fits"
    ltcube = dirs['LOCAL_DATA'] + "ltcube_july09_pass8.fits"
    srcmap = dirs['LOCAL_DATA'] + "srcMap" + roi + ".fits"
    irfs = "CALDB"
    scdata = dirs['LOCAL_DATA'] + "ft2_july09_pass8.fits"
    dispbins = 0
    srcmodel = dirs['LOCAL_DATA'] + 'srcModel' + roi + '.xml'
    srcmodelin = dirs['LOCAL_DATA'] + 'srcModel' + roi + '_tHb1.xml-initial'

    # Read data and fit
    obs = ba.BinnedObs(irfs=irfs, expCube=ltcube, srcMaps=srcmap, binnedExpMap=bexpmap)
    like = ba.BinnedAnalysis(obs, srcmodelin, optimizer="Minuit", delete_local_fixed=True)
    like.fit()

    # Freeze isotropic and refit
    indx = like.par_index("eg_v02","Normalization")
    like.freeze(indx)
    print("Freeze isotropic and refit")
    like.fit()