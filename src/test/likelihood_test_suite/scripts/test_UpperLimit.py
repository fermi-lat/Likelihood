from BinnedAnalysis import *
import IntegralUpperLimit

NOMINAL_MEMORY = 87
NOMINAL_TIME = 0.5

def test_UpperLimit(dirs):
    # test data is in dirs['LOCAL_DATA']
    srcModel = dirs['LOCAL_DATA'] + "srcModel_3_2.xml-initial"
    expCube = dirs['LOCAL_DATA'] + "ltcube_july09_pass8.fits"
    map = dirs['LOCAL_DATA'] + "srcMap_3_2.fits"
    expMap = dirs['LOCAL_DATA'] + "binned_exposure_2.fits"

    #Test code goes here

    myconf = BinnedConfig(edisp_bins=-1)
    obs = BinnedObs( irfs="CALDB", expCube=expCube, srcMaps=map, binnedExpMap=expMap)
    like = BinnedAnalysis(obs,srcModel,optimizer="NewMinuit", config=myconf, delete_local_fixed=True)
    like.fit(covar=True)
    srcname = "PSRJ0357+32"
    emin = 100.0
    emax = 300.0
    BoundCompute = IntegralUpperLimit.calc_int( like, srcname, cl=0.95, emin=emin, emax=emax, verbosity=50 )
    if BoundCompute[0] > 0:
        print("IntegralUpperLimit succeeded:",BoundCompute[1]['ul_value'])
    else:
        print("IntegralUpperLimit failed.")
