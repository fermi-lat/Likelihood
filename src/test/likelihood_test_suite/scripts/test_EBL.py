from BinnedAnalysis import *
from LikelihoodState import LikelihoodState

NOMINAL_MEMORY = 115
NOMINAL_TIME = 3.2

def test_EBL(dirs):
    # test data is in dirs['LOCAL_DATA']

    #Test code goes here
    roi = "_1055"
    suffix = "_E15"
    srcname = "P88Y3676"
    ltdir = dirs['LOCAL_DATA']
    refdir = dirs['LOCAL_DATA']
    refdir = dirs['LOCAL_DATA']
    ltname = ltdir + "ltcube_8years_zmax105.fits"
    srcmap = refdir + "srcMap" + roi + suffix + ".fits"
    bexpmap = refdir + "binned_exposure" + suffix + ".fits"
    srcmodel = dirs['LOCAL_DATA'] + "srcModel" + roi + "_testEBL.xml"

    analysis= binnedAnalysis(cmap=srcmap, bexpmap=bexpmap, expcube=ltname,
                        irfs="CALDB", optimizer="Minuit", srcmdl=srcmodel)
    analysis.fit(covar=True)

    saved_state = LikelihoodState(analysis)
    Ts = analysis.Ts(srcname)
    for srcNameEach in analysis.sourceNames():
        saved_state.restore(srcNameEach)

    print ("Freeze index")
    indx = analysis.par_index(srcname,"Index")
    analysis.freeze(indx)

    print ("Refit")
    analysis.fit(covar=True)
    Ts = analysis.Ts(srcname, reoptimize=True)
    print ("TS =", Ts)