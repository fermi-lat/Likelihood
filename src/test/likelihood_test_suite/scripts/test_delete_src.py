from BinnedAnalysis import *

NOMINAL_MEMORY = 86
NOMINAL_TIME = 0.5

def test_delete_src(dirs):
    # test data is in dirs['LOCAL_DATA']
    like = binnedAnalysis(cmap=dirs['LOCAL_DATA'] + "srcMap_1_E1.fits", 
                          expcube=dirs['LOCAL_DATA'] + "LiveTimeCube_1.fits",
                          bexpmap=dirs['LOCAL_DATA'] + "binned_exposure_E1.fits", 
                          irfs="CALDB",
                          optimizer="Minuit", 
                          srcmdl=dirs['LOCAL_DATA'] + "srcModel_1_E1.xml-initial")
    print("Full model")
    print(like.model)
    like.fit(covar=True)

    srcname = "EXTIC443"
    TabFreeze = [7, 8, 11]
    like.freeze(TabFreeze)
    print("TS("+srcname+") =", like.Ts(srcname, reoptimize=True))
    like.thaw(TabFreeze)
    like.deleteSource("LAT0005")
    logLike = -like()
    print(logLike)
