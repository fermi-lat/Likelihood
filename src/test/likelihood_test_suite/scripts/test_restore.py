# This is the test from the pylike-19 directory.
# for some reason the original runs fine on the command line
# but when convert to this functional format, crashes with
# the wrong energy bins on the `summed_like.fit()` line.

from SummedLikelihood import *
from BinnedAnalysis import *
from LikelihoodState import LikelihoodState

NOMINAL_MEMORY = 120
NOMINAL_TIME = 1.25

def test_restore(dirs):
    # test data is in dirs['LOCAL_DATA'] + 

    # Build summed like model and apply fit with original spectrum type :
    summed_like = SummedLikelihood()
    myconf = BinnedConfig(edisp_bins=-1)
    obs = BinnedObs( irfs="CALDB", 
                    expCube=dirs['LOCAL_DATA'] + "LiveTimeCube_4.fits", 
                    srcMaps=dirs['LOCAL_DATA'] + "srcMap_4_E1.fits", 
                    binnedExpMap=dirs['LOCAL_DATA'] + "binned_exposure_E1-restore.fits") 
    like = BinnedAnalysis(obs,dirs['LOCAL_DATA'] + "srcModel_4.xml-initial",
                          optimizer="Minuit", 
                          config=myconf)
    summed_like.addComponent(like)
    obs = BinnedObs(irfs="CALDB", 
                    expCube=dirs['LOCAL_DATA'] + "LiveTimeCube_4.fits", 
                    srcMaps=dirs['LOCAL_DATA'] + "srcMap_4_E3.fits", 
                    binnedExpMap=dirs['LOCAL_DATA'] + "binned_exposure_E3.fits") 
    like = BinnedAnalysis(obs,dirs['LOCAL_DATA'] + "srcModel_4.xml-initial",
                          optimizer="Minuit", 
                          config=myconf)
    summed_like.addComponent(like)
    srcName = 'PSRJ0835-4510'
    summed_like.fit()
    print ("summed_like.Ts(",srcName,"):",summed_like.Ts(srcName))
    SpectrumType = like.model[srcName].src.spectrum().genericName()
    saved_state = LikelihoodState(summed_like)
    spmodel = summed_like.model[srcName].funcs['Spectrum']
    print ("spmodel:",spmodel)
    


    TabIndx = []

    indx = summed_like.par_index(srcName, 'alpha')
    TabIndx.append(indx)
    summed_like[indx] = 5
    
    indx = summed_like.par_index(srcName, 'beta')
    TabIndx.append(indx)
    summed_like[indx].setBounds(-5,5)
    
    spmodel = summed_like.model[srcName].funcs['Spectrum']

    print ("spmodel:",spmodel)
    
    summed_like.freeze(TabIndx)
    summed_like.fit()
    summed_like.thaw(TabIndx)

    # --

    print ("summed_like.Ts(",srcName,") for powerLaw :",summed_like.Ts( srcName ))

    # Go back to previous model :
    summed_like.setSpectrum( srcName, SpectrumType)
    for srcNameEach in summed_like.sourceNames():
        saved_state.restore(srcNameEach)

    print ("\nsummed_like.model:",summed_like.model)
    print ("summed_like.Ts(",srcName,") for",SpectrumType,":",summed_like.Ts(srcName))
    print ("summed_like.NpredValue(",srcName,"):",summed_like.NpredValue(srcName))