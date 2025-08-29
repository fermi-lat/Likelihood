import time, os
from BinnedAnalysis import *

NOMINAL_MEMORY = 106
NOMINAL_TIME = 3.0

def test_freeze_thaw(dirs):
    # test data is in dirs['LOCAL_DATA']
    roi = "_1312"
    comp = "_E10"
    srcname = "P88Y0217"
    srcsfreeze2 = ["P88Y0186","P88Y0220"]
    bexpmap = dirs['LOCAL_DATA'] + "binned_exposure" + comp + ".fits"
    ltcube = dirs['LOCAL_DATA'] + "ltcube_12years_zmax105.fits"
    srcmap = dirs['LOCAL_DATA'] + "srcMap" + roi + comp + ".fits"
    irfs = "CALDB"
    dispbins = -2
    srcmodel = dirs['LOCAL_DATA'] + 'srcModel' + roi + '.xml-initial'

    myconf = BinnedConfig(edisp_bins=dispbins,save_all_srcmaps=False) # modified
    obs = BinnedObs(irfs=irfs, expCube=ltcube, srcMaps=srcmap, binnedExpMap=bexpmap)
    like = BinnedAnalysis(obs, srcmodel, optimizer="Minuit", config=myconf, delete_local_fixed=True) #modified

    like.fit()

    # Get TS values
    TSVals = {}
    for src in like.sourceNames():
        spmodel = like.model[src].funcs['Spectrum']
        normPar = spmodel.normPar().getName()
        indx = like.par_index(src,normPar)
        if like[indx].isFree():
            TSVals[src] = like.Ts(src)

    # Freeze local sources
    FreezeParams = []
    TabIndx = []
    for src2 in srcsfreeze2:
        for param in like.model[src2].funcs['Spectrum'].paramNames:
            indx = like.par_index(src2, param)
            if like[indx].isFree():
                FreezeParams.append((src2,param))
                TabIndx.append(indx)

    like.freeze(TabIndx)

    # Modify source model, originally LogParabola, and refit
    corres = {"Prefactor":"norm", "Index":"alpha", "Scale":"Eb"}
    bounds = {"Prefactor":(0,100), "Index":(0,5), "Scale":(0,1e5)}
    spmodel = like.model[srcname].funcs['Spectrum']
    parvals = {}
    parscal = {}
    for parname in ["norm","alpha","Eb"]:
        sppar = spmodel.params[parname]
        parvals[parname] = sppar.getValue()
        parscal[parname] = sppar.getScale()

    parscal["alpha"] *= -1   # Opposite sign convention for alpha and Index

    like.setSpectrum(srcname, "PowerLaw")
    for parname in ["Prefactor", "Index", "Scale"]:
        indx = like.par_index(srcname, parname)
        like[indx] = 1
        like[indx].setBounds(bounds[parname][0],bounds[parname][1])
        like[indx].setScale(parscal[corres[parname]])
        like[indx] = parvals[corres[parname]]

    like.fit()

    # Thaw parameters and refit
    TabIndx = []
    for FreezeParam in FreezeParams:
        TabIndx.append(like.par_index(FreezeParam[0],FreezeParam[1]))

    like.thaw(TabIndx)
    like.fit()
