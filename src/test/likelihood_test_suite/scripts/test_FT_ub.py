from UnbinnedAnalysis import *

NOMINAL_MEMORY = 1200
NOMINAL_TIME = 6.5

def test_FT_ub(dirs):

    sroi = "_1001"
    suffix = sroi + "_E19.fits"
    ltcube = dirs['LOCAL_DATA'] + "ltcube_10years_zmax105.fits"
    evfile = dirs['LOCAL_DATA'] + "MergedEVENTS" + suffix
    expmap = dirs['LOCAL_DATA'] + "expMap" + suffix
    irfs = "P8R3_SOURCE_V2"
    scfile = dirs['LOCAL_DATA'] + "ft2_10years_reduced.fits"
    srcmodel = dirs['LOCAL_DATA'] + 'srcModel' + sroi + '.xml-initub'

    like = unbinnedAnalysis(evfile=evfile, scfile=scfile, expmap=expmap,  \
        expcube=ltcube, irfs=irfs, optimizer="Minuit", srcmdl=srcmodel)
    for srcname in like.sourceNames():
        print (srcname, like.NpredValue(srcname))
    like.fit()
    for srcname in like.sourceNames():
        print (srcname, like.NpredValue(srcname), like.Ts(srcname))
