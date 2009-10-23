# -*- python -*-
# $Id: SConscript,v 1.83 2009/10/19 19:17:34 jchiang Exp $
# Authors: James Chiang <jchiang@slac.stanford.edu>, Pat Nolan <pln@razzle.stanford.edu>
# Version: Likelihood-15-08-00

Import('baseEnv', 'listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('LikelihoodLib', depsOnly = 1)
LikelihoodLib = libEnv.StaticLibrary('Likelihood', 
                                     listFiles(['src/*.c', 'src/*.cxx',
                                                'src/dmfit/*.cxx', 'src/dmfit/*.c']))

progEnv.Tool('LikelihoodLib')

testEnv = progEnv.Clone()
testEnv.Tool('addLibrary', library = baseEnv['cppunitLibs'])
test_LikelihoodBin = testEnv.Program('test_Likelihood',
                                     listFiles(['src/test/*.cxx']))

gtlikeBin = progEnv.Program('gtlike',listFiles(['src/likelihood/*.cxx']))

gtexpmapBin = progEnv.Program('gtexpmap', listFiles(['src/expMap/*.cxx']))

gttsmapBin = progEnv.Program('gttsmap', listFiles(['src/TsMap/*.cxx']))

gtltcubeBin = progEnv.Program('gtltcube', listFiles(['src/makeExposureCube/*.cxx']))

gtdiffrspBin = progEnv.Program('gtdiffrsp', listFiles(['src/diffuseResponses/*.cxx']))

gtsrcmapsBin = progEnv.Program('gtsrcmaps', listFiles(['src/gtsrcmaps/*.cxx']))

gtpsfBin = progEnv.Program('gtpsf', listFiles(['src/meanPsf/*.cxx']))

gtbkgBin = progEnv.Program('gtbkg', listFiles(['src/backfile/*.cxx']))

gtmodelBin = progEnv.Program('gtmodel', listFiles(['src/gtmodelmap/*.cxx']))

gtltsumBin = progEnv.Program('gtltsum', listFiles(['src/gtaddlivetime/*.cxx']))

gtfindsrcBin = progEnv.Program('gtfindsrc', listFiles(['src/gtfindsrc/*.cxx']))

progEnv.Tool('registerTargets', package = 'Likelihood', 
             staticLibraryCxts = [[LikelihoodLib,libEnv]], 
             binaryCxts = [[gtlikeBin,progEnv], [gtexpmapBin,progEnv], [gttsmapBin,progEnv],
                           [gtltcubeBin,progEnv], [gtdiffrspBin,progEnv], [gtsrcmapsBin,progEnv],
                           [gtpsfBin,progEnv], [gtbkgBin,progEnv], [gtmodelBin,progEnv],
                           [gtltsumBin,progEnv], [gtfindsrcBin,progEnv]],
             testAppCxts = [[test_LikelihoodBin, testEnv]],
             includes = listFiles(['Likelihood/*.h']), 
             pfiles = listFiles(['pfiles/*.par']),
             data = listFiles(['data/*'], recursive = True),
             xml = listFiles(['xml/*'], recursive = True))
