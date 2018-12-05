# -*- python -*-
# $Id: SConscript,v 1.384 2017/10/13 23:30:22 echarles Exp $
# Authors: James Chiang <jchiang@slac.stanford.edu>, Eric Charles <echarles@slac.stanford.edu>, Matthew Wood <mdwood@slac.stanford.edu> 
# Version: Likelihood-20-13-03

import sys
Import('baseEnv', 'listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('addLinkDeps', package='Likelihood', toBuild='shared')
LikelihoodLib = libEnv.SharedLibrary('Likelihood', 
                                     listFiles(['src/*.c', 'src/*.cxx',
                                                'src/dmfit/*.cxx', 
                                                'src/dmfit/*.c']))

progEnv.Tool('LikelihoodLib')

testEnv = progEnv.Clone()
testEnv.Tool('addLibrary', library = baseEnv['cppunitLibs'])

if sys.platform == 'darwin':
    testEnv.Append(CPPDEFINES = 'DARWIN')
test_LikelihoodBin = testEnv.Program('test_Likelihood',
                                     listFiles(['src/test/*.cxx']))

gtlikeBin = progEnv.Program('gtlike',listFiles(['src/likelihood/*.cxx']))

gtexpmapBin = progEnv.Program('gtexpmap', listFiles(['src/expMap/*.cxx']))

gttsmapBin = progEnv.Program('gttsmap', listFiles(['src/TsMap/*.cxx']))

gttscubeBin = progEnv.Program('gttscube', listFiles(['src/TsCube/TsCube.cxx']))

gthealcubeBin = progEnv.Program('gthealcube', listFiles(['src/TsCube/HealCube.cxx']))

gtltcubeBin = progEnv.Program('gtltcube', listFiles(['src/makeExposureCube/*.cxx']))

gtexpcube2Bin = progEnv.Program('gtexpcube2', listFiles(['src/gtexpcube2/*.cxx']))

gtdiffrspBin = progEnv.Program('gtdiffrsp', listFiles(['src/diffuseResponses/*.cxx']))

gtsrcmapsBin = progEnv.Program('gtsrcmaps', listFiles(['src/gtsrcmaps/*.cxx']))

gteblBin = progEnv.Program('gtebl', listFiles(['src/gtebl/*.cxx']))

gtsrcprobBin = progEnv.Program('gtsrcprob', listFiles(['src/gtsrcprob/*.cxx']))

gtpsfBin = progEnv.Program('gtpsf', listFiles(['src/meanPsf/*.cxx']))

gtbkgBin = progEnv.Program('gtbkg', listFiles(['src/backfile/*.cxx']))

gtmodelBin = progEnv.Program('gtmodel', listFiles(['src/gtmodelmap/*.cxx']))

gtltsumBin = progEnv.Program('gtltsum', listFiles(['src/gtaddlivetime/*.cxx']))

gtfindsrcBin = progEnv.Program('gtfindsrc', listFiles(['src/gtfindsrc/*.cxx']))

gteffbkgBin = progEnv.Program('gteffbkg', listFiles(['src/gteffbkg/*.cxx']))

gtalphabkgBin = progEnv.Program('gtalphabkg', listFiles(['src/gtalphabkg/*.cxx']))

gtwtsmapBin = progEnv.Program('gtwtsmap', listFiles(['src/gtwtsmap/*.cxx']))

gtdrmBin = progEnv.Program('gtdrm', listFiles(['src/gtdrm/*.cxx']))

progEnv.Tool('registerTargets', package = 'Likelihood', 
             libraryCxts = [[LikelihoodLib, libEnv]], 
             binaryCxts = [[gtlikeBin, progEnv], 
                           [gtexpmapBin, progEnv],
                           [gttsmapBin, progEnv],
                           [gttscubeBin, progEnv],
                           [gthealcubeBin, progEnv],
                           [gtltcubeBin, progEnv],
                           [gtexpcube2Bin, progEnv],
                           [gtdiffrspBin, progEnv],
                           [gtsrcmapsBin, progEnv],
                           [gtsrcprobBin, progEnv],
                           [gtpsfBin, progEnv],
                           [gteblBin, progEnv],
                           [gtbkgBin, progEnv],
                           [gtmodelBin, progEnv],
                           [gtltsumBin, progEnv],
                           [gtfindsrcBin, progEnv],
                           [gteffbkgBin, progEnv],
                           [gtalphabkgBin, progEnv],
                           [gtwtsmapBin, progEnv],
                           [gtdrmBin, progEnv]],
             testAppCxts = [[test_LikelihoodBin, testEnv]],
             includes = listFiles(['Likelihood/*.h']), 
             pfiles = listFiles(['pfiles/*.par']),
             data = listFiles(['data/*'], recursive = True),
             xml = listFiles(['xml/*'], recursive = True))
