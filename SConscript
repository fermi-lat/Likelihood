# -*- python -*-
#
# $Id: SConscript,v 1.3 2007/12/11 19:04:14 jchiang Exp $

import glob, os

Import('baseEnv', 'listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

LikelihoodLib = libEnv.StaticLibrary('Likelihood', 
                                     listFiles(['src/*.c', 'src/*.cxx']))

progEnv.Tool('LikelihoodLib')

test_LikelihoodBin = progEnv.Program('test_Likelihood',
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

progEnv.Tool('registerObjects', package = 'Likelihood', 
             libraries = [LikelihoodLib], 
             binaries = [gtlikeBin, gtexpmapBin, gttsmapBin, gtltcubeBin,
                         gtdiffrspBin, gtsrcmapsBin, gtpsfBin, gtbkgBin, 
                         gtmodelBin, gtltsumBin, gtfindsrcBin],
             testApps = [test_LikelihoodBin],
             includes = listFiles(['Likelihood/*.h']), 
             pfiles = listFiles(['pfiles/*.par']))
