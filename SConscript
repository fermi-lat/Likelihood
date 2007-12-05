import glob,os

Import('baseEnv', 'listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

LikelihoodLib = libEnv.StaticLibrary('Likelihood', listFiles(['src/*.c', 'src/*.cxx']))

progEnv.Tool('LikelihoodLib')
gtlikelihoodBin = progEnv.Program('gtlikelihood',listFiles(['src/gtexpmap/*.cxx']))
gtexpmapBin = progEnv.Program('gtexpmap', listFiles(['src/gtexpmap/*.cxx']))
gttsmapBin = progEnv.Program('gttsmap', listFiles(['src/TsMap/*.cxx']))
gtlivetimecubeBin = progEnv.Program('gtlivetimecube', listFiles(['src/makeExposureCube/*.cxx']))
gtdiffrespBin = progEnv.Program('gtdiffresp', listFiles(['src/diffuseResponses/*.cxx']))
gtcntsmapBin = progEnv.Program('gtcntsmap', listFiles(['src/gtcntsmap/*.cxx']))
gtsrcmapsBin = progEnv.Program('gtsrcmaps', listFiles(['src/gtsrcmaps/*.cxx']))
gtpsfBin = progEnv.Program('gtpsf', listFiles(['src/meanPsf/*.cxx']))
gtbackfileBin = progEnv.Program('gtbackfile', listFiles(['src/backfile/*.cxx']))
gtmodelmapBin = progEnv.Program('gtmodelmap', listFiles(['src/gtmodelmap/*.cxx']))
gtaddlivetimeBin = progEnv.Program('gtaddlivetime', listFiles(['src/gtaddlivetime/*.cxx']))

progEnv.Tool('registerObjects', package = 'Likelihood', libraries = [LikelihoodLib], includes = listFiles(['Likelihood/*.h']), pfiles = listFiles(['pfiles/*.par']),
             binaries = [gtlikelihoodBin, gtexpmapBin, gttsmapBin, gtlivetimecubeBin, gtdiffrespBin, gtcntsmapBin, gtsrcmapsBin, gtpsfBin, gtbackfileBin, gtmodelmapBin,
                         gtaddlivetimeBin])
