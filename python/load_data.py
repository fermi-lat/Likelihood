import Likelihood

Likelihood.RoiCuts_setCuts('RoiCuts.xml')
Likelihood.ScData_readData('test_scData_0000.fits')
Likelihood.LogLike_loadResponseFunctions("FRONT/BACK")
