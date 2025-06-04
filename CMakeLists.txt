##### Library ######
add_library(
  Likelihood SHARED
  src/Accumulator.cxx src/AppHelpers.cxx src/BandFunction.cxx src/BinnedConfig.cxx
  src/BinnedCountsCache.cxx src/BinnedExposure.cxx src/BinnedExposureBase.cxx
  src/BinnedHealpixExposure.cxx src/BinnedLikelihood.cxx src/BrokenPowerLaw2.cxx
  src/BrokenPowerLaw3.cxx src/BrokenPowerLawExpCutoff.cxx src/Composite2.cxx
  src/CompositeLikelihood.cxx src/CompositeSource.cxx src/Convolve.cxx
  src/ConvolveHealpix.cxx src/CountsMap.cxx src/CountsMapBase.cxx
  src/CountsMapHealpix.cxx src/CountsSpectra.cxx src/DiffRespIntegrand.cxx
  src/DiffRespIntegrand2.cxx src/DiffRespNames.cxx src/DiffuseSource.cxx src/Drm.cxx
  src/EblAtten.cxx src/EnergyBand.cxx src/EquinoxRotation.cxx src/Event.cxx
  src/EventContainer.cxx src/ExpCutoff.cxx src/ExpCutoffSEDPeak.cxx
  src/ExposureCube.cxx src/ExposureMap.cxx src/FileFunction.cxx src/FileUtils.cxx
  src/FitScanner.cxx src/FitUtils.cxx src/FluxBuilder.cxx src/GaussianError.cxx
  src/HealpixProjMap.cxx src/HistND.cxx src/LikeExposure.cxx src/LogGaussian.cxx
  src/LogLike.cxx src/LogNormal.cxx src/LogNormalMuDist.cxx src/LogParabola.cxx
  src/MapBase.cxx src/MapCubeFunction2.cxx src/MeanPsf.cxx src/ModelMap.cxx
  src/MultipleBrokenPowerLaw.cxx src/Npred.cxx src/OneSourceFunc.cxx src/OptEM.cxx
  src/PiecewisePowerLaw.cxx src/Pixel.cxx src/PointSource.cxx src/PowerLaw2.cxx
  src/PowerLawSuperExpCutoff.cxx src/PowerLawSuperExpCutoff2.cxx
  src/PowerLawSuperExpCutoff3.cxx src/PowerLawSuperExpCutoff4.cxx src/ProjMap.cxx
  src/PSFUtils.cxx src/RadialDisk.cxx src/RadialGaussian.cxx src/RadialProfile.cxx
  src/ResponseFunctions.cxx src/RoiCuts.cxx src/ScaleFactor.cxx src/ScanUtils.cxx
  src/ScData.cxx src/SkyDirFunction.cxx src/SmoothBrokenPowerLaw.cxx
  src/SmoothDoubleBrokenPowerLaw.cxx src/Snapshot.cxx src/Source.cxx
  src/SourceFactory.cxx src/SourceMap.cxx src/SourceMapCache.cxx src/SourceModel.cxx
  src/SourceModelBuilder.cxx src/SparseVector.cxx src/SpatialFunction.cxx
  src/SpatialMap.cxx src/SummedLikelihood.cxx src/TiedParameter.cxx src/TrapQuad.cxx
  src/WcsMap.cxx src/WcsMap2.cxx src/WcsMapLibrary.cxx src/WeightMap.cxx
  src/XmlBuilder.cxx src/XmlParser.cxx
  src/dmfit/dmfit_comm.c src/dmfit/dmfit_func.c src/dmfit/dmfit_load.c
  src/dmfit/DMFitFunction.cxx
)

target_link_libraries(
  Likelihood
  
  PUBLIC
  st_app
  astro
  dataSubselector
  optimizers
  evtbin
  tip
  healpix
  map_tools
  GSL::gsl
  eblAtten
  Python3::Python


  PRIVATE
  xmlBase
  st_facilities
  hoops
  st_graph
  irfLoader
  cfitsio::cfitsio
  FFTW::fftw3
  OpenMP::OpenMP_CXX
)
target_include_directories(
  Likelihood PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

add_compile_definitions(WITHOUT_NUMPY)

target_compile_definitions(Likelihood PRIVATE CONDA_FFTW)


###### Executables ######
add_executable(
  gtlike
  src/likelihood/EasyPlot.cxx
  src/likelihood/likelihood.cxx
  src/likelihood/MathUtil.cxx
)
add_executable(gtexpmap src/expMap/expMap.cxx)
add_executable(gttsmap src/TsMap/TsMap.cxx)
add_executable(gttscube src/TsCube/TsCube.cxx)
add_executable(gthealcube src/TsCube/HealCube.cxx)
add_executable(gtltcube src/makeExposureCube/makeExposureCube.cxx)
add_executable(gtexpcube2 src/gtexpcube2/gtexpcube2.cxx)
add_executable(gtdiffrsp src/diffuseResponses/diffuseResponses.cxx)
add_executable(gtsrcmaps src/gtsrcmaps/gtsrcmaps.cxx)
add_executable(gtebl src/gtebl/gtebl.cxx)
add_executable(gtsrcprob src/gtsrcprob/gtsrcprob.cxx)
add_executable(gtpsf src/meanPsf/meanPsf.cxx)
add_executable(gtbkg src/backfile/backfile.cxx)
add_executable(gtmodel src/gtmodelmap/gtmodelmap.cxx)
add_executable(gtltsum src/gtaddlivetime/gtaddlivetime.cxx)
add_executable(gtfindsrc src/gtfindsrc/gtfindsrc.cxx)
add_executable(gteffbkg src/gteffbkg/gteffbkg.cxx)
add_executable(gtalphabkg src/gtalphabkg/gtalphabkg.cxx)
add_executable(gtwtsmap src/gtwtsmap/gtwtsmap.cxx)
add_executable(gtdrm src/gtdrm/gtdrm.cxx)

target_link_libraries(gtlike PRIVATE Likelihood st_graph)
target_link_libraries(gtexpmap PRIVATE Likelihood st_stream)
target_link_libraries(gttsmap PRIVATE Likelihood)
target_link_libraries(gttscube PRIVATE Likelihood)
target_link_libraries(gthealcube PRIVATE Likelihood)
target_link_libraries(gtltcube PRIVATE Likelihood)
target_link_libraries(gtexpcube2 PRIVATE Likelihood)
target_link_libraries(gtdiffrsp PRIVATE Likelihood)
target_link_libraries(gtsrcmaps PRIVATE Likelihood)
target_link_libraries(gtebl PRIVATE Likelihood)
target_link_libraries(gtsrcprob PRIVATE Likelihood)
target_link_libraries(gtpsf PRIVATE Likelihood)
target_link_libraries(gtbkg PRIVATE Likelihood)
target_link_libraries(gtmodel PRIVATE Likelihood)
target_link_libraries(gtltsum PRIVATE Likelihood)
target_link_libraries(gtfindsrc PRIVATE Likelihood)
target_link_libraries(gteffbkg PRIVATE Likelihood)
target_link_libraries(gtalphabkg PRIVATE Likelihood)
target_link_libraries(gtwtsmap PRIVATE Likelihood)
target_link_libraries(gtdrm PRIVATE Likelihood)

add_custom_target(Likelihood-executables)
add_dependencies(
  Likelihood-executables gtlike gtexpmap gttsmap gttscube gthealcube gtltcube
  gtexpcube2 gtdiffrsp gtsrcmaps gtebl gtsrcprob gtpsf gtbkg gtmodel gtltsum gtfindsrc
  gteffbkg gtalphabkg gtwtsmap gtdrm
)

###### Tests ######

add_executable(
  test_Likelihood
  src/test/test.cxx
  src/test/SourceData.cxx
  src/test/XmlDiff.cxx
)
target_link_libraries(test_Likelihood PRIVATE Likelihood CppUnit::CppUnit)

if(APPLE)
  target_compile_definitions(test_Likelihood PRIVATE DARWIN DARWIN_F2C_FAILURE)
else()
  target_compile_definitions(test_Likelihood PRIVATE TRAP_FPE)
endif()

###############################################################
# Installation
###############################################################
install(DIRECTORY Likelihood DESTINATION ${FERMI_INSTALL_INCLUDEDIR})
install(DIRECTORY data/ DESTINATION ${FERMI_INSTALL_REFDATADIR}/Likelihood)
install(DIRECTORY pfiles/ DESTINATION ${FERMI_INSTALL_PFILESDIR})
install(DIRECTORY xml/ DESTINATION ${FERMI_INSTALL_XMLDIR}/Likelihood)

install(
  TARGETS Likelihood gtlike gtexpmap gttsmap gttscube gthealcube gtltcube gtexpcube2
  gtdiffrsp gtsrcmaps gtebl gtsrcprob gtpsf gtbkg gtmodel gtltsum gtfindsrc gteffbkg
  gtalphabkg gtwtsmap gtdrm test_Likelihood
  EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
