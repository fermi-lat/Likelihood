/**
   @page tutorial Likelihood Tutorial

   @section intro Introduction

   In order to illustrate how to use the likelihood software, this
   narrative is provided and gives a step-by-step description for
   performing a relatively simple analysis.  We will assume that you
   already have event (FT1) and spacecraft (FT2) data files in hand.

   Here are each of the steps:

   - @ref makeSubselections Since there is a computational overhead
          for each event associated with each diffuse component, it is
          useful to filter out any events that are not within the
          extraction region used for the analysis.
   - @ref makeCountsMaps These simple FITS images let us see what 
          we've got and help to pick out obvious candidate sources.
   - @ref defineROI This is the extraction region in photon direction, 
          energy, and arrival time.
   - @ref makeExposureMap This is needed for analyzing diffuse sources
          and so is required for almost any sort of analysis.
   - @ref sourceModelFile The various model parameters are read by
          likelihood from xml files.  The 
          <a href="http://www.slac.stanford.edu/exp/glast/ground/software/RM/documentation/ScienceTools/ScienceTools-v0r5/modeldef/v1r0p0/userguide.html">modelDef</a> tool will make
          this process easier.
   - @ref runLikelihood , the likelihood application itself.
   - @ref checkWithObsSim Here we create counts maps using @b obsSim
          for comparison with the data.
   - @ref makeTsMaps This is used for point source localization and for
          finding weaker sources after the stronger sources have
          been modeled.

   @section makeSubselections Make Subselections from the Event Data.
   For this tutorial, we will use data generated using @b obsSim for a
   model comprising the Third EGRET (3EG) catalog sources, Galactic diffuse,
   and extragalactic diffuse emission for a simulation time of one
   day.  Since we used @b obsSim, we can put data from each component
   into a separate FT1 file:
   @verbatim
   noric13[jchiang] ls -l *events*fits
   -rw-rw-r--    1 jchiang  glast-data  2894400 Oct 11 15:22 eg_diffuse_events_0000.fits
   -rw-rw-r--    1 jchiang  glast-data  7243200 Oct 11 15:26 galdiffuse_events_0000.fits
   -rw-rw-r--    1 jchiang  glast-data  1088640 Oct 11 15:20 ptsrcs_events_0000.fits
   noric13[jchiang] 
   @endverbatim
   We will consider data in the Virgo region within a 25 degree acceptance cone of
   the blazar 3C 279.  Here we apply @b dataSubselector to the point source data:
   @verbatim
   noric13[jchiang] ~/ST_new/dataSubselector/v2r1/rh9_gcc32/dataSubselector.exe
   Input FT1 file [Crab_front_events_0000.fits] : ptsrcs_events_0000.fits
   Output FT1 file [Crab_front_events_filtered.fits] : ptsrcs_events_filtered.fits
   RA for new search center (degrees) <0 - 360> [83.57] : 193.98
   Dec for new search center (degrees) <-90 - 90> [22.01] : -5.82
   radius of new search region (degrees) <0 - 180> [20] : 25
   Longitude coordinate lower limit (degrees) <0 - 360> [0] : 
   Longitude coordinate upper limit (degrees) <0 - 360> [360] : 
   Latitude coordinate lower limit (degrees) <-90 - 90> [-90] : 
   Latitude coordinate upper limit (degrees) <-90 - 90> [90] : 
   Coordiate system <CEL|GAL> [CEL] : 
   start time (MET in s) [0] : 
   end time (MET in s) [0] : 
   lower energy limit (MeV) [20] : 
   upper energy limit (MeV) [0] : 3e5    
   first conversion layer <0 - 15> [0] : 
   last conversion layer <0 - 15> [11] : 15
   minimum theta value (degrees) <0 - 180> [0] : 
   maximum theta value (degrees) <0 - 180> [0] : 
   minimum phi value (degrees) <0 - 360> [0] : 
   maximum phi value (degrees) <0 - 360> [0] : 
   minimum probablilty that event is a gamma ray <0 - 1> [0] : 
   maximum probablilty that event is a gamma ray <0 - 1> [0] : 
   minimum zenith angle value (degrees) <0 - 180> [0] : 
   maximum zenith angle value (degrees) <0 - 180> [0] : 
   select only events that passed background cut? [yes] : 
   select only events that passed PSF cut? [yes] : 
   select only events that passed energy resolution cut? [yes] : 
   Done.
   @endverbatim
   The @b dataSubselector interface is a bit unwieldy.  The default
   entries of "0" for most of the parameters means that no selection
   will be made.  This tool, its interface in particular, is a prime
   target for refactoring.

   @section makeCountsMaps Make Counts Maps from the Event Files 
   The next step is to create a simple counts map file to visualize
   the extraction region we wish to fit.  There are presently three
   tools available to perform this task: @b evtbin in the @b evtbin
   package, @b count_map in the @b map_tools package, and @b gtcntsmap
   in the @b Likelihood package.  For this step, since the event data
   are contained in three separate FITS files, @b gtcntsmap is the
   most convenient to use.  We start by creating a file which is a
   list of the files to be binned, and then we run the tool:
   @verbatim
   noric13[jchiang] ls -1 *filtered.fits > eventFiles
   noric13[jchiang] cat eventFiles
   eg_diffuse_events_filtered.fits
   galdiffuse_events_filtered.fits
   ptsrcs_events_filtered.fits
   noric13[jchiang]  ~/ST_new/Likelihood/v5r1p1/rh9_gcc32/gtcntsmap.exe
   event file [test_events_0000.fits] : eventFiles
   spacecraft file [test_scData_0000.fits] : ptsrcs_scData_0000.fits
   minimum energy (MeV) <20 - 3e5> [30] : 
   maximum energy (MeV) <20 - 3e5> [2e5] : 
   number of energy bounds <2 - 40> [21] : 2
   Right Ascension (or l) of map center (degrees) <0 - 360> [86.4] : 193.98
   Declination (or b) of map center (degrees) <-90 - 90> [28.9] : -5.82   
   number of longitude pixels [160] : 200
   number of latitude pixels [160] : 200
   pixel size (degrees) <1e-2 - 2> [0.25] : 
   use Galactic coordinates [no] : 
   output file name [counts_map.fits] : Virgo_map.fits
   noric13[jchiang] ds9 Virgo_map.fits
   @endverbatim
   The last command lauches the visualization tool ds9 and produces this
   display:

   @image html Virgo_map.png Virgo_map.fits

   The green circles indicate the positions of the 3EG point sources
   within 25 degrees of the location of 3C 279.

   @section defineROI Define the Region-of-Interest The
   Region-of-Interest (ROI) comprises a set of selections on photon
   arrival time, energy, and direction.  These selections are made in
   addition to any that are made using the query to the <a
   href="http://glast.gsfc.nasa.gov/cgi-bin/ssc/U1/D1WebQuery.cgi">
   GSSC database</a> or via the @b dataSubselector user-level
   selection tool</a>.  Unfortunately, because of an oversight in the
   Likelihood implementation, the ROI selections must @em include any
   that are made by these other tools.  More precisely, the data-space
   defined by the ROI must lie entirely within the intersection of the
   data-spaces defined by these other tools.  For example, if the
   user-level selection tool selects events within a time interval
   (t1, t2), then any time intervals defined in the ROI file must lie
   entirely within (t1, t2).  In future releases, this oversight will
   be corrected.
   \n\n
   Here's the ROI file we will use:
   @verbatim
   <?xml version='1.0' standalone='no'?>
   <Region-of-Interest title="Virgo Region">
      <timeInterval start="0"
                    stop="1"
                    unit="days"/>
      <energies emin="30." 
                emax="2e5"
                unit="MeV"/>
      <acceptanceCone longitude="193.98."
                      latitude="-5.82"
                      radius="20"
                      coordsys="J2000"/>
   </Region-of-Interest>
   @endverbatim
   Any number of time intervals can be specified and can be given in
   units of "seconds" or "days", referenced to the start of the
   simulation.  Data will be selected from the intersection of
   intervals given.  In future, the GTI extension of the FITS event
   file(s) will provide this same capability. A similar consideration
   applies to energy and acceptance cone cuts made, for example, using
   @b dataSubselector.  The "coordsys" attribute can either be
   "Galactic" or "J2000"; the "radius" attribute gives the
   half-opening angle of the acceptance cone; and units are in
   degrees.

   @section makeExposureMap Make an Exposure Map for a Given ROI
   With the ROI specified, we are now ready to create an exposure map.
   The type of exposure map used by likelihood differs significantly
   from the usual notion of exposure maps, which are essentially
   integrals of effective area over time.  The exposure calculation
   that likelihood uses consists of an integral of the total response
   over the entire ROI data-space:
   \f[ \epsilon(E, \hat{p}) =
        \int_{\rm ROI} dE^\prime d\hat{p}^\prime dt 
        R(E^\prime, \hat{p}^\prime; E, \hat{p}, t),
   \f]
   where primed quantities indicate measured energies, \f$E^\prime\f$
   and measured directions, \f$\hat{p}^\prime\f$.
   This exposure function can then be used to compute the expected
   numbers of events from a given source, \f$S_i(E, \hat{p})\f$:
   \f[
   N_{\rm pred} = \int dE d\hat{p} S_i(E, \hat{p}) \epsilon(E, \hat{p})
   \f]
   Since the exposure calculation involves an integral over the ROI,
   separate exposure maps must be made for every distinct ROI file,
   if, for example, one wants to subdivide an observation to look for
   secular flux variations from a particular source or sources.

   There are two separate tools for generating exposure maps.  The
   first is @b makeExposureCube.  This tool creates an all-sky map of
   the integrated livetime as a function of inclination with respect
   to the LAT z-axis.  It is otherwise identical to the @b
   exposure_cube application in the @b map_tools package, except that
   makeExposureCube also applies the time-interval cuts specified in
   the ROI file.  
   @verbatim
   noric13[jchiang] ~/ST_new/Likelihood/v5r1p1/rh9_gcc32/makeExposureCube.exe 
   Spacecraft data file [single_src_scData_0000.fits] : ptsrcs_scData_0000.fits
   Output file [!expcube_1_day.fits] : 
   Step size in cos(theta) <0. - 1.> [0.25] : 0.05
   Pixel size (degrees) [10] : 1.
   ROI file [anticenter_Roi.xml] : RoiCuts.xml
   Creating a exposure hypercube, size 1231200=360 x 180 x 19
   Working on file ptsrcs_scData_0000.fits
   .....................!
   noric13[jchiang] 
   @endverbatim
   Since @b makeExposureCube produces an all-sky map, the output FITS
   file of this tool can be used to generating exposure maps for other
   regions-of-interest that have the same time interval selections.
   Although the @expMap application (see below) can generate exposure
   maps for Likelihood without an exposure hypercube map, using one
   affords a substantial time savings.
   \n\n
   We create the exposure map using the @b expMap tool:
   @verbatim
   noric13[jchiang] ~/ST_new/Likelihood/v5r1p1/rh9_gcc32/expMap.exe 
   ROI cuts file [../data/anticenter_Roi.xml] : RoiCuts.xml
   Spacecraft file [../data/single_src_scData_0000.fits] : ptsrcs_scData_0000.fits
   Response functions to use <FRONT/BACK|FRONT|BACK|GLAST25> [GLAST25] : FRONT/BACK
   Use energy dispersion? [no] : 
   Radius of the source region (in degrees) [30] : 
   Number of longitude points <2 - 1000> [10] : 120
   Number of latitude points <2 - 1000> [10] : 120
   Number of energies <2 - 100> [5] : 20
   Exposure hypercube file [none] : expcube_1_day.fits
   Exposure file [test1.fits] : expMap.fits
   Loaded exposure map from a FITS file expcube_1_day.fits, size is 1231200
   Computing the ExposureMap....................!
   noric13[jchiang] 
   @endverbatim
   Note that we have chosen a 30 degree radius "source region", while
   the ROI acceptance cone radius is 20 degrees.  This is necessary to
   ensure that photons from sources outside the ROI are accounted for,
   at least partially, owing to the size of the instrument
   point-spread function. Half-degree pixels are a nominal choice;
   smaller pixels should result in a more accurate evaluation of the
   diffuse source fluxes but will also make the exposure map
   calculation itself lengthier.  The number of energies specifies the
   number of logarithmically spaced intervals bounded by the energy
   range given in the ROI cuts xml file.  Here is one image plane of
   the exposure map we just created:

   @image html expMap.png expMap.fits

   @section sourceModelFile Create a Source Model XML File
   Like the ROI file, the @b likelihoodApp reads the source model
   from an XML file. The <a href="http://www.slac.stanford.edu/exp/glast/ground/software/RM/documentation/ScienceTools/ScienceTools-v0r5/modeldef/v1r0p0/userguide.html">modelDef</a> application will provide a
   convenient means of defining sources and creating a file of the 
   proper format.  Here is the file we'll be using for our analysis:
   @verbatim
<source_library title="source library">
  <source name="Crab_fit" type="PointSource">
    <spectrum type="PowerLaw">
      <parameter max="1000" min="0.001" free="1" name="Prefactor" scale="1e-09" value="10" />
      <parameter max="-1" min="-3.5" free="1" name="Index" scale="1" value="-2" />
      <parameter max="200" min="50" free="0" name="Scale" scale="1" value="100" />
    </spectrum>
    <spatialModel type="SkyDirFunction">
      <parameter max="3.40282e+38" min="-3.40282e+38" free="0" name="RA" scale="1" value="83.57" />
      <parameter max="3.40282e+38" min="-3.40282e+38" free="0" name="DEC" scale="1" value="22.01" />
    </spatialModel>
  </source>
  <source name="Extragalactic Diffuse Emission" type="DiffuseSource">
    <spectrum type="PowerLaw">
      <parameter max="100" min="1e-05" free="1" name="Prefactor" scale="1e-07" value="1" />
      <parameter max="-1" min="-3.5" free="1" name="Index" scale="1" value="-2" />
      <parameter max="200" min="50" free="0" name="Scale" scale="1" value="100" />
    </spectrum>
    <spatialModel type="ConstantValue">
      <parameter max="10" min="0" free="0" name="Value" scale="1" value="1" />
    </spatialModel>
  </source>
  <source name="Geminga_fit" type="PointSource">
    <spectrum type="PowerLaw">
      <parameter max="1000" min="0.001" free="1" name="Prefactor" scale="1e-09" value="10" />
      <parameter max="-1" min="-3.5" free="1" name="Index" scale="1" value="-2" />
      <parameter max="200" min="50" free="0" name="Scale" scale="1" value="100" />
    </spectrum>
    <spatialModel type="SkyDirFunction">
      <parameter max="3.40282e+38" min="-3.40282e+38" free="0" name="RA" scale="1" value="98.49" />
      <parameter max="3.40282e+38" min="-3.40282e+38" free="0" name="DEC" scale="1" value="17.86" />
    </spatialModel>
  </source>
  <source name="PKS0528p134" type="PointSource">
    <spectrum type="PowerLaw">
      <parameter max="1000" min="0.001" free="1" name="Prefactor" scale="1e-09" value="10" />
      <parameter max="-1" min="-3.5" free="1" name="Index" scale="1" value="-2" />
      <parameter max="200" min="50" free="0" name="Scale" scale="1" value="100" />
    </spectrum>
    <spatialModel type="SkyDirFunction">
      <parameter max="3.40282e+38" min="-3.40282e+38" free="0" name="RA" scale="1" value="82.74" />
      <parameter max="3.40282e+38" min="-3.40282e+38" free="0" name="DEC" scale="1" value="13.38" />
    </spatialModel>
  </source>
</source_library>
   @endverbatim
   We'll call this source model file "anticenter_model.xml".
   \n\n
   We won't discuss the format of this file in great detail, but we
   will note some of the more salient features.  There are two kinds
   of sources that one can define, "PointSource" and "DiffuseSource".
   In turn, each type of source comprises two components, a "spectrum"
   and a "spatialModel".  Presently, one can choose from among various
   sorts of spectral types, the two most relevant ones being
   "PowerLaw" and "BrokenPowerLaw".  Additional spectral models will
   be available eventually, and it is intended that the user will be
   able to create custom spectral models that will be loadable
   dynamically by the application.
   \n\n
   There are three spatialModels presently available:
   "SkyDirFunction", "ConstantValue", and "SpatialMap".  The first one
   describes a direction on the sky and is used only for PointSources.
   The latter two are used by DiffuseSources.  ConstantValue is
   precisely what it appears to be --- it provides a constant value
   regardless of what argument value it takes.  It is used in this
   context to model the isotropic diffuse emission.  However, as a
   function, it is fairly general and could even be used in a spectral
   model in principle.  SpatialMap is specific to DiffuseSources and
   uses a FITS image file as a template for determining the
   distribution of photons on the sky.  We won't use a SpatialMap
   source in this example, but here is the XML source definition that
   one would use to describe the EGRET diffuse model:
   \n\n
   @verbatim
  <source name="Galactic Diffuse Emission" type="DiffuseSource">
    <spectrum type="PowerLaw">
      <parameter max="1000" min="0.001" free="1" name="Prefactor" scale="0.001" value="6.31" />
      <parameter max="-1" min="-3.5" free="0" name="Index" scale="1" value="-2.1" />
      <parameter max="200" min="50" free="0" name="Scale" scale="1" value="100" />
    </spectrum>
    <spatialModel file="$(LIKELIHOODROOT)/src/test/Data/gas.cel" type="SpatialMap">
      <parameter max="1000" min="0.001" free="0" name="Prefactor" scale="1" value="1" />
    </spatialModel>
  </source>
   @endverbatim
   The EGRET diffuse model is given in the FITS file gas.cel,
   which describes the Galactic diffuse emission as the photon
   distribution template.
   \n\n
   Both spectrum models and spatialModels are described by functions
   and various model parameters.  Each parameter is described by a
   specific set of attributes.  The actual value of a given parameter
   that is used in the calculation is the "value" attribute multiplied
   by the "scale" attribute.  The value attribute is what the
   optimizers see.  Using the scale attribute is necessary to ensure
   that the parameters describing the objective function,
   -log(likelihood) for this application, all have values lying
   roughly within a couple orders of magnitude of each other.
   \n\n
   The units for the spectral models are \f${\rm cm}^{-2} {\rm s}^{-1}
   {\rm MeV}^{-1}\f$ for PointSources and \f${\rm cm}^{-2} {\rm
   s}^{-1} {\rm MeV}^{-1} {\rm sr}^{-1}\f$ for DiffuseSources. 
   The "Prefactor" values in the PowerLaw models are the function
   values evaluated at the "Scale" values, e.g.,
   \f[
   {\rm PowerLaw}(x) = {\rm Prefactor} \times 
                       \left(\frac{x}{\rm Scale}\right)^{\rm Index}.
   \f]
   Each parameter has a range of valid values that can be specified.
   This is important for the Prefactor parameters since negative
   values are not allowed for source fluxes.  Lastly, there is a
   "free" attribute that determines whether the parameter will be
   allowed to be fixed or free in the fitting process.  Presently, the
   free flag attributes are disabled for the SkyDirFunction "RA" and
   "DEC" parameters since fitting for PointSource locations has not
   yet been implemented.

   @section runLikelihood Run likelihoodApp

   We are now ready to run the likelihood application.  One small
   detail remains.  Depending on the size of the dataset and a hidden
   parameter in the <a
   href="http://glast.stanford.edu/cgi-bin/cvsweb/observationSim/data/obsSim.par?cvsroot=CVS_SLAC">obsSim.par</a>
   file, @b obsSim may write the data to a series of files.  One may
   either merge these files using the @b fmerge FTOOL, or, at the
   prompt for the spacecraft or event files, one may instead provide
   the name of an ascii file containing the names of the series of
   FITS files.  On linux, one can easily create such a file by doing:
   @verbatim
   > ls -1 test_events*.fits > event_files
   > cat event_files
   test_events_0000.fits
   test_events_0001.fits
   >
   @endverbatim
   Having done this, if desired one can use @b fmerge to create a
   single event file:
   @verbatim
   > fmerge event_files test_events.fits
   @endverbatim
   For the spacecraft files, it is absolutely imperative that they be
   specified in time-order since that is the order in which the
   application ingests the data.  If the files are not in time-order,
   the application will complain and terminate.  For the DC1 data, we
   will provide the spacecraft data file.
   \n\n
   Finally, we are ready to run the application:
   @verbatim
> time likelihoodApp.exe
ROI cuts file [RoiCuts.xml] : 
Spacecraft file [test_scData_0000.fits] : 
Exposure file [anticenter_expMap.fits] : 
Response functions to use <FRONT/BACK|FRONT|BACK> [FRONT/BACK] : 
Source model file [anticenter_model.xml] : 
Event file [test_events_0000.fits] : 
Optimizer <LBFGS|MINUIT|DRMNGB> [DRMNGB] : 
Optimizer verbosity [0] : 
Fit tolerance [0.0001] : 
Source model output file [anticenter_model.xml] : 
Use OptEM? [no] : 
flux-style output file name [flux_model.xml] : 
Allow for refitting? [yes] : 
   @endverbatim
   Most of the entries prompted for are fairly obvious.  In addition
   to the various xml and FITS files, the user is prompted for a
   choice of IRFs, optimizer, and some output file names.  If "Source
   model output file" is set to the same name as given to "Source
   model file" ("anticenter_model.xml" in this case), then that file
   will be overwritten with the latest fit.  This is extremely
   convenient as it allows the user to edit that file while the
   application is waiting at the "Refit? [y]" prompt so that
   parameters can be adjusted and set free or fixed.  This would be
   similar to the use of the "newpar", "freeze", and "thaw" commands
   of <a href="http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/index.html">
   XSPEC</a>.  A prototype GUI is available as unsupported software
   during DC1 to ease this process.
   \n\n
   The "flux-style output file name" specifies the destination of xml
   definitions of sources that can be used with @b obsSim (or @b Gleam)
   to simulate an observation of the fitted source model.  We will use
   this output below for a comparison to the data.
   \n\n
   At this point, the application reads in the spacecraft and event
   data and computes exposures for each point source.  There is some
   initial delay for the exposure for the first point source as the
   application reads in the IRFs and precomputes some integrals of the
   PSF over the ROI.
   @verbatim
LogLike::getEvents:
Out of 18631 events in file test_events_0000.fits,
 1550 were accepted, and 17081 were rejected.

Creating source named Crab_fit
Computing exposure at (83.57, 22.01).....................!
Creating source named Extragalactic Diffuse Emission
Creating source named Geminga_fit
Computing exposure at (98.49, 17.86).....................!
Creating source named PKS0528p134
Computing exposure at (82.74, 13.38).....................!
adding source Crab_fit
adding source Extragalactic Diffuse Emission
adding source Geminga_fit
adding source PKS0528p134
Computing Event responses for the DiffuseSources.....................!
Drmngb return code: 4
Crab_fit:
Prefactor: 24.8706 +/- 2.60812
Index: -2.07227 +/- 0.0647763
Scale: 100
Npred: 233.069

Extragalactic Diffuse Emission:
Prefactor: 1.09014 +/- 0.066217
Index: -1.96779 +/- 0.0345433
Scale: 100
Npred: 763.044

Geminga_fit:
Prefactor: 23.9397 +/- 1.96167
Index: -1.64543 +/- 0.0349701
Scale: 100
Npred: 399.933

PKS0528p134:
Prefactor: 13.2709 +/- 2.19185
Index: -2.23041 +/- 0.112889
Scale: 100
Npred: 105.259

Writing fitted model to anticenter_model.xml
Refit? [y] 
   @endverbatim
   Here we see the results of the fit.  We won't compare in detail
   these numbers to the input values right now; however, for
   reference, the input values of the spectral indices to each of the
   sources are Crab: -2.19, Geminga: -1.66, PKS 0528+134: -2.46, and
   extragalactic diffuse: -2.1.  As an exercise, we will edit the
   model file, anticenter_model.xml, and fix the extragalactic diffuse 
   spectral index to the input value.  Here is what we would put:
   @verbatim
  <source name="Extragalactic Diffuse Emission" type="DiffuseSource">
    <spectrum type="PowerLaw">
      <parameter max="100" min="1e-05" free="1" name="Prefactor" scale="1e-07" value="1.09014" />
      <parameter max="-1" min="-3.5" free="0" name="Index" scale="1" value="-2.1" />
      <parameter max="200" min="50" free="0" name="Scale" scale="1" value="100" />
    </spectrum>
   @endverbatim
   Note that after having re-read this file the Prefactor parameter
   was over-written by the application to the fitted value.  For the
   Index parameter, we have set the free flag and value attributes.
   At this point, we hit "<return>" to instruct the application to fit
   again.
   @verbatim
Drmngb return code: 4
Crab_fit:
Prefactor: 24.892 +/- 2.64623
Index: -2.06938 +/- 0.0651139
Scale: 100
Npred: 234.076

Extragalactic Diffuse Emission:
Prefactor: 1.28741 +/- 0.0503393
Index: -2.1
Scale: 100
Npred: 770.978

Geminga_fit:
Prefactor: 23.9601 +/- 1.97848
Index: -1.64437 +/- 0.035074
Scale: 100
Npred: 401.04

PKS0528p134:
Prefactor: 13.2894 +/- 2.24432
Index: -2.22737 +/- 0.114074
Scale: 100
Npred: 105.737

Writing fitted model to anticenter_model.xml
Refit? [y] 
n
Writing flux-style xml model file to flux_model.xml
97.230u 0.260s 2:35.85 62.5%    0+0k 0+0io 4174pf+0w
   @endverbatim
   Fixing the extragalactic diffuse spectral index doesn't really
   affect the other sources in the model, but it does in fact improve
   the Prefactor normalization:  the input model value is 1.32.

   @section checkWithObsSim Use obsSim to Check the Model Fit 
   In practice, of course, we will not know what the true parameter
   values are for the sources we are fitting.  One way of comparing
   the fitted model to the data is to create a simulated data set and
   compare that to the real data.  For this purpose, we can use the
   flux_model.xml file created by the application:
   @verbatim
<source_library title="Likelihood_model">
  <source flux="0.0232931" name="Crab_fit">
    <spectrum escale="MeV">
      <particle name="gamma">
        <power_law emax="316230" emin="100" gamma="2.06938" />
      </particle>
      <celestial_dir ra="83.57" dec="22.01" />
    </spectrum>
  </source>
  <source flux="0.117125" name="Extragalactic_Diffuse_Emission">
    <spectrum escale="MeV">
      <particle name="gamma">
        <power_law emax="316230" emin="100" gamma="2.1" />
      </particle>
      <solid_angle maxcos="1.0" mincos="-0.4" />
    </spectrum>
  </source>
  <source flux="0.0369993" name="Geminga_fit">
    <spectrum escale="MeV">
      <particle name="gamma">
        <power_law emax="316230" emin="100" gamma="1.64437" />
      </particle>
      <celestial_dir ra="98.49" dec="17.86" />
    </spectrum>
  </source>
  <source flux="0.0108376" name="PKS0528p134">
    <spectrum escale="MeV">
      <particle name="gamma">
        <power_law emax="316230" emin="100" gamma="2.22737" />
      </particle>
      <celestial_dir ra="82.74" dec="13.38" />
    </spectrum>
  </source>
  <source name="all_in_flux_model.xml">
    <nestedSource sourceRef="Crab_fit" />
    <nestedSource sourceRef="Extragalactic_Diffuse_Emission" />
    <nestedSource sourceRef="Geminga_fit" />
    <nestedSource sourceRef="PKS0528p134" />
  </source>
</source_library>
   @endverbatim
   We need to hand edit the obsSim.par file and give it this file name, or
   alternatively, one could set the second entry following "XML_source_file"
   to "ql" as we've done here:
   @verbatim
XML_source_file,f,ql,"flux_model.xml",,,"File of flux-style source definitions"
Source_list,f,ql,"source_names.dat",,,"File containing source names"
Number_of_events,r,ql,86400,1,4e7,Number of events (or simulation time in seconds)
Use_as_sim_time,b,ql,yes,,,Use number of events as simulation time?
Response_functions,s,ql,"FRONT/BACK",FRONT/BACK|FRONT|BACK,,"Response functions to use"
Output_file_prefix,s,ql,"all_3EG",,,"Prefix for output files"
Maximum_effective_area,r,hl,1.21,,,Maximum effective area value
Start_time,r,hl,0,,,Simulation start time (seconds)
Pointing_history_file,s,ql,"none",,,"Pointing history file"
Maximum_number_of_rows,i,hl,200000,,,Maximum number of rows in FITS files
   @endverbatim
   We then need to edit the source_names.dat file to contain:
   @verbatim
   all_in_flux_model.xml
   @endverbatim
   run @b obsSim,
   @verbatim
   > obsSim.exe 
   File of flux-style source definitions [flux_model.xml] : 
   File containing source names [source_names.dat] : 
   Number of events (or simulation time in seconds) <1 - 4e7> [86400] : 
   Use number of events as simulation time? [yes] : 
   Response functions to use <FRONT/BACK|FRONT|BACK> [FRONT/BACK] : 
   Prefix for output files [model_fit] : 
   Pointing history file [none] : 
   Generating events for a simulation time of 86400 seconds....
   Done.
   @endverbatim
   create a counts map having the same geometry as anticenter_1day.fits,
   @verbatim
   > fcopy "model_fit_events_0000.fits[1][bin ra=65:105,dec=2:42]" model_fit_1day.fits
   @endverbatim
   and compare using ds9,

   @image html anticenter_1day.png

   @section makeTsMaps Make Test-Statistic Maps

   For serious work, one would like to find sources near the detection
   limit of the instrument.  The procedure used in the analysis of
   EGRET data was to model the strongest, most obvious sources (with
   some theoretical prejudice as to the true source positions, e.g.,
   assuming that most variable high Galactic latitude sources are
   blazars which can be localized by radio, optical, or X-ray
   observations), and then create "Test-statistic maps" to search for
   unmodeled point sources.  These TS maps are created by moving a
   putative point source through a grid of locations on the sky and
   maximizing -log(likelihood) at each grid point, with the other,
   stronger, and presumably well-identified sources included in each
   of fit.  New, fainter sources are then identified at local maxima
   of the TS map.

   Since we know there are only three point sources in data set, as an
   exercise we will comment-out the weakest source of the three from the
   model file anticenter_model.xml:
   @verbatim
   <!--
  <source name="PKS0528p134" type="PointSource">
    <spectrum type="PowerLaw">
      <parameter max="1000" min="0.001" free="1" name="Prefactor" scale="1e-09" value="13.2894" />
      <parameter max="-1" min="-3.5" free="1" name="Index" scale="1" value="-2.22737" />
      <parameter max="200" min="50" free="0" name="Scale" scale="1" value="100" />
    </spectrum>
    <spatialModel type="SkyDirFunction">
      <parameter max="3.40282e+38" min="-3.40282e+38" free="0" name="RA" scale="1" value="82.74" />
      <parameter max="3.40282e+38" min="-3.40282e+38" free="0" name="DEC" scale="1" value="13.38" />
    </spatialModel>
  </source>
-->
</source_library>
   @endverbatim
   We then create a TS map (with the same dimensions as our counts
   maps) to see how well PKS 0528+134 is identified using this method:
   @verbatim
   > time TsMap.exe
   ROI cuts file [RoiCuts.xml] : 
   Spacecraft file [all_sky_1day_scData_0000.fits] : test_scData_0000.fits
   Exposure file [anticenter_expMap.fits] : 
   Response functions to use <FRONT/BACK|COMBINED> [FRONT/BACK] : 
   Source model file [anticenter_model.xml] : 
   Event file [all_sky_1day_events_0000.fits] : test_events_0000.fits
   Optimizer <LBFGS|MINUIT|DRMNGB> [DRMNGB] : 
   Optimizer verbosity [0] : 
   Fit tolerance [0.001] : 
   Coordinate system <CEL|GAL> [CEL] : 
   Longitude maximum <-360 - 360> [105] : 
   Longitude minmum <-360 - 360> [65] : 
   Number of longitude points <2 - 200> [40] : 
   Latitude maximum <-90 - 90> [42] : 
   Latitude minimum <-90 - 90> [2] : 
   Number of latitude points <2 - 200> [40] : 
   TS map file name [TsMap.fits] : anticenter_TsMap.fits
   Creating source named Crab_fit
   Computing exposure at (83.57, 22.01).....................!
   Creating source named Extragalactic Diffuse Emission
   Creating source named Geminga_fit
   Computing exposure at (98.49, 17.86).....................!
   adding source Crab_fit
   adding source Extragalactic Diffuse Emission
   adding source Geminga_fit
   LogLike::getEvents:
   Out of 18631 events in file test_events_0000.fits,
    1550 were accepted, and 17081 were rejected.

   <...lots of output skipped...>

   Computing exposure at (105, 39.9487).....................!
   Computing exposure at (105, 40.9744).....................!
   Computing exposure at (105, 42).....................!
   4084.230u 25.040s 1:13:56.05 92.6%      0+0k 0+0io 6381pf+0w
   @endverbatim
   Note that no PKS0529p134 source is created or added in the setup to
   doing the fits over the grid of sky locations.  Here is the
   resulting TsMap file compared to the original counts map:

   @image html TsMap_compare.png
   \n
   The green circle in the left image appears at the location of PKS
   0528+134, (ra, dec) = (82.74, 13.38) degrees.

*/

/** @page All_3EG
   @section all3EG all_3EG_sources + diffuse_100mev
   @image html allsky_3EG.png allsky_3EG.fits
   
*/
