/**
   @page tutorial Likelihood Tutorial

   @section intro Introduction

   In order to illustrate how to use the Likelihood software, this
   narrative gives a step-by-step description for performing a
   relatively simple analysis.  We will assume that you already have
   event (FT1) and spacecraft (FT2) data files in hand.

   Here are each of the steps:

   - @ref makeSubselections Since there is computational overhead
          for each event associated with each diffuse component, it is
          useful to filter out any events that are not within the
          extraction region used for the analysis.
   - @ref makeCountsMaps These simple FITS images let us see what 
          we've got and help to pick out obvious candidate sources.
   - @ref defineROI This is the extraction region in terms of apparent 
          photon direction, energy, and arrival time.
   - @ref makeExposureMap This is needed for analyzing diffuse sources
          and so is required for almost any sort of analysis.
   - @ref sourceModelFile The source model xml file contains the
          various sources and their model parameters to be fit
          using the @b likelihood tool.
   - @ref computeDiffuseResps Each event must have a separate response
          precomputed for each diffuse component in the source model.
   - @ref runLikelihood
   - @ref makeTsMaps This is used for point source localization and for
          finding weaker sources after the stronger sources have
          been modeled.

   @section makeSubselections Make Subselections from the Event Data.
   For this tutorial, we will use data generated using @b obsSim for a
   model comprising the Third EGRET (3EG) catalog sources and Galactic
   diffuse and extragalactic diffuse emission for a simulation time of
   one day.  Since we have run @b obsSim to create the event data, we
   have opted to write events from each component into separate FT1
   files:
   @verbatim
   noric13[jchiang] ls -l *events*fits
   -rw-rw-r--    1 jchiang  glast-data  2894400 Oct 11 15:22 eg_diffuse_events_0000.fits
   -rw-rw-r--    1 jchiang  glast-data  7243200 Oct 11 15:26 galdiffuse_events_0000.fits
   -rw-rw-r--    1 jchiang  glast-data  1088640 Oct 11 15:20 ptsrcs_events_0000.fits
   noric13[jchiang] 
   @endverbatim
   We will consider data in the Virgo region within a 20 degree
   acceptance cone of the blazar 3C 279.  Here we apply @b
   dataSubselector to the point source data (selecting events within 25
   degrees of 3C 279 to ensure that we don't miss any near the edge):
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
   entries of "0" for most of the parameters means that no cuts will
   be applied.  This tool, its interface in particular, is a prime
   target for refactoring.

   @section makeCountsMaps Make Counts Maps from the Event Files.  The
   next step is to create a simple counts map file to visualize the
   extraction region we wish to fit.  There are presently three tools
   available to perform this task: @b evtbin in the @b evtbin package,
   @b count_map in the @b map_tools package, and @b gtcntsmap in the
   Likelihood package.  For this step, since the event data are
   contained in three separate FITS files, @b gtcntsmap is the most
   convenient to use.  We proceed by creating a file which is a list
   of the files to be binned and then run the tool on that file list:
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

   For the spacecraft data file, we have entered the one generated by
   @b obsSim for the 3EG sources.  Identical files were created for
   the diffuse component runs, and we could also have used either of
   those.  The last command launches the visualization tool <a
   href="http://hea-www.harvard.edu/RD/ds9/">ds9</a> and produces this
   display:

   @image html Virgo_map.png Virgo_map.fits

   Here I have plotted as the green circles the locations of the 3EG
   point sources within 25 degrees of 3C 279.

   @section defineROI Define the Region-of-Interest. 

   The Region-of-Interest (ROI) comprises a set of selections on
   photon arrival time, energy, and direction.  These selections are
   made in addition to any that are made using the query to the <a
   href="http://glast.gsfc.nasa.gov/cgi-bin/ssc/U1/D1WebQuery.cgi">
   GSSC database</a> or via the @b dataSubselector user-level
   selection tool.  Unfortunately, because of an oversight in the
   Likelihood implementation, the ROI selections must @em include any
   that are made by these other tools.  More precisely, the data-space
   defined by the ROI must lie entirely within the intersection of the
   data-spaces defined by these other tools.  For example, if the
   user-level selection tool selects events within a time interval
   (t1, t2), then any time intervals defined in the ROI file must lie
   entirely within (t1, t2).  Similar considerations apply to
   energy and acceptance cone cuts made using @b dataSubselector.  In
   future releases, this will be corrected.
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
   file(s) will provide this same capability. The "coordsys" attribute
   can either be "Galactic" or "J2000"; the "radius" attribute gives
   the half-opening angle of the acceptance cone; and units are in
   degrees.

   @section makeExposureMap Make an Exposure Map for a Given ROI.
   With the ROI specified, we are now ready to create an exposure map.
   The type of exposure map used by Likelihood differs significantly
   from the usual notion of exposure maps, which are essentially
   integrals of effective area over time.  The exposure calculation
   that Likelihood uses consists of an integral of the total response
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
   @b makeExposureCube also applies the time-interval cuts specified in
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
   file of this tool can be used to generating exposure maps for
   regions-of-interest in other parts of the sky that have the same
   time interval selections.  Although the @b expMap application (see
   below) can generate exposure maps for Likelihood without an
   exposure cube map, using one affords a substantial time
   savings.
   \n\n
   Creating the exposure map using the @b expMap tool, we have
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
   ensure that photons from sources outside the ROI are accounted for
   owing to the size of the instrument point-spread function.
   Half-degree pixels are a nominal choice; smaller pixels should
   result in a more accurate evaluation of the diffuse source fluxes
   but will also make the exposure map calculation itself lengthier.
   The number of energies specifies the number of logarithmically
   spaced intervals bounded by the energy range given in the ROI cuts
   xml file.  Here is one image plane of the exposure map we just
   created:

   @image html expMap.png expMap.fits

   @section sourceModelFile Create a Source Model XML File.

   Just as it does for ROI information, the @b likelihood tool reads
   the source model from an XML file.  Given the dearth of bright
   sources in the extraction region we have selected, our source model
   file will be fairly simple, consisting only of the Galactic and
   extragalactic diffuse emission, and the point sources 3C 279 and 3C 273:
   @verbatim
   <?xml version="1.0" ?>
   <source_library title="source library">
     <source name="Extragalactic Diffuse" type="DiffuseSource">
       <spectrum type="PowerLaw">
         <parameter free="1" max="100.0" min="1e-05" name="Prefactor" scale="1e-07" value="1.6"/>
         <parameter free="0" max="-1.0" min="-3.5" name="Index" scale="1.0" value="-2.1"/>
         <parameter free="0" max="200.0" min="50.0" name="Scale" scale="1.0" value="100.0"/>
       </spectrum>
       <spatialModel type="ConstantValue">
         <parameter free="0" max="10.0" min="0.0" name="Value" scale="1.0" value="1.0"/>
       </spatialModel>
     </source>
     <source name="Galactic Diffuse" type="DiffuseSource">
       <spectrum type="PowerLaw">
         <parameter free="1" max="1000.0" min="0.001" name="Prefactor" scale="0.001" value="11.0"/>
         <parameter free="0" max="-1.0" min="-3.5" name="Index" scale="1.0" value="-2.1"/>
         <parameter free="0" max="200.0" min="50.0" name="Scale" scale="1.0" value="100.0"/>
       </spectrum>
       <spatialModel file="$(LIKELIHOODROOT)/src/test/Data/gas.cel" type="SpatialMap">
         <parameter free="0" max="1000.0" min="0.001" name="Prefactor" scale="1.0" value="1.0"/>
          </spatialModel>
     </source>
     <source name="3C 273" type="PointSource">
       <spectrum type="PowerLaw">
         <parameter free="1" max="1000.0" min="0.001" name="Prefactor" scale="1e-09" value="10"/>
         <parameter free="1" max="-1.0" min="-5.0" name="Index" scale="1.0" value="-2.1"/>
         <parameter free="0" max="2000.0" min="30.0" name="Scale" scale="1.0" value="100.0"/>
       </spectrum>
       <spatialModel type="SkyDirFunction">
         <parameter free="0" max="360" min="-360" name="RA" scale="1.0" value="187.25"/>
         <parameter free="0" max="90" min="-90" name="DEC" scale="1.0" value="2.17"/>
       </spatialModel>
     </source>
     <source name="3C 279" type="PointSource">
       <spectrum type="PowerLaw">
         <parameter free="1" max="1000.0" min="0.001" name="Prefactor" scale="1e-09" value="10"/>
         <parameter free="1" max="-1.0" min="-5.0" name="Index" scale="1.0" value="-2"/>
         <parameter free="0" max="2000.0" min="30.0" name="Scale" scale="1.0" value="100.0"/>
       </spectrum>
       <spatialModel type="SkyDirFunction">
         <parameter free="0" max="360" min="-360" name="RA" scale="1.0" value="193.98"/>
         <parameter free="0" max="90" min="-90" name="DEC" scale="1.0" value="-5.82"/>
       </spatialModel>
     </source>
   </source_library>
   @endverbatim
   We'll call this source model file "Virgo_model.xml".
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
   distribution of photons on the sky.  The EGRET diffuse model is
   given in the FITS file gas.cel, which describes the Galactic
   diffuse emission, as the photon distribution template.
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

   @section computeDiffuseResps Compute the Diffuse Source Responses.
   If these quantities are not precomputed using the 
   @b diffuseResponses tool, then @b likelihood will compute them at
   runtime.  However, if multiple fits and/or sessions with the 
   @b likelihood tool are anticipated, it is probably wise to precompute
   these quantities.  The source model xml file must contain all of
   the diffuse sources to be fit.  The @b diffuseResponses tool will
   add columns to the FT1 file for each diffuse source:

   @verbatim
   salathe[jchiang] ~/ST/Likelihood/v5r0/rh9_gcc32/diffuseResponses.exe 
   Spacecraft file [oneday_scData_0000.fits] : ptsrcs_scData_0000.fits
   Source model file [diffuse_model.xml] : Virgo_model.xml
   Event file [oneday_events_0000.fits] : ptsrcs_events_filtered.fits
   ROI cuts file [RoiCuts.xml] : 
   adding source Extragalactic Diffuse
   adding source Galactic Diffuse
   adding source 3C 273
   adding source 3C 279
   .....................!
   salathe[jchiang] 
   @endverbatim

   @section runLikelihood Run likelihood.

   We are now ready to run the @b likelihood application:
   @verbatim
salathe[jchiang] ~/ST/Likelihood/v5r0/rh9_gcc32/likelihood.exe 
Response functions to use <FRONT/BACK|FRONT|BACK|GLAST25> [FRONT/BACK] : 
Use energy dispersion? [no] : 
Source model file [Virgo_model.xml] : 
Source model output file [none] : 
flux-style output file name [flux_model.xml] : 
Statistic to use <BINNED|UNBINNED|OPTEM> [UNBINNED] : 
Optimizer <LBFGS|MINUIT|DRMNGB> [MINUIT] : 
Optimizer verbosity [1] : 
Fit tolerance [0.0001] : 
Allow for refitting? [yes] : 
ROI cuts file [RoiCuts.xml] : 
Spacecraft file [ptsrcs_scData_0000.fits] : 
Event file [eventFiles] : 
Exposure file [expMap.fits] : 
Counts map file [srcMaps.fits] : none
Exposure hypercube file [expcube_1_day.fits] : none
   @endverbatim
   Most of the entries prompted for are fairly obvious.  In addition
   to the various xml and FITS files, the user is prompted for a
   choice of IRFs, the type of statistic to use, the optimizer, 
   and some output file names.  

   Four response function combinations are available:

   - @b FRONT DC1 IRFs for events converting in the "front" part of the LAT
           (layers 0--11).
   - @b BACK  DC1 IRFs for "back"-converting events (layers 12-15)
   - @b FRONT/BACK Both the DC1 FRONT and BACK IRFs
   - @b GLAST25 AO era IRFs for the FRONT and BACK of the LAT

   @b obsSim has the same option of generating events for these
   combinations of IRFs and the choice used for @b likelihood must be
   the same as that chosen for @b obsSim.  The @b dataSubselector tool
   also allows one to filter on conversion layer and so can allow one
   to select FRONT or BACK only datasets.

   Three statistics are available:

   - @b UNBINNED This is a standard unbinned analysis, described in this
     tutorial.  If this option is chosen then parameters for the ROI
     cuts file, Spacecraft file, Event file, and Exposure file must be
     given.

   - @b BINNED This is a binned analysis, which is still in development.
     In this case, the Counts map file and Exposure hypercube file
     must be given.  The ROI cuts file Spacecraft file, Event file and
     Exposure file are ignored in this case.

   - @b OPTEM This uses an expectation-maximization algorithm to
     partition the events, effectively, among the individual sources
     so that the fitting procedure can be performed for each source in
     isolation rather than for all the parameters at once.  The cost
     is that the algorithm must loop over all the sources several
     times until it converges.  This is also an unbinned analysis so
     the parameters relevant for this algorithm are the same as for
     @b UNBINNED.

   There are three optimizers from which to choose.  @b Minuit will
   generally give the best performance.

   If "Source model output file" is set to the same name as given to
   "Source model file" ("Virgo_model.xml" in this case), then that
   file will be overwritten with the latest fit.  This is 
   convenient as it allows the user to edit that file while the
   application is waiting at the "Refit? [y]" prompt so that
   parameters can be adjusted and set free or fixed.  This would be
   similar to the use of the "newpar", "freeze", and "thaw" commands
   of <a
   href="http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/index.html">
   XSPEC</a>.
   \n\n
   The "flux-style output file name" specifies the destination of xml
   definitions of sources that can be used with @b obsSim (or @b Gleam)
   to simulate an observation of the fitted source model.
   \n\n
   The application proceeds by reading in the spacecraft and event
   data, and if necessary, computing event responses for each diffuse
   source.
   @verbatim
LogLike::getEvents:
Out of 1292 events in file eg_diffuse_events_filtered.fits,
 814 were accepted, and 478 were rejected.

LogLike::getEvents:
Out of 759 events in file galdiffuse_events_filtered.fits,
 481 were accepted, and 278 were rejected.

LogLike::getEvents:
Out of 305 events in file ptsrcs_events_filtered.fits,
 257 were accepted, and 48 were rejected.

Computing Event responses for the DiffuseSources.....................!
 **********
 **    1 **SET PRINT    0.000    
 **********
 **********
 **    2 **SET NOWARN 
 **********

 PARAMETER DEFINITIONS:
    NO.   NAME         VALUE      STEP SIZE      LIMITS
     1 'Prefactor '    1.6000       1.0000        0.10000E-04   100.00    
     2 'Prefactor '    11.000       1.0000        0.10000E-02   1000.0    
     3 'Prefactor '    2.4330       1.0000        0.10000E-02   1000.0    
     4 'Index     '   -2.5800       1.0000        -5.0000      -1.0000    
     5 'Prefactor '    7.1230       1.0000        0.10000E-02   1000.0    
     6 'Index     '   -1.9600       1.0000        -5.0000      -1.0000    
 **********
 **    3 **SET ERR   0.5000    
 **********
 **********
 **    4 **SET GRAD    1.000    
 **********
 **********
 **    5 **MIGRAD    200.0       3173.    
 **********

 MIGRAD MINIMIZATION HAS CONVERGED.

 MIGRAD WILL VERIFY CONVERGENCE AND ERROR MATRIX.

 FCN=   15861.44     FROM MIGRAD    STATUS=CONVERGED     84 CALLS       85 TOTAL
                     EDM=  0.27E+01    STRATEGY= 1      ERR MATRIX NOT POS-DEF

  EXT PARAMETER                APPROXIMATE        STEP         FIRST   
  NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE 
   1  Prefactor     1.6899       0.15110       0.14753E-01    12.658    
   2  Prefactor     9.2402       0.92370       0.95978E-02    23.084    
   3  Prefactor     4.7127        17.650       0.20766E-01   0.27761    
   4  Index        -2.8179        3.0898       0.37782       0.90273    
   5  Prefactor     11.038        3.1939       0.11921E-01   -2.8765    
   6  Index        -2.2091       0.16818       0.11625       -1.4825    
                               ERR DEF= 0.500    
 **********
 **    6 **HESSE 
 **********

 FCN=   15861.44     FROM HESSE     STATUS=OK            42 CALLS      127 TOTAL
                     EDM=  0.15E+00    STRATEGY= 1      ERROR MATRIX ACCURATE 

  EXT PARAMETER                                INTERNAL      INTERNAL  
  NO.   NAME        VALUE          ERROR       STEP SIZE       VALUE   
   1  Prefactor     1.6899       0.19468       0.65213E-04   -1.3101    
   2  Prefactor     9.2402        2.2072       0.10257E-03   -1.3783    
   3  Prefactor     4.7127        1.4955       0.22439E-03   -1.4334    
   4  Index        -2.8179       0.30961       0.14985E-01   0.91172E-01
   5  Prefactor     11.038        1.9122       0.76577E-03   -1.3603    
   6  Index        -2.2091       0.11027       0.51120E-03   0.40656    
                               ERR DEF= 0.500    
Final values: 
  Prefactor  = 1.68988
  Prefactor  = 9.24022
  Prefactor  = 4.71275
  Index      = -2.81791
  Prefactor  = 11.0378
  Index      = -2.2091
Minuit fit quality: 3   estimated distance: 0.149679
Minuit parameter uncertainties:
  1  0.194686
  2  2.20744
  3  1.49564
  4  0.310871
  5  1.91235
  6  0.11034
Computing TS values for each source (4 total)
....!

Extragalactic Diffuse:
Prefactor: 1.68988 +/- 0.194686
Index: -2.1
Scale: 100
Npred: 945.06

Galactic Diffuse:
Prefactor: 9.24022 +/- 2.20744
Index: -2.1
Scale: 100
Npred: 424.658

3C 273:
Prefactor: 4.71275 +/- 1.49564
Index: -2.81791 +/- 0.310871
Scale: 100
Npred: 44.9021
TS value: 23.6949

3C 279:
Prefactor: 11.0378 +/- 1.91235
Index: -2.2091 +/- 0.11034
Scale: 100
Npred: 135.984
TS value: 215.798

-log(Likelihood): 15861.4

Refit? [y] 

   @endverbatim
   At this point we can choose to edit the source model file, setting
   various spectral parameters free or fixed, and in the latter case,
   also setting them to some nominal value.  For our refit, we fix the
   diffuse component prefactors to their nominal values as follows:
   @verbatim
  <source name="Extragalactic Diffuse" type="DiffuseSource">
    <spectrum type="PowerLaw">
      <parameter free="0" max="100.0" min="1e-05" name="Prefactor" scale="1e-07" value="1.45"/>
      <parameter free="0" max="-1.0" min="-3.5" name="Index" scale="1.0" value="-2.1"/>
      <parameter free="0" max="200.0" min="50.0" name="Scale" scale="1.0" value="100.0"/>
    </spectrum>
    <spatialModel type="ConstantValue">
      <parameter free="0" max="10.0" min="0.0" name="Value" scale="1.0" value="1.0"/>
    </spatialModel>
  </source>
  <source name="Galactic Diffuse" type="DiffuseSource">
    <spectrum type="PowerLaw">
      <parameter free="0" max="1000.0" min="0.001" name="Prefactor" scale="0.001" value="11.0"/>
      <parameter free="0" max="-1.0" min="-3.5" name="Index" scale="1.0" value="-2.1"/>
      <parameter free="0" max="200.0" min="50.0" name="Scale" scale="1.0" value="100.0"/>
    </spectrum>
    <spatialModel file="$(LIKELIHOODROOT)/src/test/Data/gas.cel" type="SpatialMap">
      <parameter free="0" max="1000.0" min="0.001" name="Prefactor" scale="1.0" value="1.0"/>
    </spatialModel>
  </source>
   @endverbatim
   We hit "return" to instruct the application to fit again. Here is the
   output:
   @verbatim
 **********
 **    1 **SET PRINT    0.000    
 **********
 **********
 **    2 **SET NOWARN 
 **********

 PARAMETER DEFINITIONS:
    NO.   NAME         VALUE      STEP SIZE      LIMITS
     1 'Prefactor '    2.4330       1.0000        0.10000E-02   1000.0    
     2 'Index     '   -2.5800       1.0000        -5.0000      -1.0000    
     3 'Prefactor '    7.1230       1.0000        0.10000E-02   1000.0    
     4 'Index     '   -1.9600       1.0000        -5.0000      -1.0000    
 **********
 **    3 **SET ERR   0.5000    
 **********
 **********
 **    4 **SET GRAD    1.000    
 **********
 **********
 **    5 **MIGRAD    200.0       3174.    
 **********

 MIGRAD MINIMIZATION HAS CONVERGED.

 MIGRAD WILL VERIFY CONVERGENCE AND ERROR MATRIX.

 FCN=   15862.83     FROM MIGRAD    STATUS=CONVERGED     55 CALLS       56 TOTAL
                     EDM=  0.21E+01    STRATEGY= 1      ERR MATRIX NOT POS-DEF

  EXT PARAMETER                APPROXIMATE        STEP         FIRST   
  NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE 
   1  Prefactor     5.3499        18.836       0.20766E-01   -1.5066    
   2  Index        -2.8663        3.0634       0.36280       0.82559    
   3  Prefactor     11.471        6.5707       0.11921E-01   -15.667    
   4  Index        -2.2263       0.33522       0.11611       -3.1938    
                               ERR DEF= 0.500    
 **********
 **    6 **HESSE 
 **********

 FCN=   15862.83     FROM HESSE     STATUS=OK            25 CALLS       81 TOTAL
                     EDM=  0.51E-01    STRATEGY= 1      ERROR MATRIX ACCURATE 

  EXT PARAMETER                                INTERNAL      INTERNAL  
  NO.   NAME        VALUE          ERROR       STEP SIZE       VALUE   
   1  Prefactor     5.3499        1.4291       0.21879E-03   -1.4244    
   2  Index        -2.8663       0.28400       0.14565E-01   0.66878E-01
   3  Prefactor     11.471        1.8645       0.76635E-03   -1.3562    
   4  Index        -2.2263       0.10781       0.51061E-03   0.39719    
                               ERR DEF= 0.500    
Final values: 
  Prefactor  = 5.34992
  Index      = -2.86634
  Prefactor  = 11.4714
  Index      = -2.22634
Minuit fit quality: 3   estimated distance: 0.0513511
Minuit parameter uncertainties:
  1  1.42914
  2  0.284963
  3  1.86459
  4  0.107871
Computing TS values for each source (4 total)
....!

Extragalactic Diffuse:
Prefactor: 1.45
Index: -2.1
Scale: 100
Npred: 810.907

Galactic Diffuse:
Prefactor: 11
Index: -2.1
Scale: 100
Npred: 505.534

3C 273:
Prefactor: 5.34992 +/- 1.42914
Index: -2.86634 +/- 0.284963
Scale: 100
Npred: 50.7384
TS value: 31.8937

3C 279:
Prefactor: 11.4714 +/- 1.86459
Index: -2.22634 +/- 0.107871
Scale: 100
Npred: 139.793
TS value: 246.244

-log(Likelihood): 15862.8

Refit? [y] 
n
Writing flux-style xml model file to flux_model.xml
salathe[jchiang] 
   @endverbatim

   @section makeTsMaps Make Test-Statistic Maps.

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

   Let's comment out 3C 273 from the source model xml file and see if we
   can find evidence for it in the data.
   @verbatim

salathe[jchiang] ~/ST/Likelihood/v5r0/rh9_gcc32/TsMap.exe
ROI cuts file [RoiCuts.xml] : 
Spacecraft file [ptsrcs_escData_0000.fits] : ptsrcs_scData_0000.fits
Exposure file [expMap.fits] : 
Response functions to use <FRONT/BACK|FRONT|BACK|GLAST25> [FRONT/BACK] : 
Source model file [Virgo_model.xml] : 
Event file [eventFiles] : 
Optimizer <LBFGS|MINUIT|DRMNGB> [MINUIT] : 
Optimizer verbosity [0] : 
Fit tolerance [0.001] : 
Coordinate system <CEL|GAL> [CEL] : 
Longitude minmum <-360 - 360> [180] : 
Longitude maximum <-360 - 360> [200] : 
Number of longitude points <2 - 200> [40] : 
Latitude minimum <-90 - 90> [-10] : 
Latitude maximum <-90 - 90> [10] : 
Number of latitude points <2 - 200> [40] : 
TS map file name [TsMap.fits] : 
LogLike::getEvents:
Out of 1292 events in file eg_diffuse_events_filtered.fits,
 814 were accepted, and 478 were rejected.

LogLike::getEvents:
Out of 759 events in file galdiffuse_events_filtered.fits,
 481 were accepted, and 278 were rejected.

LogLike::getEvents:
Out of 305 events in file ptsrcs_events_filtered.fits,
 257 were accepted, and 48 were rejected.

Computing Event responses for the DiffuseSources.....................!
   @endverbatim

   Here is the resulting TS map:

   @image html TsMap.png TsMap.fits

   The various circles indicate the locations of the 3EG sources, all
   of which were included in the simulation.  The red circle is 3C
   279, the only source that was included in the TS map model; and the
   blue circle is at the location of 3C 273.  Several other TS
   enhancements appear in the map and seem to be associated with other
   3EG sources (green circles) that were not included in the the TS
   map model.
*/
