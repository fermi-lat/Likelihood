/**
   @page tutorial Likelihood Tutorial

   @section intro Introduction

   In order to illustrate how to use the likelihood software, this
   narrative is provided and gives a step-by-step description for
   performing an analysis on a highly simplified data set.

   Here are each of the steps:
   - @ref createData (or retrieve the real thing from the web server.)
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

   @section createData Create Simulated Data Using obsSim

   In practice, the data to be analyzed will be downloaded from the 
   <a href="http://glast.gsfc.nasa.gov/cgi-bin/ssc/U1/D1WebQuery.cgi">
   GSSC server</a>.  However, here we will run the @b obsSim application to
   generate data with a known set of sources so that we can compare
   the results of our fits to the input parameters of the source model
   that was fed to the simulation.
   \n\n
   The first step is to create a list of sources.  The sources that
   are provided for simulation are given in the xml files in the <a
   href="http://glast.stanford.edu/cgi-bin/cvsweb/observationSim/xml/?cvsroot=CVS_SLAC">observationSim/xml</a>
   subdirectory.  Sources in the <a
   href="http://glast.stanford.edu/cgi-bin/cvsweb/observationSim/xml/source_library.xml?cvsroot=CVS_SLAC">source_library.xml</a>
   and <a
   href="http://glast.stanford.edu/cgi-bin/cvsweb/observationSim/xml/3EG_catalog_20-1e6MeV.xml?cvsroot=CVS_SLAC">3EG_catalog_20-1e6MeV.xml</a>
   files are provided by default.  (Note that one should only use
   photon sources from the former with @b obsSim ).  If one wishes to
   use sources from some of the other xml files in <a
   href="http://glast.stanford.edu/cgi-bin/cvsweb/observationSim/xml/?cvsroot=CVS_SLAC">observationSim/xml</a>
   or from some other source model xml file, there is a hidden
   parameter in the <a
   href="http://glast.stanford.edu/cgi-bin/cvsweb/observationSim/data/obsSim.par?cvsroot=CVS_SLAC">obsSim.par</a>
   file (<a href="http://www.slac.stanford.edu/exp/glast/ground/software/RM/documentation/ScienceTools/ScienceTools-v0r5/Likelihood/v2r6p5/userguide.html">What's a .par file?</a>) that will allow one to specify that file.  
   We will see how to
   do this when we rerun @b obsSim again below.  
   \n\n
   The list of sources should be put into a file, which we will call
   "source_names.dat".  For this simple example, we choose the
   following sources:
   @verbatim
   diffuse-100mev
   anticenter
   @endverbatim
   The first source is an isotropic diffuse component intended to
   model the extragalactic diffuse emission measured by EGRET.  The
   second source is actually a composite of the three brightest
   sources in the Galactic anti-center region.  
   An individual source is defined in the xml libraries like this:
   @verbatim
   <!-- Geminga, from the 3EG catalog -->
   <source name="Geminga" flux="0.03529">
       <spectrum escale="MeV">
           <particle name="gamma"> 
               <power_law emin="100." emax="100000." gamma="1.66"/>
           </particle>
           <celestial_dir ra="98.49" dec="17.86"/>
       </spectrum>
   </source>
   @endverbatim
   and composite sources, such as "anticenter", are defined like this:
   @verbatim
   <!-- Strong point sources in the Galactic anticenter region -->
   <source name = "anticenter">
       <nestedSource sourceRef="Crab" />
       <nestedSource sourceRef="Geminga" />
       <nestedSource sourceRef="PKS0528" />
   </source>
   @endverbatim
   It should be noted that these particular sources will only produce
   photons with energies greater than 100 MeV.  This will be very
   important for how we define the ROI later on.
   \n\n
   As a convenience, a composite source "all_3EG_sources" is provided in the <a
   href="http://glast.stanford.edu/cgi-bin/cvsweb/observationSim/xml/3EG_catalog_20-1e6MeV.xml?cvsroot=CVS_SLAC">3EG_catalog_20-1e6MeV.xml</a>
   model file.
   \n\n
   Once the desired sources are listed in the target file, we can now
   run @b obsSim to generate some data.  Be sure the <a
   href="http://glast.stanford.edu/cgi-bin/cvsweb/observationSim/data/obsSim.par?cvsroot=CVS_SLAC">obsSim.par</a> file
   is in the directory pointed to by the PFILES environment variable; or
   if PFILES is unset, be sure it is in your current working directory.
   Here's our @b obsSim session:
   @verbatim
   > time obsSim.exe 
   File containing source names [source_names.dat] : 
   Number of events (or simulation time in seconds) <1 - 4e7> [86400] : 
   Use number of events as simulation time? [yes] : 
   Response functions to use <FRONT/BACK|FRONT|BACK> [FRONT/BACK] : 
   Prefix for output files [diffuse_test] : test
   Pointing history file [none] : 
   Generating events for a simulation time of 86400 seconds....
   Done.
   23.600u 0.370s 0:32.75 73.1%    0+0k 0+0io 4561pf+0w
   @endverbatim

   @section makeCountsMaps Make Counts Maps from the Event Files
   @b obsSim creates the following FITS files:
   @verbatim
   > ls -l test*.fits
   -rw-rw-r--    1 jchiang  gl        3231360 Dec  7 13:28 test_events_0000.fits
   -rw-rw-r--    1 jchiang  gl         279360 Dec  7 13:28 test_scData_0000.fits
   @endverbatim
   The first file contains FITS binary tables of the photon events,
   having been generated for each source, then passed through the
   instrument response functions (IRFs); and the second contains the
   spacecraft data --- time stamps, orbit and attitude information, etc.
   \n\n 
   Various tools are available for displaying data in FITS files. 
   One can use ds9 directly to do the binning for you:
   @verbatim
   > ds9 "test_events_0000.fits[1][bin=ra,dec]"
   @endverbatim
   but here we use the FTOOL @b fcopy and the extended filename syntax
   of cfitsio to create a counts map file:
   @verbatim
   > fcopy "test_events_0000.fits[1][bin ra,dec]" all_sky_1day.fits
   @endverbatim
   The <a href="http://glast.gsfc.nasa.gov/ssc/dev/binned_analysis/EventBin_DC1.html">EventBin</a> application, part of the Science Tools
   distribution, provides capabilities for creating counts maps, in addition
   to spectra and light curves that are particularly well-suited for LAT 
   analysis.  In any case, we will use ds9 to view the counts maps:

   @image html allsky.png all_sky_1day.fits
   Here is a corresponding image using @ref all3EG.
   \n\n
   Since we would like to compare these data against a model and will be
   concentrating on the Galactic anti-center, let's create a more serviceable
   FITS file, concentrating on that region:
   @verbatim
   > fcopy "test_events_0000.fits[1][bin ra=65:105,dec=2:42]" anticenter_1day.fits
   @endverbatim

   @section defineROI Define the Region-of-Interest 
   The Region-of-Interest (ROI) comprises a set of selections on photon
   arrival time, energy, and direction.  These selections are made in
   addition to any that are made using the query to the 
   <a href="http://glast.gsfc.nasa.gov/cgi-bin/ssc/U1/D1WebQuery.cgi">
   GSSC database</a>
   or via the 
   <a href="http://glast.gsfc.nasa.gov/ssc/dev/databases/dataSubselector.html">
   user-level selection tool</a>.  Unfortunately, because of an
   oversight in the Likelihood implementation, the ROI
   selections must @em include any that are made by these other
   tools.  More precisely, the data-space defined by the ROI must lie
   entirely within the intersection of the data-spaces defined by these
   other tools.  For example, if the user-level selection tool selects
   events within a time interval (t1, t2), then any time intervals
   defined in the ROI file must lie entirely within (t1, t2).
   In future releases, this oversight will be corrected.
   \n\n
   Here's the ROI file we will use:
   @verbatim
   <?xml version='1.0' standalone='no'?>
   <Region-of-Interest title="Anticenter Region">
      <timeInterval start="0"
                    stop="1"
                    unit="days"/>
      <energies emin="100." 
                emax="3.1623e5"
                unit="MeV"/>
      <acceptanceCone longitude="180."
                      latitude="0."
                      radius="25"
                      coordsys="Galactic"/>
   </Region-of-Interest>
   @endverbatim
   Any number of time intervals can be specified and can be given in
   units of seconds or days, referenced to the start of the DC1
   simulation.  Data will be selected from the intersection of
   intervals given.  In future, the GTI extension of the FITS event
   and spacecraft files will provide this same capability.  We choose
   a minimum energy of 100 MeV.  This is necessary since the sources
   we have used in the simulation do not produce events with energies
   less than 100 MeV. (Note that we really should consider energy
   dispersion in our choice of emin.)  A similar consideration would
   have to be made if the user-level selection tool also applied an
   energy cut.  Lastly, an acceptance cone defines the part of the
   sky from which photons will be accepted.  The "coordsys" attribute
   can either be "Galactic" or "J2000"; the "radius" attribute gives
   the half-opening angle of the acceptance cone; and units are in
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
   \n\n
   We create the exposure map using the @b expMap tool:
   @verbatim
   > time expMap.exe
   ROI cuts file [RoiCuts.xml] : 
   Spacecraft file [test_scData_0000.fits] : 
   Response functions to use <FRONT/BACK|FRONT|BACK> [FRONT/BACK] : 
   Radius of the source region (in degrees) [35] : 
   Number of longitude points <2 - 1000> [70] : 
   Number of latitude points <2 - 1000> [70] : 
   Number of energies <2 - 100> [10] : 
   Exposure file [anticenter_expMap.fits] : 
   Computing the ExposureMap....................!
   464.390u 1.110s 8:25.68 92.0%	0+0k 0+0io 4071pf+0w
   @endverbatim
   Note that we have chosen a 35 degree radius "source region", while
   the ROI acceptance cone radius is 25 degrees.  This is necessary to
   ensure that photons from sources outside the ROI are accounted for,
   at least partially, owing to the size of the instrument
   point-spread function.  Half-degree pixels are a nominal choice;
   smaller pixels should result in a more accurate evaluation of the
   diffuse source fluxes but will also make the exposure map
   calculation itself lengthier.  The number of energies specifies the
   number of logarithmically spaced intervals bounded by the energy
   range given in the ROI cuts xml file.  Here is one image plane of the
   exposure map we just created:

   @image html expMap.png anticenter_expMap.fits

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
   At this point, we hit <return> to instruct the application to fit
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
