// Mainpage for doxygen

/** @mainpage package Likelihood

 @authors James Chiang, Pat Nolan, Karl Young, Toby Burnett, and others

 @section intro Introduction

 This package implements an extended maximum likelihood (EML)
 calculation for analyzing LAT event data.

 These data are generally assumed to be in a format consistent with
 that produced by the Level 1 Event Data Extractor, otherwise known as
 U1.  The data may alternatively have been generated by the
 observation simulator (O2).  Accordingly, use of this tool for
 analysis of Event data requires access to a complete set of
 accompanying spacecraft orbit and attitude information, obtained
 using the Pointing, Livetime History Extractor (U3) or the orbit
 simulator tool (O1), as well as access to appropriate instrument
 response function data (i.e., CALDB).

 However, the classes and methods used here are intended to be
 sufficiently general so that any properly implemented objective
 function should be able to be analyzed with this package, whether it
 is LAT-specific or not.

 @section LatStatModel The Unbinned log-Likelihood

 For LAT event analysis, the default statistical model we assume is
 the unbinned log-likelihood:

 \f[
 \log L = \sum_j \left[\log \left(\sum_i M_i(x_j; \vec{\alpha_i})\right)\right] 
        - \sum_i \left[\int dx M_i(x; \vec{\alpha_i})\right]
 \f]

 Here \f$x_j\f$ is the \f$j\f$th photon Event, as specified by
 apparent energy, direction, and arrival time. The function \f$M_i(x;
 \vec{\alpha_i})\f$ returns the flux density in units of counts per
 energy-time-area-solid angle (i.e., photon fluxes convolved through
 the instrument response) for the \f$i\f$th Source at a point \f$x\f$
 in the Event configuration space, hereafter known as the "data
 space".  Each \f$M_i\f$ is defined, in part, by a vector of parameter
 values \f$\vec{\alpha_i}\f$; collectively, the \f$\vec{\alpha_i}\f$
 vectors form the space over which the objective function is to be
 optimized.  The integral over the data space in the second term is
 the predicted number of Events expected to be seen from Source
 \f$i\f$.

 @section classes Important Classes

 Cast in this form, the problem lends itself to being described by the
 following classes.  Some of these classes now reside in the optimizers
 and latResponse packages.

   - optimizers::Function Objects of this class act as "function
   objects" in that the function call operator () is overloaded so
   that instances behave like ordinary C functions.  Several methods
   are also provided for accessing the model Parameters and
   derivatives with respect to those Parameters, either singly or in
   groups.  The behavior of this class is greatly facilitated by the
   Parameter and Arg classes.

   - optimizers::Parameter This is essentially an n-tuple containing
   model parameter information (and access methods) comprising the
   parameter value, scale factor, name, upper and lower bounds and
   whether the parameter is to be considered free or fixed in the
   fitting process.

   - optimizers::Arg This class wraps arguments to Function objects so
   that Function's derivative passing mechanisms can be inherited by
   subclasses regardless of the actual type of the underlying
   argument.  For example, in the log-likelihood, we define
   Likelihood::logSrcModel, a Function subclass that returns the
   quantity inside the square brackets of the first term on the rhs.
   Acting as a function, logSrcModel naturally wants to take an Event
   as its argument, so we wrap an Event object with EventArg.
   Similarly, the quantity in the square brackets of the second term
   on the rhs we implement as the Likelihood::Npred class.  This class
   wants to have a Source object as its argument, which we wrap with
   the SrcArg class.

   - Likelihood::Source An abstract base class for gamma-ray sources.
   It specifies four key methods (as pure virtual functions); the
   latter two methods are wrapped by the Npred class in order to give
   them Function behavior:
      - fluxDensity(...): counts per energy-time-area-solid angle
      - fluxDensityDeriv(...): derivative wrt a Parameter
      - Npred(): predicted number of photons in the ROI
      - NpredDeriv(...): derivative of Npred wrt a Parameter

   - Likelihood::Event An n-tuple containing photon event arrival
   time, apparent energy and direction, as well as spacecraft attitude
   information at the event arrival time and event-specific response
   function data for use with components of the diffuse emission
   model.

   - latResponse::[IAeff, IPsf, IEdisp] These classes provide abstract
   interfaces to the instrument response functions comprising the
   effective area, the point-spread function, and the energy
   dispersion.  Concrete implementations exist for the Glast25
   parameterizations as well as for the EGRET response functions.

   - Likelihood::RoiCuts An n-tuple Singleton class that contains the
   "region-of-interest" cuts.  These are essentially the bounds of the
   data space as a function of arrival time, apparent energy, apparent
   direction, zenith angle, etc..

   - Likelihood::ScData A Singleton object that contains the
   spacecraft data n-tuples (ScNtuple).

   - Likelihood::SourceFactory This class provides a common access
   point for retrieving and storing Sources. Sources that have been
   constructed to comprise position and spectral information can be
   stored here, then cloned for later use.  The sources can also be
   read in from and output to an xml file.

   - optimizers::Optimizer An abstract base class for the algorithms
   which maximize the desired objective functions.  Choice of
   optimization methods are encapsulated in three sub-classes which
   simply wrap existing implementations that are/were originally
   available as Fortran code: Lbfgs, Minuit, and Drmngb.  

 <hr>
 @section notes Release Notes
  release.notes

 <hr>
 @section requirements requirements
 @verbinclude requirements

 <hr> 
 @todo Energy dispersion
 @todo Generalize Npred calculation, e.g., zenith angle cuts, fit-able 
       source locations
 @todo Refactor Statistic and FITS-related classes (Table, FitsImage, etc.)
 @todo Use more realistic response function data
 @todo Analyze EGRET data
 @todo Boost.Python
 */

/**
 @page userGuide User's Guide

 @section likeApp likelihood application

 This is an FTOOLS-like interface to the Likelihood class library.
 It uses HOOPS to obtain the command-line parameters and therefore
 requires a .par file called "likelihood.par":

@verbinclude likelihood.par

 Here is a sample session in which 1000 events of simulated data,
 produced for two 3EG sources, 3C 279 and 3C 273, are analyzed:

 @verbatim
glast-guess1[jchiang] ~/SciTools/dev/Likelihood/v2r3/rh72_gcc2953/like_test.exe
ROI cuts file [$(LIKELIHOODROOT)/xml/RoiCuts.xml] : 
Spacecraft file [virgo_scFiles] : 
Exposure file [none] : 
Response functions to use <FRONT/BACK|COMBINED> [FRONT/BACK] : 
Source model file [virgo_region.xml] : 
Event file [virgo_region_events_0000.fits] : 
Optimizer <LBFGS|MINUIT|DRMNGB> [DRMNGB] : 
Optimizer verbosity [1] : 0
Fit tolerance [0.0001] : 
Source model output file [/dev/null] : 
Use OptEM? [no] : 
flux-style output file name (can be 'none') [flux_model.xml] : 
Creating source named 3C 279
Computing exposure at (193.98, -5.82).....................!
Creating source named 3C 273
Computing exposure at (187.25, 2.17).....................!
adding source 3C 273
adding source 3C 279
LogLike::getEvents:
Out of 1000 events in file virgo_region_events_0000.fits,
 854 were accepted, and 146 were rejected.

Drmngb return code: 3
3C 273:
Prefactor: 2.52004 +/- 0.186379
Index: -2.56014 +/- 0.0871032
Scale: 100
Npred: 227.436

3C 279:
Prefactor: 6.46561 +/- 0.298801
Index: -1.90003 +/- 0.0296193
Scale: 100
Npred: 626.56
Writing fitted model to /dev/null
Writing flux-style xml model file to flux_model.xml
glast-guess1[jchiang] 
@endverbatim

 - The <a
   href="http://glast.stanford.edu/cgi-bin/cvsweb/Likelihood/xml/RoiCuts.xml?cvsroot=CVS_SLAC">ROI
   cuts file</a> contains basic information about the extraction
   region, valid time ranges, energies, etc.:

@verbatim
<?xml version='1.0' standalone='no'?>
<!DOCTYPE Region-of-Interest SYSTEM "RoiCuts.dtd">
<Region-of-Interest title="Virgo Region">
   <timeInterval start="0"
                 stop="7"
                 unit="days"/>
   <energies emin="30." 
             emax="3.1623e5"
             unit="MeV"/>
   <acceptanceCone longitude="193.98"
                   latitude="-5.82"
                   radius="30"
                   coordsys="Celestial"/>
</Region-of-Interest>
@endverbatim

 - The Spacecraft file is either a single FITS file or a list of FITS
   files:

@verbatim
virgo_region_scData_0000.fits
virgo_region_scData_0001.fits
@endverbatim

 - If there are diffuse components in the source model, then one must
   provide an exposure map that has been computed using the @ref
   expMap.  If the data contain only point sources, then an exposure
   map can be omitted.

 - Presently, only Front/Back and Combined response functions for the
   Glast25 parameterizations are available.

 - The Source model file contains an xml description of the various
   sources to be modeled.  Here's the source model file used in the
   above fit:

@verbatim
<source_library title="prototype sources" function_library="$(LIKELIHOODROOT)/xml/A1_Functions.xml">
  <source name="3C 279" type="PointSource">
    <spectrum type="PowerLaw">
      <parameter max="1000" min="0.001" free="1" name="Prefactor" scale="1e-09" value="1" />
      <parameter max="-1" min="-3.5" free="1" name="Index" scale="1" value="-2" />
      <parameter max="200" min="50" free="0" name="Scale" scale="1" value="100" />
    </spectrum>
    <spatialModel type="SkyDirFunction">
      <parameter max="3.40282e+38" min="-3.40282e+38" free="0" name="RA" scale="1" value="193.98" />
      <parameter max="3.40282e+38" min="-3.40282e+38" free="0" name="DEC" scale="1" value="-5.82" />
    </spatialModel>
  </source>
  <source name="3C 273" type="PointSource">
    <spectrum type="PowerLaw">
      <parameter max="1000" min="0.001" free="1" name="Prefactor" scale="1e-09" value="1" />
      <parameter max="-1" min="-3.5" free="1" name="Index" scale="1" value="-2" />
      <parameter max="200" min="50" free="0" name="Scale" scale="1" value="100" />
    </spectrum>
    <spatialModel type="SkyDirFunction">
      <parameter max="3.40282e+38" min="-3.40282e+38" free="0" name="RA" scale="1" value="187.25" />
      <parameter max="3.40282e+38" min="-3.40282e+38" free="0" name="DEC" scale="1" value="2.17" />
    </spatialModel>
  </source>
</source_library>
@endverbatim

 - The Event file also can either be a single FITS file or a list.

 - The optimizer package provides three optimizers from which to
   choose.  Minuit and Drmngb provide error estimates based on
   covariance matrix computed from the inverse Hessian.

 - After fitting the source model can be written to an output file.
   If the same file name as the input source model file is given, then
   the file will be overwritten.

 - OptEM is an alternative method for performing the optimization
   using the expectation maximization method.

 - Finally, the source model can be written out as a flux-package
   style xml file that's suitable for use with either obsSim or Gleam.
   Here's the resulting file for the above fit:

@verbatim
<source_library title="Likelihood_model">
  <source flux="0.0105858" name="_3C_273">
    <spectrum escale="MeV">
      <particle name="gamma">
        <power_law emax="316230" emin="30" gamma="2.56014" />
      </particle>
      <celestial_dir ra="187.25" dec="2.17" />
    </spectrum>
  </source>
  <source flux="0.0212465" name="_3C_279">
    <spectrum escale="MeV">
      <particle name="gamma">
        <power_law emax="316230" emin="30" gamma="1.90003" />
      </particle>
      <celestial_dir ra="193.98" dec="-5.82" />
    </spectrum>
  </source>
  <source name="all_in_flux_model.xml">
    <nestedSource sourceRef="_3C_273" />
    <nestedSource sourceRef="_3C_279" />
  </source>
</source_library>
@endverbatim

 @section TsMap TsMap application

 @section expMap expMap application

*/


