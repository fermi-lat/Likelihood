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
   - @ref makeExposureMap This is needed for analyzing diffuse sources.
   - @ref sourceModelFile The source model XML file contains the
          various sources and their model parameters to be fit
          using the @b gtlikelihood tool.
   - @ref computeDiffuseResps Each event must have a separate response
          precomputed for each diffuse component in the source model.
   - @ref runLikelihood
   - @ref makeTsMaps These are used for point source localization and for
          finding weaker sources after the stronger sources have
          been modeled.

   @section makeSubselections Make Subselections from the Event Data.
   For this tutorial, we will use data generated using @b gtobssim for a
   model comprising the Third EGRET (3EG) catalog sources and Galactic
   diffuse and extragalactic diffuse emission for a simulation time of
   one day.  Since we have run @b gtobssim to create the event data, we
   have opted to write events from each component into separate FT1
   files:
   @verbatim
   newtoby[jchiang] ls -lrt *events_0000.fits
   -rw-r--r--    1 jchiang  jchiang   2894400 Mar  7 12:08 eg_diffuse_events_0000.fits
   -rw-r--r--    1 jchiang  jchiang   1088640 Mar  7 12:08 ptsrcs_events_0000.fits
   -rw-r--r--    1 jchiang  jchiang   7243200 Mar  7 12:08 galdiffuse_events_0000.fits
   @endverbatim

   We will consider data in the Virgo region within a 20 degree
   acceptance cone of the blazar 3C 279.  Here we apply @b gtselect to
   the point source data; the same selections need to be applied to
   the files <tt>galdiffuse_events_0000.fits</tt> and 
   <tt>eg_diffuse_events_0000.fits</tt>:
   @verbatim
   newtoby[jchiang] gtselect
   Input FT1 file [ptsrcs_events_0000.fits] : 
   Output FT1 file [ptsrcs_events_filtered.fits] : 
   RA for new search center (degrees) <0 - 360> [193.98] : 
   Dec for new search center (degrees) <-90 - 90> [-5.82] : 
   radius of new search region (degrees) <0 - 180> [25] : 20
   start time (MET in s) [0] : 
   end time (MET in s) [86400] : 
   lower energy limit (MeV) [30] : 
   upper energy limit (MeV) [200000] : 
   first conversion layer <0 - 15> [0] : 
   last conversion layer <0 - 15> [15] : 
   Done.
   newtoby[jchiang] 
   @endverbatim
   There are numerous hidden parameters that correspond to
   FT1 columns that may also be selected upon.  However, as is customary
   with FTOOLS' hidden parameters, the typical user need not worry
   about about them.

   @b gtselect writes descriptions of the data selections to a
   series of "Data Sub-Space" keywords in the <tt>EVENTS</tt>
   extension header.  These keywords are used by the exposure-related
   tools and by @b gtlikelihood for calculating various quantities,
   such as the predicted number of detected events given the source
   model.  These keywords @em must be same for all of the filtered
   event files considered in a given analysis.  @b gtlikelihood will check to
   ensure that all off the DSS keywords are the same in all of the FT1 files.
   See <a
   href="http://glast.gsfc.nasa.gov/ssc/dev/binned_analysis/dss_keywords.html">Yasushi
   Ikebe's DSS keyword page</a> and <a
   href="http://confluence.slac.stanford.edu/display/ST/Data+SubSpace+keywords?showComments=false">this
   confluence page</a> for a discussion of DSS keywords.

   I have written a small non-FTOOL that can be used to view the DSS
   keywords in a given extension, where the <tt>EVENTS</tt> extension is 
   assumed by default:

   @verbatim
   newtoby[jchiang] viewCuts
   usage: viewCuts <filename> [<extname>]
   newtoby[jchiang] viewCuts ptsrcs_events_filtered.fits
   DSTYP1: POS(RA,DEC)
   DSUNI1: deg
   DSVAL1: CIRCLE(193.98,-5.82,20)
   
   DSTYP2: TIME
   DSUNI2: s
   DSVAL2: :86400
   
   DSTYP3: ENERGY
   DSUNI3: MeV
   DSVAL3: 30:200000

   DSTYP4: CONVERSION_LAYER
   DSUNI4: dimensionless
   DSVAL4: :15
   
   DSTYP5: CALIB_VERSION[1]
   DSUNI5: dimensionless
   DSVAL5: 1:1
   
   DSTYP6: CALIB_VERSION[2]
   DSUNI6: dimensionless
   DSVAL6: 1:1
   
   DSTYP7: CALIB_VERSION[3]
   DSUNI7: dimensionless
   DSVAL7: 1:1

   newtoby[jchiang] 
   @endverbatim

   The <tt>CALIB_VERSION</tt> column selections refer to the 
   IM-based event quality flags used in DC1 and are examples of
   hidden parameter selections that the typical user need not be
   concerned with.

   @section makeCountsMaps Make Counts Maps from the Event Files.  

   Next, we create a counts map file to visualize the extraction
   region we wish to fit.  There are presently three tools available
   to perform this task: @b gtbin in the @b evtbin package, @b
   count_map in the @b map_tools package, and @b gtcntsmap in the
   Likelihood package.  Since the event data are contained in three
   separate FITS files and since the first two tools cannot read in
   more than one input event file at a time, @b gtcntsmap is the most
   convenient to use for these data.  We start by creating a file
   which is a list of the files to be binned and then run the tool on
   that file list:

   @verbatim
   newtoby[jchiang] ls -1 *filtered.fits > eventFiles
   newtoby[jchiang] cat eventFiles
   eg_diffuse_events_filtered.fits
   galdiffuse_events_filtered.fits
   ptsrcs_events_filtered.fits
   newtoby[jchiang] gtcntsmap
   Event data file [eventFiles] : eventFiles
   Spacecraft data file [ptsrcs_scData_0000.fits] : ptsrcs_scData_0000.fits
   Output file name [countsMap.fits] : Virgo_map.fits 
   Minimum energy (MeV) <20 - 3e5> [30] : 
   Maximum energy (MeV) <20 - 3e5> [200000] : 
   Number of energy bounds <2 - 40> [2] : 2
   Right Ascension of map center (degrees) <0 - 360> [193.98] : 
   Declination of map center (degrees) <-90 - 90> [-5.82] : 
   Number of longitude pixels [200] : 
   Number of latitude pixels [200] : 
   Pixel size (degrees) <1e-2 - 2> [0.5] : 0.25
   newtoby[jchiang] ds9 Virgo_map.fits
   @endverbatim
   
   For the spacecraft data file, we have entered the one generated by
   @b gtobssim for the 3EG sources.  Identical files were created for
   the diffuse component runs, and we could also have used either of
   those instead.  The last command launches the visualization tool <a
   href="http://hea-www.harvard.edu/RD/ds9/">ds9</a> and produces this
   display:

   @image html Virgo_map.png Virgo_map.fits

   Here I have plotted as the green circles the locations of the 3EG
   point sources within 20 degrees of 3C 279.

   @section makeExposureMap Make an Exposure Map 
   We are now ready to create an exposure map.
   The type of exposure map used by Likelihood differs significantly
   from the usual notion of exposure maps, which are essentially
   integrals of effective area over time.  The exposure calculation
   that Likelihood uses consists of an integral of the total response
   over the entire region-of-interest (ROI) data-space:
   \f[ \epsilon(E, \hat{p}) =
        \int_{\rm ROI} dE^\prime d\hat{p}^\prime dt 
        R(E^\prime, \hat{p}^\prime; E, \hat{p}, t),
   \f]
   where primed quantities indicate measured energies, \f$E^\prime\f$,
   and measured directions, \f$\hat{p}^\prime\f$.
   This exposure function can then be used to compute the expected
   numbers of events from a given source:
   \f[
   N_{\rm pred} = \int dE d\hat{p} S_i(E, \hat{p}) \epsilon(E, \hat{p}),
   \f]
   where \f$S_i(E, \hat{p})\f$ is the photon intensity from source 
   @em i.  Since the exposure calculation involves an integral over the
   ROI, separate exposure maps must be made for every distinct set of
   DSS cuts if, for example, one wants to subdivide an observation to
   look for secular flux variations from a particular source or
   sources.

   There are two tools needed for generating exposure maps.  The first
   is @b gtlivetimecube.  This tool creates a HealPix table, covering
   the full sky, of the integrated livetime as a function of
   inclination with respect to the LAT z-axis. It is otherwise
   identical to the @b exposure_cube application in the @b map_tools
   package, except that @b gtlivetimecube applies the time-range and
   GTI cuts specified in the DSS keywords of the filtered event files.

   @verbatim
   newtoby[jchiang] gtlivetimecube 
   Event data file [ptsrcs_events_0000.fits] : 
   Spacecraft data file [./tmpy3XDg1] : ptsrcs_scData_0000.fits
   Output file [expCubeFile.fits] : 
   Step size in cos(theta) <0. - 1.> [0.05] : 0.025
   Pixel size (degrees) [1] : 
   Working on file ptsrcs_scData_0000.fits
   .....................!
   newtoby[jchiang] 
   @endverbatim
   Since @b gtlivetimecube produces a FITS file covering the entire
   sky map, the output of this tool can be used for generating
   exposure maps for regions-of-interest in other parts of the sky
   that have the same time interval selections.  Although the @b
   gtexpmap application (see below) can generate exposure maps for
   Likelihood without a livetime file, using one affords a substantial
   time savings.
   \n\n
   Creating the exposure map using the @b gtexpmap tool, we have
   @verbatim
   newtoby[jchiang] gtexpmap
   Event data file [ptsrcs_events_filtered.fits] : 
   Spacecraft data file [ptsrcs_scData_0000.fits] : 
   Exposure hypercube file [expCubeFile.fits] : 
   output file name [expMap.fits] : 
   Response functions <DC1|G25|TEST> [DC1] : 
   Radius of the source region (in degrees) [35] : 30
   Number of longitude points <2 - 1000> [140] : 120
   Number of latitude points <2 - 1000> [140] : 120
   Number of energies <2 - 100> [20] : 
   Computing the ExposureMap....................!
   newtoby[jchiang] 
   @endverbatim
   
   Note that we have chosen a 30 degree radius "source region", while
   the acceptance cone radius specified for @b gtselect was 20
   degrees.  This is necessary to ensure that photons from sources
   outside the ROI are accounted for owing to the size of the
   instrument point-spread function.  Half-degree pixels are a nominal
   choice; smaller pixels should result in a more accurate evaluation
   of the diffuse source fluxes but will also make the exposure map
   calculation itself lengthier.  The number of energies specifies the
   number of logarithmically spaced intervals bounded by the energy
   range given in the DSS keywords.  Here is one image plane of the
   exposure map we just created:

   @image html expMap.png expMap.fits

   @section sourceModelFile Create a Source Model XML File.

   The @b gtlikelihood tool reads the source model from an XML file.
   Given the dearth of bright sources in the extraction region we have
   selected, our source model file will be fairly simple, comprising
   only the Galactic and extragalactic diffuse emission, and point
   sources to represent the blazars 3C 279 and 3C 273:
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
   \n\n

   - Two types of sources can be defined, PointSource and DiffuseSource.

   - Each type comprises a spectrum and a spatialModel component.

   - Model parameters are described by a set of attributes. The
     actual value of a given parameter that is used in the
     calculation is the value attribute multiplied by the scale
     attribute. The value attribute is what the optimizers see. Using
     the scale attribute is necessary to ensure that the parameters
     describing the objective function, -log(likelihood) for this
     application, all have values lying roughly within an
     order-of-magnitude of each other.

   - Each parameter has a range of valid values that can be specified.

   - The free attribute determines whether the parameter will be
     allowed to be fixed or free in the fitting process. Free flag
     attributes are currently disabled for spatial model parameters
     since fitting for these parameters has not been implemented
     (primarily due to the enormous overhead associated with
     computing the energy-dependent response functions for each
     source component).

   - The units for the spectral models are 
     \f${\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}\f$
     for point sources and 
     \f${\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}\,{\rm sr}^{-1}\f$
     for diffuse sources

   - Several spectral functions are available:

   - @b PowerLaw This function has the form
   \f[
   \frac{dN}{dE} = N_0 \left(\frac{E}{E_0}\right)^\gamma
   \f]
   where the parameters in the XML definition have the following
   mappings:
     - Prefactor = \f$N_0\f$
     - Index = \f$\gamma\f$
     - Scale = \f$E_0\f$

   - @b BrokenPowerLaw
   \f[
   \frac{dN}{dE} = N_0 \times\left\{\begin{array}{ll}
                                    (E/E_b)^{\gamma_1} & \mbox{if $E < E_b$}\\
                                    (E/E_b)^{\gamma_2} & \mbox{otherwise}
                                    \end{array} \right.
   \f]
   where
     - Prefactor = \f$N_0\f$
     - Index1 = \f$\gamma_1\f$
     - Index2 = \f$\gamma_2\f$
     - BreakValue = \f$E_b\f$

   - @b PowerLaw2 This function uses the integrated flux as a free
   parameter rather than the Prefactor:
   \f[
   \frac{dN}{dE} = \frac{N(\gamma+1)E^{\gamma}}
                        {E_{\rm max}^{\gamma+1} - E_{\rm min}^{\gamma+1}}
   \f]
   where 
   <ul>
     <li> Integral = \f$N\f$
     <li> Index = \f$\gamma\f$
     <li> LowerLimit = \f$E_{\rm min}\f$
     <li> UpperLimit = \f$E_{\rm max}\f$
   </ul><br>
   The UpperLimit and LowerLimit parameters are always treated as
   fixed, and as should be apparent from this definition, the flux
   given by the Integral parameter is over the range (LowerLimit,
   UpperLimit). Use of this model allows the errors on the integrated
   flux to be evaluated directly by likelihood, obviating the need to
   propagate the errors if one is using the PowerLaw form.

   - @b BrokenPowerLaw2 Similar to PowerLaw2, the integral flux
   is the free parameter rather than the Prefactor:
   \f[
   \frac{dN}{dE} = N_0(N, E_{\rm min}, E_{\rm max}, \gamma_1, \gamma_2)
                   \times\left\{\begin{array}{ll}
                                (E/E_b)^{\gamma_1} & \mbox{if $E < E_b$}\\
                                (E/E_b)^{\gamma_2} & \mbox{otherwise}
                                \end{array} \right.
   \f]
   where
   \f[
   \newcommand{\emin}{{E_{\rm min}}}
   \newcommand{\emax}{{E_{\rm max}}}
   \newcommand{\pfrac}[2]{\left(\frac{#1}{#2}\right)}
   \newcommand{\Int}{{\displaystyle\int}}
   N_0(N, E_{\rm min}, E_{\rm max}, \gamma_1, \gamma_2)
        = N \times \left\{\begin{array}{ll}
                   \left[ \Int_\emin^\emax \pfrac{E}{E_b}^{\gamma_1} dE\right]^{-1}
                   & \mbox{$\emax < E_b$}\\
                   \left[ \Int_\emin^\emax \pfrac{E}{E_b}^{\gamma_2} dE\right]^{-1}
                   & \mbox{$\emin > E_b$}\\
                   \left[ \Int_\emin^{E_b} \pfrac{E}{E_b}^{\gamma_1} dE
                   + \Int_{E_b}^\emax \pfrac{E}{E_b}^{\gamma_2} dE
                   \right]^{-1}
                   & \mbox{otherwise}\end{array}\right.
   \f]
   and 
     - Integral = \f$N\f$
     - Index1 = \f$\gamma_1\f$
     - Index2 = \f$\gamma_2\f$
     - BreakValue = \f$E_b\f$
     - LowerLimit = \f$E_{\rm min}\f$
     - UpperLimit = \f$E_{\rm max}\f$

   - @b LogParabola This is typically used for modeling Blazar spectra.
   \f[
   \newcommand{\pfrac}[2]{\left(\frac{#1}{#2}\right)}
   \frac{dN}{dE} = N_0\pfrac{E}{E_b}^{-(\alpha + \beta\log(E/E_b))}
   \f]
   where
     - norm = \f$N_0\f$
     - alpha = \f$\alpha\f$
     - beta = \f$\beta\f$
     - Eb = \f$E_b\f$

   - @b ExpCutoff An exponentially cut-off power-law used for modeling
     blazar spectra subject to absorption by the extragalactic
     background light (EBL).  This model was implemented by Luis Reyes
     (<tt>lreyes@milkyway.gsfc.nasa.gov</tt>).
   \f[
   \newcommand{\pfrac}[2]{\left(\frac{#1}{#2}\right)}
   \frac{dN}{dE} = N_0 \times \left\{\begin{array}{ll}
                              \pfrac{E}{E_0}^\gamma & \mbox{$E < E_b$}\\
                              \pfrac{E}{E_0}^\gamma
         \exp\left[ - ( (E - E_b)/p_1 + p_2\log(E/E_b) + p_3\log^2(E/E_b) )
             \right] & \mbox{otherwise} \end{array}\right.
   \f]
   where
     - Prefactor = \f$N_0\f$
     - Index = \f$\gamma\f$
     - Scale = \f$E_0\f$
     - Ebreak = \f$E_b\f$
     - P1 = \f$p_1\f$
     - P2 = \f$p_2\f$
     - P3 = \f$p_3\f$

   - @b BPLExpCutoff An exponentially cut-off broken power-law,
     implemented by Jennifer Carson (<tt>carson@slac.stanford.edu</tt>).
   \f[
   \newcommand{\pfrac}[2]{\left(\frac{#1}{#2}\right)}
   \newcommand{\eabs}{{E_{\rm abs}}}
   \frac{dN}{dE} = N_0 \times \left\{\begin{array}{ll}
          \pfrac{E}{E_b}^{\gamma_1} & \mbox{$E < E_b$ and $E < \eabs$}\\
          \pfrac{E}{E_b}^{\gamma_2} & \mbox{$E > E_b$ and $E < \eabs$}\\
          \pfrac{E}{E_b}^{\gamma_1}\exp[-(E - \eabs)/p_1] 
                          & \mbox{$E < E_b$ and $E > \eabs$}\\
          \pfrac{E}{E_b}^{\gamma_2}\exp[-(E - \eabs)/p_1] 
                          & \mbox{$E > E_b$ and $E > \eabs$}\end{array}\right.
   \f]
   where
     - Prefactor = \f$N_0\f$
     - Index1 = = \f$\gamma_1\f$
     - Index2 = \f$\gamma_2\f$
     - BreakValue = \f$E_b\f$
     - Eabs = \f$E_{\rm abs}\f$
     - P1 = \f$p_1\f$

   - @b Gaussian A Gaussian function that can be used to model an emission
     line.
   \f[
   \frac{dN}{dE} = \frac{N_0}{\sigma\sqrt{2\pi}} 
                   \exp\left[\frac{-( E - \bar{E} )^2}{2\sigma^2}\right]
   \f]
   where
     - Prefactor = \f$N_0\f$
     - Mean = \f$\bar{E}\f$
     - Sigma = \f$\sigma\f$

   - @b ConstantValue A constant-valued function, independent of energy.
   \f[
   \frac{dN}{dE} = N_0
   \f]
   where 
     - Value = \f$N_0\f$

   - @b FileFunction A function defined using an input ASCII file 
   with columns of energy and differential flux values.  The energy 
   units are assumed to be MeV and the flux values are assumed to 
   \f${\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}\f$
   for a point source and 
   \f${\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}\,{\rm sr}^{-1}\f$
   for a diffuse source.  The sole parameter is a multiplicative
   normalization.
   \f[
   \frac{dN}{dE} = N_0\left.\frac{dN}{dE}\right|_{\rm file}
   \f]
   where
     - Normalization = \f$N_0\f$

   - @b BandFunction This function is used to model GRB spectra
   \f[
   \newcommand{\pfrac}[2]{{(#1/#2)}}
   \frac{dN}{dE} = N_0 \times \left\{\begin{array}{ll}
      \pfrac{E}{E_0}^\alpha \exp\left[-(E/E_p)\right] 
         & \mbox{if $E < E_p(\alpha - \beta)$}\\
      \pfrac{E}{E_0}^\beta\left[\pfrac{E_p}{E_0} 
      (\alpha - \beta)\right]^{\alpha - \beta} 
      \exp(\beta - \alpha) & \mbox{otherwise}
   \end{array} \right.
   \f]
   where
     - norm = \f$N_0\f$
     - alpha = \f$\alpha\f$
     - beta = \f$\beta\f$
     - Ep = \f$E_p\f$
     - Scale = \f$E_0\f$

   - @b PLSuperExpCutoff For modeling pulsars, implemented by Damien Parent 
   (<tt>parent@cenbg.in2p3.fr</tt>)
   \f[
   \newcommand{\pfrac}[2]{\left(\frac{#1}{#2}\right)}
   \frac{dN}{dE} = N_0 \pfrac{E}{E_0}^{\gamma_1}
                   \exp\left(-\pfrac{E}{E_c}^{\gamma_2}\right)
   \f]
   where
     - Prefactor = \f$N_0\f$
     - Index1 = \f$\gamma_1\f$
     - Scale = \f$E_0\f$
     - Cutoff = \f$E_c\f$
     - Index2 = \f$\gamma_2\f$

   - @b SmoothBrokenPowerLaw  Implemented by Benoit Lott 
   (<tt>lott@cenbg.in2p3.fr</tt>)
   \f[
   \frac{dN}{dE} = N_0\left(\frac{E}{E_0}\right)^{\alpha_1} 
                   \displaystyle \left(1+\left(\frac{E}{E_b}\right)
                   ^\frac{\alpha_1-\alpha_2}{\beta}\right)^{-\beta}
   \f]
   where
     - Prefactor = \f$N_0\f$
     - Index1 = \f$\alpha_1\f$
     - Scale = \f$E_0\f$
     - Index2 = \f$\alpha_2\f$
     - BreakValue = \f$E_b\f$
     - Beta = \f$\beta\f$

   - Four spatial models are available:
     - SkyDirFunction describes a direction on the sky and is used
       only for point sources.
     - ConstantValue provides a constant value regardless of what
       argument value it takes. In the current context, ConstantValue
       is used to model the isotropic diffuse emission. As a function,
       however, ConstantValue is fairly general and can even be used
       in a spectral model; as it is when the spatial model is a
       MapCubeFunction.
     - SpatialMap uses a FITS image file as a template for determining
       the distribution of photons on the sky. The EGRET diffuse model
       is given in the FITS file gas.cel, which describes the
       interstellar emission.
     - MapCubeFunction used for diffuse sources that are modeled by a
       3 dimensional FITS map (two sky coordinates and energy),
       thereby allowing arbitrary spectral variation as a function of
       sky position.

   @section computeDiffuseResps Compute the Diffuse Source Responses.
   If these quantities are not precomputed using the @b gtdiffresp
   tool, then @b gtlikelihood will compute them at runtime.  However,
   if multiple fits and/or sessions with the @b gtlikelihood tool are
   anticipated, it is probably wise to precompute these quantities.
   The source model XML file must contain all of the diffuse sources
   to be fit.  The @b gtdiffresp tool will add columns to the FT1 file
   for each diffuse source:

   @verbatim
   newtoby[jchiang] gtdiffresp
   Event data file [diffuse_events_0000.fits] : ptsrcs_events_filtered.fits
   Spacecraft data file [diffuse_scData_0000.fits] : ptsrcs_scData_0000.fits
   Source model file [diffuse_model.xml] : Virgo_model.xml
   Response functions to use <DC1|G25|TEST> [DC1] : 
   adding source 3C 273
   adding source 3C 279
   adding source Extragalactic Diffuse
   adding source Galactic Diffuse
   ......................!
   @endverbatim
   Since 3C 273 and 3C 279 are point sources, diffuse responses are not
   calculated for these components.

   @section runLikelihood Run gtlikelihood.

   We are now ready to run the @b gtlikelihood application:
   @verbatim
   newtoby[jchiang] gtlikelihood
   Statistic to use <BINNED|UNBINNED|OPTEM> [UNBINNED] : 
   Spacecraft file [oneday_scData_0000.fits] : ptsrcs_scData_0000.fits
   Event file [test_events_0000.fits] : eventFiles
   Exposure file [none] : expMap.fits
   Source model file [my_source_model.xml] : Virgo_model.xml
   Source model output file [none] : Virgo_model.xml
   flux-style output file name [flux_model.xml] : 
   Response functions to use <DC1|G25|TEST> [DC1] : 
   Optimizer <LBFGS|MINUIT|DRMNGB> [MINUIT] : 
   Allow for refitting? [no] : yes
   @endverbatim
   Most of the entries prompted for are fairly obvious.  In addition
   to the various XML and FITS files, the user is prompted for a
   choice of IRFs, the type of statistic to use, the optimizer, 
   and some output file names.  

   Three response function specifications are available:

   - @b DC1 The response functions created for DC1 analysis.
   - @b G25 AO era IRFs for the FRONT and BACK of the LAT
   - @b TEST A set of idealized response functions based on an AllGamma
        data set generated post-DC1 using GlastRelease v4r2.

   @b gtobssim has the same option of generating events for these IRFs
   and the choice used for @b gtlikelihood must be the same as that
   chosen for @b gtobssim.

   Three statistics are available:

   - @b UNBINNED This is a standard unbinned analysis, described in this
     tutorial.  If this option is chosen then parameters for the
     spacecraft file, event file, and exposure file must be
     given.

   - @b BINNED This is a binned analysis, which is still in development.
     In this case, the Counts map file and Exposure hypercube file
     must be given.  The spacecraft file, event file and
     exposure file are ignored in this case.

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
   The "flux-style output file name" specifies the destination of XML
   definitions of sources that can be used with @b gtobssim (or @b Gleam)
   to simulate an observation of the fitted source model.
   \n\n
   The application proceeds by reading in the spacecraft and event
   data, and if necessary, computing event responses for each diffuse
   source.
   \n\n
   Here is the output from our fit:
   @verbatim
   
    **********
    **    1 **SET PRINT    0.000    
    **********
    **********
    **    2 **SET NOWARN 
    **********
   
    PARAMETER DEFINITIONS:
       NO.   NAME         VALUE      STEP SIZE      LIMITS
        1 'Prefactor '    10.000       1.0000        0.10000E-02   1000.0    
        2 'Index     '   -2.1000       1.0000        -5.0000      -1.0000    
        3 'Prefactor '    10.000       1.0000        0.10000E-02   1000.0    
        4 'Index     '   -2.0000       1.0000        -5.0000      -1.0000    
        5 'Prefactor '    1.6000       1.0000        0.10000E-04   100.00    
        6 'Prefactor '    11.000       1.0000        0.10000E-02   1000.0    
    **********
    **    3 **SET ERR   0.5000    
    **********
    **********
    **    4 **SET GRAD    1.000    
    **********
    **********
    **    5 **MIGRAD    200.0       3184.    
    **********
   
    MIGRAD MINIMIZATION HAS CONVERGED.
   
    MIGRAD WILL VERIFY CONVERGENCE AND ERROR MATRIX.
   
    FCN=   15883.48     FROM MIGRAD    STATUS=CONVERGED     87 CALLS       88 TOTAL
                        EDM=  0.17E+01    STRATEGY= 1      ERR MATRIX NOT POS-DEF
     EXT PARAMETER                APPROXIMATE        STEP         FIRST   
     NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE 
      1  Prefactor     4.9455        18.081       0.10063E-01   -4.3787    
      2  Index        -2.8631        3.0635       0.38323       0.30560    
      3  Prefactor     10.964        4.3290       0.10063E-01   -5.8494    
      4  Index        -2.1936       0.23045       0.11464        1.8059    
      5  Prefactor     1.6140       0.17338       0.14832E-01    5.7092    
      6  Prefactor     9.7163        1.0606       0.95978E-02    31.706    
                                  ERR DEF= 0.500    
    **********
    **    6 **HESSE 
    **********
   
    FCN=   15883.48     FROM HESSE     STATUS=OK            42 CALLS      130 TOTAL
                     EDM=  0.40E+00    STRATEGY= 1      ERROR MATRIX ACCURATE 
   
     EXT PARAMETER                                INTERNAL      INTERNAL  
     NO.   NAME        VALUE          ERROR       STEP SIZE       VALUE   
      1  Prefactor     4.9455        1.4755       0.22298E-03   -1.4300    
      2  Index        -2.8631       0.30475       0.15183E-01   0.68517E-01
      3  Prefactor     10.964        1.9151       0.76098E-03   -1.3610    
      4  Index        -2.1936       0.10962       0.50432E-03   0.41504    
      5  Prefactor     1.6140       0.18767       0.65564E-04   -1.3160    
      6  Prefactor     9.7163        2.1563       0.99601E-04   -1.3733    
                                  ERR DEF= 0.500    
   Final values: 
     Prefactor  = 4.94554
     Index      = -2.86307
     Prefactor  = 10.9642
     Index      = -2.19355
     Prefactor  = 1.61396
     Prefactor  = 9.71634
   Minuit fit quality: 3   estimated distance: 0.40036
   Minuit parameter uncertainties:
     1  1.47564
     2  0.305948
     3  1.91523
     4  0.109684
     5  0.187678
     6  2.15644
   Computing TS values for each source (4 total)
   ....!
   
   3C 273:
   Prefactor: 4.94554 +/- 1.47564
   Index: -2.86307 +/- 0.305948
   Scale: 100
   Npred: 46.9157
   ROI distance: 10.4265
   TS value: 29.3858
   
   3C 279:
   Prefactor: 10.9642 +/- 1.91523
   Index: -2.19355 +/- 0.109684
   Scale: 100
   Npred: 136.453
   ROI distance: 0
   TS value: 239.298
   
   Extragalactic Diffuse:
   Prefactor: 1.61396 +/- 0.187678
   Index: -2.1
   Scale: 100
   Npred: 917.184
   
   Galactic Diffuse:
   Prefactor: 9.71634 +/- 2.15644
   Index: -2.1
   Scale: 100
   Npred: 449.394
   
   -log(Likelihood): 15883.5

   Writing fitted model to fitted_model.xml
   Refit? [y] 
   @endverbatim
   At this point we can choose to edit the source model file, setting
   various spectral parameters free or fixed, and in the latter case,
   also setting them to some nominal value.  For our refit, we fix the
   diffuse component prefactors to their nominal values as follows:
   @verbatim
  <source name="Extragalactic Diffuse" type="DiffuseSource">
    <spectrum type="PowerLaw">
      <parameter free="0" max="100.0" min="1e-05" name="Prefactor" scale="1e-07" value="1.6"/>
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
        1 'Prefactor '    4.9455       1.0000        0.10000E-02   1000.0    
        2 'Index     '   -2.8631       1.0000        -5.0000      -1.0000    
        3 'Prefactor '    10.964       1.0000        0.10000E-02   1000.0    
        4 'Index     '   -2.1936       1.0000        -5.0000      -1.0000    
    **********
    **    3 **SET ERR   0.5000    
    **********
    **********
    **    4 **SET GRAD    1.000    
    **********
    **********
    **    5 **MIGRAD    200.0       3177.    
    **********
   
    MIGRAD MINIMIZATION HAS CONVERGED.
   
    MIGRAD WILL VERIFY CONVERGENCE AND ERROR MATRIX.
   
    FCN=   15884.86     FROM MIGRAD    STATUS=CONVERGED     42 CALLS       43 TOTAL
                        EDM=  0.10E+01    STRATEGY= 1      ERR MATRIX NOT POS-DEF
   
     EXT PARAMETER                APPROXIMATE        STEP         FIRST   
     NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE 
      1  Prefactor     4.9451        18.473       0.14331E-01    6.2063    
      2  Index        -2.8636        3.0558       0.38148       0.18429    
      3  Prefactor     10.964        5.2255       0.96134E-02    6.7041    
      4  Index        -2.2012       0.26830       0.11581      -0.84071E-01
                                  ERR DEF= 0.500    
    **********
    **    6 **HESSE 
    **********
   
    FCN=   15884.86     FROM HESSE     STATUS=OK            25 CALLS       68 TOTAL
                        EDM=  0.24E-01    STRATEGY= 1      ERROR MATRIX ACCURATE 
   
     EXT PARAMETER                                INTERNAL      INTERNAL  
     NO.   NAME        VALUE          ERROR       STEP SIZE       VALUE   
      1  Prefactor     4.9451        1.4075       0.22285E-03   -1.4301    
      2  Index        -2.8636       0.30020       0.15135E-01   0.68237E-01
      3  Prefactor     10.964        1.8496       0.76435E-03   -1.3610    
      4  Index        -2.2012       0.10865       0.50853E-03   0.41086    
                                  ERR DEF= 0.500    
   Final values: 
     Prefactor  = 4.94512
     Index      = -2.86363
     Prefactor  = 10.9636
     Index      = -2.2012
   Minuit fit quality: 3   estimated distance: 0.0238738
   Minuit parameter uncertainties:
     1  1.40756
     2  0.301342
     3  1.84966
     4  0.108712
   Computing TS values for each source (4 total)
   ....!
   
   3C 273:
   Prefactor: 4.94512 +/- 1.40756
   Index: -2.86363 +/- 0.301342
   Scale: 100
   Npred: 46.9096
   ROI distance: 10.4265
   TS value: 27.1845
   
   3C 279:
   Prefactor: 10.9636 +/- 1.84966
   Index: -2.2012 +/- 0.108712
   Scale: 100
   Npred: 135.763
   ROI distance: 0
   TS value: 231.296
   
   Extragalactic Diffuse:
   Prefactor: 1.6
   Index: -2.1
   Scale: 100
   Npred: 909.249
   
   Galactic Diffuse:
   Prefactor: 11
   Index: -2.1
   Scale: 100
   Npred: 508.765
   
   -log(Likelihood): 15884.9
   
   Writing fitted model to Virgo_model.xml
   Refit? [y] 
   n
   Writing flux-style xml model file to flux_model.xml
   newtoby[jchiang] 
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

   Let's comment out 3C 273 from the source model XML file and see if we
   can find evidence for it in the data.
   @verbatim
   newtoby[jchiang] gttsmap
   Event data file [eventFiles] : 
   Spacecraft data file [ptsrcs_scData_0000.fits] : 
   Exposure map file [expMap.fits] : 
   Exposure hypercube file [expCubeFile.fits] : 
   Source model file [Virgo_model.xml] : 
   TS map file name [TsMap.fits] : 
   Response functions to use <DC1|G25|TEST> [DC1] : 
   Optimizer <LBFGS|MINUIT|DRMNGB> [DRMNGB] : 
   Fit tolerance [0.001] : 
   RA minmum <-360 - 360> [180] : 
   RA maximum <-360 - 360> [200] : 
   Number of RA points <2 - 200> [40] : 
   Dec minimum <-90 - 90> [-10] : 
   Dec maximum <-90 - 90> [10] : 
   Number of Dec points <2 - 200> [40] : 
   Use Galactic coordinates [no] : 
   .....................!
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
