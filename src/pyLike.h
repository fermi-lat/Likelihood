/**
   @page pyLike Likelihood and Python

   @section intro Introduction

   Using <a href="http://www.swig.org/">SWIG</a> via the <a
   href="http://www.slac.stanford.edu/~jchiang/SwigPolicy/doxy-html/main.html">SwigPolicy</a>
   package, the C++ class library of Likelihood has been exposed to <a
   href="http://www.python.org/">Python</a>.  These classes by
   themselves do not constitute a suitable interactive interface.
   However, armed with these classes and with a modest effort, 
   two small python modules, <a href="http://glast.stanford.edu/cgi-bin/cvsweb/Likelihood/python/SrcAnalysis.py?hideattic=1&cvsroot=CVS_SLAC">SrcAnalysis.py</a>
   and
   <a href="http://glast.stanford.edu/cgi-bin/cvsweb/Likelihood/python/SrcModel.py?hideattic=1&cvsroot=CVS_SLAC">SrcModel.py</a>, provide something fairly reasonable.

   @section demo Interface Demo

   There are two classes exposed directly to the user in this
   implementation.  The <tt>Observation</tt> class contains and
   provides a single access point for all of the data associated with
   a specific observation. In this context, an "observation" is
   defined in terms of the extraction region in the data space of
   photon arrival time, measured energy, and direction.  An
   observation also comprises ancilliary information specific to the 
   data selections such as the exposure map and response functions to
   be used.  Its constructor looks like
@verbatim
class Observation(object):
    def __init__(self, eventFile=None, scFile=None, expMap=None, irfs='TEST'):
@endverbatim
   and has these parameters
   - <tt>eventFile</tt>: event data file name(s)
   - <tt>scFile</tt>: spacecraft data file name(s)
   - <tt>expMap</tt>: exposure map (produced by @b gtexpmap)
   - <tt>irfs</tt>: response functions, e.g., <tt>DC1</tt>, <tt>G25</tt>,
     or <tt>TEST</tt>

   One creates an <tt>Observation</tt> object like this:
@verbatim
>>> my_obs = Observation(eventFiles, 'demo_scData_0000.fits', 'demo_expMap.fits', 'TEST')
@endverbatim
   Here, <tt>eventFiles</tt> is an ascii file containing the names of
   the event files,
@verbatim
salathe[jchiang] cat eventFiles
eg_diffuse_events_0000.fits
galdiffuse_events_0000.fits
ptsrcs_events_0000.fits
salathe[jchiang] 
@endverbatim
   One could also have entered in a tuple or list if file names or use
   <tt>glob</tt> to generate the list.  For these data, the following 
   are equivalent to the above:
@verbatim
>>> my_obs = Observation(('eg_diffuse_events_0000.fits', 
                          'galdiffuse_events_0000.fits',
                          'ptsrcs_events_0000.fits'),
                         'demo_scData_0000.fits', 'demo_expMap.fits', 'TEST')
@endverbatim
or
@verbatim
>>> my_obs = Observation(glob.glob('*events*.fits'), 'demo_scData_0000.fits',
                         'demo_expMap.fits', 'TEST')
@endverbatim

   One may omit all of the arguments in the <tt>Observation</tt> class
   constructor, in which case a small GUI dialog is launched that
   allows one to browse the file system, use wild cards for specifying
   groups of files, etc..  So, entering
@verbatim
>>> my_obs = Observation()
@endverbatim
   launches this dialog:

   @image html Obs_dialog0.png
   \n
   The buttons on the left will open file dialog boxes using the value
   shown in the text entry fields as a filter.  One can include wild
   cards or comma-separated lists of files in the text entry field.  
   Setting the entries as follows will create an <tt>Observation</tt> 
   object equivalent to the previous examples:

   @image html Obs_dialog1.png

   \n\n
   The second class is <tt>SrcAnalysis</tt>,
@verbatim
class SrcAnalysis(object):
    def __init__(self, observation, srcModel=None, optimizer='Minuit'):
@endverbatim
    It has parameters
    - <tt>observation</tt>: an Observation object
    - <tt>srcModel</tt>: xml source model file
    - <tt>optimizer</tt>: <tt>'Minuit'</tt>, <tt>'Drmngb'</tt>, 
      or <tt>'Lbfgs'</tt>.

    Factoring into separate analysis and observation classes affords
    one considerable flexibility in mixing and matching observations
    and models in a single Python session or script while preserving
    computational resources, so that one can do something like this:
@verbatim
>>> analysis1 = SrcAnalysis(my_obs, "model1.xml")
>>> analysis2 = SrcAnalysis(my_obs, "model2.xml")
@endverbatim

    The two <tt>SrcAnalysis</tt> objects will access the same
    observational data in memory, but at the same time allow for
    distinct source models to be fit to those data without interfering
    with one another.  This will be useful in comparing, via a
    likelihood ratio test, for example, how well one model compares to
    another.

    Just as with the <tt>Observation</tt> class, a small GUI is provided
    that allows one to browse the file system in case one forgets the
    name of the model XML file, for example:

@verbatim
>>> analysis = SrcAnalysis(my_obs)
@endverbatim
    launches

    @image html SrcAnalysis_dialog.png
    \n\n

    The <tt>Observation</tt> and <tt>SrcAnalysis</tt> classes both have
    <tt>__repr__</tt> methods implemented to allow one to see easily what
    data these objects comprise:

@verbatim
>>> print my_obs
['galdiffuse_events_0000.fits', 'eg_diffuse_events_0000.fits', 'ptsrcs_events_0000.fits']
['demo_scData_0000.fits']
demo_expMap.fits
TEST
>>> 
>>> print analysis
demo_model.xml
['galdiffuse_events_0000.fits', 'eg_diffuse_events_0000.fits', 'ptsrcs_events_0000.fits']
['demo_scData_0000.fits']
demo_expMap.fits
TEST
Minuit
>>>
@endverbatim

For this demo, I've prepared a small script that creates the
<tt>Observation</tt> and <tt>SrcAnalysis</tt> objects with the proper
parameter values:

@verbatim
salathe[jchiang] cat demo.py
import glob
from SrcAnalysis import *

eventFiles = glob.glob('*events*.fits')
obs = Observation(eventFiles, 'demo_scData_0000.fits',
                  expMap='demo_expMap.fits', irfs='TEST')
like = SrcAnalysis(obs, 'demo_model.xml')
@endverbatim

This just needs to be imported at the Python prompt (after sourcing the 
package setup script to set the PYTHONPATH environment variable):

@verbatim
salathe[jchiang] source ../../Likelihood/v7r1p3/cmt/setup.csh
salathe[jchiang] python
Python 2.3.3 (#1, Apr 24 2004, 23:59:52) 
[GCC 3.3.2 20031022 (Red Hat Linux 3.3.2-1)] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from demo import *
@endverbatim

<tt>dir(...)</tt> shows an object's attributes:

@verbatim
>>> dir(like)
['_Nobs', '__call__', '__class__', '__delattr__', '__dict__', '__doc__', 
'__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', 
'__reduce_ex__', '__repr__', '__setattr__', '__str__', '__weakref__', '_errors', 
'_fileList', '_plot_model', '_plot_residuals', '_readData', '_readEvents', 
'_readScData', '_srcCnts', 'disp', 'energies', 'events', 'fit', 'logLike', 
'model', 'optimizer', 'plot', 'plotData', 'resids']
@endverbatim

Let's create some references to the object's methods for Xspec-like syntax:

@verbatim
>>> model = like.model
>>> plot = like.plot
>>> fit = like.fit
@endverbatim

<tt>SrcModel.__repr__()</tt> has been overloaded to provide a summary of
the source model.  The columns are

- parameter index
- parameter name
- parameter value
- error estimate (zero before fitting)
- lower bound 
- upper bound 
- parameter scaling (in parentheses)
- fixed (or not)

@verbatim
>>> model
Extragalactic Diffuse
   Spectrum: PowerLaw
0      Prefactor:  1.450e+00  0.000e+00  1.000e-05  1.000e+02 ( 1.000e-07)
1          Index: -2.100e+00  0.000e+00 -3.500e+00 -1.000e+00 ( 1.000e+00) fixed
2          Scale:  1.000e+02  0.000e+00  5.000e+01  2.000e+02 ( 1.000e+00) fixed

Galactic Diffuse
   Spectrum: PowerLaw
3      Prefactor:  1.100e+01  0.000e+00  1.000e-03  1.000e+03 ( 1.000e-03)
4          Index: -2.100e+00  0.000e+00 -3.500e+00 -1.000e+00 ( 1.000e+00) fixed
5          Scale:  1.000e+02  0.000e+00  5.000e+01  2.000e+02 ( 1.000e+00) fixed

my_3EG_J0530p1323
   Spectrum: PowerLaw
6      Prefactor:  1.365e+01  0.000e+00  1.000e-05  1.000e+03 ( 1.000e-09)
7          Index: -2.460e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
8          Scale:  1.000e+02  0.000e+00  3.000e+01  2.000e+03 ( 1.000e+00) fixed

my_3EG_J0534p2200
   Spectrum: BrokenPowerLaw
9      Prefactor:  8.000e-02  0.000e+00  1.000e-03  1.000e+03 ( 1.000e-09)
10        Index1: -2.190e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
11        Index2: -4.890e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
12    BreakValue:  1.500e+03  0.000e+00  3.000e+01  2.000e+03 ( 1.000e+00)

my_3EG_J0633p1751
   Spectrum: BrokenPowerLaw
13     Prefactor:  2.000e+00  0.000e+00  1.000e-03  1.000e+03 ( 1.000e-09)
14        Index1: -1.660e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
15        Index2: -3.100e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
16    BreakValue:  2.000e+03  0.000e+00  3.000e+01  2.000e+03 ( 1.000e+00)

@endverbatim

One can set parameter values using the index:

@verbatim
>>> model[0] = 1.595
@endverbatim

The <tt>plot()</tt> method has been bound to HippoDraw, but any suitable 
plotting package, e.g., matplotlib, Biggles, pyROOT, could be used. On SLAC
linux, HippoDraw has been installed system-wide and should be available 
from Python.

@verbatim
>>> plot()
@endverbatim

@image html demo_plot_1.png
\n\n
optimizer::Parameter member functions are automatically dispatched to
from Python (via the <tt>__getattr__</tt> method):

@verbatim
>>> model[0].setFree(0)
>>> model[16].setBounds(300, 3000)
@endverbatim

Let's create a couple more Xspec-like commands for convenience, using Python's
<tt>lambda</tt> expressions:

@verbatim
>>> thaw = lambda x: model[x].setFree(1)
>>> freeze = lambda x: model[x].setFree(0)
>>> freeze(3)
>>> model
Extragalactic Diffuse
   Spectrum: PowerLaw
0      Prefactor:  1.595e+00  0.000e+00  1.000e-05  1.000e+02 ( 1.000e-07) fixed
1          Index: -2.100e+00  0.000e+00 -3.500e+00 -1.000e+00 ( 1.000e+00) fixed
2          Scale:  1.000e+02  0.000e+00  5.000e+01  2.000e+02 ( 1.000e+00) fixed

Galactic Diffuse
   Spectrum: PowerLaw
3      Prefactor:  1.100e+01  0.000e+00  1.000e-03  1.000e+03 ( 1.000e-03) fixed
4          Index: -2.100e+00  0.000e+00 -3.500e+00 -1.000e+00 ( 1.000e+00) fixed
5          Scale:  1.000e+02  0.000e+00  5.000e+01  2.000e+02 ( 1.000e+00) fixed

my_3EG_J0530p1323
   Spectrum: PowerLaw
6      Prefactor:  1.365e+01  0.000e+00  1.000e-05  1.000e+03 ( 1.000e-09)
7          Index: -2.460e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
8          Scale:  1.000e+02  0.000e+00  3.000e+01  2.000e+03 ( 1.000e+00) fixed

my_3EG_J0534p2200
   Spectrum: BrokenPowerLaw
9      Prefactor:  8.000e-02  0.000e+00  1.000e-03  1.000e+03 ( 1.000e-09)
10        Index1: -2.190e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
11        Index2: -4.890e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
12    BreakValue:  1.500e+03  0.000e+00  3.000e+01  2.000e+03 ( 1.000e+00)

my_3EG_J0633p1751
   Spectrum: BrokenPowerLaw
13     Prefactor:  2.000e+00  0.000e+00  1.000e-03  1.000e+03 ( 1.000e-09)
14        Index1: -1.660e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
15        Index2: -3.100e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
16    BreakValue:  2.000e+03  0.000e+00  3.000e+02  3.000e+03 ( 1.000e+00)

@endverbatim

Now, perform a fit, suppressing Minuit's screen output; the resulting value 
<tt>-log(Likelihood)</tt> is returned:

@verbatim
>>> fit(verbosity=0)
59334.689815397171
@endverbatim

Overplot the new result in red, the default (note the IDL-like syntax):

@verbatim
>>> plot(oplot=True)
@endverbatim

@image html demo_plot_2.png
\n\n

A specific source can be accessed by its full name or by a fragment thereof:

@verbatim
>>> print model["633"]
my_3EG_J0633p1751
   Spectrum: BrokenPowerLaw
13     Prefactor:  1.416e-01  8.264e-02  1.000e-03  1.000e+03 ( 1.000e-09)
14        Index1: -1.722e+00  6.136e-02 -5.000e+00 -1.000e+00 ( 1.000e+00)
15        Index2: -3.143e+00  8.971e-01 -5.000e+00 -1.000e+00 ( 1.000e+00)
16    BreakValue:  2.000e+03  6.771e+02  3.000e+02  3.000e+03 ( 1.000e+00)
@endverbatim

The <tt>logLike</tt> attribute is a Likelihood::LogLike object, and so its
member functions are exposed as Python methods:

@verbatim
>>> like.logLike.writeXml("fitted_model.xml")
>>> 
@endverbatim
*/

