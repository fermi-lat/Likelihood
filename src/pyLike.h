/**
   @page pyLike Likelihood and Python

   @section intro Introduction

   Using <a href="http://www.swig.org/">SWIG</a> via the <a
   href="http://www.slac.stanford.edu/~jchiang/SwigPolicy/doxy-html/main.html">SwigPolicy</a>
   package, I've exposed the C++ classes of Likelihood to <a
   href="http://www.python.org/">Python</a>.  These classes by
   themselves do not constitute a suitable interactive interface.
   However, armed with these classes and with a modest effort, I've
   written a couple of small python modules, SrcAnalysis.py and
   SrcModel.py, that provide something fairly reasonable.

   @section demo Interface Demo

   There is one just one class exposed directly to the user in this
   implementation.  Here is the signature of its constructor:

@verbatim
class SrcAnalysis(object):
    def __init__(self, srcModel, eventFile, scFile, expMap=None, irfs='TEST',
                 optimizer='Minuit'):
@endverbatim
    
    The parameters:
    - <tt>srcModel</tt>: xml source model file
    - <tt>eventFile</tt>: event data file name(s)
    - <tt>scFile</tt>: spacecraft data file name(s)
    - <tt>expMap</tt>: exposure map (produced by @b gtexpmap)
    - <tt>irfs</tt>: response functions, e.g., <tt>'DC1'</tt>, <tt>'G25'</tt>,
      or <tt>'TEST'</tt>
    - <tt>optimizer</tt>: <tt>'Minuit'</tt>, <tt>'Drmngb'</tt>, 
      or <tt>'Lbfgs'</tt>.

For this demo, I've prepared a small script that creates the SrcAnalysis
object with the proper parameter values:

@verbatim
salathe[jchiang] cat demo.py
import glob
from SrcAnalysis import SrcAnalysis

eventFiles = glob.glob('*events*.fits')
like = SrcAnalysis('demo_model.xml',
                   eventFiles,
                   'demo_scData_0000.fits',
                   expMap='demo_expMap.fits',
                   irfs='TEST')
@endverbatim

This just needs to be imported at the Python prompt:

@verbatim
salathe[jchiang] python
Python 2.3.3 (#1, Apr 24 2004, 23:59:52) 
[GCC 3.3.2 20031022 (Red Hat Linux 3.3.2-1)] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from demo import *
@endverbatim

Here are the event files found by <tt>glob</tt>:

@verbatim
>>> eventFiles
['galdiffuse_events_0000.fits', 'eg_diffuse_events_0000.fits', 'ptsrcs_events_0000.fits']
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
plotting package, e.g., matplotlib, Biggles, pyROOT, could be used.

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
<tt>lambda</tt> functions:

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

Overplot the new result in red (note the IDL-like syntax):

@verbatim
>>> plot(oplot=True, color='red')
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

