/**
   @page pyLike Likelihood and Python

   @section intro Introduction

   Using <a href="http://www.swig.org/">SWIG</a> via the <a
   href="http://www.slac.stanford.edu/~jchiang/SwigPolicy/doxy-html/main.html">SwigPolicy</a>
   package, I've exposed the C++ classes of Likelihood to <a
   href="http://www.python.org/">Python</a>.  These classes by
   themselves do not constitute a suitable interactive interface.
   However, with a modest amount effort, I've written a couple of
   small python modules using these classes, SrcAnalysis.py and
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
    - <tt>optimizer</tt>: <tt>'Minuit'</tt>, <tt>'Dmnrgb'</tt>, 
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

Create some references to the object's methods for Xspec-like syntax:

@verbatim
>>> model = like.model; plot = like.plot; fit = like.fit
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
   Spectrum: PowerLaw
9      Prefactor:  2.700e+01  0.000e+00  1.000e-05  1.000e+03 ( 1.000e-09)
10         Index: -2.190e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
11         Scale:  1.000e+02  0.000e+00  3.000e+01  2.000e+03 ( 1.000e+00) fixed

my_3EG_J0633p1751
   Spectrum: PowerLaw
12     Prefactor:  2.329e+01  0.000e+00  1.000e-05  1.000e+03 ( 1.000e-09)
13         Index: -1.660e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
14         Scale:  1.000e+02  0.000e+00  3.000e+01  2.000e+03 ( 1.000e+00) fixed
@endverbatim

One can set parameter values using the numeric index:

@verbatim
>>> model[0] = 1.595
@endverbatim

The <tt>plot()</tt> method has been bound to HippoDraw, but any plotting
package can be used.

@verbatim
>>> plot()
@endverbatim

@image html demo_plot_1.png
\n\n
optimizer::Parameter member functions are automatically dispatched to
from Python (via the <tt>__getattr__</tt> method):

@verbatim
>>> model[1].setFree(1)
>>> model[3].setFree(0)
>>> model[4].setFree(1)
>>> model
Extragalactic Diffuse
   Spectrum: PowerLaw
0      Prefactor:  1.595e+00  0.000e+00  1.000e-05  1.000e+02 ( 1.000e-07)
1          Index: -2.100e+00  0.000e+00 -3.500e+00 -1.000e+00 ( 1.000e+00)
2          Scale:  1.000e+02  0.000e+00  5.000e+01  2.000e+02 ( 1.000e+00) fixed

Galactic Diffuse
   Spectrum: PowerLaw
3      Prefactor:  1.100e+01  0.000e+00  1.000e-03  1.000e+03 ( 1.000e-03) fixed
4          Index: -2.100e+00  0.000e+00 -3.500e+00 -1.000e+00 ( 1.000e+00)
5          Scale:  1.000e+02  0.000e+00  5.000e+01  2.000e+02 ( 1.000e+00) fixed

my_3EG_J0530p1323
   Spectrum: PowerLaw
6      Prefactor:  1.365e+01  0.000e+00  1.000e-05  1.000e+03 ( 1.000e-09)
7          Index: -2.460e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
8          Scale:  1.000e+02  0.000e+00  3.000e+01  2.000e+03 ( 1.000e+00) fixed

my_3EG_J0534p2200
   Spectrum: PowerLaw
9      Prefactor:  2.700e+01  0.000e+00  1.000e-05  1.000e+03 ( 1.000e-09)
10         Index: -2.190e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
11         Scale:  1.000e+02  0.000e+00  3.000e+01  2.000e+03 ( 1.000e+00) fixed

my_3EG_J0633p1751
   Spectrum: PowerLaw
12     Prefactor:  2.329e+01  0.000e+00  1.000e-05  1.000e+03 ( 1.000e-09)
13         Index: -1.660e+00  0.000e+00 -5.000e+00 -1.000e+00 ( 1.000e+00)
14         Scale:  1.000e+02  0.000e+00  3.000e+01  2.000e+03 ( 1.000e+00) fixed
@endverbatim

Now, perform a fit, suppressing Minuit's screen output; the resulting value 
<tt>-log(Likelihood)</tt> is returned:

@verbatim
>>> fit(verbosity=0)
59350.111358513532
@endverbatim

Overplot the new result in red (note the IDL-like syntax):

@verbatim
>>> plot(oplot=True, color='red')
@endverbatim

@image html demo_plot_2.png
\n\n
@verbatim
>>> model
Extragalactic Diffuse
   Spectrum: PowerLaw
0      Prefactor:  1.513e+00  1.415e-01  1.000e-05  1.000e+02 ( 1.000e-07)
1          Index: -1.922e+00  6.930e-02 -3.500e+00 -1.000e+00 ( 1.000e+00)
2          Scale:  1.000e+02  0.000e+00  5.000e+01  2.000e+02 ( 1.000e+00) fixed

Galactic Diffuse
   Spectrum: PowerLaw
3      Prefactor:  1.100e+01  0.000e+00  1.000e-03  1.000e+03 ( 1.000e-03) fixed
4          Index: -2.135e+00  2.585e-02 -3.500e+00 -1.000e+00 ( 1.000e+00)
5          Scale:  1.000e+02  0.000e+00  5.000e+01  2.000e+02 ( 1.000e+00) fixed

my_3EG_J0530p1323
   Spectrum: PowerLaw
6      Prefactor:  1.432e+01  2.200e+00  1.000e-05  1.000e+03 ( 1.000e-09)
7          Index: -2.545e+00  1.339e-01 -5.000e+00 -1.000e+00 ( 1.000e+00)
8          Scale:  1.000e+02  0.000e+00  3.000e+01  2.000e+03 ( 1.000e+00) fixed

my_3EG_J0534p2200
   Spectrum: PowerLaw
9      Prefactor:  2.665e+01  2.687e+00  1.000e-05  1.000e+03 ( 1.000e-09)
10         Index: -2.290e+00  7.180e-02 -5.000e+00 -1.000e+00 ( 1.000e+00)
11         Scale:  1.000e+02  0.000e+00  3.000e+01  2.000e+03 ( 1.000e+00) fixed

my_3EG_J0633p1751
   Spectrum: PowerLaw
12     Prefactor:  2.984e+01  2.488e+00  1.000e-05  1.000e+03 ( 1.000e-09)
13         Index: -1.925e+00  4.423e-02 -5.000e+00 -1.000e+00 ( 1.000e+00)
14         Scale:  1.000e+02  0.000e+00  3.000e+01  2.000e+03 ( 1.000e+00) fixed
@endverbatim

Individual sources can be accessed:

@verbatim
>>> model["Extra"]
Extragalactic Diffuse
   Spectrum: PowerLaw
0      Prefactor:  1.513e+00  1.415e-01  1.000e-05  1.000e+02 ( 1.000e-07)
1          Index: -1.922e+00  6.930e-02 -3.500e+00 -1.000e+00 ( 1.000e+00)
2          Scale:  1.000e+02  0.000e+00  5.000e+01  2.000e+02 ( 1.000e+00) fixed
@endverbatim

The <tt>logLike</tt> attribute @em is a Likelihood::LogLike object and so its
member functions are exposed as Python methods:

@verbatim
>>> like.logLike.writeXml("fitted_model.xml")
>>> 
@endverbatim

For completeness, here's the xml file that was created:

@verbatim
salathe[jchiang] cat fitted_model.xml
<?xml version="1.0" standalone="no"?>
<!DOCTYPE source_library SYSTEM "$(LIKELIHOODROOT)/xml/A1_Sources.dtd" >
<source_library title="source library">
  <source name="Extragalactic Diffuse" type="DiffuseSource">
    <spectrum type="PowerLaw">
      <parameter error="0.141459" free="1" max="100" min="1e-05" name="Prefactor" scale="1e-07" value="1.5127" />
      <parameter error="0.0692986" free="1" max="-1" min="-3.5" name="Index" scale="1" value="-1.92173" />
      <parameter free="0" max="200" min="50" name="Scale" scale="1" value="100" />
    </spectrum>
    <spatialModel type="ConstantValue">
      <parameter free="0" max="10" min="0" name="Value" scale="1" value="1" />
    </spatialModel>
  </source>
  <source name="Galactic Diffuse" type="DiffuseSource">
    <spectrum type="PowerLaw">
      <parameter free="0" max="1000" min="0.001" name="Prefactor" scale="0.001" value="11" />
      <parameter error="0.025852" free="1" max="-1" min="-3.5" name="Index" scale="1" value="-2.13475" />
      <parameter free="0" max="200" min="50" name="Scale" scale="1" value="100" />
    </spectrum>
    <spatialModel file="$(LIKELIHOODROOT)/src/test/Data/gas.cel" type="SpatialMap">
      <parameter free="0" max="1000" min="0.001" name="Prefactor" scale="1" value="1" />
    </spatialModel>
  </source>
  <source name="my_3EG_J0530p1323" type="PointSource">
    <spectrum type="PowerLaw">
      <parameter error="2.19972" free="1" max="1000" min="1e-05" name="Prefactor" scale="1e-09" value="14.3212" />
      <parameter error="0.133909" free="1" max="-1" min="-5" name="Index" scale="1" value="-2.54496" />
      <parameter free="0" max="2000" min="30" name="Scale" scale="1" value="100" />
    </spectrum>
    <spatialModel type="SkyDirFunction">
      <parameter free="0" max="360" min="-360" name="RA" scale="1" value="82.74" />
      <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="13.38" />
    </spatialModel>
  </source>
  <source name="my_3EG_J0534p2200" type="PointSource">
    <spectrum type="PowerLaw">
      <parameter error="2.68722" free="1" max="1000" min="1e-05" name="Prefactor" scale="1e-09" value="26.6549" />
      <parameter error="0.071799" free="1" max="-1" min="-5" name="Index" scale="1" value="-2.29019" />
      <parameter free="0" max="2000" min="30" name="Scale" scale="1" value="100" />
    </spectrum>
    <spatialModel type="SkyDirFunction">
      <parameter free="0" max="360" min="-360" name="RA" scale="1" value="83.57" />
      <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="22.01" />
    </spatialModel>
  </source>
  <source name="my_3EG_J0633p1751" type="PointSource">
    <spectrum type="PowerLaw">
      <parameter error="2.48764" free="1" max="1000" min="1e-05" name="Prefactor" scale="1e-09" value="29.8384" />
      <parameter error="0.0442259" free="1" max="-1" min="-5" name="Index" scale="1" value="-1.9255" />
      <parameter free="0" max="2000" min="30" name="Scale" scale="1" value="100" />
    </spectrum>
    <spatialModel type="SkyDirFunction">
      <parameter free="0" max="360" min="-360" name="RA" scale="1" value="98.49" />
      <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="17.86" />
    </spatialModel>
  </source>
</source_library>
salathe[jchiang] 
@endverbatim

*/
