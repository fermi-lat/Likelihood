/**
 * @file EasyPlot.cxx
 * @brief Implementation of a friendly user interface to st_graph.
 * @author J. Chiang
 *
 * $Header$
 */

#include <iostream>
#include <stdexcept>

#include "st_graph/Engine.h"
#include "st_graph/IEventReceiver.h"
#include "st_graph/IFrame.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

#include "EasyPlot.h"

typedef std::vector<double> vector;
typedef st_graph::ValueSpreadSequence<vector::const_iterator> valuesWithErrors;
typedef st_graph::LowerBoundSequence<vector::const_iterator> lowerBounds;
typedef st_graph::PointSequence<vector::const_iterator> values;

EasyPlot::EasyPlot(const std::string & title,
                   unsigned int xsize, unsigned int ysize)
   : m_mainFrame(0), m_plotFrame(0) {
   st_graph::Engine & engine(st_graph::Engine::instance());
   m_mainFrame = engine.createMainFrame(0, xsize, ysize);
   m_plotFrame = engine.createPlotFrame(m_mainFrame, title, xsize, ysize);
}

EasyPlot::~EasyPlot() throw() {
   try {
      for (unsigned int i = 0; i < m_plots.size(); i++) {
         delete m_plots[i];
      }
      delete m_plotFrame;
      delete m_mainFrame;
   } catch (std::exception & eObj) {
      std::cerr << eObj.what() << std::endl;
   } catch (...) {
   }
}

void EasyPlot::scatter(const std::vector<double> & x,
                       const std::vector<double> & y,
                       const std::vector<double> & xerr,
                       const std::vector<double> & yerr) {
   st_graph::Engine & engine(st_graph::Engine::instance());
   st_graph::IPlot * plot = 
      engine.createPlot(m_plotFrame, "scat", 
                        valuesWithErrors(x.begin(), x.end(), xerr.begin()),
                        valuesWithErrors(y.begin(), y.end(), yerr.begin()));
   m_plots.push_back(plot);
}

void EasyPlot::scatter(const std::vector<double> & x,
                       const std::vector<double> & y,
                       const std::vector<double> & yerr) {
// It would be nice if one could work in plot frame coordinates or if
// one could specify the plotting symbols.  Instead we use data
// coordinates and are forced to guess-timate the proper symbol size
// assuming the plot x-scale is set using the x-axis min and max.
// These guess-timates will likely worsen as other data sets are added
// and the plot x-scale changes.
   if (x.size() == 0) {
      return;
   }
   std::vector<double> xerr;
   scatterPlotErrorBars(x, xerr);
   scatter(x, y, xerr, yerr);
}

void EasyPlot::scatter(const std::vector<double> & x,
                       const std::vector<double> & y) {
   if (x.size() == 0) {
      return;
   }
   std::vector<double> xerr;
   scatterPlotErrorBars(x, xerr);
   std::vector<double> yerr;
   scatterPlotErrorBars(y, yerr);
   scatter(x, y, xerr, yerr);
}

void EasyPlot::linePlot(const std::vector<double> & x,
                        const std::vector<double> & y) {
   if (x.size() == 0) {
      return;
   }
   st_graph::Engine & engine(st_graph::Engine::instance());
   st_graph::IPlot * plot = 
      engine.createPlot(m_plotFrame, "scat", values(x.begin(), x.end()),
                        values(y.begin(), y.end()));
   m_plots.push_back(plot);
}
void EasyPlot::histogram(const std::vector<double> &x,
                         const std::vector<double> &y,
                         const std::vector<double> &xwidth) {
   st_graph::Engine & engine(st_graph::Engine::instance());
   std::vector<double> xx;
   for (unsigned int i = 0; i < x.size(); i++) {
      xx.push_back(x[i] - xwidth[i]/2.);
   }
   st_graph::IPlot * plot = 
      engine.createPlot(m_plotFrame, "hist", 
                        lowerBounds(xx.begin(), xx.end()),
                        values(y.begin(), y.end()));
   m_plots.push_back(plot);
}

void EasyPlot::histogram(const std::vector<double> & x, 
                         const std::vector<double> & y) {
// This routine computes Voronoi cells and resets the x values to the
// cell midpoints since Sequence cannot handle asymmetric "spreads".

// Compute mid-points in x:
   std::vector<double> xmids;
   for (unsigned int i = 0; i < x.size()-1; i++) {
      xmids.push_back((x[i+1] + x[i])/2.);
   }

// Recompute abscissa points.
   std::vector<double> xx;
   xx.push_back(x[0]);
   for (unsigned int i = 0; i < xmids.size()-1; i++) {
      xx.push_back((xmids[i+1] + xmids[i])/2.);
   }
   xx.push_back(x.back());

// Cell widths are computed from the original x grid.
   std::vector<double> xerr;
   xerr.reserve(x.size());
   xerr.push_back(x[1] - x[0]);
   for (unsigned int i = 0; i < x.size()-2; i++) {
      xerr.push_back((x[i+2] - x[i])/2.);
   }
   xerr.push_back(x.back() - *(x.end()-2));

   histogram(xx, y, xerr);
}

void EasyPlot::scatterPlotErrorBars(const std::vector<double> & x,
                                    std::vector<double> & xerr,
                                    unsigned int nbins) const {
   unsigned int npts(x.size());
   double xmax(x[0]);
   double xmin(x[0]);
   for (unsigned int i = 0; i < npts; i++) {
      if (x[i] < xmin) xmin = x[i];
      if (x[i] > xmin) xmax = x[i];
   }
   double xstep = (xmax - xmin)/nbins;
   xerr.reserve(npts);
   for (unsigned int i = 0; i < npts; i++) {
      xerr.push_back(xstep);
   }
   return;
}

void EasyPlot::run() {
   st_graph::Engine & engine(st_graph::Engine::instance());
   engine.run();
}
