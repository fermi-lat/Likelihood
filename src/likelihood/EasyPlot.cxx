/**
 * @file EasyPlot.cxx
 * @brief Implementation of friendly user interface.
 * @author J. Chiang
 *
 * $Header$
 */

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "st_graph/Engine.h"
#include "st_graph/IEventReceiver.h"
#include "st_graph/IFrame.h"
#include "st_graph/IPlot.h"
#include "st_graph/PlotHist.h"
#include "st_graph/ValueSet.h"

#include "EasyPlot.h"

EasyPlot::EasyPlot(unsigned int xsize, unsigned int ysize)
   : m_mainFrame(0), m_plotFrame(0) {
   st_graph::Engine & engine(st_graph::Engine::instance());
   m_mainFrame = engine.createMainFrame(0, xsize, ysize);
   m_plotFrame = engine.createPlotFrame(m_mainFrame, xsize, ysize);
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
   st_graph::IPlot * plot = engine.createPlot(m_plotFrame, "Scatter", "foo",
                                              st_graph::ValueSet(x, xerr),
                                              st_graph::ValueSet(y, yerr));
   m_plots.push_back(plot);
}

void EasyPlot::scatter(const std::vector<double> & x,
                       const std::vector<double> & y,
                       const std::vector<double> & yerr) {
// It would be nice if one could work in plot frame coordinates or if
// one could specify the plotting symbols.  Instead we use data
// coordinates and are forced to guesstimate the proper symbol size
// assuming the plot x-scale is set using the x-axis min and max.
// This will likely break as other data sets are added and the plot 
// x-scale changes.
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

void EasyPlot::histogram(const std::vector<double> &x,
                         const std::vector<double> &y,
                         const std::vector<double> &xerr) {
   st_graph::Engine & engine(st_graph::Engine::instance());
   std::vector<double> yerr(y.size(), 0);
   st_graph::IPlot * plot = engine.createPlot(m_plotFrame, "hist", "foo",
                                              st_graph::ValueSet(x, xerr),
                                              st_graph::ValueSet(y, yerr));
   m_plots.push_back(plot);
}

void EasyPlot::histogram(const std::vector<double> & x, 
                         const std::vector<double> & y) {

// Compute mid-points in x:
   std::vector<double> xmids;
   for (unsigned int i = 0; i < x.size()-1; i++) {
      xmids.push_back((x[i+1] + x[i])/2.);
   }

// Recompute abscissa points
   std::vector<double> xx;
   xx.push_back(x[0]);
   for (unsigned int i = 0; i < xmids.size()-1; i++) {
      xx.push_back((xmids[i+1] + xmids[i])/2.);
   }
   xx.push_back(x.back());

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
