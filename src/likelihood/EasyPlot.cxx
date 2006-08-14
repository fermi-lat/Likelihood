/**
 * @file EasyPlot.cxx
 * @brief Implementation of a friendly user interface to st_graph.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/likelihood/EasyPlot.cxx,v 1.10 2006/07/03 17:05:32 jchiang Exp $
 */
#include <iostream>
#include <stdexcept>

#include "st_graph/Axis.h"
#include "st_graph/Engine.h"
#include "st_graph/IEventReceiver.h"
#include "st_graph/IFrame.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

#include "EasyPlot.h"

typedef std::vector<double> vctr;
typedef st_graph::ValueSpreadSequence<vctr::const_iterator> valuesWithErrors;
typedef st_graph::LowerBoundSequence<vctr::const_iterator> lowerBounds;
typedef st_graph::PointSequence<vctr::const_iterator> values;

EasyPlot::EasyPlot(st_graph::IFrame * mainFrame, const std::string & title,
                   bool logX, bool logY,
                   const std::string &xAxisTitle, const std::string &yAxisTitle,
                   unsigned int xsize, unsigned int ysize)
   : m_mainFrame(mainFrame), m_plotFrame(0), m_logX(logX), m_logY(logY), m_xAxisTitle(xAxisTitle), m_yAxisTitle(yAxisTitle) {
   st_graph::Engine & engine(st_graph::Engine::instance());
   m_plotFrame = engine.createPlotFrame(m_mainFrame, title, xsize, ysize);
}

EasyPlot::~EasyPlot() throw() {
   try {
      for (unsigned int i = 0; i < m_plots.size(); i++) {
         delete m_plots[i];
      }
      delete m_plotFrame;
   } catch (std::exception & eObj) {
      std::cerr << eObj.what() << std::endl;
   } catch (...) {
   }
}

void EasyPlot::scatter(const std::vector<double> & x,
                       const std::vector<double> & y,
                       const std::vector<double> & xerr,
                       const std::vector<double> & yerr,
                       int color,
                       const std::string & line_style) {
   st_graph::Engine & engine(st_graph::Engine::instance());
   st_graph::IPlot * plot = 
      engine.createPlot(m_plotFrame, "scat", 
                        valuesWithErrors(x.begin(), x.end(), xerr.begin()),
                        valuesWithErrors(y.begin(), y.end(), yerr.begin()));
   setScale(plot);
   plot->setLineColor(color);
   plot->setLineStyle(line_style);
   m_plots.push_back(plot);
}

void EasyPlot::scatter(const std::vector<double> & x,
                       const std::vector<double> & y,
                       const std::vector<double> & yerr,
                       int color,
                       const std::string & line_style) {
   st_graph::Engine & engine(st_graph::Engine::instance());
   st_graph::IPlot * plot = 
      engine.createPlot(m_plotFrame, "scat", 
                        values(x.begin(), x.end()),
                        valuesWithErrors(y.begin(), y.end(), yerr.begin()));
   setScale(plot);
   plot->setLineColor(color);
   plot->setLineStyle(line_style);
   m_plots.push_back(plot);
}

void EasyPlot::scatter(const std::vector<double> & x,
                       const std::vector<double> & y,
                       int color,
                       const std::string & line_style) {
   if (x.size() == 0) {
      return;
   }
   std::vector<double> xerr;
   scatterPlotErrorBars(x, xerr);
   std::vector<double> yerr;
   scatterPlotErrorBars(y, yerr);
   scatter(x, y, xerr, yerr, color, line_style);
}

void EasyPlot::linePlot(const std::vector<double> & x,
                        const std::vector<double> & y,
                        int color,
                        const std::string & line_style) {
   if (x.size() == 0) {
      return;
   }
   st_graph::Engine & engine(st_graph::Engine::instance());
   st_graph::IPlot * plot = 
      engine.createPlot(m_plotFrame, "scat", values(x.begin(), x.end()),
                        values(y.begin(), y.end()));
   setScale(plot);
   plot->setLineColor(color);
   plot->setLineStyle(line_style);
   plot->setCurveType("curve");
   m_plots.push_back(plot);
}
void EasyPlot::histogram(const std::vector<double> &x,
                         const std::vector<double> &y,
                         const std::vector<double> &xwidth,
                         int color,
                         const std::string & line_style) {
   st_graph::Engine & engine(st_graph::Engine::instance());
#if 0
   std::vector<double> xx;
   for (unsigned int i = 0; i < x.size(); i++) {
      xx.push_back(x[i] - xwidth[i]/2.);
   }
#endif
   st_graph::IPlot * plot = 
      engine.createPlot(m_plotFrame, "hist", 
                        lowerBounds(x.begin(), x.end()),
                        values(y.begin(), y.end()));
   plot->setLineColor(color);
   plot->setLineStyle(line_style);
   setScale(plot);
   m_plots.push_back(plot);
}

void EasyPlot::histogram(const std::vector<double> & x, 
                         const std::vector<double> & y,
                         int color,
                         const std::string & line_style) {
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

   histogram(xx, y, xerr, color, line_style);
}

st_graph::IFrame * EasyPlot::getPlotFrame() { return m_plotFrame; }

void EasyPlot::scatterPlotErrorBars(const std::vector<double> & x,
                                    std::vector<double> & xerr,
                                    unsigned int nbins) const {
   (void)(nbins);
   unsigned int npts(x.size());
   double xmax(x[0]);
   double xmin(x[0]);
   for (unsigned int i = 0; i < npts; i++) {
      if (x[i] < xmin) xmin = x[i];
      if (x[i] > xmin) xmax = x[i];
   }
//   double xstep = (xmax - xmin)/nbins;
   double xstep = (xmax - xmin)/npts;
   xerr.reserve(npts);
   for (unsigned int i = 0; i < npts; i++) {
      xerr.push_back(xstep);
   }
   return;
}

void EasyPlot::setScale(st_graph::IPlot * plot) {
  std::vector<st_graph::Axis> & axes(plot->getAxes());
  axes.at(0).setScaleMode(m_logX ? st_graph::Axis::eLog : st_graph::Axis::eLinear);
  axes.at(0).setTitle(m_xAxisTitle);
  axes.at(1).setScaleMode(m_logY ? st_graph::Axis::eLog : st_graph::Axis::eLinear);
  axes.at(1).setTitle(m_yAxisTitle);
}

void EasyPlot::run() {
   st_graph::Engine & engine(st_graph::Engine::instance());
   engine.run();
}
