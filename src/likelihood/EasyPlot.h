/**
 * @file EasyPlot.h
 * @brief Friendly user interface to st_graph plotting.
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/likelihood/EasyPlot.h,v 1.4 2006/06/27 15:56:05 peachey Exp $
 */

#include "st_graph/IPlot.h"

namespace st_graph {
   class IFrame;
}

/**
 * @class EasyPlot
 * @brief Simple interface to st_graph plotting.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/likelihood/EasyPlot.h,v 1.4 2006/06/27 15:56:05 peachey Exp $
 */

class EasyPlot {

public:

   EasyPlot(st_graph::IFrame * mainFrame, const std::string &title,
            bool logX = false, bool logY = false,
            const std::string &xAxisTitle = "", const std::string &yAxisTitle = "",
            unsigned int xsize=400, unsigned int ysize=400);

   ~EasyPlot() throw();

   void scatter(const std::vector<double> & x, 
                const std::vector<double> & y,
                const std::vector<double> & xerr, 
                const std::vector<double> & yerr,
                int color = st_graph::Color::eBlack,
                const std::string & line_style = "none");

   void scatter(const std::vector<double> & x, 
                const std::vector<double> & y,
                const std::vector<double> & yerr,
                int color = st_graph::Color::eBlack,
                const std::string & line_style = "none");

   void scatter(const std::vector<double> & x, 
                const std::vector<double> & y,
                int color = st_graph::Color::eBlack,
                const std::string & line_style = "none");

   void linePlot(const std::vector<double> & x,
                 const std::vector<double> & y,
                 int color = st_graph::Color::eBlack,
                 const std::string & line_style = "solid");

   void histogram(const std::vector<double> & x,
                  const std::vector<double> & y,
                  int color = st_graph::Color::eBlack,
                  const std::string & line_style = "solid");

   void histogram(const std::vector<double> & x,
                  const std::vector<double> & y,
                  const std::vector<double> & xerr,
                  int color = st_graph::Color::eBlack,
                  const std::string & line_style = "solid");

   st_graph::IFrame * getPlotFrame();

   static void run();

private:

   st_graph::IFrame * m_mainFrame;
   st_graph::IFrame * m_plotFrame;
   std::vector<st_graph::IPlot *> m_plots;
   bool m_logX;
   bool m_logY;
   std::string m_xAxisTitle;
   std::string m_yAxisTitle;

   void scatterPlotErrorBars(const std::vector<double> & x,
                             std::vector<double> & xerr,
                             unsigned int nbins=100) const;

   void setScale(st_graph::IPlot * plot);
};
