/**
 * @file EasyPlot.h
 * @brief Friendly user interface to st_graph plotting.
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/likelihood/EasyPlot.h,v 1.3 2004/11/28 06:58:23 jchiang Exp $
 */

namespace st_graph {
   class IFrame;
   class IPlot;
}

/**
 * @class EasyPlot
 * @brief Simple interface to st_graph plotting.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/likelihood/EasyPlot.h,v 1.3 2004/11/28 06:58:23 jchiang Exp $
 */

class EasyPlot {

public:

   EasyPlot(st_graph::IFrame * mainFrame, const std::string &title,
            bool logX = false, bool logY = false,
            unsigned int xsize=400, unsigned int ysize=400);

   ~EasyPlot() throw();

   void scatter(const std::vector<double> & x, 
                const std::vector<double> & y,
                const std::vector<double> & xerr, 
                const std::vector<double> & yerr);

   void scatter(const std::vector<double> & x, 
                const std::vector<double> & y,
                const std::vector<double> & yerr);

   void scatter(const std::vector<double> & x, 
                const std::vector<double> & y);

   void linePlot(const std::vector<double> & x,
                 const std::vector<double> & y);

   void histogram(const std::vector<double> & x,
                  const std::vector<double> & y);

   void histogram(const std::vector<double> & x,
                  const std::vector<double> & y,
                  const std::vector<double> & xerr);

   st_graph::IFrame * getPlotFrame();

   static void run();

private:

   st_graph::IFrame * m_mainFrame;
   st_graph::IFrame * m_plotFrame;
   std::vector<st_graph::IPlot *> m_plots;
   bool m_logX;
   bool m_logY;

   void scatterPlotErrorBars(const std::vector<double> & x,
                             std::vector<double> & xerr,
                             unsigned int nbins=100) const;

   void setScale(st_graph::IPlot * plot);
};
