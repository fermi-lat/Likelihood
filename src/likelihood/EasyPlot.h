/**
 * @file EasyPlot.h
 * @brief Friendly user interface to st_graph plotting.
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header$
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
 * $Header$
 */

class EasyPlot {

public:

   EasyPlot(const std::string &title,
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

   static void run();

private:

   st_graph::IFrame * m_mainFrame;
   st_graph::IFrame * m_plotFrame;
   std::vector<st_graph::IPlot *> m_plots;

   void scatterPlotErrorBars(const std::vector<double> & x,
                             std::vector<double> & xerr,
                             unsigned int nbins=100) const;

};
