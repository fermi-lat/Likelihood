/**
 * @file EasyPlot.h
 * @brief Friendly user interface to st_graph plotting.
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/likelihood/EasyPlot.h,v 1.1 2004/09/21 14:51:20 jchiang Exp $
 */

namespace st_graph {
   class IFrame;
   class IPlot;
}

class EasyPlot {

public:

   EasyPlot(unsigned int xsize=400, unsigned int ysize=400);

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

   void scatterPlotSymbolSizes(const std::vector<double> & x,
                               std::vector<double> & xerr,
                               unsigned int nbins=100) const;

};
