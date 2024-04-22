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

namespace EasyPlot {
class Color {

public:

  enum Color_e { eWhite, eBlack, eRed, eGreen, eBlue, eYellow, eMagenta, eCyan, eNumberOfColors };

  static int nextColor(int current_color);

  static char getColor(int current_color);

};
 
class MPLPlot {

public:

   static void scatter(const std::vector<double> & x, 
		       const std::vector<double> & y,
		       const std::vector<double> & xerr, 
		       const std::vector<double> & yerr,
		       char color = 'k',
		       const std::string & line_style = "none");

   static void scatter(const std::vector<double> & x, 
		       const std::vector<double> & y,
		       const std::vector<double> & yerr,
		       char color = 'k',
		       const std::string & line_style = "none");

   static void scatter(const std::vector<double> & x, 
		       const std::vector<double> & y,
		       char color = 'k',
		       const std::string & line_style = "none");

   static void linePlot(const std::vector<double> & x,
			const std::vector<double> & y,
		        char color = 'k',
			const std::string & line_style = "solid");

   static void logLog(const std::vector<double> & x,
		      const std::vector<double> & y,
		      char color = 'k');

   static void semilogx(const std::vector<double> & x,
		        const std::vector<double> & y,
		        char color,
			const std::string & line_style,
			const std::string & plot_title);
  
   static void ylim(const double left, const double right);

   static void xlim(const double left, const double right);

   static void subplot(long nrows, long ncols, long plot_number);

   static void figure(long number);
  
   static void showPlot();
   // void histogram(const std::vector<double> & x,
   //                const std::vector<double> & y,
   //                int color = st_graph::Color::eBlack,
   //                const std::string & line_style = "solid");

   // static void histogram(const std::vector<long> & x,
   // 			 const std::vector<double> & y,
   // 			 const std::vector<double> & xwidth,
   // 			 int color = st_graph::Color::eBlack,
   // 			 const std::string & line_style = "solid");

   st_graph::IFrame * getPlotFrame();

   static void run();

private:

  // st_graph::IFrame * m_mainFrame;
  //st_graph::IFrame * m_plotFrame;
  //std::vector<st_graph::IPlot *> m_plots;
   bool m_logX;
   bool m_logY;
   std::string m_xAxisTitle;
   std::string m_yAxisTitle;

   // void scatterPlotErrorBars(const std::vector<double> & x,
   //                           std::vector<double> & xerr,
   //                           unsigned int nbins=100) const;

   //void setScale(st_graph::IPlot * plot);
};

} // namespace EasyPlot
