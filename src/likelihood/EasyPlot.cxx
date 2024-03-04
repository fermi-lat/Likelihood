#include "matplotlibcpp.h"
#include "EasyPlot.h" 
#include <vector>

namespace plt = matplotlibcpp;

/* Patterning file after the plotting definitions in EasyPlot.cxx
 * Need to make sure that references to that file get repointed to 
 * this one
 */

namespace EasyPlot {
  // Scatter Type I
  void MPLPlot::scatter(const std::vector<double> & x,
	       const std::vector<double> & y,
	       const std::vector<double> & xerr,
	       const std::vector<double> & yerr,
	       const int color,
	       const std::string & line_style)
  {
    // Keyword Map
    const std::map<std::string, std::string> & keywords = {
      {"fmt", "o"},
      {"color", std::to_string(color)},
      {"linestyle", line_style}
    };
    plt::errorbar(x, y, yerr, xerr, keywords);
    plt::show();
  }

  // Scatter Type II
  void MPLPlot::scatter(const std::vector<double> & x,
	       const std::vector<double> & y,
	       const std::vector<double> & yerr,
	       const int color,
	       const std::string & line_style)
  {
    // Keyword Map
    const std::map<std::string, std::string> & keywords = {
      {"fmt", "o"},
      {"color", std::to_string(color)},
      {"linestyle", line_style}
    };
    plt::errorbar(x, y, yerr, keywords);
    plt::show();
  }

  // Scatter Type III
  void MPLPlot::scatter(const std::vector<double> & x,
	       const std::vector<double> & y,
	       const int color,
	       const std::string & line_style)
  {
    // Keyword Map
    const std::map<std::string, std::string> & keywords = {
      {"fmt", "o"},
      {"color", std::to_string(color)},
      {"linestyle", line_style}
    };
    const double s=1.0; // The marker size in points**2                    
    plt::scatter(x, y, s, keywords);
    plt::show();
  }

  // Line Type I
  void MPLPlot::linePlot(const std::vector<double> & x,
		const std::vector<double> & y,
		int color,
		const std::string & line_style)
  {
    // Keyword Map
    const std::map<std::string, std::string> & keywords = {
      {"color", std::to_string(color)},
      {"linestyle", line_style}
    };
    plt::plot(x, y, keywords);
    plt::show();
  }

  // // Histogram I
  // void MPLPlot::histogram(const std::vector<long> &x,
  // 		 const std::vector<double> &y,
  // 		 const std::vector<double> &xwidth,
  // 		 int color,
  // 		 const std::string & line_style)
  // {
  //   plt::hist(y); // Do we need width?
  //   plt::save("./output.png");
  // }

  // Scatter + Error Bars
  //void scatterPlotErrorBars(const std::vector<double> & x,
  //                          std::vector<double> & xerr,
  //                          unsigned int nbins) const
  //{

} // namespace EasyPlot

