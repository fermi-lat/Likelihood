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

    // Not sure how to handle err at this time
    // Not sure how to handle 'line style'
    const double s=1.0; // Marker size in points**2
    const std::map<std::string, std::string> & keywords = {}; // keyword map

    std::vector<const int> colorvec;
    for (int i=0; i<x.size(); i++)
      colorvec.push_back(color);
    
    plt::scatter_colored(x, y, colorvec, s, keywords);
    //plt::scatter(x, y, s, keywords); // Ignore colorvector for now

    plt::show();
    plt::save("./output.png");
  }

  // Scatter Type II
  void MPLPlot::scatter(const std::vector<double> & x,
	       const std::vector<double> & y,
	       const std::vector<double> & yerr,
	       const int color,
	       const std::string & line_style)
  {
    // Not sure how to handle err at this time                                                                                                     
    // Not sure how to handle 'line style'                                                                                                         
    const double s=1.0; // Marker size in points**2                                                                                                
    const std::map<std::string, std::string> & keywords = {}; // keyword map

    std::vector<const int> colorvec;
    for (int i=0; i<x.size(); i++)
      colorvec.push_back(color);
  
    plt::scatter_colored(x, y, colorvec, s, keywords);
    //plt::scatter(x, y, s, keywords); // Ignore colorvector for now
    plt::show();
    plt::save("./output.png");
  }

  // Scatter Type III
  void MPLPlot::scatter(const std::vector<double> & x,
	       const std::vector<double> & y,
	       const int color,
	       const std::string & line_style)
  {
    // Not sure how to handle err at this time                                                                                                     
    // Not sure how to handle 'line style'                                                                                                         
    const double s=1.0; // Marker size in points**2                                                                                                
    const std::map<std::string, std::string> & keywords = {}; // keyword map

    std::vector<const int> colorvec;
    for (int i=0; i<x.size(); i++)
      colorvec.push_back(color);
    
    plt::scatter_colored(x, y, colorvec, s, keywords);
    //plt::scatter(x, y, s, keywords); // Ignore colorvector for now
    plt::save("./output.png");
  }

  // Line Type I
  void MPLPlot::linePlot(const std::vector<double> & x,
		const std::vector<double> & y,
		int color,
		const std::string & line_style)
  {
    plt::plot(x, y, line_style);
    plt::save("./output.png");
  }

  // Histogram I
  void MPLPlot::histogram(const std::vector<long> &x,
		 const std::vector<double> &y,
		 const std::vector<double> &xwidth,
		 int color,
		 const std::string & line_style)
  {
    plt::hist(y); // Do we need width?
    plt::save("./output.png");
  }

  // Scatter + Error Bars
  //void scatterPlotErrorBars(const std::vector<double> & x,
  //                          std::vector<double> & xerr,
  //                          unsigned int nbins) const
  //{

} // namespace EasyPlot

