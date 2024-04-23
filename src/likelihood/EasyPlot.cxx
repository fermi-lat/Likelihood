#include "matplotlibcpp.h"
#include "EasyPlot.h" 
#include <vector>
#include <string>

namespace plt = matplotlibcpp;

/* Patterning file after the plotting definitions in EasyPlot.cxx
 * Need to make sure that references to that file get repointed to 
 * this one
 */

namespace EasyPlot {
  
  using StringMap = std::map<std::string, std::string>;
  //using FuncMap = std::map<std::string,
    
  // Helper function to get char value to properly handle colors
  int Color::nextColor(int current_color) {
    // Go to the next color char representation in the sequence.
    ++current_color;

    // Skip White and Yellow because they don't show up well/at all
    if (eYellow == current_color)
      current_color = ++current_color;
    if (eNumberOfColors == current_color) // This wraps when the last color is reache
      current_color = eBlack;

    return current_color;
  }

  char Color::getColor(int current_color) {
    // Map the index value back to the mpl char shorthand
    char colors[8] = { 'w' /*White*/,
		       'k' /*Black*/,
		       'r' /*Red*/,
		       'g' /*Green*/,
		       'b' /*Blue*/,
		       'y' /*Yellow*/,
		       'm' /*Magenta*/,
		       'c' /*Cyan*/
    };
    return colors[current_color];
  }

  //only put defaults in header file
  StringMap MPLPlot::plotKeywords(const std::string & fmt,
				  char color,
				  const std::string & linestyle)
  {
    std::string c(1,color); 
    StringMap kwmap = {
      {"color", c} // Default is black, represented by 'k'
    };
    if (!fmt.empty()) kwmap.insert({"fmt", fmt});
    if (!linestyle.empty()) kwmap.insert({"linestyle", linestyle});
    return kwmap;
  }
  
  void MPLPlot::plotParams(const std::string & xlabel,
			   const std::string & ylabel,
			   const std::string & xscale,
			   const std::string & yscale,
			   const std::string & title)
  {
    if (!xlabel.empty()) plt::xlabel(xlabel); // Default is ''
    if (!ylabel.empty()) plt::ylabel(ylabel); // Default is ''
    if (!title.empty()) plt::title(title); // Default is ''
    plt::xscale(xscale); // Default is linear
    plt::yscale(yscale); // Default is linear
  }

  // Scatter Type I
  void MPLPlot::scatter(const std::vector<double> & x,
	       const std::vector<double> & y,
	       const std::vector<double> & xerr,
	       const std::vector<double> & yerr,
	       const char color,
	       const std::string & line_style)
  {
    // Keyword Map - scatter format (-fmt='o') and plot color (-color=str(color))
    const StringMap keywords = plotKeywords(".", color, line_style);
    plt::errorbar(x, y, yerr, xerr, keywords);
  }

  // Scatter Type II
  void MPLPlot::scatter(const std::vector<double> & x,
	       const std::vector<double> & y,
	       const std::vector<double> & yerr,
	       const char color,
	       const std::string & line_style)
  {
    // Keyword Map - scatter format (-fmt='o') and plot color (-color=str(color))
    const StringMap keywords = plotKeywords(".",color,line_style);
    plt::errorbar(x, y, yerr, keywords);
  }

  // Scatter Type III
  void MPLPlot::scatter(const std::vector<double> & x,
	       const std::vector<double> & y,
	       const char color,
	       const std::string & line_style)
  {
    // Keyword Map - scatter format (-fmt='o') and plot color (-color=str(color))
    const StringMap keywords = plotKeywords(".",color,line_style);
    const double s=1.0; // The marker size in points**2
    plt::scatter(x, y, s, keywords);
  }

  // Line Type I
  void MPLPlot::linePlot(const std::vector<double> & x,
		const std::vector<double> & y,
		const char color,
		const std::string & line_style)
  {
    // Keyword Map - scatter format (-fmt='o') and plot color (-color=str(color))
    const StringMap keywords = plotKeywords("",color,line_style);
    plt::plot(x, y, keywords);
  }

  // void MPLPlot::logLog(const std::vector<double> & x,
  // 		       const std::vector<double> & y,
  // 		       char color)
  // {
  //   // Keyword Map - scatter format (-fmt='o') and plot color (-color=str(color))                                                                                           
  //   const StringMap keywords = plotKeywords("",color,"");
  //   plt::loglog(x,y,keywords);
  // }

  // void MPLPlot::semilogx(const std::vector<double> & x,
  // 			 const std::vector<double> & y,
  // 			 char color,
  // 			 const std::string & line_style)
  // {
  //   // Keyword Map - scatter format (-fmt='o') and plot color (-color=str(color))
  //   const StringMap keywords = plotKeywords("",color,line_style);
  //   plt::named_semilogx(plt_title, x, y, keywords);
  // }
  
  void MPLPlot::showPlot()
  {
    plt::show();
  }
  
  // Histogram I - Don't need Histograms unless otherwise indicated
  // void MPLPlot::histogram(const std::vector<long> &x,
  // 		 const std::vector<double> &y,
  // 		 const std::vector<double> &xwidth,
  // 		 int color,
  // 		 const std::string & line_style)
  // {
  //   plt::hist(y); // Do we need width?
  // }

  // Formatting/Utility function passthrough
  void MPLPlot::ylim(const double left, const double right)
  {
    plt::ylim(left,right);
  }

  void MPLPlot::xlim(const double left, const double right)
  {
    plt::xlim(left,right);
  }

  void MPLPlot::subplot(long nrows, long ncols, long plot_number)
  {
    plt::subplot(nrows, ncols, plot_number);
  }

  void MPLPlot::figure(long number)
  {
    plt::figure(number);
  }
} // namespace EasyPlot
