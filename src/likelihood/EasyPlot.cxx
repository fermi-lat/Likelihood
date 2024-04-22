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
      

  // Scatter Type I
  void MPLPlot::scatter(const std::vector<double> & x,
	       const std::vector<double> & y,
	       const std::vector<double> & xerr,
	       const std::vector<double> & yerr,
	       const char color,
	       const std::string & line_style)
  {
    // Keyword Map - scatter format (-fmt='o') and plot color (-color=str(color))
    std::string c(1,color);
    const std::map<std::string, std::string> & keywords = {
      {"fmt", "o"},
      {"color", c},
      {"linestyle", line_style}
    };
    
    plt::errorbar(x, y, yerr, xerr, keywords);
    //plt::show();  
  }

  // Scatter Type II
  void MPLPlot::scatter(const std::vector<double> & x,
	       const std::vector<double> & y,
	       const std::vector<double> & yerr,
	       const char color,
	       const std::string & line_style)
  {
    // Keyword Map - scatter format (-fmt='o') and plot color (-color=str(color))
    std::string c(1,color);
    std::string xscale = "log";
    std::string yscale = "log";
    const std::map<std::string, std::string> & keywords = {
      {"fmt", "o"},
      {"color", c},
      {"linestyle", line_style}
    };
    plt::errorbar(x, y, yerr, keywords);
    plt::xscale(xscale);
    plt::yscale(yscale);
    //plt::show();
  }

  // Scatter Type III
  void MPLPlot::scatter(const std::vector<double> & x,
	       const std::vector<double> & y,
	       const char color,
	       const std::string & line_style)
  {
    // Keyword Map - scatter format (-fmt='o') and plot color (-color=str(color))
    std::string c(1,color);
    const std::map<std::string, std::string> & keywords = {
      {"fmt", "o"},
      {"color", c},
      {"linestyle", line_style}
      //{"xscale", "log"}
    };
    const double s=1.0; // The marker size in points**2
    plt::scatter(x, y, s, keywords);
    //plt::show();
  }

  // Line Type I
  void MPLPlot::linePlot(const std::vector<double> & x,
		const std::vector<double> & y,
		char color,
		const std::string & line_style)
  {
    // Keyword Map - scatter format (-fmt='o') and plot color (-color=str(color))
    std::string c(1,color);
    const std::map<std::string, std::string> & keywords = {
      {"color", c},
      {"linestyle", line_style}
    };
    plt::plot(x, y, keywords);
    //plt::show();
  }

  void MPLPlot::logLog(const std::vector<double> & x,
		       const std::vector<double> & y,
		       char color)
  {
    std::string name = "plotname";
    // Keyword Map - scatter format (-fmt='o') and plot color (-color=str(color))                                                                                           
    std::string c(1, color);
    const std::map<std::string, std::string> & keywords = {
      {"color", c}
    };
    plt::named_loglog(name,x,y,keywords);
    //plt::show();
  }

  void MPLPlot::semilogx(const std::vector<double> & x,
			 const std::vector<double> & y,
			 char color,
			 const std::string & line_style,
			 const std::string & plt_title)
  {
    std::string c(1,color);
    const std::map<std::string, std::string> & keywords = {
      {"marker", "o"},
      {"color", c},
      {"linestyle", line_style}
    };
    plt::named_semilogx(plt_title, x, y, keywords);
  }
  
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
