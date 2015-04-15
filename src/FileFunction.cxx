/** 
 * @file FileFunction.cxx
 * @brief Implementation for the FileFunction Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/FileFunction.cxx,v 1.11 2015/03/21 05:38:03 jchiang Exp $
 */

#include <cmath>

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "facilities/Util.h"

#include "st_facilities/Util.h"

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/FileFunction.h"

namespace Likelihood {

FileFunction::FileFunction(double Normalization) 
   : optimizers::Function("FileFunction", 1, "Normalization"),
     m_filename("") {
   addParam("Normalization", Normalization, true);
}

double FileFunction::value(const optimizers::Arg & xarg) const {
   double x(std::log(dynamic_cast<const optimizers::dArg &>(xarg).getValue()));
   double norm(m_parameter[0].getTrueValue());
   if (x == m_x.front()) {
      return norm*std::exp(m_y.front());
   } else if (x == m_x.back()) {
      return norm*std::exp(m_y.back());
   }
   return norm*interpolateFlux(x);
}

double FileFunction::
derivByParamImp(const optimizers::Arg & xarg,
                const std::string & paramName) const {
   if (paramName != "Normalization") {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "FileFunction::derivByParam");
   }
   double x(std::log(dynamic_cast<const optimizers::dArg &>(xarg).getValue()));
   double scale(m_parameter[0].getScale());
   return scale*interpolateFlux(x);
}

double FileFunction::interpolateFlux(double logEnergy) const {
   double flux(0);
   try {
      flux = std::exp(st_facilities::Util::interpolate(m_x, m_y, logEnergy));
   } catch (std::range_error & eObj) {
      std::ostringstream message;
      message << "Requested energy, " << std::exp(logEnergy) 
              << ", lies outside the range of the input file, "
              << std::exp(m_x.front()) << ", " << std::exp(m_x.back());
      throw std::range_error(message.str());
   }
   if (flux != flux) {
      std::runtime_error("FileFunction::interpolateFlux: Nan encountered");
   }
   return flux;
}

void FileFunction::readFunction(const std::string & filename) {
   m_filename = filename;
   std::string input_file(filename);
   facilities::Util::expandEnvVar(&input_file);
   st_facilities::Util::file_ok(input_file);
   std::vector<std::string> lines;
//   bool removeWindowsCRs(true);
   bool removeWindowsCRs(false);
   st_facilities::Util::readLines(input_file, lines, "#", removeWindowsCRs);
   m_x.clear();
   m_y.clear();
   for (size_t i = 0; i < lines.size(); i++) {
      std::vector<std::string> tokens;
      facilities::Util::stringTokenize(lines.at(i), " \t", tokens);
      double xval(facilities::Util::stringToDouble(tokens.at(0)));
      double yval(facilities::Util::stringToDouble(tokens.at(1)));
      m_x.push_back(std::log(xval));
      m_y.push_back(std::log(yval));
   }
}   

void FileFunction::setSpectrum(const std::vector<double> & energy,
                               const std::vector<double> & dnde) {
   if (energy.size() != m_x.size() || dnde.size() != m_y.size()) {
      throw std::runtime_error("FileFunction::setSpectrum: inconsistent "
                               "array sizes for input spectra.");
   }
   for (size_t k(0); k < m_x.size(); k++) {
      m_x[k] = std::log(energy[k]);
      m_y[k] = std::log(dnde[k]);
   }
}

const std::vector<double> & FileFunction::log_energy() const {
   return m_x;
}

const std::vector<double> & FileFunction::log_dnde() const {
   return m_y;
}

} // namespace Likelihood
