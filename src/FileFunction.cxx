/** 
 * @file FileFunction.cxx
 * @brief Implementation for the FileFunction Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/FileFunction.cxx,v 1.5 2009/01/17 21:42:12 jchiang Exp $
 */

#include <cmath>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "facilities/Util.h"

#include "st_facilities/Util.h"

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/FileFunction.h"

namespace Likelihood {

FileFunction::FileFunction(double Normalization) : m_filename("") {
   setMaxNumParams(1);

   addParam("Normalization", Normalization, true);

// Set FuncType and ArgType for use with CompositeFunction hierarchy.
   m_funcType = Addend;
   m_argType = "dArg";

   m_genericName = "FileFunction";
   m_normParName = "Normalization";
}

double FileFunction::value(optimizers::Arg & xarg) const {
   double x(std::log(dynamic_cast<optimizers::dArg &>(xarg).getValue()));
   double norm = m_parameter[0].getTrueValue();
   return norm*std::exp(st_facilities::Util::interpolate(m_x, m_y, x));
}

double FileFunction::
derivByParam(optimizers::Arg & xarg, const std::string & paramName) const {
   if (paramName != "Normalization") {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "FileFunction::derivByParam");
   }
   double x(std::log(dynamic_cast<optimizers::dArg &>(xarg).getValue()));
   double scale(m_parameter[0].getScale());
   return scale*std::exp(st_facilities::Util::interpolate(m_x, m_y, x));
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
   return flux;
}

void FileFunction::readFunction(const std::string & filename) {
   m_filename = filename;
   std::string input_file(filename);
   facilities::Util::expandEnvVar(&input_file);
   st_facilities::Util::file_ok(input_file);
   std::vector<std::string> lines;
   bool removeWindowsCRs(true);
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

double FileFunction::
integral(optimizers::Arg &, optimizers::Arg &) const {
   throw std::runtime_error("FileFunction::integral is disabled.");
}

} // namespace Likelihood
