/** 
 * @file FileFunction.cxx
 * @brief Implementation for the FileFunction Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/FileFunction.cxx,v 1.5 2005/10/19 22:13:54 jchiang Exp $
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
}

double FileFunction::value(optimizers::Arg & xarg) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();

   double norm = m_parameter[0].getTrueValue();

   return norm*st_facilities::Util::interpolate(m_x, m_y, x);
}

double FileFunction::
derivByParam(optimizers::Arg & xarg, const std::string & paramName) const {
   if (paramName != "Normalization") {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "FileFunction::derivByParam");
   }

   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();
   
   return st_facilities::Util::interpolate(m_x, m_y, x);
}

void FileFunction::readFunction(const std::string & filename) {
   m_filename = filename;
   facilities::Util::expandEnvVar(&m_filename);
   st_facilities::Util::file_ok(m_filename);
   std::vector<std::string> lines;
   bool removeWindowsCRs(true);
   st_facilities::Util::readLines(m_filename, lines, "#", removeWindowsCRs);
   m_x.clear();
   m_y.clear();
   for (size_t i = 0; i < lines.size(); i++) {
      std::vector<std::string> tokens;
      facilities::Util::stringTokenize(lines.at(i), " \t", tokens);
      m_x.push_back(facilities::Util::stringToDouble(tokens.at(0)));
      m_y.push_back(facilities::Util::stringToDouble(tokens.at(1)));
   }
}   

double FileFunction::
integral(optimizers::Arg &, optimizers::Arg &) const {
   throw std::runtime_error("FileFunction::integral is disabled.");
}

} // namespace Likelihood
