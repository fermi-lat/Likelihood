/**
 * @file DiffRespNames.cxx
 * @brief Class to keep track of diffuse response column names.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/DiffRespNames.cxx,v 1.1 2008/11/29 15:53:49 jchiang Exp $
 */

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "facilities/Util.h"

#include "Likelihood/DiffRespNames.h"

namespace Likelihood {

DiffRespNames::DiffRespNames(const std::string & rootName) 
   : m_rootName(rootName) {}

size_t DiffRespNames::size() const {
   return m_colnames.size();
}

const std::vector<std::string> & DiffRespNames::colnames() const {
   return m_colnames;
}

const std::string & DiffRespNames::operator[](size_t indx) const {
   return m_colnames.at(indx);
}

const std::string & DiffRespNames::operator[](const std::string & key) const {
   if (key.substr(0, m_rootName.size()) != m_rootName) {
      throw DiffRespNameError("invalid root name for "
                              "diffuse response column");
   }
   size_t indx = facilities::Util::atoi(key.substr(m_rootName.size()));
   return operator[](indx);
}

void DiffRespNames::addColumn(const std::string & diffRspName) {
   if (!std::count(m_colnames.begin(), m_colnames.end(), diffRspName)) {
      m_colnames.push_back(diffRspName);
   }
}

std::string DiffRespNames::key(const std::string & colname) const {
   std::vector<std::string>::const_iterator it = 
      std::find(m_colnames.begin(), m_colnames.end(), colname);
   if (it == m_colnames.end()) {
      throw DiffRespNameError("diffuse response named "
                              + colname + " not found");
   }
   size_t indx = it - m_colnames.begin();
   std::ostringstream key_string;
   key_string << m_rootName << indx;
   return key_string.str();
}

} // namespace Likelihood
