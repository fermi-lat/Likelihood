/**
 * @file SourceData.h
 * @brief Small class to encapsulate data expected to be read from
 * xml test files by the Likelihood test program.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef my_LikelihoodTests_SourceData_h
#define my_LikelihoodTests_SourceData_h

#include <cassert>

#include <map>
#include <string>
#include <vector>

#include "optimizers/Parameter.h"

using optimizers::Parameter;

class SourceData {

public:

   SourceData();

   ~SourceData() {}

   void getSrcNames(std::vector<std::string> &srcNames) {
      srcNames.clear();
      for (std::map<std::string, unsigned int>::iterator it = m_srcId.begin();
           it != m_srcId.end(); it++) {
         srcNames.push_back(it->first);
      }
   }

   const std::string & srcType(const std::string & srcName) {
      assert(m_srcId.count(srcName));
      return m_srcTypes[m_srcId[srcName]];
   }

   const std::string & spatialModel(const std::string & srcName) {
      assert(m_srcId.count(srcName));
      return m_spatialModels[m_srcId[srcName]];
   }

   const std::vector<std::string> & paramNames() const {
      return m_paramNames;
   }

   const std::vector<Parameter> & parameters(const std::string & srcName) {
      assert(m_srcId.count(srcName));
      return m_parameters[m_srcId[srcName]];
   }

   const Parameter & paramObject(const std::string & srcName,
                                 const std::string & paramName) {
      assert(m_srcId.count(srcName));
      const std::vector<Parameter> & params = parameters(srcName);
      for (unsigned int i = 0; i < params.size(); i++) {
         if (params[i].getName() == paramName) {
            return params[i];
         }
      }
      bool paramFound(false);
      assert(paramFound);
      return params[0];
   }

private:

   std::map<std::string, unsigned int> m_srcId;
   std::vector<std::string> m_srcTypes;
   std::vector<std::string> m_spatialModels;

   std::vector<std::string> m_paramNames;
   std::vector< std::vector<Parameter> > m_parameters;
   void setParameters();

};

#endif // my_LikelihoodTests_SourceData_h
