/** @file SourceModel.h
 * @brief Declaration of SourceModel class
 * $Header:
 */

#ifndef SourceModel_h
#define SourceModel_h

#include <vector>
#include <string>

#include "Likelihood/Function.h"
#include "Likelihood/Source.h"

namespace Likelihood {

/** 
 * @class SourceModel
 *
 * @brief This base class provides derived classes with methods that
 * allow for Sources to be added and deleted and for the Parameters
 * and derivatives describing those Sources to be accessed.
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceModel.h,v 1.4 2003/03/04 17:45:32 jchiang Exp $ */

class SourceModel : public Function {
    
public:
   
   SourceModel(){setMaxNumParams(0);};
   SourceModel(const SourceModel &rhs);
   virtual ~SourceModel(){};

   //! setParam method to include function and source name checking
   void setParam(const Parameter &param, const std::string &funcName,
                 const std::string &srcName);

   //! group parameter access (note name mangling for inheritance 
   //! from Function)
   virtual std::vector<double>::const_iterator setParamValues_(
      std::vector<double>::const_iterator);
   virtual std::vector<double>::const_iterator setFreeParamValues_(
      std::vector<double>::const_iterator);

   Parameter* getParam(const std::string &paramName, 
                       const std::string &funcName,
                       const std::string &srcName) const;

   //! add and delete sources by name
   void addSource(Source *src);
   void deleteSource(const std::string &srcName);

   //! delete all the sources
   void deleteAllSources() {
      m_sources.clear();
      m_parameter.clear();
   }

   unsigned int getNumSrcs() const {return m_sources.size();}
   void getSrcNames(std::vector<std::string> &) const;

   // this is a bit convoluted, but necessary for some derived classes 
   // (e.g., logLike_gauss)
   double evaluate_at(Arg &) const;
   virtual double value(Arg &x) const {return evaluate_at(x);};

protected:

   static std::vector<Source *> m_sources;

   //! method to sync the m_parameter vector with those of the 
   //! m_sources' Functions
   void m_syncParams();

   //! disable this since parameters may no longer have unique names
   double derivByParam(Arg &, const std::string &) const {return 0;}

   void fetchDerivs(Arg &x, std::vector<double> &derivs, bool getFree) const;

};

} // namespace Likelihood

#endif // SourceModel_h
