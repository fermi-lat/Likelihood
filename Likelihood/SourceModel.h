/** @file SourceModel.h
 * @brief Declaration of SourceModel class
 * $Header:
 */

#ifndef SourceModel_h
#define SourceModel_h

#include <vector>
#include <string>

#include "../Likelihood/Function.h"
#include "../Likelihood/Source.h"

namespace Likelihood {

/** 
 * @class SourceModel
 *
 * @brief Container for source models.
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/Likelihood/SourceModel.h,v 1.2 2003/02/23 22:28:59 jchiang Exp $ */

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

   unsigned int getNumSrcs() const {return m_sources.size();}
   void getSrcNames(std::vector<std::string> &) const;

   // this is a bit convoluted, but necessary for derived classes 
   // (e.g., Statistic)
   double evaluate_at(double) const;
   virtual double value(double x) const {return evaluate_at(x);};
   virtual double operator()(double x) const {return value(x);};

   virtual void getDerivs(double, std::vector<double>&) const;
   virtual void getFreeDerivs(double, std::vector<double>&) const;

protected:

   std::vector<Source *> m_sources;

   //! method to sync the m_parameter vector with those of the 
   //! m_sources' Functions
   void m_syncParams();

};

} // namespace Likelihood

#endif // SourceModel_h
