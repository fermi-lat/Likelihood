/** @file SourceModel.h
 * @brief Declaration of SourceModel class
 * $Header:
 */

#ifndef SourceModel_h
#define SourceModel_h

#include <vector>
#include <string>

#include "Function.h"
#include "Source.h"

namespace Likelihood {

/** 
 * @class SourceModel
 *
 * @brief Container for source models.
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/Likelihood/SourceModel.h,v 1.1 2003/02/19 01:34:33 jchiang Exp $ */

class SourceModel : public Function {
    
public:
   
   SourceModel(){setMaxNumParams(0);};
   SourceModel(const SourceModel &rhs);
   virtual ~SourceModel();

   //! overloaded setParam to include function name checking
   void setParam(const Parameter &param, const std::string &fName);

   //! overloaded setParamValues to ensure synching with 
   //! m_function parameters
   void setParamValues(const std::vector<double> &paramVec);

   Parameter* getParam(const std::string &paramName, 
		       const std::string &fName) const;

   //! add and delete sources by name
   void addSource(Function *func, const std::string &fName);
   void addSource(Source *src, const std::string &srcName);
   void deleteSource(const std::string &fName);

   //! function access
   Function * getFunc(const std::string &fName);
   unsigned int getNumSrcs() const {return m_functions.size();};
   void getSrcNames(std::vector<std::string> &) const;
   void getNumSrcParams(std::vector<int> &) const;

// this is a bit convoluted, but necessary for derived classes (e.g., Statistic)
   double evaluate_at(double) const;
   virtual double value(double x) const {return evaluate_at(x);};
   virtual double operator()(double x) const {return value(x);};
   virtual void getDerivs(double, std::vector<double>&) const;

protected:

   std::vector<Function *> m_functions;
   std::vector<Source *> m_sources;

private:
   
   //! update the parameters in m_functions with those contained
   //! in m_parameters
   void m_syncParams();

   //! need to keep track of which parameters go with which function,
   //! so we set up a vector of iterators to Parameter

   std::vector<std::vector<Parameter>::iterator> m_paramIterators;

};

} // namespace Likelihood

#endif // SourceModel_h
