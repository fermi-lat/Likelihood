/** @file SourceModel.h
 * @brief Declaration of SourceModel class
 * $Header:
 */

#ifndef SourceModel_h
#define SourceModel_h

#include <vector>
#include <string>

#include "Function.h"

namespace Likelihood {

/** 
 * @class SourceModel
 *
 * @brief Container for source models.
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/SourceModel.h,v 1.1.1.1 2003/01/30 23:23:03 burnett Exp $ */

class SourceModel : public Function {
    
public:
   
   SourceModel(){setMaxNumParams(0);};
   ~SourceModel();

   //! overloaded setParam to include function name checking
   void setParam(const Parameter param, const std::string fName);

   //! overloaded setParamValues to ensure synching with 
   //! m_function parameters
   void setParamValues(const std::vector<double> paramVec);

   Parameter* getParam(const std::string paramName, 
		       const std::string fName) const;

   //! add and delete sources by name
   void addSource(Function *func, const std::string fName);
   void deleteSource(const std::string fName);

   //! function access
   Function * getFunc(const std::string fName);
   int getNumSrcs() const {return m_functions.size();};
   std::vector<std::string> getSrcNames() const;
   std::vector<int> getNumSrcParams() const;

   virtual double value(const double) const;
   virtual double operator()(const double x) const {return value(x);};
   virtual std::vector<double> getDerivs(const double) const;

private:
   
   //! update the parameters in m_functions with those contained
   //! in m_parameters
   void m_syncParams();

   std::vector<Function *> m_functions;

   //! need to keep track of which parameters go with which function,
   //! so we set up a vector of iterators to Parameter

   std::vector<std::vector<Parameter>::iterator> m_paramIterators;

};

} // namespace Likelihood

#endif // SourceModel_h
