/** @file logSrcModel.h
 * @brief Declaration of logSrcModel class
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef logSrcModel_h
#define logSrcModel_h

#include <vector>
#include <string>

#include "Likelihood/SourceModel.h"
#include "Likelihood/EventArg.h"

namespace Likelihood {

/** 
 * @class logSrcModel
 *
 * @brief A SourceModel subclass that returns as its Function
 * value(Arg &), the log of the sum of source flux densities for Arg
 * cast as an EventArg.  It forms an additive Function component to
 * the "data sum" in the log-likelihood statistic.  Since it is a
 * Function, its fetchDerivs() method can be applied transparently
 * using the get[Free]Derivs() methods inherited from the Function
 * base class.
 *
 * @authors J. Chiang
 *    
 * $Header$
 */

class logSrcModel : public SourceModel {
    
public:
   
   logSrcModel(){setMaxNumParams(0);}
   logSrcModel(const logSrcModel &rhs);
   virtual ~logSrcModel(){}

   double value(Arg &xarg) const;
   double derivByParam(Arg&, std::string &) const {return 0;}

   // would be nice if this wasn't necessary...
   void mySyncParams() {syncParams();}

protected:

   void fetchDerivs(Arg &x, std::vector<double> &derivs, 
                    bool getFree) const;

};

} // namespace Likelihood

#endif // logSrcModel_h
