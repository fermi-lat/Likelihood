/** @file logSrcModel.h
 * @brief Declaration of logSrcModel class
 * $Header:
 */

#ifndef logSrcModel_h
#define logSrcModel_h

#include <vector>
#include <string>

#include "../Likelihood/SourceModel.h"
#include "../Likelihood/EventArg.h"

namespace Likelihood {

/** 
 * @class logSrcModel
 *
 * @brief A SourceModel subclass that returns as its Function
 * value(Arg &), the quantity log(evaluate_at(EventArg &)).  Thus it
 * forms an additive Function component to the "data sum" in the
 * log-likelihood statistic.  Since it is a Function, the Statistic
 * fetchDerivs() method (which is inherited from SourceModel) can be
 * applied transparently.
 *
 * @authors J. Chiang
 *    
 * $Header: */

class logSrcModel : public SourceModel {
    
public:
   
   logSrcModel(){setMaxNumParams(0);};
   logSrcModel(const logSrcModel &rhs);
   virtual ~logSrcModel(){};

   virtual double value(Arg &xarg) const;
   virtual double derivByParam(Arg&, std::string &) const {return 0;}

protected:

   virtual void fetchDerivs(Arg &x, std::vector<double> &derivs, 
                            bool getFree) const;

};

} // namespace Likelihood

#endif // logSrcModel_h
