/** 
 * @file logSrcModel.h
 * @brief Declaration of logSrcModel class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/logSrcModel.h,v 1.9 2003/08/06 20:52:04 jchiang Exp $
 */

#ifndef Likelihood_logSrcModel_h
#define Likelihood_logSrcModel_h

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/logSrcModel.h,v 1.9 2003/08/06 20:52:04 jchiang Exp $
 */

class logSrcModel : public SourceModel {
    
public:
   
   logSrcModel() {setMaxNumParams(0); m_genericName = "logSrcModel";}
   virtual ~logSrcModel() {}

   double value(optimizers::Arg &xarg) const;
   double derivByParam(optimizers::Arg&, std::string &) const {return 0;}

   // would be nice if this wasn't necessary...
   void mySyncParams() {syncParams();}

protected:

   void fetchDerivs(optimizers::Arg &x, std::vector<double> &derivs, 
                    bool getFree) const;

};

} // namespace Likelihood

#endif // Likelihood_logSrcModel_h
