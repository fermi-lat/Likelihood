/**
 * @file LikelihoodBase.h
 * @brief Base class for likelihood functions.
 * @author J. Chiang <jchiang@slac.stanford.edu>
 * 
 * $Header$
 */

#ifndef Likelihood_LikelihoodBase_h
#define Likelihood_LikelihoodBase_h

#include <string>
#include <vector>

#include "optimizers/Statistic.h"



namespace Likelihood {

   class Source;

/**
 * @class LikelihoodBase
 */

class LikelihoodBase : public optimizers::Statistic {

   virtual ~LikelihoodBase();
   
   virtual double value() const = 0;
   
   virtual void getFreeDerivs(std::vector<double> & derivs) const = 0;

   /// Member functions also needed by pyLikelihood classes.

   /// Used by AnalysisBase
   virtual void syncParams();
   virtual void addSource(Source * src, bool fromClone);
   virtual void setFreeParamValues(const std::vector<double> & params);
   virtual Source * getSource(const std::string & srcName);
   virtual void getFreeParamValues(std::vector<double> & params) const;
   virtual Source * deleteSource(const std::string & srcName);
   virtual void saveCurrentFit();
   virtual void restoreBestFit();
   virtual void syncSrcParams(const std::string & srcName);
   virtual void getSrcNames(std::vector<std::string> & names) const;
   void readXml(const std::string & xmlModel);
   virtual void writeXml(const std::string & xmlfile) const;
   
   virtual double NpredValue(const std::string & srcName) const = 0;

protected:

   LikelihoodBase();

   LikelihoodBase(const LikelihoodBase & other) 
      : optimizers::Statistic(other) {}

   // These are pure virtual in optimizers::Function, so we
   // provide inaccessible concrete implementations.
   double value(optimizers::Arg & ) const {
      return value();
   }

   double derivByParam(optimizers::Arg &, const std::string &) const {
      return 0;
   }

   Function * clone() const {
      return 0;
   }

};

} //namespace Likelihood

#endif // Likelihood_LikelihoodBase_h

