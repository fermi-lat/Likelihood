/**
 * @file SummedLikelihood.h
 *
 * @brief Statistic object that comprises LogLike objects that have
 * identical source models.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SummedLikelihood.h,v 1.5 2008/09/26 00:27:15 jchiang Exp $
 */

#ifndef Likelihood_SummedLikelihood_h
#define Likelihood_SummedLikelihood_h

#include <map>
#include <string>
#include <vector>

#include "optimizers/Statistic.h"

#include "Likelihood/LogLike.h"

namespace Likelihood {

/*
 * @class SummedLikelihood
 */

class SummedLikelihood : public optimizers::Statistic {

public:

   SummedLikelihood() : optimizers::Statistic() {}

   virtual ~SummedLikelihood() throw() {}

   void addComponent(LogLike & component);

   virtual double value() const;
   virtual void getFreeParams(std::vector<optimizers::Parameter> &params) const;
   virtual void setFreeParamValues(const std::vector<double> & values);
   virtual unsigned int getNumFreeParams() const;
   virtual void getFreeDerivs(std::vector<double> & derivs) const;

   void syncParams();

   double NpredValue(const std::string &) const {return 0;}

protected:

   double value(optimizers::Arg&) const {return value();}

   double derivByParam(optimizers::Arg&, const std::string&) const {return 0;}

   optimizers::Function * clone() const {return 0;}

private:

   typedef std::vector<LogLike *> ComponentVector_t;
   typedef ComponentVector_t::iterator ComponentIterator_t;
   typedef ComponentVector_t::const_iterator ComponentConstIterator_t;

   ComponentVector_t m_components;

   std::string m_normParName;
   std::string m_commonFuncName;

};

} // namespace Likelihood

#endif // Likelihood_SummedLikelihood_h
