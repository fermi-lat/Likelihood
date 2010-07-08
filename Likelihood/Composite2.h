/**
 * @file Composite2.h
 *
 * @brief Composite log-likelihood that allows for tying together of
 * arbitrary combinations of parameters.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Composite2.h,v 1.6 2010/05/17 21:17:49 jchiang Exp $
 */

#ifndef Likelihood_Composite2_h
#define Likelihood_Composite2_h

#include <map>
#include <string>
#include <vector>

#include "optimizers/Statistic.h"

#include "Likelihood/LogLike.h"
#include "Likelihood/TiedParameter.h"

namespace Likelihood {

/*
 * @class Composite2
 */

class Composite2 : public optimizers::Statistic {

public:

   Composite2() : optimizers::Statistic() {}

   virtual ~Composite2() throw() {}

   void addComponent(LogLike & component);

   void tieParameters(const TiedParameter::ParVector_t & pars);

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

   /// Using the LogLike pointer as the map key and associate
   /// the list of tied parameter indices for that component.
   typedef std::map<LogLike *, std::vector<size_t> > ComponentMap_t;
   typedef ComponentMap_t::iterator ComponentIterator_t;
   typedef ComponentMap_t::const_iterator ComponentConstIterator_t;

   ComponentMap_t m_components;

   std::vector<TiedParameter *> m_tiedPars;

};

} // namespace Likelihood

#endif // Likelihood_Composite2_h
