/**
 * @file CompositeLikelihood.h
 *
 * @brief log-likelihood comprised of a number of LogLike objects that
 * have the same spectral parameters.  Used for Dark Matter group
 * analysis of dwarf spheroidal satellites.
 * 
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_CompositeLikelihood_h
#define Likelihood_CompositeLikelihood_h

#include <map>
#include <string>

#include "optimizers/Statistic.h"

namespace Likelihood {

/*
 * @class CompositeLikelihood
 */

class CompositeLikelihood : public optimizers::Statistic {

public:

   CompositeLikelihood();

   virtual ~CompositeLikelihood() throw();

   void addComponent(const std::string & srcName, LogLike & component);

   virtual double value() const;
   virtual void getFreeParams(std::vector<Parameter> & params) const;
   virtual void setFreeParamValues(std::vector<double> & values);
   virtual unsigned int getNumFreeParams() const;
   virtual void getFreeDerivs(std::vector<double> & derivs) const;

   void syncParams();

   double NpredValue(const std::string & srcName) const;

private:

   typedef std::map<std::string, LogLike *> ComponentMap_t;
   typedef ComponentMap_t::iterator ComponentIterator_t;
   typedef ComponentMap_t::const_iterator ComponentConstIterator_t;

   ComponentMap_t m_components;

};

} // namespace Likelihood

#endif // Likelihood_CompositeLikelihood_h
