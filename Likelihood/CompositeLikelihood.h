/**
 * @file CompositeLikelihood.h
 *
 * @brief log-likelihood comprised of a number of LogLike objects that
 * have the same spectral parameters.  Used for Dark Matter group
 * analysis of dwarf spheroidal satellites.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/CompositeLikelihood.h,v 1.7 2011/03/02 11:02:08 cohen Exp $
 */

#ifndef Likelihood_CompositeLikelihood_h
#define Likelihood_CompositeLikelihood_h

#include <map>
#include <string>
#include <vector>

#include "optimizers/Statistic.h"

#include "Likelihood/LogLike.h"

namespace Likelihood {

/*
 * @class CompositeLikelihood
 */

class CompositeLikelihood : public optimizers::Statistic {

public:

   CompositeLikelihood() : optimizers::Statistic() {}

   virtual ~CompositeLikelihood() throw() {}

   void addComponent(const std::string & srcName, LogLike & component);

   virtual double value() const;
   virtual void getFreeParams(std::vector<optimizers::Parameter> &params) const;
   virtual void setFreeParamValues(const std::vector<double> & values);
   virtual unsigned int getNumFreeParams() const;
   virtual void getFreeDerivs(std::vector<double> & derivs) const;
   virtual void fetchParamValues(std::vector<double> & values, 
                                 bool getFree) const;
   void syncParams();

   double NpredValue(const std::string &) const {return 0;}

protected:

   double value(optimizers::Arg&) const {return value();}

   double derivByParam(optimizers::Arg&, const std::string&) const {return 0;}

   optimizers::Function * clone() const {return 0;}

private:

   /// Using the LogLike pointer as the map key allows for the same
   /// source name to be used in different components.
   typedef std::map<LogLike *, std::string> ComponentMap_t;
   typedef ComponentMap_t::iterator ComponentIterator_t;
   typedef ComponentMap_t::const_iterator ComponentConstIterator_t;

   ComponentMap_t m_components;

   std::string m_normParName;
   std::string m_commonFuncName;

};

} // namespace Likelihood

#endif // Likelihood_CompositeLikelihood_h
