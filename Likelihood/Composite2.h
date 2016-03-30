/**
 * @file Composite2.h
 *
 * @brief Composite log-likelihood that allows for tying together of
 * arbitrary combinations of parameters.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Composite2.h,v 1.9 2015/03/21 05:38:02 jchiang Exp $
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

   virtual ~Composite2() throw();

   void addComponent(LogLike & component);

   void tieParameters(const TiedParameter::ParVector_t & pars);

   virtual double value() const;
   virtual void getFreeParams(std::vector<optimizers::Parameter> &params) const;
   virtual void fetchParamValues(std::vector<double> &values, bool getFree) const;
   virtual void setFreeParamValues(const std::vector<double> & values);
   virtual unsigned int getNumFreeParams() const;
   virtual void getFreeDerivs(std::vector<double> & derivs) const;

   /// Set the errors for the free parameters.
   void setErrors(const std::vector<double> & errors);

   /// Find the parameter index (in the optimizer context, i.e., the
   /// index number among the free parameters) in the composite model
   /// for the specified parameter.  This is for use with the
   /// minosError function in the python interface.
   ///
   /// @return The desired index.  -1 will be returned if the
   /// parameter is fixed, or if the par_index does not exist in the 
   /// specified LogLike object.
   ///
   /// @param like The LogLike object with the desired parameter.
   /// @param par_index The "absolute" index (i.e., including the
   /// fixed pars) of the parameter in the LogLike object.
   ///
   int findIndex(const LogLike & like, size_t par_index) const;

   void syncParams();

   double NpredValue(const std::string &, bool weighted=false) const {return 0;}

   TiedParameter & getTiedParam(const LogLike & like, size_t i);
   void setTiedParamValue(const LogLike & like, size_t i, double value);

protected:

   double value(const optimizers::Arg &) const {
      return value();
   }

   double derivByParamImp(const optimizers::Arg &, const std::string&) const {
      return 0;
   }

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
