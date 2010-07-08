/**
 * @file TiedParameter.h
 * @brief Handle parameters that are "tied together" in the Xspec
 * sense.  This is a subclass of optimizers::Parameter.
 * 
 * @author J. Chiang
 * 
 * $Header$
 */

#ifndef Likelihood_TiedParameter_h
#define Likelihood_TiedParameter_h

#include <string>
#include <vector>

#include "optimizers/Parameter.h"

namespace Likelihood {

class TiedParameter : public optimizers::Parameter {

public:

   TiedParameter(const std::string & name, double value,
                 double minValue, double maxValue,
                 bool isFree, double error);

   TiedParameter(const TiedParameter & other);

   virtual ~TiedParameter() throw();

   TiedParameter & operator=(const TiedParameter & rhs);
   
   void addParam(optimizers::Parameter & par);

   void removeParam(optimizers::Parameter & par);

   bool has_member(const optimizers::Parameter & par) const;

   virtual void setName(const std::string name);

   virtual void setValue(double value);

   virtual void setScale(double scale);

   virtual void setTrueValue(double trueValue);

   virtual void setBounds(double minValue, double maxValue);

   virtual void setFree(bool free);

   virtual void setAlwaysFixed(bool flag);

   virtual void setError(double error);

private:

   typedef std::vector<optimizers::Parameter *> ParVector_t;
   ParVector_t m_pars;

};

} // namespace Likelihood

#endif // Likelihood_TiedParameter_h

