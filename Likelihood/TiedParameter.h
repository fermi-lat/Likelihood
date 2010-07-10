/**
 * @file TiedParameter.h
 * @brief Handle parameters that are "tied together" in the Xspec
 * sense.  This is a subclass of optimizers::Parameter.
 * 
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/TiedParameter.h,v 1.2 2010/07/08 23:10:36 jchiang Exp $
 */

#ifndef Likelihood_TiedParameter_h
#define Likelihood_TiedParameter_h

#include <string>
#include <vector>
#include <utility>

#include "optimizers/Parameter.h"

namespace Likelihood {

class LogLike;

class TiedParameter : public optimizers::Parameter {

public:

   TiedParameter();

   TiedParameter(const TiedParameter & other);

   virtual ~TiedParameter() throw();

   TiedParameter & operator=(const TiedParameter & rhs);
   
   void addParam(LogLike & like, size_t i);

   void removeParam(LogLike & like, size_t i);

   bool has_member(const LogLike & like, size_t i) const;

   virtual void setName(const std::string name);

   virtual void setValue(double value);

   virtual void setScale(double scale);

   virtual void setTrueValue(double trueValue);

   virtual void setBounds(double minValue, double maxValue);

   virtual void setFree(bool free);

   virtual void setAlwaysFixed(bool flag);

   virtual void setError(double error);

   typedef std::vector<std::pair<LogLike *, size_t> > ParVector_t;
   typedef ParVector_t::iterator ParVectorIterator_t;
   typedef ParVector_t::const_iterator ParVectorConstIterator_t;

   const ParVector_t & pars() const {
      return m_pars;
   }

   bool containsIndex(const LogLike & like, size_t par_index) const;

private:

   ParVector_t m_pars;

};

} // namespace Likelihood

#endif // Likelihood_TiedParameter_h

