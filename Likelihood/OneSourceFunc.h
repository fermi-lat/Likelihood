/**
 * @file OneSourceFunc.h
 * @brief Extended likelihood function for one source
 *
 * @author P. Nolan

 * $Header:
 */

#include "Event.h"
#include "Source.h"
#include <vector>
#include "optimizers/Function.h"
#include "optimizers/Arg.h"
#include "optimizers/Exception.h"
#include "optimizers/ParameterNotFound.h"

namespace Likelihood {
  
  class OneSourceFunc: public optimizers::Function {
    
  public:

    OneSourceFunc(Source * src,  // The Source of interest
		  std::vector<Event>& events, 
		  std::vector<double> * weights = NULL);

    virtual double value(optimizers::Arg& arg) const;
    virtual double derivByParam(optimizers::Arg& x, 
				const std::string& paramName) const;
    virtual std::vector<double>::const_iterator 
      setFreeParamValues_(std::vector<double>::const_iterator it);
    virtual std::vector<double>::const_iterator
      setParamValues_(std::vector<double>::const_iterator);
    virtual Function *clone() const {return new OneSourceFunc(*this);}
    virtual void setParams(std::vector<optimizers::Parameter>&)
      throw(optimizers::Exception, optimizers::ParameterNotFound);

  protected:

    void syncParams(void);
    
  private:
    
    Source * m_src;
    std::vector<Event>& m_events;
    std::vector<double> * m_weights;
    
  };  // class OneSourceFunc
} // namespace Likelihood
