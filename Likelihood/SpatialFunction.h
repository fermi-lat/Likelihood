/** 
 * @file SpatialFunction.h
 * @brief Declaration for the SpatialFunction Function class
 * @author M. Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialFunction.h,v 1.27 2015/03/21 05:38:03 jchiang Exp $
 *
 */

#ifndef Likelihood_SpatialFunction_h
#define Likelihood_SpatialFunction_h

#include <utility>

#include "optimizers/Function.h"

#include "Likelihood/MapBase.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/SkyDirFunction.h"

namespace astro {
   class SkyDir;
}

namespace Likelihood {

#ifndef SWIG
using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
using XERCES_CPP_NAMESPACE_QUALIFIER DOMNode;
#endif // SWIG

class Event;
class ResponseFunctions;

/** 
 * @class SpatialFunction
 *
 * @brief Base class for spatial functions.
 *
 */
    
class SpatialFunction : public optimizers::Function {

public:
  
   SpatialFunction(const std::string& name, unsigned nparam);

   SpatialFunction(const std::string& name, unsigned nparam, double ra, double dec);
                         
   SpatialFunction(const SpatialFunction &);

   SpatialFunction & operator=(const SpatialFunction &);

   virtual ~SpatialFunction();

   virtual double value(const astro::SkyDir &) const = 0;
   virtual double value(const astro::SkyDir &, double energy, const MeanPsf& psf) const = 0;
   virtual double value(double delta, double energy, const MeanPsf& psf) const = 0;

   const astro::SkyDir& dir() const;

   void setCenter(double ra, double dec);

   virtual void update();

   virtual void setParam(const std::string &paramName, double paramValue) {
     optimizers::Function::setParam(paramName, paramValue);
     update();
   }

#ifndef SWIG
   /// Set the Parameters from a Function DOM_Element.
   virtual void setParams(const DOMElement * elt);
#endif // SWIG

protected:

   virtual double value(const optimizers::Arg & dir) const = 0;

   virtual double derivByParamImp(const optimizers::Arg & dir,
                                  const std::string & parName) const = 0;

private:

   void init(double ra = 0, double dec = 0);
   astro::SkyDir m_dir;
};

} // namespace Likelihood

#endif // Likelihood_SpatialFunction_h
