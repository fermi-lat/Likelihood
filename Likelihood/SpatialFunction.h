/** 
 * @file SpatialFunction.h
 * @brief Declaration for the SpatialFunction Function class
 * @author M. Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialFunction.h,v 1.1 2015/10/17 17:19:14 mdwood Exp $
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
 * @class ResponseFunctor
 *
 * @brief Base class functor for evaluating IRF response in spatial
 * convolutions.
 *
 */

class ResponseFunctor {

public:

   virtual ~ResponseFunctor() { }
   virtual double operator()(double energy, double separation) const = 0;
};

/** 
 * @class BinnedResponseFunctor
 *
 * @brief Functor for spatial response.  This class is used to
 * evaluate the PSF when generating source maps for binned analysis.
 *
 */

class BinnedResponseFunctor : public ResponseFunctor {

public:
   BinnedResponseFunctor(const Likelihood::MeanPsf & psf):
   ResponseFunctor(), m_psf(psf) { }

   double operator()(double energy, double separation) const {
       return m_psf(energy,separation);
   }

private:
  const Likelihood::MeanPsf& m_psf;
};

/** 
 * @class UnbinnedResponseFunctor
 *
 * @brief Functor for unbinned (diffuse) response.  This class is used
 * to evaluate the diffuse response for unbinned analysis.
 *
 */

class UnbinnedResponseFunctor : public ResponseFunctor {

public:

   UnbinnedResponseFunctor(const Likelihood::ResponseFunctions & respFuncs, 
			   double theta, double phi, int type):
   ResponseFunctor(), m_respFuncs(respFuncs), 
     m_theta(theta), m_phi(phi), m_type(type) { }

   double operator()(double energy, double separation) const {
     double psf_val = m_respFuncs.psf(m_type).value(separation,energy,m_theta,m_phi);
     double aeff_val = m_respFuncs.aeff(m_type).value(energy,m_theta,m_phi);
     return psf_val*aeff_val;
   }
	     
private:
   const Likelihood::ResponseFunctions& m_respFuncs;
   double                   m_theta;
   double                   m_phi;
   int                      m_type;
};

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

   virtual double spatialResponse(const astro::SkyDir &, double energy, const MeanPsf& psf) const = 0;
   virtual double spatialResponse(double delta, double energy, const MeanPsf& psf) const = 0;

   virtual double diffuseResponse(const Event & evt,
				  const ResponseFunctions & respFuncs) const;

   virtual double diffuseResponse(const ResponseFunctor& fn, double energy,
				  double separation) const = 0;

   virtual double getDiffRespLimits(const astro::SkyDir &, 
				    double & mumin, double & mumax,
				    double & phimin, double & phimax) const = 0;

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
