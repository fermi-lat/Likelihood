/** 
 * @file SourceModel.h
 * @brief Declaration of SourceModel class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceModel.h,v 1.24 2003/09/29 15:32:04 jchiang Exp $
 */

#ifndef Likelihood_SourceModel_h
#define Likelihood_SourceModel_h

#include <vector>
#include <string>

#include "optimizers/Function.h"
#include "Likelihood/Source.h"

namespace optimizers {
   class FunctionFactory;
}

namespace Likelihood {

/** 
 * @class SourceModel
 *
 * @brief This base class provides derived classes with methods that
 * allow for Sources to be added and deleted and for the Parameters
 * and derivatives describing those Sources to be accessed.
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceModel.h,v 1.24 2003/09/29 15:32:04 jchiang Exp $ 
 */

class SourceModel : public optimizers::Function {
    
public:
   
   SourceModel() {
      setMaxNumParams(0); 
      m_genericName = "SourceModel";
      s_refCount++;
   }
   SourceModel(const SourceModel &rhs) : optimizers::Function(rhs) 
      {s_refCount++;}

   virtual ~SourceModel();

   //! setParam method to include function and source name checking
   virtual void setParam(const optimizers::Parameter &param, 
                         const std::string &funcName,
                         const std::string &srcName) 
      throw(optimizers::Exception);

   //! group parameter access (note name mangling for inheritance 
   //! from Function)
   virtual std::vector<double>::const_iterator setParamValues_(
      std::vector<double>::const_iterator);
   virtual std::vector<double>::const_iterator setFreeParamValues_(
      std::vector<double>::const_iterator);

   virtual double getParamValue(const std::string &paramName, 
                                const std::string &funcName,
                                const std::string &srcName) const
      throw(optimizers::Exception, optimizers::ParameterNotFound) {
      return getParam(paramName, funcName, srcName).getValue();
   }
   
   virtual optimizers::Parameter getParam(const std::string &paramName, 
                              const std::string &funcName,
                              const std::string &srcName) const
      throw(optimizers::Exception, optimizers::ParameterNotFound);

   virtual void setParamBounds(const std::string &paramName,
                               const std::string &funcName,
                               const std::string &srcName,
                               double lower, double upper)
      throw(optimizers::ParameterNotFound, optimizers::OutOfBounds);

   virtual void setParamScale(const std::string &paramName,
                              const std::string &funcName,
                              const std::string &srcName,
                              double scale) 
      throw(optimizers::ParameterNotFound);

   virtual void setParamTrueValue(const std::string &paramName,
                                  const std::string &funcName,
                                  const std::string &srcName,
                                  double paramValue) 
      throw(optimizers::ParameterNotFound, optimizers::OutOfBounds);

   virtual void setParams(std::vector<optimizers::Parameter> &params) 
      throw(optimizers::Exception, optimizers::ParameterNotFound) 
      {setParams_(params, false);}

   virtual void setFreeParams(std::vector<optimizers::Parameter> &params) 
      throw(optimizers::Exception, optimizers::ParameterNotFound) 
      {setParams_(params, true);}

   //! add and delete sources by name
   void addSource(Source *src);
   void deleteSource(const std::string &srcName) throw(optimizers::Exception);

   //! delete all the sources
   void deleteAllSources();

   //! return a Source pointer by name
   Source * getSource(const std::string &srcName);

   unsigned int getNumSrcs() const {return s_sources.size();}
   void getSrcNames(std::vector<std::string> &) const;

   // this is a bit convoluted, but necessary for some derived classes 
   // (e.g., logLike_gauss)
   double evaluate_at(optimizers::Arg &) const;
   virtual double value(optimizers::Arg &x) const {return evaluate_at(x);}

   /// Create the source model by reading an XML file.
   void readXml(const std::string &xmlFile,
                optimizers::FunctionFactory &funcFactory);

   /// Write an XML file for the current source model.
   void writeXml(const std::string &xmlFile,
                 const std::string &functionLibrary);

protected:

   static int s_refCount;

   static std::vector<Source *> s_sources;

   //! method to sync the m_parameter vector with those of the 
   //! s_sources' Functions
   void syncParams();

   //! disable this since parameters may no longer have unique names
   double derivByParam(optimizers::Arg &, const std::string &) 
      const {return 0;}

   void fetchDerivs(optimizers::Arg &x, std::vector<double> &derivs, 
                    bool getFree) const;

   void setParams_(std::vector<optimizers::Parameter> &, bool)
      throw(optimizers::Exception, optimizers::ParameterNotFound);

};

} // namespace Likelihood

#endif // Likelihood_SourceModel_h
