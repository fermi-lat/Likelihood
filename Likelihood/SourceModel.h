/** 
 * @file SourceModel.h
 * @brief Declaration of SourceModel class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceModel.h,v 1.45 2004/09/16 04:38:59 jchiang Exp $
 */

#ifndef Likelihood_SourceModel_h
#define Likelihood_SourceModel_h

#include <map>
#include <vector>
#include <string>

#include "map_tools/Exposure.h"

#include "optimizers/Statistic.h"

#include "Likelihood/Pixel.h"
#include "Likelihood/Source.h"

namespace optimizers {
   class FunctionFactory;
}

namespace Likelihood {

   class CountsMap;

/** 
 * @class SourceModel
 *
 * @brief This base class provides derived classes with methods that
 * allow for Sources to be added and deleted and for the Parameters
 * and derivatives describing those Sources to be accessed.
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceModel.h,v 1.45 2004/09/16 04:38:59 jchiang Exp $ 
 */

class SourceModel : public optimizers::Statistic {

public:
   
   SourceModel(bool verbose=false) : m_verbose(verbose) {
      setMaxNumParams(0); 
      m_genericName = "SourceModel";
      s_refCount++;
   }

   SourceModel(const SourceModel &rhs) : optimizers::Statistic(rhs) 
      {s_refCount++;}

   virtual ~SourceModel();

   /// setParam method to include function and source name checking
   virtual void setParam(const optimizers::Parameter &param, 
                         const std::string &funcName,
                         const std::string &srcName) 
      throw(optimizers::Exception);

   /// group parameter access (note name mangling for inheritance 
   /// from Function)
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

   /// This needs to be re-implemented, delegating to the base class
   /// method, since all member functions with the same name get
   /// hidden by a local declaration, even if the signatures differ.
   virtual void getFreeDerivs(optimizers::Arg &x, 
                              std::vector<double> &derivs) const {
      Function::getFreeDerivs(x, derivs);
   }

   /// Add a source.
   void addSource(Source *src);

   /// Delete a source by name and return a copy.
   Source * deleteSource(const std::string &srcName) 
      throw(optimizers::Exception);

   /// delete all the sources
   void deleteAllSources();

   /// return a Source pointer by name
   Source * getSource(const std::string &srcName);

   /// @return reference to the Source map.
   const std::map<std::string, Source *> & sources() {return s_sources;}

   unsigned int getNumSrcs() const {return s_sources.size();}
   void getSrcNames(std::vector<std::string> &) const;
   bool hasSrcNamed(const std::string & srcName) const;

   virtual double value(optimizers::Arg &x) const;

   /// Create the source model by reading an XML file.
   virtual void readXml(std::string xmlFile,
                        optimizers::FunctionFactory &funcFactory,
                        bool requireExposure=true);

   /// Re-read an XML file, updating only the Parameters in the
   /// source model.
   virtual void reReadXml(std::string xmlFile);

   /// Write an XML file for the current source model.
   virtual void writeXml(std::string xmlFile,
                         const std::string &functionLibrary="",
                         const std::string &srcLibTitle="source library");

   /// Write a flux-style xml file for the current source model.
   virtual void write_fluxXml(std::string xmlFile);

   /// Create a counts map based on the current model.
   CountsMap * createCountsMap(const CountsMap & dataMap) const;

protected:

   static int s_refCount;

   static std::map<std::string, Source *> s_sources;

   /// method to sync the m_parameter vector with those of the 
   /// s_sources' Functions
   void syncParams();

   /// disable this since parameters may no longer have unique names
   double derivByParam(optimizers::Arg &, const std::string &) 
      const {return 0;}

   void fetchDerivs(optimizers::Arg &x, std::vector<double> &derivs, 
                    bool getFree) const;

   void setParams_(std::vector<optimizers::Parameter> &, bool)
      throw(optimizers::Exception, optimizers::ParameterNotFound);

   /// Although these member functions are required by being a
   /// Statistic subclass, they are not needed for any practical use
   /// of SourceModel objects themselves, so we implement them here in
   /// the protected area.  SourceModel subclasses that need them,
   /// e.g., LogLike, will have to re-implement them.
   virtual double value() const {
      return 0;
   }

   virtual void getFreeDerivs(std::vector<double> &derivs) const {
      derivs.clear();
   }

protected:

   static void getPixels(const CountsMap & countsMap,
                         std::vector<Pixel> & pixels);

   void computeModelMap(const std::vector<Pixel> & pixels,
                        const std::vector<double> & energies,
                        std::vector<double> & modelMap) const;

private:

   bool m_verbose;

   static double computeSolidAngle(std::vector<double>::const_iterator lon,
                                   std::vector<double>::const_iterator lat,
                                   const astro::SkyProj & proj);

   static void getPixels(const CountsMap & countsMap, 
                         std::vector<astro::SkyDir> & pixelDirs,
                         std::vector<double> & pixelSolidAngles);

};

} // namespace Likelihood

#endif // Likelihood_SourceModel_h
