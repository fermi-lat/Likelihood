/** 
 * @file SourceModel.h
 * @brief Declaration of SourceModel class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceModel.h,v 1.70 2012/06/14 02:01:24 jchiang Exp $
 */

#ifndef Likelihood_SourceModel_h
#define Likelihood_SourceModel_h

#include <map>
#include <vector>
#include <stdexcept>
#include <string>

#include "optimizers/Statistic.h"

#include "Likelihood/Pixel.h"
#include "Likelihood/Observation.h"
#include "Likelihood/Source.h"

namespace st_stream {
   class StreamFormatter;
}

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
 */

class SourceModel : public optimizers::Statistic {

public:
   
   SourceModel(const Observation & observation, bool verbose=false);

   SourceModel(const SourceModel & rhs);

   virtual ~SourceModel();

   /// setParam method to include function and source name checking
   virtual void setParam(const optimizers::Parameter &param, 
                         const std::string &funcName,
                         const std::string &srcName);

   /// group parameter access (note name mangling for inheritance 
   /// from Function)
   virtual std::vector<double>::const_iterator setParamValues_(
      std::vector<double>::const_iterator);
   virtual std::vector<double>::const_iterator setFreeParamValues_(
      std::vector<double>::const_iterator);

   virtual double getParamValue(const std::string &paramName, 
                                const std::string &funcName,
                                const std::string &srcName) const {
      return getParam(paramName, funcName, srcName).getValue();
   }
   
   virtual optimizers::Parameter getParam(const std::string &paramName, 
                                          const std::string &funcName,
                                          const std::string &srcName) const;

   virtual void setParamTrueValue(const std::string &paramName,
                                  const std::string &funcName,
                                  const std::string &srcName,
                                  double paramValue);

   virtual void setParams(std::vector<optimizers::Parameter> &params) {
      setParams_(params, false);
   }

   virtual void setFreeParams(std::vector<optimizers::Parameter> &params) {
      setParams_(params, true);
   }

   /// This needs to be re-implemented, delegating to the base class
   /// method, since all member functions with the same name get
   /// hidden by a local declaration, even if the signatures differ.
   virtual void getFreeDerivs(optimizers::Arg &x, 
                              std::vector<double> &derivs) const {
      Function::getFreeDerivs(x, derivs);
   }

   /// Add a source.
   virtual void addSource(Source *src, bool fromClone=true);

   /// Delete a source by name and return a copy.
   virtual Source * deleteSource(const std::string &srcName);

   /// delete all the sources
   void deleteAllSources();

   /// return a Source pointer by name
   Source * getSource(const std::string &srcName);

   /// @return reference to the desired source
   const Source & source(const std::string & srcName) const;

   /// @return reference to the Source map.
   const std::map<std::string, Source *> & sources() const {
      return m_sources;
   }

   unsigned int getNumSrcs() const {
      return m_sources.size();
   }

   void getSrcNames(std::vector<std::string> &) const;

   bool hasSrcNamed(const std::string & srcName) const;

   virtual double value(optimizers::Arg &x) const;

   /// Create the source model by reading an XML file.
   virtual void readXml(std::string xmlFile,
                        optimizers::FunctionFactory & funcFactory,
                        bool requireExposure=true,
                        bool addPointSources=true,
                        bool loadMaps=true);

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
   virtual CountsMap * createCountsMap(const CountsMap & dataMap) const;

   virtual CountsMap * createCountsMap() const {
      throw std::runtime_error("SourceModel::createCountsMap needs to be "
                             + std::string("reimplemented in this subclass."));
   }

   virtual const Observation & observation() const {
      return m_observation;
   }

   /// method to sync the m_parameter vector with those of the 
   /// m_sources' Functions
   virtual void syncParams();

   const std::vector<optimizers::Parameter> & parameters() const {
      return m_parameter;
   }

   std::vector<optimizers::Parameter> & parameters() {
      return m_parameter;
   }

   /// Default methods to set energy bounds of the analysis.
   /// These should be re-implement in derived classes 
   virtual void set_ebounds(double emin, double emax) {
      throw std::runtime_error("SourceModel::set_ebounds not implemented.");
   }

   virtual void unset_ebounds() {
      throw std::runtime_error("SourceModel::unset_ebounds not implemented.");
   }

protected:

   const Observation & m_observation;

   virtual SourceModel * clone() const {
      return new SourceModel(*this);
   }

   std::map<std::string, Source *> m_sources;

   void fetchDerivs(optimizers::Arg &x, std::vector<double> &derivs, 
                    bool getFree) const;

   void setParams_(std::vector<optimizers::Parameter> &, bool);

   std::vector<Source *> m_freeSrcs;

   bool m_useNewImp;

   void findFreeSrcs();

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

   /// disable this since parameters may no longer have unique names
   double derivByParamImp(optimizers::Arg &, const std::string &) const {
      return 0;
   }

private:

   bool m_verbose;

   st_stream::StreamFormatter * m_formatter;

   void computeModelMap(const std::vector<Pixel> & pixels,
                        const std::vector<double> & energies,
                        std::vector<float> & modelMap) const;

};

} // namespace Likelihood

#endif // Likelihood_SourceModel_h
