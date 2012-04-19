/** 
 * @file SourceModel.cxx
 * @brief SourceModel class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/SourceModel.cxx,v 1.95 2012/02/07 00:24:28 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "xmlBase/Dom.h"
#include "xmlBase/XmlParser.h"

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "astro/SkyDir.h"

#include "optimizers/Arg.h"
#include "optimizers/FunctionFactory.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/CountsMap.h"
#include "Likelihood/Exception.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SourceFactory.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/FluxBuilder.h"
#include "Likelihood/SourceModelBuilder.h"
#include "Likelihood/SourceModel.h"
#include "Likelihood/XmlParser.h"

namespace Likelihood {

using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

SourceModel::SourceModel(const Observation & observation, bool verbose) 
   : m_observation(observation), m_useNewImp(true), m_verbose(verbose), 
     m_formatter(new st_stream::StreamFormatter("SourceModel", "", 2)) {
   setMaxNumParams(0); 
   m_genericName = "SourceModel";
   char * useOldImp(::getenv("USE_OLD_LOGLIKE"));
   if (useOldImp) {
      m_useNewImp = false;
   }
}

SourceModel::SourceModel(const SourceModel &rhs) : optimizers::Statistic(rhs),
   m_observation(rhs.m_observation), m_verbose(rhs.m_verbose) {
   delete m_formatter;
   m_formatter = new st_stream::StreamFormatter("SourceModel", "", 2);
}

SourceModel::~SourceModel() {
   std::map<std::string, Source *>::iterator it = m_sources.begin();
   for ( ; it != m_sources.end(); ++it) {
      delete it->second;
   }
   m_sources.clear();
   delete m_formatter;
}

void SourceModel::setParam(const optimizers::Parameter & param, 
                           const std::string & funcName,
                           const std::string & srcName) {
   if (m_sources.count(srcName)) {
      // Source::FuncMap srcFuncs = (*m_sources[srcName]).getSrcFuncs();
      // if (srcFuncs.count(funcName)) {
      //    srcFuncs[funcName]->setParam(param);
      //    syncParams();
      //    return;
      // }
      if ("Spectrum" == funcName) {
         m_sources[srcName]->spectrum().setParam(param);
         syncParams();
         return;
      }
   }
   std::string errorMessage = "SourceModel::setParam:\n Function " 
      + funcName + " for Source " + srcName + " was not found.";
   throw optimizers::Exception(errorMessage);
}
 
std::vector<double>::const_iterator 
SourceModel::setParamValues_(std::vector<double>::const_iterator it) {
   std::map<std::string, Source *>::iterator srcIt = m_sources.begin();
   for ( ; srcIt != m_sources.end(); ++srcIt) {
      // Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      // Source::FuncMap::iterator func_it = srcFuncs.begin();
      // for ( ; func_it != srcFuncs.end(); func_it++) {
      //    it = (*func_it).second->setParamValues_(it);
      // }
      it = srcIt->second->spectrum().setParamValues_(it);
   }
   syncParams();
   return it;
}

std::vector<double>::const_iterator 
SourceModel::setFreeParamValues_(std::vector<double>::const_iterator it) {
   std::map<std::string, Source *>::iterator srcIt = m_sources.begin();
   for (size_t i = 0 ; srcIt != m_sources.end(); ++srcIt, i++) {
      // Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      // Source::FuncMap::iterator func_it = srcFuncs.begin();
      // for ( ; func_it != srcFuncs.end(); func_it++) {
      //    it = func_it->second->setFreeParamValues_(it);
      // }
      it = srcIt->second->spectrum().setFreeParamValues_(it);
   }
   syncParams();
   return it;
}

optimizers::Parameter SourceModel::getParam(const std::string &paramName,
                                            const std::string &funcName,
                                            const std::string &srcName) const {
   if (m_sources.count(srcName)) {
      std::vector<optimizers::Parameter> params;
      const Source * my_source = m_sources.find(srcName)->second;
      // const Source::FuncMap & srcFuncs = my_source->getSrcFuncs();
      // if (srcFuncs.count(funcName)) {    //check for funcName
      //    try {
      //       const optimizers::Function * my_function =
      //          srcFuncs.find(funcName)->second;
      //       my_function->getParams(params);
      //    } catch (optimizers::Exception &eObj) {
      //       m_formatter->err() << eObj.what() << std::endl;
      //       throw;
      //    }
      if ("Spectrum" == funcName) {
         try {
            my_source->spectrum().getParams(params);
         } catch (optimizers::Exception & eObj) {
            m_formatter->err() << eObj.what() << std::endl;
            throw;
         }
         for (unsigned int j = 0; j < params.size(); j++) {
            if (paramName == params[j].getName()) {
               return params[j];
            }
         }
         throw optimizers::ParameterNotFound(paramName, funcName, 
                                             "SourceModel::getParam");
      }
      std::string errorMessage = "SourceModel::getParam:\n Function "
         + funcName + " was not found in Source " 
         + srcName + ".";
      throw optimizers::Exception(errorMessage);
   }
   std::string errorMessage = "SourceModel::getParam:\nSource "
      + srcName + " was not found.";
   throw optimizers::Exception(errorMessage);
}

void SourceModel::setParamTrueValue(const std::string &paramName,
                                    const std::string &funcName,
                                    const std::string &srcName,
                                    double paramValue) {
   optimizers::Parameter my_param = getParam(paramName, funcName, srcName);
   my_param.setTrueValue(paramValue);
   setParam(my_param, funcName, srcName);
   syncParams();
}

void SourceModel::setParams_(std::vector<optimizers::Parameter> &params, 
                             bool setFree) {
// Ensure the number of Parameters matches.
   unsigned int numParams;
   if (setFree) {
      numParams = getNumFreeParams();
   } else {
      numParams = getNumParams();
   }
   if (params.size() != numParams) {
      std::string errorMessage = std::string("SourceModel::setParams:\n") 
         + "Inconsistent number of Parameters.";
      throw optimizers::Exception(errorMessage);
   }

// Assume ordering of Parameters in params matches that given by the
// ordering of the Sources and their Functions.
   int k = 0;  // params' index
   std::map<std::string, Source *>::iterator srcIt = m_sources.begin();
   for ( ; srcIt != m_sources.end(); ++srcIt) {
      // Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      // Source::FuncMap::iterator func_it = srcFuncs.begin();
      // for (; func_it != srcFuncs.end(); func_it++) {
      //    unsigned int numParams;
      //    if (setFree) {
      //       numParams = func_it->second->getNumFreeParams();
      //    } else { 
      //       numParams = func_it->second->getNumParams();
      //    }
      //    for (unsigned int j = 0; j < numParams; j++) {
      //       func_it->second->setParam(params[k]);
      //       k++;
      //    }
      // }
      size_t numParams;
      if (setFree) {
         numParams = srcIt->second->spectrum().getNumFreeParams();
      } else {
         numParams = srcIt->second->spectrum().getNumParams();
      }
      for (size_t j(0); j < numParams; j++) {
         srcIt->second->spectrum().setParam(params[k]);
         k++;
      }
   }
   syncParams();
}

void SourceModel::addSource(Source *src, bool fromClone) {
   if (!m_sources.count(src->getName())) {
      m_sources[src->getName()] = fromClone ? src->clone() : src;
      m_sources[src->getName()]->setObservation(&m_observation);
      syncParams();
   } else {
      throw Exception("Likelihood::SourceModel:\nSource named " 
                      + src->getName() + " alread exists.");
   }
}
 
Source * SourceModel::deleteSource(const std::string &srcName) {
   std::map<std::string, Source *>::iterator it = m_sources.find(srcName);
   if (it != m_sources.end()) {
      Source * mySource = it->second;
      m_sources.erase(it);
      syncParams();
      return mySource;
   }
   std::string errorMessage = "SourceModel::deleteSource:\n" 
      + srcName + " was not found.";
   throw optimizers::Exception(errorMessage);
}

void SourceModel::deleteAllSources() {
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source * source = deleteSource(srcNames[i]);
      delete source;
   }
   m_parameter.clear();
}

Source * SourceModel::getSource(const std::string & srcName) {
   if (m_sources.count(srcName)) {
      return m_sources[srcName];
   }
   return 0;
}

const Source & SourceModel::source(const std::string & srcName) const {
   std::map<std::string, Source *>::const_iterator my_src =
      m_sources.find(srcName);
   if (my_src == m_sources.end()) {
      throw std::runtime_error("SourceModel::source: Source " + 
                               srcName + " not found.");
   }
   return *(my_src->second);
}

void SourceModel::getSrcNames(std::vector<std::string> &names) const {
   names.clear();
   std::map<std::string, Source *>::const_iterator it = m_sources.begin();
   for ( ; it != m_sources.end(); ++it) {
      names.push_back(it->first);
   }
}

double SourceModel::value(optimizers::Arg &x) const {
   double my_val = 0.;
   std::map<std::string, Source *>::const_iterator srcIt = m_sources.begin();
   for ( ; srcIt != m_sources.end(); ++srcIt) {
      // Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      // Source::FuncMap::iterator func_it = srcFuncs.begin();
      // for (; func_it != srcFuncs.end(); func_it++) {
      //    my_val += (*func_it).second->value(x);
      // }
      my_val += srcIt->second->spectrum().value(x);
   }
   return my_val;
}

void SourceModel::syncParams() { // remake parameter vector from scratch 
   m_parameter.clear();

   std::map<std::string, Source *>::iterator srcIt = m_sources.begin();
   for ( ; srcIt != m_sources.end(); ++srcIt) {
      // Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      // Source::FuncMap::iterator func_it = srcFuncs.begin();
      // for (; func_it != srcFuncs.end(); func_it++) {
      //    std::vector<optimizers::Parameter> params;
      //    func_it->second->getParams(params);
      //    for (size_t ip = 0; ip < params.size(); ip++) {
      //       m_parameter.push_back(params[ip]);
      //    }
      // }
      std::vector<optimizers::Parameter> params;
      srcIt->second->spectrum().getParams(params);
      for (size_t ip(0); ip < params.size(); ip++) {
         m_parameter.push_back(params.at(ip));
      }
   }
   if (m_useNewImp) {
      findFreeSrcs();
   }
}

void SourceModel::fetchDerivs(optimizers::Arg &x,
                              std::vector<double> &derivs, 
                              bool getFree) const {
   derivs.clear();

   std::map<std::string, Source *>::const_iterator srcIt = m_sources.begin();
   for ( ; srcIt != m_sources.end(); ++srcIt) {
      // Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      // Source::FuncMap::iterator func_it = srcFuncs.begin();
      // for (; func_it != srcFuncs.end(); func_it++) {
      //    std::vector<double> my_derivs;
      //    if (getFree) {
      //       (*func_it).second->getFreeDerivs(x, my_derivs);
      //    } else {
      //       (*func_it).second->getDerivs(x, my_derivs);
      //    }
      //    for (unsigned int j = 0; j < my_derivs.size(); j++) 
      //       derivs.push_back(my_derivs[j]);
      // }
      std::vector<double> my_derivs;
      if (getFree) {
         srcIt->second->spectrum().getFreeDerivs(x, my_derivs);
      } else {
         srcIt->second->spectrum().getDerivs(x, my_derivs);
      }
      for (size_t j(0); j < my_derivs.size(); j++) {
         derivs.push_back(my_derivs.at(j));
      }
   }
}

void SourceModel::readXml(std::string xmlFile,
                          optimizers::FunctionFactory &funcFactory,
                          bool requireExposure,
                          bool addPointSources,
                          bool loadMaps) {

// Expand any environment variables in the xmlFile name.
   facilities::Util::expandEnvVar(&xmlFile);

// Create a SourceFactory to read in the xml file.
   SourceFactory srcFactory(m_observation);
   try {
      srcFactory.readXml(xmlFile, funcFactory, requireExposure,
                         addPointSources, loadMaps);
   } catch (xmlBase::DomException & eObj) {
      m_formatter->err() << eObj.what() << std::endl;
      std::ostringstream message;
      message << "\nError reading in the xml model file.\n"
              << "Please check that you are using the correct xml "
              << "format for this tool." << std::endl;
      throw Exception(message.str());
   }

// Loop over the sources that are now contained in srcFactory and add
// each one to the source model (removing it from the srcFactory to avoid
// making a copy).
   std::vector<std::string> srcNames;
   srcFactory.fetchSrcNames(srcNames);

   std::vector<std::string>::iterator nameIt = srcNames.begin();
   for ( ; nameIt != srcNames.end(); nameIt++) {
      Source * src = srcFactory.releaseSource(*nameIt);
      if (m_verbose) {
         m_formatter->info() << "adding source " << *nameIt << std::endl;
      }
      addSource(src, false);
   }
   syncParams();
}

void SourceModel::reReadXml(std::string xmlFile) {

   facilities::Util::expandEnvVar(&xmlFile);

   xmlBase::XmlParser * parser = XmlParser_instance();

   DOMDocument * doc = parser->parse(xmlFile.c_str());

   if (doc == 0) {
      std::string errorMessage = "SourceFactory::readXml:\nInput xml file, "
         + xmlFile + ", not parsed successfully.";
      throw Exception(errorMessage);
   }

   DOMElement * source_library = doc->getDocumentElement();
   if (!xmlBase::Dom::checkTagName(source_library, "source_library")) {
      throw Exception("SourceModel::reReadXml:\nsource_library not found");
   }

// Loop through source DOMElements and Source objects in parallel.
   std::vector<DOMElement *> srcs;
   xmlBase::Dom::getChildrenByTagName(source_library, "source", srcs);

   for (unsigned int j = 0; j < srcs.size(); j++) {
      std::string srcName = xmlBase::Dom::getAttribute(srcs[j], "name");
      Source * my_source(0);
      if (m_sources.count(srcName)) {
         my_source = m_sources[srcName];
      } else {
         continue;
      }
      std::string srcType = xmlBase::Dom::getAttribute(srcs[j], "type");
// Get spectrum and spatialModel elements
      std::vector<DOMElement *> child;
      xmlBase::Dom::getChildrenByTagName(srcs[j], "spectrum", child);
      DOMElement * spectrum = child[0];

      xmlBase::Dom::getChildrenByTagName(srcs[j], "spatialModel", child);
      DOMElement * spatialModel = child[0];

      my_source->getSrcFuncs()["Spectrum"]->setParams(spectrum);
      if (srcType == "PointSource") {
// Extract (RA, Dec) from the parameter elements and use
// PointSource::setDir(...) method to ensure that exposure gets
// recalculated.
         double ra(0), dec(0);
         std::vector<DOMElement *> params;
         xmlBase::Dom::getChildrenByTagName(spatialModel, "parameter", params);

         std::vector<DOMElement *>::const_iterator paramIt = params.begin();
         for ( ; paramIt != params.end(); paramIt++) {
            std::string name = xmlBase::Dom::getAttribute(*paramIt, "name");
            if (name == "RA")
               ra = std::atof(xmlBase::Dom::getAttribute(*paramIt,
                                                         "value").c_str());
                  
            if (name == "DEC")
               dec = std::atof(xmlBase::Dom::getAttribute(*paramIt,
                                                          "value").c_str());
                  
         }
         astro::SkyDir newDir(ra, dec);
         double tol = 1e-4;
         PointSource * ptSrc = dynamic_cast<PointSource *>(my_source);
         if (newDir.difference(ptSrc->getDir() ) > tol) {
// Reset the direction, re-computing the PointSource exposure.
            ptSrc->setDir(newDir);
         }
      } else if (srcType == "DiffuseSource") {
         my_source->getSrcFuncs()["SpatialDist"]->setParams(spatialModel);
      } else {
         std::ostringstream message;
         message << "SourceModel::reReadXml: "
                 << "Unknown Source type: " << srcType;
         throw std::runtime_error(message.str());
      }
   }
   syncParams();
   delete doc;
}

void SourceModel::writeXml(std::string xmlFile,
                           const std::string & functionLibrary,
                           const std::string & srcLibTitle) {

   SourceModelBuilder builder(functionLibrary, srcLibTitle);

   std::map<std::string, Source *>::iterator srcIt = m_sources.begin();
   for ( ; srcIt != m_sources.end(); srcIt++) {
      builder.addSource(*(srcIt->second));
   }

   builder.write(xmlFile);
}

void SourceModel::write_fluxXml(std::string xmlFile) {

   std::pair<double, double> ebounds = m_observation.roiCuts().getEnergyCuts();
   FluxBuilder builder(ebounds.first, ebounds.second);

   std::map<std::string, Source *>::iterator srcIt = m_sources.begin();
   for ( ; srcIt != m_sources.end(); srcIt++) {
      builder.addSource(*(srcIt->second));
   }

   builder.write(xmlFile);
}

bool SourceModel::hasSrcNamed(const std::string & srcName) const {
   std::vector<std::string> names;
   getSrcNames(names);
   std::vector<std::string>::iterator it 
      = std::find(names.begin(), names.end(), srcName);
   if (it == names.end()) {
      return false;
   }
   return true;
}

CountsMap * SourceModel::createCountsMap(const CountsMap & dataMap) const {
   const std::vector<Pixel> & pixels(dataMap.pixels());

   std::vector<double> energies;
   dataMap.getAxisVector(2, energies);

   std::vector<float> map;
   computeModelMap(pixels, energies, map);

   CountsMap * modelMap = new CountsMap(dataMap);
   modelMap->setImage(map);
   return modelMap;
}

void SourceModel::computeModelMap(const std::vector<Pixel> & pixels,
                                  const std::vector<double> & energies,
                                  std::vector<float> & modelMap) const {
   modelMap.clear();
   modelMap.reserve(pixels.size()*(energies.size()-1));
   for (unsigned int k = 0; k < energies.size()-1; k++) {
      for (unsigned int j = 0; j < pixels.size(); j++) {
         modelMap.push_back(
            pixels[j].modelCounts(energies[k], energies[k+1],
                                  *(const_cast<SourceModel *>(this))));
      }
   }
}

void SourceModel::findFreeSrcs() {
   m_freeSrcs.clear();
   std::map<std::string, Source *>::const_iterator src(m_sources.begin());
   for ( ; src != m_sources.end(); ++src) {
      // Source::FuncMap & srcFuncs(src->second->getSrcFuncs());
      // Source::FuncMap::const_iterator func(srcFuncs.begin());
      // bool haveFreePars(false);
      // for ( ; func != srcFuncs.end(); ++func) {
      //    if (func->second->getNumFreeParams() > 0) {
      //       haveFreePars = true;
      //    }
      // }
      bool haveFreePars = (src->second->spectrum().getNumFreeParams() > 0);
      if (haveFreePars) {
         m_freeSrcs.push_back(src->second);
      }
   }
}

} // namespace Likelihood
