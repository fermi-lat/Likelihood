/** 
 * @file SourceModel.cxx
 * @brief SourceModel class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceModel.cxx,v 1.63 2004/12/22 16:57:34 jchiang Exp $
 */

#include <cassert>
#include <cmath>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "fitsio.h"

#include "xmlBase/Dom.h"
#include "xmlBase/XmlParser.h"

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "optimizers/Arg.h"
#include "optimizers/FunctionFactory.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/CountsMap.h"
#include "Likelihood/Exception.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SourceFactory.h"
#include "Likelihood/TrapQuad.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/FluxBuilder.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SourceModelBuilder.h"
#include "Likelihood/SourceModel.h"

#include "Verbosity.h"

namespace {
   void fitsReportError(FILE *stream, int status) {
      fits_report_error(stream, status);
      if (status != 0) {
         throw std::runtime_error("writeExposureFile: cfitsio error.");
      }
   }
} // unnamed namespace

namespace Likelihood {

XERCES_CPP_NAMESPACE_USE

int SourceModel::s_refCount = 0;
std::map<std::string, Source *> SourceModel::s_sources;

SourceModel::~SourceModel() {
   s_refCount--;
   if (s_refCount == 0) {
      std::map<std::string, Source *>::iterator it = s_sources.begin();
      for (; it != s_sources.end(); it++)
         delete (*it).second;
   }
   s_sources.clear();
}

void SourceModel::setParam(const optimizers::Parameter &param, 
                           const std::string &funcName,
                           const std::string &srcName) {
   if (s_sources.count(srcName)) {
      Source::FuncMap srcFuncs = (*s_sources[srcName]).getSrcFuncs();
      if (srcFuncs.count(funcName)) {
//          srcFuncs[funcName]->setParam(param.getName(), 
//                                       param.getValue(),
//                                       param.isFree());
         srcFuncs[funcName]->setParam(param);
         syncParams();
         return;
      }
   }
   std::string errorMessage = "SourceModel::setParam:\n Function " 
      + funcName + " for Source " + srcName + " was not found.";
   throw optimizers::Exception(errorMessage);
}
 
std::vector<double>::const_iterator SourceModel::setParamValues_(
   std::vector<double>::const_iterator it) { 
   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) 
         it = (*func_it).second->setParamValues_(it);
   }
   syncParams();
   return it;
}

std::vector<double>::const_iterator SourceModel::setFreeParamValues_(
   std::vector<double>::const_iterator it) { 
   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) 
         it = (*func_it).second->setFreeParamValues_(it);
   }
   syncParams();
   return it;
}

optimizers::Parameter SourceModel::getParam(const std::string &paramName,
                                            const std::string &funcName,
                                            const std::string &srcName) const {
   if (s_sources.count(srcName)) {
      std::vector<optimizers::Parameter> params;
      Source::FuncMap srcFuncs = s_sources[srcName]->getSrcFuncs();
      if (srcFuncs.count(funcName)) {    //check for funcName
         try {
            srcFuncs[funcName]->getParams(params);
         } catch (optimizers::Exception &eObj) {
            std::cerr << eObj.what() << std::endl;
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
      std::string errorMessage = "SourceModel::getParam:\n Function"
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
   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         unsigned int numParams;
         if (setFree) {
            numParams = func_it->second->getNumFreeParams();
         } else { 
            numParams = func_it->second->getNumParams();
         }
         for (unsigned int j = 0; j < numParams; j++) {
            func_it->second->setParam(params[k]);
            k++;
         }
      }
   }
   syncParams();
}

void SourceModel::addSource(Source *src) {
   if (!s_sources.count(src->getName())) {
      s_sources[src->getName()] = src->clone();
      syncParams();
   } else {
      throw Exception("Likelihood::SourceModel:\nSource named " 
                      + src->getName() + " alread exists.");
   }
}
 
Source * SourceModel::deleteSource(const std::string &srcName) {
   std::map<std::string, Source *>::iterator it = s_sources.find(srcName);
   if (it != s_sources.end()) {
      Source * mySource = it->second;
      s_sources.erase(it);
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

Source * SourceModel::getSource(const std::string &srcName) {
   if (s_sources.count(srcName)) {
      return s_sources[srcName];
   }
   return 0;
}

void SourceModel::getSrcNames(std::vector<std::string> &names) const {
   names.clear();
   std::map<std::string, Source *>::iterator it = s_sources.begin();
   for ( ; it != s_sources.end(); ++it) {
      names.push_back(it->first);
   }
}

double SourceModel::value(optimizers::Arg &x) const {
   double my_val = 0.;
   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         my_val += (*func_it).second->value(x);
      }
   }
   return my_val;
}

void SourceModel::syncParams() { // remake parameter vector from scratch 
   m_parameter.clear();

   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         std::vector<optimizers::Parameter> params;
         (*func_it).second->getParams(params);
         for (unsigned int ip = 0; ip < params.size(); ip++)
            m_parameter.push_back(params[ip]);
      }
   }
}

void SourceModel::fetchDerivs(optimizers::Arg &x,
                              std::vector<double> &derivs, 
                              bool getFree) const {
   derivs.clear();

   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         std::vector<double> my_derivs;
         if (getFree) {
            (*func_it).second->getFreeDerivs(x, my_derivs);
         } else {
            (*func_it).second->getDerivs(x, my_derivs);
         }
         for (unsigned int j = 0; j < my_derivs.size(); j++) 
            derivs.push_back(my_derivs[j]);
      }
   }
}

void SourceModel::readXml(std::string xmlFile,
                          optimizers::FunctionFactory &funcFactory,
                          bool requireExposure) {

// Expand any environment variables in the xmlFile name.
   facilities::Util::expandEnvVar(&xmlFile);

// Create a SourceFactory to read in the xml file.
   SourceFactory srcFactory;
   try {
      srcFactory.readXml(xmlFile, funcFactory, requireExposure);
   } catch (xmlBase::DomException & eObj) {
      if (print_output()) {
         std::cout << "SourceModel::readXml:\n DomException: " 
                   << eObj.what() << std::endl;
      }
   }

// Loop over the sources that are now contained in srcFactory and add
// each one to the source model.
   std::vector<std::string> srcNames;
   srcFactory.fetchSrcNames(srcNames);

   std::vector<std::string>::iterator nameIt = srcNames.begin();
   for ( ; nameIt != srcNames.end(); nameIt++) {
      Source *src = srcFactory.create(*nameIt);
      if (print_output() && m_verbose) {
         std::cout << "adding source " << *nameIt << std::endl;
      }
      addSource(src);
   }
   syncParams();
}

void SourceModel::reReadXml(std::string xmlFile) {

   facilities::Util::expandEnvVar(&xmlFile);

   xmlBase::XmlParser *parser = new xmlBase::XmlParser();

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
      if (s_sources.count(srcName)) {
         my_source = s_sources[srcName];
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
               ra = 
                  std::atof(xmlBase::Dom::getAttribute(*paramIt, "value").c_str());
            if (name == "DEC")
               dec = 
                  std::atof(xmlBase::Dom::getAttribute(*paramIt, "value").c_str());
         }
         astro::SkyDir newDir(ra, dec);
         double tol = 1e-4;
         if (newDir.difference(
                dynamic_cast<PointSource *>(my_source)->getDir() ) > tol) {
// Reset the direction, re-computing the PointSource exposure.
            my_source->setDir(newDir);
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
   delete parser;
}

void SourceModel::writeXml(std::string xmlFile,
                           const std::string &functionLibrary,
                           const std::string &srcLibTitle) {

   SourceModelBuilder builder(functionLibrary, srcLibTitle);

   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); srcIt++) {
      builder.addSource(*(srcIt->second));
   }

   builder.write(xmlFile);
}

void SourceModel::write_fluxXml(std::string xmlFile) {

   FluxBuilder builder;

   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); srcIt++) {
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

   if (ExposureCube::instance() == 0) {
      std::runtime_error("SourceModel::createCountsMap:\n"
                         + std::string("Exposure cube not available."));
   }

   std::vector<Pixel> pixels;
   dataMap.getPixels(pixels);

   std::vector<double> energies;
   dataMap.getAxisVector(2, energies);

   std::vector<double> map;
   computeModelMap(pixels, energies, map);

   CountsMap * modelMap = new CountsMap(dataMap);
   modelMap->setImage(map);
   return modelMap;
}

void SourceModel::computeModelMap(const std::vector<Pixel> & pixels,
                                  const std::vector<double> & energies,
                                  std::vector<double> & modelMap) const {
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

} // namespace Likelihood
