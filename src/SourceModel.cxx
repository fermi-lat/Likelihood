/** 
 * @file SourceModel.cpp
 * @brief SourceModel class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceModel.cxx,v 1.100 2016/10/13 02:03:06 echarles Exp $
 */

#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <fstream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// RapidXML headers
#include "xmlBase/rapidxml.hpp"

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "astro/SkyDir.h"

#include "optimizers/Arg.h"
#include "optimizers/FunctionFactory.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/CountsMapBase.h"
#include "Likelihood/CompositeSource.h"
#include "Likelihood/Exception.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/SourceMap.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SourceFactory.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/FluxBuilder.h"
#include "Likelihood/SourceModelBuilder.h"
#include "Likelihood/SourceModel.h"
#include "Likelihood/XmlParser.h"

namespace Likelihood {

namespace {

// RAII wrapper for RapidXML document that owns its buffer
class XmlDocument {
public:
   XmlDocument() = default;
   ~XmlDocument() = default;
   
   // Non-copyable
   XmlDocument(const XmlDocument&) = delete;
   auto operator=(const XmlDocument&) -> XmlDocument& = delete;
   
   // Movable
   XmlDocument(XmlDocument&&) = default;
   auto operator=(XmlDocument&&) -> XmlDocument& = default;
   
   auto parseFile(const std::string& filename) -> bool {
      std::ifstream file(filename);
      if (!file.is_open()) {
         return false;
      }
      
      // Read entire file into buffer
      std::ostringstream ss;
      ss << file.rdbuf();
      m_buffer = ss.str();
      
      try {
         // RapidXML modifies the buffer, so we need a non-const copy
         m_doc.parse<rapidxml::parse_default>(&m_buffer[0]);
         return true;
      } catch (const rapidxml::parse_error&) {
         return false;
      }
   }
   
   auto first_node(const char* name = nullptr) -> rapidxml::xml_node<>* {
      return m_doc.first_node(name);
   }
   
   auto get() -> rapidxml::xml_document<>* { return &m_doc; }

private:
   std::string m_buffer;
   rapidxml::xml_document<> m_doc;
};

// Helper function to collect child nodes by tag name
void collectChildNodes(rapidxml::xml_node<>* parent, 
                       const char* tagName,
                       std::vector<rapidxml::xml_node<>*>& children) {
   children.clear();
   if (parent == nullptr) return;
   
   for (rapidxml::xml_node<>* child = parent->first_node(tagName);
        child != nullptr;
        child = child->next_sibling(tagName)) {
      children.push_back(child);
   }
}

// Helper function to get string attribute value
auto getStringAttribute(rapidxml::xml_node<>* node, 
                                              const char* attrName) -> std::optional<std::string> {
   if (node == nullptr) return std::nullopt;
   
   rapidxml::xml_attribute<>* attr = node->first_attribute(attrName);
   if (attr == nullptr) return std::nullopt;
   
   return std::string(attr->value(), attr->value_size());
}

// Helper function to get double attribute value
auto getDoubleAttribute(rapidxml::xml_node<>* node,
                                         const char* attrName) -> std::optional<double> {
   auto strVal = getStringAttribute(node, attrName);
   if (!strVal) return std::nullopt;
   
   try {
      return std::stod(*strVal);
   } catch (...) {
      return std::nullopt;
   }
}

// Helper function to parse a Parameter from XML node
auto parseParameterNode(rapidxml::xml_node<>* paramNode) -> std::optional<optimizers::Parameter> {
   if (paramNode == nullptr) return std::nullopt;
   
   auto nameOpt = getStringAttribute(paramNode, "name");
   if (!nameOpt) return std::nullopt;
   
   auto valueOpt = getDoubleAttribute(paramNode, "value");
   if (!valueOpt) return std::nullopt;
   
   optimizers::Parameter param(*nameOpt, *valueOpt);
   
   // Set optional attributes
   if (auto scaleOpt = getDoubleAttribute(paramNode, "scale")) {
      param.setScale(*scaleOpt);
   }
   
   if (auto minOpt = getDoubleAttribute(paramNode, "min")) {
      if (auto maxOpt = getDoubleAttribute(paramNode, "max")) {
         param.setBounds(*minOpt, *maxOpt);
      }
   }
   
   if (auto freeOpt = getStringAttribute(paramNode, "free")) {
      param.setFree(*freeOpt == "1" || *freeOpt == "true");
   }
   
   return param;
}

} // anonymous namespace

SourceModel::SourceModel(const Observation & observation, bool verbose) 
   : optimizers::Statistic("SourceModel", 0),
     m_observation(observation), 
     m_useNewImp(true), 
     m_verbose(verbose), 
     m_formatter(new st_stream::StreamFormatter("SourceModel", "", 2)) {
   if (char* useOldImp = std::getenv("USE_OLD_LOGLIKE"); useOldImp != nullptr) {
      m_useNewImp = false;
   }
}

SourceModel::SourceModel(const SourceModel &rhs) 
   : optimizers::Statistic(rhs),
     m_observation(rhs.m_observation), 
     m_useNewImp(rhs.m_useNewImp),
     m_verbose(rhs.m_verbose),
     m_formatter(new st_stream::StreamFormatter("SourceModel", "", 2)) {
   for (const auto& [name, source] : rhs.m_sources) {
      m_sources[name] = source->clone();
   }
   findFreeSrcs();   
}

SourceModel::~SourceModel() {
   for (auto& [name, source] : m_sources) {
      delete source;
   }
   m_sources.clear();
   delete m_formatter;
}

void SourceModel::addSource(Source *src, bool fromClone, SourceMap* /* srcMap */, bool /* loadMap */) {
   if (!m_sources.count(src->getName())) {
      m_sources[src->getName()] = fromClone ? src->clone() : src;
      m_sources[src->getName()]->setObservation(&m_observation);
      syncParams();
   } else {
      throw Exception("Likelihood::SourceModel:\nSource named " 
                      + src->getName() + " already exists.");
   }
}
 
auto SourceModel::deleteSource(const std::string &srcName) -> Source * {
   if (auto it = m_sources.find(srcName); it != m_sources.end()) {
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
   for (const auto& name : srcNames) {
      Source * source = deleteSource(name);
      delete source;
   }
   m_parameter.clear();
}

auto SourceModel::getSource(const std::string & srcName) -> Source * {
   if (auto it = m_sources.find(srcName); it != m_sources.end()) {
      return it->second;
   }
   return nullptr;
}

auto SourceModel::source(const std::string & srcName) const -> const Source & {
   if (auto my_src = m_sources.find(srcName); my_src != m_sources.end()) {
      return *(my_src->second);
   }
   throw std::runtime_error("SourceModel::source: Source " + 
                            srcName + " not found.");
}

void SourceModel::getSources(const std::vector<std::string>& srcNames, 
                             std::vector<const Source*>& srcs) const {
   srcs.clear();
   srcs.reserve(srcNames.size());
   for (const auto& name : srcNames) {
      const Source& src = source(name);
      srcs.push_back(&src);
   }
}

void SourceModel::getSrcNames(std::vector<std::string> &names) const {
   names.clear();
   names.reserve(m_sources.size());
   for (const auto& [name, source] : m_sources) {
      names.push_back(name);
   }
}

auto SourceModel::splitCompositeSource(const std::string& compName,
                                                        std::vector<std::string>& srcNames) -> optimizers::Function* {
   Source* src = deleteSource(compName);
   if (src == nullptr) {
      throw std::runtime_error("SourceModel::splitCompositeSource did not find source: " + compName);
   }
   auto* cmp = dynamic_cast<CompositeSource*>(src);
   if (cmp == nullptr) {
      throw std::runtime_error("SourceModel::splitCompositeSource source is not composite: " + compName);
   }
   srcNames.clear();
   cmp->getSrcNames(srcNames);
   for (const auto& name : srcNames) {
      cmp->give_source(*this, name);
   }
   optimizers::Function* specFunc = cmp->spectrum().clone();
   delete cmp;
   return specFunc;
}
  
auto SourceModel::steal_source(SourceModel& other,
                                  const std::string& srcName,
                                  SourceMap* srcMap) -> Source* {
   Source* src = other.deleteSource(srcName);
   if (src == nullptr) {
      throw std::runtime_error("SourceModel::steal_source could not find source: " + srcName);
   }
   addSource(src, false, srcMap);
   return src;
}

auto SourceModel::give_source(SourceModel& other,
                                 const std::string& srcName,
                                 SourceMap* srcMap) -> Source* {
   Source* src = deleteSource(srcName);
   if (src == nullptr) {
      throw std::runtime_error("SourceModel::give_source could not find source: " + srcName);
   }
   other.addSource(src, false, srcMap);
   return src;
}

auto SourceModel::mergeSources(const std::string& compName,
                                           const std::vector<std::string>& srcNames,
                                           const std::string& specFuncName) -> CompositeSource* {
   auto cmp = new CompositeSource(m_observation, compName, specFuncName);
   initialize_composite(*cmp);
   for (const auto& name : srcNames) {
      cmp->steal_source(*this, name);
   }
   addSource(cmp, false);
   return cmp;
}

auto SourceModel::value(const optimizers::Arg &) const -> double {
   return 0;
}

void SourceModel::setParam(const optimizers::Parameter & param, 
                           const std::string & funcName,
                           const std::string & srcName) {
   if (m_sources.count(srcName)) {
      if (funcName == "Spectrum") {
         m_sources[srcName]->spectrum().setParam(param);
         syncParams();
         return;
      }
   }
   std::string errorMessage = "SourceModel::setParam:\n Function " 
      + funcName + " for Source " + srcName + " was not found.";
   throw optimizers::Exception(errorMessage);
}
 
auto 
SourceModel::setParamValues_(std::vector<double>::const_iterator it) -> std::vector<double>::const_iterator {
   for (auto& [name, source] : m_sources) {
      it = source->spectrum().setParamValues_(it);
   }
   syncParams();
   return it;
}

auto 
SourceModel::setFreeParamValues_(std::vector<double>::const_iterator it) -> std::vector<double>::const_iterator {
   for (auto& [name, source] : m_sources) {
      it = source->spectrum().setFreeParamValues_(it);
   }
   syncParams();
   return it;
}

auto SourceModel::getParam(const std::string &paramName,
                                            const std::string &funcName,
                                            const std::string &srcName) const -> optimizers::Parameter {
   if (auto srcIt = m_sources.find(srcName); srcIt != m_sources.end()) {
      if (funcName == "Spectrum") {
         return srcIt->second->spectrum().getParam(paramName);
      }
   }
   throw optimizers::ParameterNotFound(paramName, funcName, 
                                       "SourceModel::getParam");
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
   const size_t numParams = setFree ? getNumFreeParams() : getNumParams();
   
   if (params.size() != numParams) {
      throw optimizers::Exception(
         "SourceModel::setParams:\nInconsistent number of Parameters.");
   }

   // Assume ordering of Parameters in params matches that given by the
   // ordering of the Sources and their Functions.
   size_t k = 0;  // params' index
   for (auto& [name, source] : m_sources) {
      const size_t srcNumParams = setFree 
         ? source->spectrum().getNumFreeParams()
         : source->spectrum().getNumParams();
      
      for (size_t j = 0; j < srcNumParams; ++j) {
         source->spectrum().setParam(params[k]);
         ++k;
      }
   }
   syncParams();
}

void SourceModel::syncParams() {
   m_parameter.clear();
   for (const auto& [name, source] : m_sources) {
      std::vector<optimizers::Parameter> srcParams;
      source->spectrum().getParams(srcParams);
      for (const auto& param : srcParams) {
         m_parameter.push_back(param);
      }
   }
   findFreeSrcs();
}

void SourceModel::fetchDerivs(optimizers::Arg &x, std::vector<double> &derivs, 
                              bool getFree) const {
   derivs.clear();
   for (const auto& [name, source] : m_sources) {
      std::vector<double> my_derivs;
      if (getFree) {
         source->spectrum().getFreeDerivs(x, my_derivs);
      } else {
         source->spectrum().getDerivs(x, my_derivs);
      }
      for (const auto& deriv : my_derivs) {
         derivs.push_back(deriv);
      }
   }
}

void SourceModel::readXml(const std::string& xmlFile,
                          optimizers::FunctionFactory &funcFactory,
                          bool requireExposure,
                          bool addPointSources,
                          bool loadMaps) {
   // Expand any environment variables in the xmlFile name.
   std::string expandedFile = xmlFile;
   facilities::Util::expandEnvVar(&expandedFile);

   XmlDocument doc;
   if (!doc.parseFile(expandedFile)) {
      std::string errorMessage = "SourceModel::readXml:\nInput xml file, "
         + expandedFile + ", not parsed successfully.";
      throw Exception(errorMessage);
   }
   
   rapidxml::xml_node<>* source_library = doc.first_node("source_library");
   
   if (source_library == nullptr) {
      throw Exception("SourceModel::readXml:\nsource_library not found");
   }
   
   readXml(source_library, expandedFile, funcFactory, requireExposure, 
           addPointSources, loadMaps);
}

void SourceModel::readXml(rapidxml::xml_node<>* srcLibrary,
                          const std::string& xmlFile,
                          optimizers::FunctionFactory & funcFactory,
                          bool requireExposure,
                          bool addPointSources,
                          bool loadMaps) {
   // Create a SourceFactory to read in the xml file.
   SourceFactory srcFactory(m_observation);
   try {
      srcFactory.readXml(srcLibrary, xmlFile, funcFactory, requireExposure,
                         addPointSources, loadMaps);
   } catch (std::exception & eObj) {
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

   for (const auto& name : srcNames) {
      Source * src = srcFactory.releaseSource(name);
      if (m_verbose) {
         m_formatter->info() << "adding source " << name << std::endl;
      }
      addSource(src, false);
   }
   syncParams();
}

void SourceModel::reReadXml(const std::string& xmlFile) {
   std::string expandedFile = xmlFile;
   facilities::Util::expandEnvVar(&expandedFile);

   XmlDocument doc;
   if (!doc.parseFile(expandedFile)) {
      std::string errorMessage = "SourceModel::reReadXml:\nInput xml file, "
         + expandedFile + ", not parsed successfully.";
      throw Exception(errorMessage);
   }
   
   rapidxml::xml_node<>* source_library = doc.first_node("source_library");
   
   if (source_library == nullptr) {
      throw Exception("SourceModel::reReadXml:\nsource_library not found");
   }
   
   reReadXml(source_library);
}

void SourceModel::reReadXml(rapidxml::xml_node<>* source_library) {
   // Loop through source xml nodes and Source objects in parallel.
   std::vector<rapidxml::xml_node<>*> srcs;
   collectChildNodes(source_library, "source", srcs);

   for (auto* srcNode : srcs) {
      auto srcNameOpt = getStringAttribute(srcNode, "name");
      if (!srcNameOpt) {
         continue;
      }
      std::string srcName = *srcNameOpt;
      
      // Find corresponding Source in model
      Source* my_source = getSource(srcName);
      if (my_source == nullptr) {
         continue;  // Source not in model, skip
      }

      // Get source type
      auto srcTypeOpt = getStringAttribute(srcNode, "type");
      if (!srcTypeOpt) {
         continue;
      }
      std::string srcType = *srcTypeOpt;

      // Get spectrum node and update parameters
      rapidxml::xml_node<>* spectrumNode = srcNode->first_node("spectrum");
      if (spectrumNode != nullptr) {
         std::vector<optimizers::Parameter> spectrumModel;
         std::vector<rapidxml::xml_node<>*> specParams;
         collectChildNodes(spectrumNode, "parameter", specParams);
         
         for (auto* paramNode : specParams) {
            auto paramOpt = parseParameterNode(paramNode);
            if (paramOpt) {
               spectrumModel.push_back(*paramOpt);
            }
         }
         
         if (!spectrumModel.empty()) {
            my_source->spectrum().setParams(spectrumModel);
         }
      }

      // Handle source type specific processing
      if (srcType == "PointSource") {
         rapidxml::xml_node<>* spatialModelNode = srcNode->first_node("spatialModel");
         if (spatialModelNode != nullptr) {
            double ra = 0.0, dec = 0.0;
            
            std::vector<rapidxml::xml_node<>*> params;
            collectChildNodes(spatialModelNode, "parameter", params);
            
            for (auto* paramNode : params) {
               auto nameOpt = getStringAttribute(paramNode, "name");
               auto valueOpt = getDoubleAttribute(paramNode, "value");
               
               if (nameOpt && valueOpt) {
                  if (*nameOpt == "RA") {
                     ra = *valueOpt;
                  } else if (*nameOpt == "DEC") {
                     dec = *valueOpt;
                  }
               }
            }
            
            astro::SkyDir newDir(ra, dec);
            constexpr double tol = 1e-4;
            auto* ptSrc = dynamic_cast<PointSource*>(my_source);
            if (ptSrc != nullptr && newDir.difference(ptSrc->getDir()) > tol) {
               // Reset the direction, re-computing the PointSource exposure.
               ptSrc->setDir(newDir);
            }
         }
      } else if (srcType == "DiffuseSource") {
         rapidxml::xml_node<>* spatialModelNode = srcNode->first_node("spatialModel");
         if (spatialModelNode != nullptr) {
            std::vector<optimizers::Parameter> spatialModel;
            std::vector<rapidxml::xml_node<>*> spatialParams;
            collectChildNodes(spatialModelNode, "parameter", spatialParams);
            
            for (auto* paramNode : spatialParams) {
               auto paramOpt = parseParameterNode(paramNode);
               if (paramOpt) {
                  spatialModel.push_back(*paramOpt);
               }
            }
            
            if (!spatialModel.empty()) {
               my_source->getSrcFuncs()["SpatialDist"]->setParams(spatialModel);
            }
         }
      } else {
         std::ostringstream message;
         message << "SourceModel::reReadXml: "
                 << "Unknown Source type: " << srcType;
         throw std::runtime_error(message.str());
      }
   }
   syncParams();
}

void SourceModel::writeXml(const std::string& xmlFile,
                           const std::string & functionLibrary,
                           const std::string & srcLibTitle) {
   SourceModelBuilder builder(functionLibrary, srcLibTitle);
   writeXml(builder);
   builder.write(xmlFile);
}

void SourceModel::write_fluxXml(const std::string& xmlFile) {
   auto ebounds = m_observation.roiCuts().getEnergyCuts();
   FluxBuilder builder(ebounds.first, ebounds.second);
   write_fluxXml(builder);
   builder.write(xmlFile);
}

void SourceModel::writeXml(SourceModelBuilder& builder) {
   builder.addSourceModel(*this);
}

void SourceModel::write_fluxXml(FluxBuilder& builder) {
   builder.addSourceModel(*this);   
}

auto SourceModel::hasSrcNamed(const std::string & srcName) const -> bool {
   return m_sources.find(srcName) != m_sources.end();
}

auto SourceModel::createCountsMap(const CountsMapBase & dataMap) const -> CountsMapBase * {
   const std::vector<Pixel> & pixels = dataMap.pixels();

   std::vector<double> energies;
   dataMap.getEnergies(energies);

   std::vector<float> map;
   computeModelMap(pixels, energies, map);

   CountsMapBase * modelMap = dataMap.clone();
   modelMap->setImage(map);
   return modelMap;
}

void SourceModel::computeModelMap(const std::vector<Pixel> & pixels,
                                  const std::vector<double> & energies,
                                  std::vector<float> & modelMap) const {
   modelMap.clear();
   modelMap.reserve(pixels.size() * (energies.size() - 1));
   
   for (size_t k = 0; k < energies.size() - 1; ++k) {
      for (size_t j = 0; j < pixels.size(); ++j) {
         modelMap.push_back(
            pixels[j].modelCounts(energies[k], energies[k+1],
                                  *const_cast<SourceModel*>(this)));
      }
   }
}

void SourceModel::findFreeSrcs() {
   m_freeSrcs.clear();
   for (const auto& [name, source] : m_sources) {
      if (source->spectrum().getNumFreeParams() > 0) {
         m_freeSrcs.push_back(source);
      }
   }
}

} // namespace Likelihood
