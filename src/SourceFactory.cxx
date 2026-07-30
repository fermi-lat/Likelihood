/**
 * @file SourceFactory.cpp
 * @brief Implementation for the SourceFactory class, which applies the
 * Prototype pattern to return clones of various gamma-ray Sources.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceFactory.cxx,v 1.87 2016/10/28 00:44:25 echarles Exp $
 */

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

// RapidXML-based XML framework (replaces Xerces-C)
#include "xmlBase/rapidxml.hpp"
#include "xmlBase/safe_xml_parser.hpp"

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "optimizers/Exception.h"
#include "optimizers/FunctionFactory.h"

#include "Likelihood/CompositeSource.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Event.h"
#include "Likelihood/Exception.h"
#include "Likelihood/FileFunction.h"
#include "Likelihood/DMFitFunction.h"
#include "Likelihood/MapBase.h"
#include "Likelihood/MultipleBrokenPowerLaw.h"
#include "Likelihood/PiecewisePowerLaw.h"
#include "Likelihood/Observation.h"
#include "Likelihood/RadialDisk.h"
#include "Likelihood/RadialGaussian.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/RadialProfile.h"
#include "Likelihood/ScaleFactor.h"
#include "Likelihood/SourceFactory.h"

namespace Likelihood {

  // ============================================================================
  // Construction / Destruction
  // ============================================================================

  SourceFactory::SourceFactory(const Observation& observation, bool verbose)
    : m_verbose(verbose)
    , m_requireExposure(true)
    , m_observation(observation)
    , m_formatter(new st_stream::StreamFormatter("SourceFactory", "", 2))
  {
  }

  SourceFactory::~SourceFactory() {
    for (auto& [name, src] : m_prototypes) {
      delete src;
    }
    m_prototypes.clear();
    delete m_formatter;
  }

  // ============================================================================
  // Source Management
  // ============================================================================

  Source* SourceFactory::create(const std::string& name) {
    auto it = m_prototypes.find(name);
    if (it == m_prototypes.end()) {
      std::ostringstream message;
      message << "SourceFactory::create:\n"
	      << "Cannot create Source named '" << name << "'";
      throw Exception(message.str());
    }
    return it->second->clone();
  }

  Source* SourceFactory::releaseSource(const std::string& name) {
    auto it = m_prototypes.find(name);
    if (it == m_prototypes.end()) {
      std::ostringstream message;
      message << "SourceFactory::releaseSource:\n"
	      << "Source named '" << name << "' not found";
      throw Exception(message.str());
    }
    Source* src = it->second;
    m_prototypes.erase(it);
    return src;
  }

  void SourceFactory::addSource(const std::string& name, Source* src, bool fromClone) {
    if (m_prototypes.count(name)) {
      std::ostringstream message;
      message << "SourceFactory::addSource:\n"
	      << "A Source named '" << name << "' already exists.";
      throw Exception(message.str());
    }

    if (fromClone) {
      m_prototypes[name] = src->clone();
    } else {
      m_prototypes[name] = src;
    }
  }

  void SourceFactory::replaceSource(Source* src, bool fromClone) {
    const std::string& name = src->getName();
    auto it = m_prototypes.find(name);
    if (it != m_prototypes.end()) {
      delete it->second;
      m_prototypes.erase(it);
    }
    addSource(name, src, fromClone);
  }

  void SourceFactory::fetchSrcNames(std::vector<std::string>& srcNames) const {
    srcNames.clear();
    srcNames.reserve(m_prototypes.size());
    for (const auto& [name, src] : m_prototypes) {
      srcNames.push_back(name);
    }
  }

  std::vector<std::string> SourceFactory::fetchSrcNames() const {
    std::vector<std::string> srcNames;
    srcNames.reserve(m_prototypes.size());
    for (const auto& [name, src] : m_prototypes) {
      srcNames.push_back(name);
    }
    return srcNames;
  }

  // ============================================================================
  // RapidXML Helper Methods
  // ============================================================================

  std::string SourceFactory::getAttributeValue(const rapidxml::xml_node<>* node,
					       const char* attrName,
					       const std::string& defaultValue) {
    if (!node) return defaultValue;

    auto* attr = node->first_attribute(attrName);
    if (attr && attr->value()) {
      return std::string(attr->value(), attr->value_size());
    }
    return defaultValue;
  }

  double SourceFactory::getAttributeValueAsDouble(const rapidxml::xml_node<>* node,
						  const char* attrName,
						  double defaultValue) {
    std::string value = getAttributeValue(node, attrName, "");
    if (value.empty()) return defaultValue;

    try {
      return std::stod(value);
    } catch (const std::exception&) {
      return defaultValue;
    }
  }

  std::vector<rapidxml::xml_node<>*> SourceFactory::collectChildren(
								    const rapidxml::xml_node<>* parent,
								    const char* childName) {

    std::vector<rapidxml::xml_node<>*> children;
    if (!parent) return children;

    for (auto* child = parent->first_node(childName);
	 child != nullptr;
	 child = child->next_sibling(childName)) {
      children.push_back(child);
    }
    return children;
  }

  // ============================================================================
  // XML Reading
  // ============================================================================

  void SourceFactory::readXml(const std::string& xmlFile,
			      optimizers::FunctionFactory& funcFactory,
			      bool requireExposure,
			      bool addPointSources,
			      bool loadMaps) {
    // Expand environment variables in the file path
    std::string expandedPath = xmlFile;
    facilities::Util::expandEnvVar(&expandedPath);

    // Load and parse the XML file
    std::ifstream file(expandedPath, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
      throw Exception("SourceFactory::readXml: Cannot open file: " + xmlFile);
    }

    auto size = file.tellg();
    file.seekg(0, std::ios::beg);

    m_xmlBuffer.resize(static_cast<size_t>(size) + 1);
    if (!file.read(m_xmlBuffer.data(), size)) {
      throw Exception("SourceFactory::readXml: Failed to read file: " + xmlFile);
    }
    m_xmlBuffer[static_cast<size_t>(size)] = '\0';

    // Parse the XML
    try {
      m_xmlDoc.clear();
      m_xmlDoc.parse<rapidxml::parse_default>(m_xmlBuffer.data());
    } catch (const rapidxml::parse_error& e) {
      std::ostringstream message;
      message << "SourceFactory::readXml: XML parse error in file " << xmlFile
	      << ": " << e.what();
      throw Exception(message.str());
    }

    // Get the root element (source_library)
    auto* source_library = m_xmlDoc.first_node("source_library");
    if (!source_library) {
      throw Exception("SourceFactory::readXml: No source_library element in " + xmlFile);
    }

    // Call the overloaded version with the parsed node
    readXml(source_library, xmlFile, funcFactory, requireExposure, addPointSources, loadMaps);
  }

  void SourceFactory::readXml(rapidxml::xml_node<>* source_library,
			      const std::string& xmlFile,
			      optimizers::FunctionFactory& funcFactory,
			      bool requireExposure,
			      bool addPointSources,
			      bool loadMaps) {
    // Prepare the FunctionFactory using the xml file specified in the source_library tag
    m_requireExposure = requireExposure;
    std::string function_library = getAttributeValue(source_library, "function_library", "");

    if (!function_library.empty() && function_library.find("xml") != std::string::npos) {
      facilities::Util::expandEnvVar(&function_library);
      try {
	funcFactory.readXml(function_library);
      } catch (optimizers::Exception& eObj) {
	m_formatter->err() << eObj.what() << std::endl;
	throw;
      }
    }

    std::vector<Source*> sources;
    makeSources(xmlFile, source_library, sources, funcFactory,
		requireExposure, addPointSources, loadMaps);

    for (auto* src : sources) {
      addSource(src->getName(), src);
      delete src;
    }
  }

  void SourceFactory::makeSources(const std::string& xmlFile,
				  const rapidxml::xml_node<>* source_library,
				  std::vector<Source*>& sources,
				  optimizers::FunctionFactory& funcFactory,
				  bool requireExposure,
				  bool addPointSources,
				  bool loadMaps) {
    // Loop through source elements
    auto srcs = collectChildren(source_library, "source");

    for (auto* srcNode : srcs) {
      Source* src = nullptr;

      // Get the type and name of this source
      std::string srcType = getAttributeValue(srcNode, "type", "");
      std::string srcName = getAttributeValue(srcNode, "name", "");

      m_currentSrcName = srcName;

      m_formatter->info(3) << "Creating source named " << srcName << std::endl;

      // Retrieve the spectrum element
      auto spectrumChildren = collectChildren(srcNode, "spectrum");
      if (spectrumChildren.size() != 1) {
	std::ostringstream message;
	message << "Error parsing xml model file: \n"
		<< xmlFile << "\n"
		<< "for source " << srcName << "\n"
		<< "Missing spectral model component.";
	throw Exception(message.str());
      }
      auto* spectrum = spectrumChildren[0];

      // Process based on source type
      if (srcType == "CompositeSource") {
	auto sourceLibChildren = collectChildren(srcNode, "source_library");
	if (sourceLibChildren.size() != 1) {
	  std::ostringstream message;
	  message << "Error parsing xml model file: \n"
		  << xmlFile << "\n"
		  << "for source " << srcName << ".\n"
		  << "Missing source_library component.\n"
		  << "Please check that you are using the correct xml format.";
	  throw Exception(message.str());
	}
	auto* nestedSourceLibrary = sourceLibChildren[0];
	src = makeCompositeSource(xmlFile, spectrum, nestedSourceLibrary, funcFactory,
				  requireExposure, addPointSources, loadMaps);
      } else {
	auto spatialModelChildren = collectChildren(srcNode, "spatialModel");
	if (spatialModelChildren.size() != 1) {
	  std::ostringstream message;
	  message << "Error parsing xml model file: \n"
		  << xmlFile << "\n"
		  << "for source " << srcName << ".\n"
		  << "Missing spatial model component.\n"
		  << "Please check that you are using the correct xml format.";
	  throw Exception(message.str());
	}
	auto* spatialModel = spatialModelChildren[0];

	if (addPointSources && srcType == "PointSource") {
	  src = makePointSource(spectrum, spatialModel, funcFactory);
	} else if (srcType == "DiffuseSource") {
	  src = makeDiffuseSource(spectrum, spatialModel, funcFactory, loadMaps);
	}

	if (src != nullptr) {
	  // Determine if psf can be applied
	  std::string apply_psf = getAttributeValue(spatialModel, "apply_psf", "");
	  if (apply_psf != "true" && apply_psf != "false" && !apply_psf.empty()) {
	    throw std::runtime_error("Invalid value for apply_psf attribute in xml definition of " + src->getName());
	  }
	  src->set_psf_flag(apply_psf != "false");

	  // Determine if exposure can be applied
	  std::string apply_exposure = getAttributeValue(spatialModel, "apply_exposure", "");
	  if (apply_exposure != "true" && apply_exposure != "false" && !apply_exposure.empty()) {
	    throw std::runtime_error("Invalid value for apply_exposure attribute in xml definition of " + src->getName());
	  }
	  src->set_exposure_flag(apply_exposure != "false");
	}
      }

      if (src != nullptr) {
	src->setName(srcName);
	sources.push_back(src);
      }
    }
  }

  // ============================================================================
  // Source Creation Methods
  // ============================================================================

  Source* SourceFactory::makePointSource(const rapidxml::xml_node<>* spectrum,
					 const rapidxml::xml_node<>* spatialModel,
					 optimizers::FunctionFactory& funcFactory) {
    std::string funcType = getAttributeValue(spatialModel, "type", "");

    if (funcType != "SkyDirFunction") {
      std::ostringstream message;
      message << "SourceFactory::readXml:\n"
	      << "Trying to create a PointSource with a spatialModel of type "
	      << funcType << ".";
      throw Exception(message.str());
    }

    // Extract (RA, Dec) from the parameter elements
    double ra = 0, dec = 0;
    auto params = collectChildren(spatialModel, "parameter");

    for (auto* paramNode : params) {
      std::string name = getAttributeValue(paramNode, "name", "");
      if (name == "RA") {
	ra = getAttributeValueAsDouble(paramNode, "value", 0.0);
      }
      if (name == "DEC") {
	dec = getAttributeValueAsDouble(paramNode, "value", 0.0);
      }
    }

    checkRoiDist(ra, dec);

    auto* src = new PointSource(ra, dec, m_observation, m_requireExposure);
    setSpectrum(src, spectrum, funcFactory);

    return src;
  }

  Source* SourceFactory::makeDiffuseSource(const rapidxml::xml_node<>* spectrum,
                                          const rapidxml::xml_node<>* spatialModel,
                                          optimizers::FunctionFactory& funcFactory,
                                          bool loadMap) {
   // Get the spatial model type and create the corresponding function
   std::string type = getAttributeValue(spatialModel, "type", "");
   
   optimizers::Function* spatialDist = funcFactory.create(type);
   
   // Set spatial distribution parameters from XML
   auto params = collectChildren(spatialModel, "parameter");
   for (auto* paramNode : params) {
      std::string name = getAttributeValue(paramNode, "name", "");
      spatialDist->parameter(name).extractDomData(paramNode);
   }
   
   bool mapBasedIntegral = false;
   
   if (type == "SpatialMap" || type == "MapCubeFunction") {
      std::string fitsFile = getAttributeValue(spatialModel, "file", "");
      facilities::Util::expandEnvVar(&fitsFile);
      dynamic_cast<MapBase*>(spatialDist)->readFitsFile(fitsFile, "", loadMap);
      std::string map_based_integral = getAttributeValue(spatialModel, "map_based_integral", "");
      mapBasedIntegral = (map_based_integral == "true");
   } else if (type == "RadialProfile") {
      std::string tpl_file = getAttributeValue(spatialModel, "file", "");
      facilities::Util::expandEnvVar(&tpl_file);
      dynamic_cast<RadialProfile*>(spatialDist)->readTemplateFile(tpl_file);
   } else if (type == "RadialGaussian" || type == "RadialDisk") {
      dynamic_cast<SpatialFunction*>(spatialDist)->setParams(spatialModel);
      //spatialDist->setParams(spatialModel);
   }
   
   // Create the DiffuseSource with the spatial distribution function
   DiffuseSource* src = new DiffuseSource(spatialDist, m_observation, 
                                          m_requireExposure, mapBasedIntegral);
   
   // Set the spectrum
   setSpectrum(src, spectrum, funcFactory);
   
   return src;
}

  Source* SourceFactory::makeCompositeSource(const std::string& xmlFile,
					     const rapidxml::xml_node<>* spectrum,
					     rapidxml::xml_node<>* source_library,
					     optimizers::FunctionFactory& funcFactory,
					     bool requireExposure,
					     bool addPointSources,
					     bool loadMap) {
    std::vector<Source*> sources;
    makeSources(xmlFile, source_library, sources, funcFactory,
		requireExposure, addPointSources, loadMap);

    auto* comp_src = new CompositeSource(m_observation);
    SourceMap* srcMap = nullptr;

    for (auto* src : sources) {
      addSource(src->getName(), src);
      comp_src->addSource(src, srcMap, true);
      releaseSource(src->getName());
      delete src;
    }

    setSpectrum(comp_src, spectrum, funcFactory);
    return comp_src;
  }

  void SourceFactory::setSpectrum(Source* src,
				  const rapidxml::xml_node<>* spectrum,
				  optimizers::FunctionFactory& funcFactory) {
    std::string type = getAttributeValue(spectrum, "type", "");

    optimizers::Function* spec = funcFactory.create(type);

    // Fetch the parameter elements
    auto params = collectChildren(spectrum, "parameter");

    // Handle special function types
    if (type == "MultipleBPL") {
      addParamsToMultipleBPL(spec, params, src);
    }
    if (type == "PiecewisePowerLaw") {
      addParamsToPiecewisePL(spec, params, src);
    }

    // Set parameter values from XML
    for (auto* paramNode : params) {
      std::string name = getAttributeValue(paramNode, "name", "");
      spec->parameter(name).extractDomData(paramNode);
    }

    // Handle FileFunction types
    if (type == "FileFunction") {
      std::string filename = getAttributeValue(spectrum, "file", "");
      dynamic_cast<FileFunction*>(spec)->readFunction(filename);
    }
    if (type == "ScaleFactor::FileFunction") {
      std::string filename = getAttributeValue(spectrum, "file", "");
      dynamic_cast<FileFunction*>(dynamic_cast<ScaleFactor*>(spec)->spectrum())->readFunction(filename);
    }
    if (type == "DMFitFunction") {
      std::string filename = getAttributeValue(spectrum, "file", "");
      dynamic_cast<DMFitFunction*>(spec)->readFunction(filename);
    }

    // Check for spectral scaling
    std::string scaling_file = getAttributeValue(spectrum, "scaling_file", "");
    if (!scaling_file.empty()) {
      FileFunction scalingFunction;
      scalingFunction.readFunction(scaling_file);
      spec->setScalingFunction(scalingFunction);
    }

    src->setSpectrum(spec);

    // Determine if energy dispersion can be applied
    std::string apply_edisp = getAttributeValue(spectrum, "apply_edisp", "");
    if (apply_edisp != "true" && apply_edisp != "false" && !apply_edisp.empty()) {
      throw std::runtime_error("Invalid value for apply_edisp attribute in xml definition of " + src->getName());
    }
    src->set_edisp_flag(apply_edisp != "false");

    delete spec;
  }

  // ============================================================================
  // Parameter Handling for Special Functions
  // ============================================================================

  void SourceFactory::addParamsToMultipleBPL(optimizers::Function* spec,
					     const std::vector<rapidxml::xml_node<>*>& params,
					     const Source* src) const {
    double normalization = 1.0;
    std::vector<double> photonIndexes;
    std::vector<double> breakEnergies;

    for (auto* paramNode : params) {
      std::string name = getAttributeValue(paramNode, "name", "");

      if (name == "Normalization") {
	normalization = getAttributeValueAsDouble(paramNode, "value", 1.0);
      }
      if (name.find("Index") == 0) {
	size_t idx = static_cast<size_t>(std::atoi(name.substr(5).c_str()));
	if (idx >= photonIndexes.size()) {
	  photonIndexes.resize(idx + 1);
	}
	photonIndexes[idx] = getAttributeValueAsDouble(paramNode, "value", 0.0);
      }
      if (name.find("Break") == 0) {
	size_t idx = static_cast<size_t>(std::atoi(name.substr(5).c_str()));
	if (idx >= breakEnergies.size()) {
	  breakEnergies.resize(idx + 1);
	}
	double value = getAttributeValueAsDouble(paramNode, "value", 0.0);
	double scale = getAttributeValueAsDouble(paramNode, "scale", 1.0);
	breakEnergies[idx] = value * scale;
      }
    }

    dynamic_cast<MultipleBrokenPowerLaw*>(spec)->addParams(normalization, photonIndexes, breakEnergies);
  }

  void SourceFactory::addParamsToPiecewisePL(optimizers::Function* spec,
                                            const std::vector<rapidxml::xml_node<>*>& params,
                                            const Source* src) const {
   // Extract dNdE# and Energy# parameters
   std::vector<double> dNdEs;
   std::vector<double> energies;
   
   for (auto* paramNode : params) {
      std::string name = getAttributeValue(paramNode, "name", "");
      
      if (name.substr(0, 4) == "dNdE") {
         size_t idx = static_cast<size_t>(std::atoi(name.substr(4).c_str()));
         if (idx >= dNdEs.size()) {
            dNdEs.resize(idx + 1);
         }
         // Placeholder value - actual value will be set by extractDomData later
         dNdEs[idx] = 0.0;
      }
      if (name.substr(0, 6) == "Energy") {
         size_t idx = static_cast<size_t>(std::atoi(name.substr(6).c_str()));
         if (idx >= energies.size()) {
            energies.resize(idx + 1);
         }
         double value = getAttributeValueAsDouble(paramNode, "value", 0.0);
         double scale = getAttributeValueAsDouble(paramNode, "scale", 1.0);
         energies[idx] = value * scale;
      }
   }
   
   // Use default index values as in original code
   double indexL = -2.0;
   double indexH = -3.0;
   
   dynamic_cast<PiecewisePowerLaw*>(spec)->addParams(indexL, indexH, dNdEs, energies);
}
  
  void SourceFactory::checkRoiDist(double ra, double dec) const {
    std::vector<double> roiPars(m_observation.roiCuts().roiCone());
    astro::SkyDir roiDir(roiPars.at(0), roiPars.at(1));
    double radius = roiPars.at(2);
    astro::SkyDir srcDir(ra, dec);
    double sep = srcDir.difference(roiDir) * 180.0 / M_PI;

    if (sep > radius + 10) {
      m_formatter->warn() << "WARNING: Point source " << m_currentSrcName
			  << " lies " << sep
			  << " degrees from the ROI center at RA, Dec = "
			  << roiPars.at(0) << ", " << roiPars.at(1)
			  << ' ' << radius << std::endl;
    }
  }

} // namespace Likelihood~
