/** 
 * @file SourceFactory.cxx
 * @brief Implementation for the SourceFactory class, which applies the
 * Prototype pattern to return clones of various gamma-ray Sources.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceFactory.cxx,v 1.33 2004/06/05 15:22:15 jchiang Exp $
 */

#include <xercesc/util/XercesDefs.hpp>

#include "xml/Dom.h"
#include "xml/XmlParser.h"

#include "facilities/Util.h"

#include "optimizers/Exception.h"
#include "optimizers/FunctionFactory.h"

#include "Likelihood/Exception.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/SpectrumFactory.h"
#include "Likelihood/SourceFactory.h"

namespace {
   std::string rootPath() {
      char *root = ::getenv("LIKELIHOODROOT");
      if (!root) {
         return std::string("..");
      } else {
         return std::string(root);
      }
   }
} // unnamed namespace

namespace Likelihood {

XERCES_CPP_NAMESPACE_USE

SourceFactory::SourceFactory(bool verbose) : m_verbose(verbose) {
}

SourceFactory::~SourceFactory() {
   std::map<std::string, Source *>::iterator it = m_prototypes.begin();
   for (; it != m_prototypes.end(); it++)
      delete it->second;
}

Source *SourceFactory::create(const std::string &name) throw(Exception) {
   if (!m_prototypes.count(name)) {
      std::string errorMessage 
         = "SourceFactory::create:\nCannot create Source named " + name + ".";
      throw Exception(errorMessage);
   }
   return m_prototypes[name]->clone();
}

void SourceFactory::addSource(const std::string &name, Source* src, 
                              bool fromClone) 
   throw(Exception) {
   if (!m_prototypes.count(name)) {
      if (fromClone) {
         m_prototypes[name] = src->clone();
      } else {
         m_prototypes[name] = src;
      }
   } else {
      std::string errorMessage = "SourceFactory::addSource:\n A Source named "
         + name + " already exists.";
      throw Exception(errorMessage);
   }
}

void SourceFactory::replaceSource(Source* src, bool fromClone) {
   if (m_prototypes.count(src->getName())) {
      if (fromClone) {
         m_prototypes[src->getName()] = src->clone();
      } else {
         m_prototypes[src->getName()] = src;
      }
   } else {
      std::cerr << "SourceFactory::replaceSource: A Source named "
                << src->getName() << " does not yet exist.\n"
                << "Adding it instead. "
                << std::endl;
      addSource(src->getName(), src, fromClone);
   }
}

void SourceFactory::readXml(const std::string &xmlFile,
                            optimizers::FunctionFactory &funcFactory,
                            bool requireExposure)
   throw(Exception) {
   m_requireExposure = requireExposure;

   xml::XmlParser * parser = new xml::XmlParser();

   DOMDocument * doc = parser->parse(xmlFile.c_str());

   if (doc == 0) { // xml file not parsed successfully
      std::string errorMessage = "SourceFactory::readXml:\nInput xml file, "
         + xmlFile + " not parsed successfully.";
      throw Exception(errorMessage);
   }

// Direct Xerces API call...still available in Xerces 2.4.0:
   DOMElement * source_library = doc->getDocumentElement();
   if (!xml::Dom::checkTagName(source_library, "source_library")) {
      throw Exception("SourceFactory::readXml:\nsource_library not found in "
         + xmlFile);
   }

// Prepare the FunctionFactory object using the xml file specified in
// the source_library tag.
   std::string function_library 
      = xml::Dom::getAttribute(source_library, "function_library");
   if (function_library.find("xml") != std::string::npos) {
      facilities::Util::expandEnvVar(&function_library);
      try {
         funcFactory.readXml(function_library);
      } catch(optimizers::Exception &eObj) {
         std::cout << eObj.what() << std::endl;
      }
   }

// Loop through source elements, adding each as a Source object prototype.
   std::vector<DOMElement *> srcs;
   xml::Dom::getChildrenByTagName(source_library, "source", srcs);
   std::vector<DOMElement *>::const_iterator srcIt = srcs.begin();
   for ( ; srcIt != srcs.end(); srcIt++) {

// Get the type of this source which is either PointSource or
// DiffuseSource (CompositeSource pending)...
      std::string srcType = xml::Dom::getAttribute(*srcIt, "type");
// and its name.
      std::string srcName = xml::Dom::getAttribute(*srcIt, "name");

      if (m_verbose) std::cout << "Creating source named "
                               << srcName << std::endl;

// Retrieve the spectrum and spatialModel elements (there should only
// be one of each).
      std::vector<DOMElement *> child;

      DOMElement * spectrum;
      try {
         xml::Dom::getChildrenByTagName(*srcIt, "spectrum", child);
         spectrum = child[0];
      } catch (optimizers::Exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      }

      xml::Dom::getChildrenByTagName(*srcIt, "spatialModel", child);
      DOMElement * spatialModel = child[0];

// The processing logic for the spatialModel depends on the source
// type, so we must consider each case individually:
      Source *src = 0;
      if (srcType == "PointSource") {
         src = makePointSource(spectrum, spatialModel, funcFactory);
      } else if (srcType == "DiffuseSource") {
         src = makeDiffuseSource(spectrum, spatialModel, funcFactory);
      }

// Add the source to the vector of prototypes.
      if (src != 0) {
         src->setName(srcName);
         addSource(srcName, src);
         delete src;
      }
   }
   delete parser;
}

void SourceFactory::fetchSrcNames(std::vector<std::string> &srcNames) {
   if (!srcNames.empty()) srcNames.clear();
   std::map<std::string, Source *>::const_iterator it = m_prototypes.begin();
   for (; it != m_prototypes.end(); it++)
      srcNames.push_back(it->first);
}

Source * SourceFactory::makePointSource(const DOMElement * spectrum, 
                                        const DOMElement * spatialModel,
                                        optimizers::FunctionFactory 
                                        &funcFactory) {
   std::string funcType = xml::Dom::getAttribute(spatialModel, "type");
   if (funcType != "SkyDirFunction") {
      std::string errorMessage = std::string("SourceFactory::readXml:\n") 
         + "Trying to create a PointSource with a spatialModel of type "
         + funcType + ".";
      throw Exception(errorMessage);
   }
// For v1r0 and prior versions, setting the direction of the 
// PointSource object forces the exposure to be calculated and
// so requires the ROI cuts and spacecraft data to have been
// specified.  *If* this is desired behavior, then some checks
// should be made and perhaps exceptions thrown if the ROI and
// spacecraft info are not available.
//
// Extract the (RA, Dec) from the parameter elements.
   double ra(0), dec(0);
   std::vector<DOMElement *> params;
   xml::Dom::getChildrenByTagName(spatialModel, "parameter", params);
   std::vector<DOMElement *>::const_iterator paramIt = params.begin();
   for ( ; paramIt != params.end(); paramIt++) {
      std::string name = xml::Dom::getAttribute(*paramIt, "name");
      if (name == "RA") 
         ra = ::atof( xml::Dom::getAttribute(*paramIt, "value").c_str() );
      if (name == "DEC") 
         dec = ::atof( xml::Dom::getAttribute(*paramIt, "value").c_str() );
   }

   Source *src = new PointSource();
//    dynamic_cast<PointSource *>(src)->setDir(ra, dec);
   bool updateExposure(true);
   src->setDir(ra, dec, updateExposure, m_verbose);

   try {
      setSpectrum(src, spectrum, funcFactory);
      return src;
   } catch (optimizers::Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   } catch (Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   } catch (std::exception &eObj) {
      std::cout << eObj.what() << std::endl;
   } catch (...) {
      std::cerr << "Unexpected exception from SourceFactory::setSpectrum" 
                << std::endl;
   }
   return 0;
}

Source * SourceFactory::makeDiffuseSource(const DOMElement * spectrum, 
                                          const DOMElement * spatialModel,
                                          optimizers::FunctionFactory 
                                          &funcFactory) {
   std::string type = xml::Dom::getAttribute(spatialModel, "type");
   optimizers::Function *spatialDist = funcFactory.create(type);
   std::vector<DOMElement *> params;
   xml::Dom::getChildrenByTagName(spatialModel, "parameter", params);
   std::vector<DOMElement *>::const_iterator paramIt = params.begin();
   for ( ; paramIt != params.end(); paramIt++) {
      optimizers::Parameter parameter;
      parameter.extractDomData(*paramIt);
      spatialDist->setParam(parameter);
   }
   if (type == "SpatialMap") {
      std::string fitsFile 
         = xml::Dom::getAttribute(spatialModel, "file");
      dynamic_cast<SpatialMap *>(spatialDist)->readFitsFile(fitsFile);
   }
   Source *src;
   try {
      src = new DiffuseSource(spatialDist, m_requireExposure);
      setSpectrum(src, spectrum, funcFactory);
      return src;
   } catch (optimizers::Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   } catch (Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   }
   return 0;
}

void SourceFactory::setSpectrum(Source *src, const DOMElement * spectrum, 
                                optimizers::FunctionFactory &funcFactory) {
   std::string type = xml::Dom::getAttribute(spectrum, "type");
   optimizers::Function *spec = funcFactory.create(type);

// Fetch the parameter elements (if any).
   std::vector<DOMElement *> params;
   xml::Dom::getChildrenByTagName(spectrum, "parameter", params);
   if (params.size() > 0) {
      std::vector<DOMElement *>::const_iterator paramIt = params.begin();
      for ( ; paramIt != params.end(); paramIt++) {
         optimizers::Parameter parameter;
         parameter.extractDomData(*paramIt);
         spec->setParam(parameter);
      }
   }  
   src->setSpectrum(spec);
}

} // namespace Likelihood
