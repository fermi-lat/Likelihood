/** 
 * @file SourceFactory.cxx
 * @brief Implementation for the SourceFactory class, which applies the
 * Prototype pattern to return clones of various gamma-ray Sources.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceFactory.cxx,v 1.26 2003/10/01 17:03:22 jchiang Exp $
 */

#include <sstream>

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include <xercesc/dom/DOM_Element.hpp>
#include <xercesc/dom/DOM_NodeList.hpp>

#include "facilities/Util.h"

#include "optimizers/Exception.h"
#include "optimizers/Dom.h"
#include "optimizers/FunctionFactory.h"

#include "Likelihood/Exception.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/ConstantValue.h"
#include "Likelihood/SpectrumFactory.h"
#include "Likelihood/SourceFactory.h"

namespace {

std::string rootPath() {
   std::string root_path;
   char *root = ::getenv("LIKELIHOODROOT");
   if (!root) {
      return std::string("..");
   } else {
      return std::string(root);
   }
}

} // unnamed namespace

namespace Likelihood {

SourceFactory::SourceFactory() {
}

SourceFactory::~SourceFactory() {
   std::map<std::string, Source *>::iterator it = m_prototypes.begin();
   for (; it != m_prototypes.end(); it++)
      delete it->second;
}

Source *SourceFactory::create(const std::string &name) throw(Exception) {
   if (!m_prototypes.count(name)) {
      std::ostringstream errorMessage;
      errorMessage << "SourceFactory::create: "
                   << "Cannot create Source named "
                   << name << ".\n";
      throw Exception(errorMessage.str());
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
      std::ostringstream errorMessage;
      errorMessage << "SourceFactory::addSource: A Source named "
                   << name << " already exists.\n";
      throw Exception(errorMessage.str());
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
                            optimizers::FunctionFactory &funcFactory)
   throw(Exception) {

   xml::XmlParser *parser = new xml::XmlParser();

   DOM_Document doc = parser->parse(xmlFile.c_str());

   if (doc == 0) { // xml file not parsed successfully
      std::ostringstream errorMessage;
      errorMessage << "SourceFactory::readXml: "
                   << "Input xml file, " << xmlFile 
                   << " not parsed successfully.\n";
      throw Exception(errorMessage.str());
   }

   DOM_Element source_library = doc.getDocumentElement();
   optimizers::Dom::checkTag(source_library, "source_library", 
                             "SourceFactory::readXml");

// Prepare the FunctionFactory object using the xml file specified in
// the source_library tag.
   std::string function_library 
      = xml::Dom::getAttribute(source_library, "function_library");
   facilities::Util::expandEnvVar(&function_library);
   try {
      funcFactory.readXml(function_library);
   } catch(optimizers::Exception &eObj) {
      std::cout << eObj.what() << std::endl;
// // Rethrow as a Likelihood::Exception.
//       throw Exception(eObj.what());
   }

// Loop through source elements, adding each as a Source object prototype.
   std::vector<DOM_Element> srcs;
   optimizers::Dom::getElements(source_library, "source", srcs);
   std::vector<DOM_Element>::const_iterator srcIt = srcs.begin();
   for ( ; srcIt != srcs.end(); srcIt++) {

// Get the type of this source which is either PointSource or
// DiffuseSource (CompositeSource pending)...
      std::string srcType = xml::Dom::getAttribute(*srcIt, "type");
// and its name.
      std::string srcName = xml::Dom::getAttribute(*srcIt, "name");

      std::cout << "Creating source named "
                << srcName << std::endl;

// Retrieve the spectrum and spatialModel elements (there should only
// be one of each).
      std::vector<DOM_Element> child;

      DOM_Element spectrum;
      try {
         optimizers::Dom::getElements(*srcIt, "spectrum", child);
         spectrum = child[0];
      } catch (optimizers::Exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      }

      optimizers::Dom::getElements(*srcIt, "spatialModel", child);
      DOM_Element spatialModel = child[0];

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

Source * SourceFactory::makePointSource(const DOM_Element &spectrum, 
                                        const DOM_Element &spatialModel,
                                        optimizers::FunctionFactory 
                                        &funcFactory) {
   std::string funcType = xml::Dom::getAttribute(spatialModel, "type");
   if (funcType != "SkyDirFunction") {
      std::ostringstream errorMessage;
      errorMessage << "SourceFactory::readXml: "
                   << "Trying to create a PointSource with "
                   << "a spatialModel of type "
                   << funcType << "." << std::endl;
      throw Exception(errorMessage.str());
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
   std::vector<DOM_Element> params;
   optimizers::Dom::getElements(spatialModel, "parameter", params);
   std::vector<DOM_Element>::const_iterator paramIt = params.begin();
   for ( ; paramIt != params.end(); paramIt++) {
      std::string name = xml::Dom::getAttribute(*paramIt, "name");
      if (name == "RA") 
         ra = ::atof( xml::Dom::getAttribute(*paramIt, "value").c_str() );
      if (name == "DEC") 
         dec = ::atof( xml::Dom::getAttribute(*paramIt, "value").c_str() );
   }

   Source *src = new PointSource();
//    dynamic_cast<PointSource *>(src)->setDir(ra, dec);
   src->setDir(ra, dec);

   try {
      setSpectrum(src, spectrum, funcFactory);
      return src;
   } catch (optimizers::Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   } catch (Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   } catch (...) {
      std::cerr << "other exception from setSpectrum" << std::endl;
   }
   return 0;
}

Source * SourceFactory::makeDiffuseSource(const DOM_Element &spectrum, 
                                          const DOM_Element &spatialModel,
                                          optimizers::FunctionFactory 
                                          &funcFactory) {
   std::string type = xml::Dom::getAttribute(spatialModel, "type");
   optimizers::Function *spatialDist = funcFactory.create(type);
   std::vector<DOM_Element> params;
   optimizers::Dom::getElements(spatialModel, "parameter", params);
   std::vector<DOM_Element>::const_iterator paramIt = params.begin();
   for ( ; paramIt != params.end(); paramIt++) {
      optimizers::Parameter parameter;
      parameter.extractDomData(*paramIt);
      spatialDist->setParam(parameter);
   }
   if (type == "SpatialMap") {
      std::string fitsFile 
         = xml::Dom::getAttribute(spatialModel, "file");
      facilities::Util::expandEnvVar(&fitsFile);
      dynamic_cast<SpatialMap *>(spatialDist)->readFitsFile(fitsFile);
   }
   Source *src;
   try {
      src = new DiffuseSource(spatialDist);
      setSpectrum(src, spectrum, funcFactory);
      return src;
   } catch (optimizers::Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   } catch (Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   }
   return 0;
}

void SourceFactory::setSpectrum(Source *src, const DOM_Element &spectrum, 
                                optimizers::FunctionFactory &funcFactory) {
   std::string type = xml::Dom::getAttribute(spectrum, "type");
   optimizers::Function *spec = funcFactory.create(type);

// Fetch the parameter elements (if any).
   std::vector<DOM_Element> params;
   optimizers::Dom::getElements(spectrum, "parameter", params);
   if (params.size() > 0) {
      std::vector<DOM_Element>::const_iterator paramIt = params.begin();
      for ( ; paramIt != params.end(); paramIt++) {
         optimizers::Parameter parameter;
//         optimizers::Dom::readParamData(*paramIt, parameter);
         parameter.extractDomData(*paramIt);
         spec->setParam(parameter);
      }
   }  
   src->setSpectrum(spec);
}

} // namespace Likelihood
