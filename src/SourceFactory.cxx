/** 
 * @file SourceFactory.cxx
 * @brief Implementation for the SourceFactory class, which applies the
 * Prototype pattern to return clones of various gamma-ray Sources.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceFactory.cxx,v 1.72 2011/11/22 01:50:01 jchiang Exp $
 */

#include <xercesc/util/XercesDefs.hpp>

#include "xmlBase/Dom.h"
#include "xmlBase/XmlParser.h"

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "optimizers/Exception.h"
#include "optimizers/FunctionFactory.h"

#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Event.h"
#include "Likelihood/Exception.h"
#include "Likelihood/FileFunction.h"
#include "Likelihood/DMFitFunction.h"
#include "Likelihood/MapBase.h"
#include "Likelihood/Observation.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/RadialProfile.h"
#include "Likelihood/SourceFactory.h"
#include "Likelihood/XmlParser.h"

namespace Likelihood {

//XERCES_CPP_NAMESPACE_USE
using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

SourceFactory::SourceFactory(const Observation & observation, bool verbose) 
   : m_verbose(verbose), m_observation(observation), 
     m_formatter(new st_stream::StreamFormatter("SourceFactory", "", 2)) {
}

SourceFactory::~SourceFactory() {
   std::map<std::string, Source *>::iterator it = m_prototypes.begin();
   for (; it != m_prototypes.end(); it++) {
      delete it->second;
   }
   delete m_formatter;
}

Source * SourceFactory::create(const std::string &name) {
   if (!m_prototypes.count(name)) {
      std::string errorMessage 
         = "SourceFactory::create:\nCannot create Source named " + name + ".";
      throw Exception(errorMessage);
   }
   return m_prototypes[name]->clone();
}

Source * SourceFactory::releaseSource(const std::string &name) {
   if (!m_prototypes.count(name)) {
      std::string errorMessage 
         = "SourceFactory::releaseSource:\nNo Source named " + name + ".";
      throw Exception(errorMessage);
   }
   Source * src = m_prototypes[name];
   m_prototypes.erase(name);
   return src;
}

void SourceFactory::addSource(const std::string &name, Source* src, 
                              bool fromClone) {
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
      m_formatter->info() << "SourceFactory::replaceSource: A Source named "
                          << src->getName() << " does not yet exist.\n"
                          << "Adding it instead. "
                          << std::endl;
      addSource(src->getName(), src, fromClone);
   }
}

void SourceFactory::readXml(const std::string & xmlFile,
                            optimizers::FunctionFactory & funcFactory,
                            bool requireExposure,
                            bool addPointSources, 
                            bool loadMaps) {
   m_requireExposure = requireExposure;

   xmlBase::XmlParser * parser = XmlParser_instance();

   DOMDocument * doc = parser->parse(xmlFile.c_str());

   if (doc == 0) { // xml file not parsed successfully
      std::string errorMessage = "SourceFactory::readXml:\nInput xml file, "
         + xmlFile + " not parsed successfully.";
      throw Exception(errorMessage);
   }

// Direct Xerces API call...still available in Xerces 2.6.0:
   DOMElement * source_library = doc->getDocumentElement();
   if (!xmlBase::Dom::checkTagName(source_library, "source_library")) {
      throw Exception("SourceFactory::readXml:\nsource_library not found in "
         + xmlFile);
   }

// Prepare the FunctionFactory object using the xml file specified in
// the source_library tag.
   std::string function_library 
      = xmlBase::Dom::getAttribute(source_library, "function_library");
   if (function_library.find("xml") != std::string::npos) {
      facilities::Util::expandEnvVar(&function_library);
      try {
         funcFactory.readXml(function_library);
      } catch (optimizers::Exception &eObj) {
         m_formatter->err() << eObj.what() << std::endl;
         throw;
      }
   }

// Loop through source elements, adding each as a Source object prototype.
   std::vector<DOMElement *> srcs;
   xmlBase::Dom::getChildrenByTagName(source_library, "source", srcs);
   std::vector<DOMElement *>::const_iterator srcIt = srcs.begin();
   for ( ; srcIt != srcs.end(); srcIt++) {

// Get the type of this source which is either PointSource or
// DiffuseSource (CompositeSource pending)...
      std::string srcType = xmlBase::Dom::getAttribute(*srcIt, "type");
// and its name.
      std::string srcName = xmlBase::Dom::getAttribute(*srcIt, "name");

      m_currentSrcName = srcName;

      m_formatter->info(3) << "Creating source named "
                           << srcName << std::endl;

// Retrieve the spectrum and spatialModel elements (there should only
// be one of each).
      std::vector<DOMElement *> child;

      DOMElement * spectrum;
      try {
         xmlBase::Dom::getChildrenByTagName(*srcIt, "spectrum", child);
         if (child.size() != 1) {
            std::ostringstream message;
            message << "Error parsing xml model file: \n"
                    << xmlFile << "\n"
                    << "for source " << srcName << "\n"
                    << "Missing spectral model component.";
            throw Exception(message.str());
         }
         spectrum = child[0];
      } catch (optimizers::Exception &eObj) {
         m_formatter->err() << eObj.what() << std::endl;
         throw;
      }
      
      xmlBase::Dom::getChildrenByTagName(*srcIt, "spatialModel", child);
      if (child.size() != 1) {
         std::ostringstream message;
         message << "Error parsing xml model file: \n"
                 << xmlFile << "\n"
                 << "for source " << srcName << ".\n"
                 << "Missing spatial model component.\n"
                 << "Please check that you are using the correct xml format.";
         throw Exception(message.str());
      }
      DOMElement * spatialModel = child[0];

// The processing logic for the spatialModel depends on the source
// type, so we must consider each case individually:
      Source * src = 0;
      if (addPointSources && srcType == "PointSource") {
         src = makePointSource(spectrum, spatialModel, funcFactory);
      } else if (srcType == "DiffuseSource") {
         src = makeDiffuseSource(spectrum, spatialModel, funcFactory,
                                 loadMaps);
      }

// Add the source to the vector of prototypes.
      if (src != 0) {
         src->setName(srcName);
         addSource(srcName, src);
         delete src;
      }
   }

   delete doc;
}

void SourceFactory::fetchSrcNames(std::vector<std::string> &srcNames) {
   if (!srcNames.empty()) srcNames.clear();
   std::map<std::string, Source *>::const_iterator it = m_prototypes.begin();
   for (; it != m_prototypes.end(); it++)
      srcNames.push_back(it->first);
}

Source * SourceFactory::
makePointSource(const DOMElement * spectrum, 
                const DOMElement * spatialModel,
                optimizers::FunctionFactory & funcFactory) {
   std::string funcType = xmlBase::Dom::getAttribute(spatialModel, "type");
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
   xmlBase::Dom::getChildrenByTagName(spatialModel, "parameter", params);
   std::vector<DOMElement *>::const_iterator paramIt = params.begin();
   for ( ; paramIt != params.end(); paramIt++) {
      std::string name = xmlBase::Dom::getAttribute(*paramIt, "name");
      if (name == "RA") {
         ra = ::atof( xmlBase::Dom::getAttribute(*paramIt, "value").c_str() );
      }
      if (name == "DEC") {
         dec = ::atof( xmlBase::Dom::getAttribute(*paramIt, "value").c_str() );
      }
   }

   Source * src(0);

   checkRoiDist(ra, dec);

   if (m_requireExposure) {
      src = new PointSource(ra, dec, m_observation, m_verbose);
   } else { // for BinnedLikelihood, skip the exposure calculation
      src = new PointSource();
      dynamic_cast<PointSource *>(src)->setDir(ra, dec, m_requireExposure,
                                               m_verbose);
   }

   try {
      setSpectrum(src, spectrum, funcFactory);
      return src;
   } catch (std::exception &eObj) {
      m_formatter->err() << eObj.what() << std::endl;
      throw;
   } catch (...) {
      m_formatter->err() << "Unexpected exception from "
                         << "SourceFactory::setSpectrum" 
                         << std::endl;
      throw;
   }
   return 0;
}

Source * SourceFactory::
makeDiffuseSource(const DOMElement * spectrum, 
                  const DOMElement * spatialModel,
                  optimizers::FunctionFactory & funcFactory,
                  bool loadMap) {
   std::string type = xmlBase::Dom::getAttribute(spatialModel, "type");
   optimizers::Function * spatialDist = funcFactory.create(type);
   std::vector<DOMElement *> params;
   xmlBase::Dom::getChildrenByTagName(spatialModel, "parameter", params);
   std::vector<DOMElement *>::const_iterator paramIt = params.begin();
   for ( ; paramIt != params.end(); paramIt++) {
      std::string name = xmlBase::Dom::getAttribute(*paramIt, "name");
      spatialDist->parameter(name).extractDomData(*paramIt);
   }
   if (type == "SpatialMap" || type == "MapCubeFunction") {
      std::string fitsFile 
         = xmlBase::Dom::getAttribute(spatialModel, "file");
      dynamic_cast<MapBase *>(spatialDist)->readFitsFile(fitsFile, "", loadMap);
   } else if (type == "RadialProfile") {
      std::string tpl_file(xmlBase::Dom::getAttribute(spatialModel, "file"));
      dynamic_cast<RadialProfile *>(spatialDist)->readTemplateFile(tpl_file);
   }
   Source * src;
   try {
      src = new DiffuseSource(spatialDist, m_observation, m_requireExposure);
      setSpectrum(src, spectrum, funcFactory);
      delete spatialDist;
      return src;
   } catch (std::exception &eObj) {
      m_formatter->err() << eObj.what() << std::endl;
      throw;
   } catch (...) {
      m_formatter->err() << "Unexpected exception from "
                         << "SourceFactory::setSpectrum" 
                         << std::endl;
      throw;
   }
   return 0;
}

void SourceFactory::setSpectrum(Source * src, const DOMElement * spectrum, 
                                optimizers::FunctionFactory & funcFactory) {
   std::string type = xmlBase::Dom::getAttribute(spectrum, "type");
   optimizers::Function * spec = funcFactory.create(type);

// Fetch the parameter elements (if any).
   std::vector<DOMElement *> params;
   xmlBase::Dom::getChildrenByTagName(spectrum, "parameter", params);
   if (params.size() > 0) {
      std::vector<DOMElement *>::const_iterator paramIt = params.begin();
      for ( ; paramIt != params.end(); paramIt++) {
         std::string name = xmlBase::Dom::getAttribute(*paramIt, "name");
         spec->parameter(name).extractDomData(*paramIt);
      }
   }

// If FileFunction, read in the data:
   if (type == "FileFunction") {
      std::string filename = xmlBase::Dom::getAttribute(spectrum, "file");
      dynamic_cast<FileFunction *>(spec)->readFunction(filename);
   }
   if (type == "DMFitFunction") {
      std::string filename = xmlBase::Dom::getAttribute(spectrum, "file");
      dynamic_cast<DMFitFunction *>(spec)->readFunction(filename);
   }

   src->setSpectrum(spec);
   delete spec;
}

void SourceFactory::checkRoiDist(double ra, double dec) const {
   std::vector<double> roiPars(m_observation.roiCuts().roiCone());
   astro::SkyDir roiDir(roiPars.at(0), roiPars.at(1));
   double radius(roiPars.at(2));
   astro::SkyDir srcDir(ra, dec);
   double sep(srcDir.difference(roiDir)*180./M_PI);
   if (sep > radius + 10) {
      m_formatter->warn() << "WARNING: Point source " << m_currentSrcName
                          << " lies " << sep 
                          << " degrees from the ROI center at RA, Dec = " 
                          << roiPars.at(0) << ", " << roiPars.at(1)
                          << std::endl;
   }
}

} // namespace Likelihood
