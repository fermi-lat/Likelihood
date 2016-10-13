/**
 * @file SourceModelBuilder.cxx
 * @brief Implementation for class to provide methods to write xml
 * files for the Likelihood package source models.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceModelBuilder.cxx,v 1.21 2015/03/03 21:41:31 jchiang Exp $
 */

#include <fstream>
#include <sstream>

#include <xercesc/dom/DOM.hpp>

#include "facilities/Util.h"

#include "optimizers/Dom.h"

#include "Likelihood/CompositeSource.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/DMFitFunction.h"
#include "Likelihood/FileFunction.h"
#include "Likelihood/MapCubeFunction2.h"
#include "Likelihood/RadialProfile.h"
#include "Likelihood/ScaleFactor.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceModel.h"
#include "Likelihood/SourceModelBuilder.h"

namespace Likelihood {

using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

SourceModelBuilder::SourceModelBuilder(const std::string &functionLibrary,
                                       const std::string &srcLibTitle) 
  : XmlBuilder(){
   m_srcLib = optimizers::Dom::createElement(m_doc, "source_library");
   if (functionLibrary != "") {
      xmlBase::Dom::addAttribute(m_srcLib, "function_library",
                                 functionLibrary);
   }
   xmlBase::Dom::addAttribute(m_srcLib, "title", srcLibTitle);
}


SourceModelBuilder::~SourceModelBuilder() {}

void SourceModelBuilder::addSourceModel(const SourceModel& srcModel) {   
   append_source_model(m_srcLib, srcModel);
}

void SourceModelBuilder::addSource(const Source & src) {
   append_source(m_srcLib, src);
}

void SourceModelBuilder::write(std::string xmlFile) {
   facilities::Util::expandEnvVar(&xmlFile);
   std::ofstream outFile(xmlFile.c_str());
   outFile << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
   xmlBase::Dom::prettyPrintElement(m_srcLib, outFile, std::string(""));
}

DOMElement * SourceModelBuilder::likelihoodSource(const Source & src) {
   DOMElement * srcElt = optimizers::Dom::createElement(m_doc, "source");
   xmlBase::Dom::addAttribute(srcElt, "name", src.getName());
   optimizers::Dom::appendChild(srcElt, spectralPart(src));
   switch ( src.srcType() ) {
   case Source::Composite:
     addComposite(srcElt, src);
     break;
   case Source::Point:
   case Source::Diffuse:
     addSpatialPart(srcElt, src);
     break;
   default:
     throw std::runtime_error("SourceModelBuilder::likelihoodSource: unknown source type.");     
   }
   
   return srcElt;
}

DOMElement * SourceModelBuilder::spectralPart(const Source & src) {
   const Source::FuncMap & srcFuncs = src.getSrcFuncs();

   DOMElement * specElt = optimizers::Dom::createElement(m_doc, "spectrum");
   Source::FuncMap::const_iterator itrFind = srcFuncs.find("Spectrum");
   if ( itrFind == srcFuncs.end() ) {
      throw std::runtime_error("SourceModelBuilder::spectralPart: no spectrum for source " + src.getName() );
   }
   const optimizers::Function* spectrum = itrFind->second;
   
   xmlBase::Dom::addAttribute(specElt, "type",spectrum->genericName());

   if (!src.use_edisp()) {
      // Explicitly set the apply_edisp attribute to "false".
      xmlBase::Dom::addAttribute(specElt, "apply_edisp", "false");
   }

   const FileFunction * fileFunc(0);
   const ScaleFactor * scaleFactor(dynamic_cast<const ScaleFactor *>(spectrum));
   if (scaleFactor != 0) {
      fileFunc = dynamic_cast<const FileFunction *>(scaleFactor->spectrum());
   } else {
      fileFunc = dynamic_cast<const FileFunction *>(spectrum);
   }
   if (fileFunc != 0) {
      xmlBase::Dom::addAttribute(specElt, "file", fileFunc->filename());
   }

   const DMFitFunction * dmFitFunc 
      = dynamic_cast<const DMFitFunction *>(spectrum);
   if (dmFitFunc != 0) {
      xmlBase::Dom::addAttribute(specElt, "file", dmFitFunc->filename());
   }

   // If the source spectrum has a scaling function, add the filename
   // as the scaling_file attribute.
   const optimizers::Function * scalingFunc = spectrum->scalingFunction();
   if (scalingFunc) {
      xmlBase::Dom::addAttribute(specElt, "scaling_file", 
                                 dynamic_cast<const FileFunction *>(scalingFunc)->filename());
   }

   optimizers::Function* nc_spectrum = const_cast<optimizers::Function*>(spectrum);
   nc_spectrum->appendParamDomElements(m_doc, specElt);
   return specElt;
}

void SourceModelBuilder::addSpatialPart(DOMElement * srcElt, const Source & src) {
   const Source::FuncMap & srcFuncs = src.getSrcFuncs();

   DOMElement * spatialElt 
      = optimizers::Dom::createElement(m_doc, "spatialModel");
   
   Source::FuncMap::const_iterator find_pos = srcFuncs.find("Position");
   Source::FuncMap::const_iterator find_spatial = srcFuncs.find("SpatialDist");

   if (find_pos != srcFuncs.end()) {
      xmlBase::Dom::addAttribute(srcElt, "type", "PointSource");
      xmlBase::Dom::addAttribute(spatialElt, "type", "SkyDirFunction");
      optimizers::Function* nc_pos = const_cast<optimizers::Function*>(find_pos->second);
      nc_pos->appendParamDomElements(m_doc, spatialElt);
      optimizers::Dom::appendChild(srcElt, spatialElt);
   } else if (find_spatial != srcFuncs.end()) {
      xmlBase::Dom::addAttribute(srcElt, "type", "DiffuseSource");
      const optimizers::Function* spatial = find_spatial->second;
      std::string type = spatial->genericName();
      xmlBase::Dom::addAttribute(spatialElt, "type", type);
      if (type == "SpatialMap") {
         std::string file = 
            dynamic_cast<const SpatialMap *>(spatial)->fitsFile();
         xmlBase::Dom::addAttribute(spatialElt, "file", file);
      } else if (type == "MapCubeFunction") {
         std::string file = 
           dynamic_cast<const MapCubeFunction2*>(spatial)->fitsFile();
         xmlBase::Dom::addAttribute(spatialElt, "file", file);
      } else if (type == "RadialProfile") {
         std::string file = 
            dynamic_cast<const RadialProfile *>(spatial)->templateFile();
         xmlBase::Dom::addAttribute(spatialElt, "file", file);
      }
      const DiffuseSource * diffuseSource(dynamic_cast<const DiffuseSource *>(&src));
      if (diffuseSource->mapBasedIntegral()) {
         xmlBase::Dom::addAttribute(spatialElt, "map_based_integral", "true");
      }
      optimizers::Function* nc_spatial = const_cast<optimizers::Function*>(spatial);

      nc_spatial->appendParamDomElements(m_doc, spatialElt);
      optimizers::Dom::appendChild(srcElt, spatialElt);
   }
}


void SourceModelBuilder::addComposite(DomElement * srcElt, const Source & src) {
  xmlBase::Dom::addAttribute(srcElt, "type", "CompositeSource");
  const CompositeSource* comp = dynamic_cast<const CompositeSource*>(&src);
  if ( comp == 0 ) {
    throw std::runtime_error("SourceModelBuilder::addComposite: source not a composite");
  }
  const std::string& compXmlFile = comp->xmlFile();
  DOMElement * compElt = optimizers::Dom::createElement(m_doc, "source_library");
  if ( compXmlFile.empty() ) {
    append_source_model(compElt,comp->sourceModel());    
  } else {
    xmlBase::Dom::addAttribute(compElt, "xmlFile", compXmlFile);    
  }
  optimizers::Dom::appendChild(srcElt, compElt);
}

void SourceModelBuilder::append_source(DomElement * parent, const Source & src) {
   optimizers::Dom::appendChild(parent, likelihoodSource(src));
} 

void SourceModelBuilder::append_source_model(DomElement * parent, const SourceModel& srcModel) {
   std::map<std::string, Source *>::const_iterator srcIt = srcModel.sources().begin();
   for ( ; srcIt != srcModel.sources().end(); srcIt++ ) {
      append_source(parent,*(srcIt->second));
   }
}

} // Likelihood
