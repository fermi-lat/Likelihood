/**
 * @file SourceModelBuilder.cxx
 * @brief Implementation for class to provide methods to write xml
 * files for the Likelihood package source models.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/SourceModelBuilder.cxx,v 1.16 2011/11/28 09:58:02 cohen Exp $
 */

#include <fstream>
#include <sstream>

#include <xercesc/dom/DOM.hpp>

#include "facilities/Util.h"

#include "optimizers/Dom.h"

#include "Likelihood/DiffuseSource.h"
#include "Likelihood/DMFitFunction.h"
#include "Likelihood/FileFunction.h"
#include "Likelihood/MapCubeFunction2.h"
#include "Likelihood/RadialProfile.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceModelBuilder.h"

namespace Likelihood {

using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

SourceModelBuilder::SourceModelBuilder(const std::string &functionLibrary,
                                       const std::string &srcLibTitle) 
   : XmlBuilder() {
   m_srcLib = optimizers::Dom::createElement(m_doc, "source_library");
   if (functionLibrary != "") {
      xmlBase::Dom::addAttribute(m_srcLib, "function_library",
                                 functionLibrary);
   }
   xmlBase::Dom::addAttribute(m_srcLib, "title", srcLibTitle);
}

SourceModelBuilder::~SourceModelBuilder() {}

void SourceModelBuilder::addSource(Source & src) {
   optimizers::Dom::appendChild(m_srcLib, likelihoodSource(src));
}

void SourceModelBuilder::write(std::string xmlFile) {
   facilities::Util::expandEnvVar(&xmlFile);
   std::ofstream outFile(xmlFile.c_str());
   outFile << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
   xmlBase::Dom::prettyPrintElement(m_srcLib, outFile, std::string(""));
}

DOMElement * SourceModelBuilder::likelihoodSource(Source & src) {
   DOMElement * srcElt = optimizers::Dom::createElement(m_doc, "source");
   xmlBase::Dom::addAttribute(srcElt, "name", src.getName());
   optimizers::Dom::appendChild(srcElt, spectralPart(src));
   addSpatialPart(srcElt, src);
   return srcElt;
}

DOMElement * SourceModelBuilder::spectralPart(Source & src) {
   Source::FuncMap & srcFuncs = src.getSrcFuncs();

   DOMElement * specElt = optimizers::Dom::createElement(m_doc, "spectrum");
   xmlBase::Dom::addAttribute(specElt, "type",
                              srcFuncs["Spectrum"]->genericName());

   FileFunction * fileFunc(dynamic_cast<FileFunction *>(srcFuncs["Spectrum"]));
   if (fileFunc != 0) {
      xmlBase::Dom::addAttribute(specElt, "file", fileFunc->filename());
   }

   DMFitFunction * dmFitFunc 
      = dynamic_cast<DMFitFunction *>(srcFuncs["Spectrum"]);
   if (dmFitFunc != 0) {
      xmlBase::Dom::addAttribute(specElt, "file", dmFitFunc->filename());
   }

   srcFuncs["Spectrum"]->appendParamDomElements(m_doc, specElt);
   return specElt;
}

void SourceModelBuilder::addSpatialPart(DOMElement * srcElt, Source & src) {
   Source::FuncMap & srcFuncs = src.getSrcFuncs();

   DOMElement * spatialElt 
      = optimizers::Dom::createElement(m_doc, "spatialModel");
   if (srcFuncs.count("Position")) {
      xmlBase::Dom::addAttribute(srcElt, "type", "PointSource");
      xmlBase::Dom::addAttribute(spatialElt, "type", "SkyDirFunction");
      srcFuncs["Position"]->appendParamDomElements(m_doc, spatialElt);
      optimizers::Dom::appendChild(srcElt, spatialElt);
   } else if (srcFuncs.count("SpatialDist")) {
      xmlBase::Dom::addAttribute(srcElt, "type", "DiffuseSource");
      std::string type = srcFuncs["SpatialDist"]->genericName();
      xmlBase::Dom::addAttribute(spatialElt, "type", type);
      if (type == "SpatialMap") {
         std::string file = 
            dynamic_cast<SpatialMap *>(srcFuncs["SpatialDist"])->fitsFile();
         xmlBase::Dom::addAttribute(spatialElt, "file", file);
      } else if (type == "MapCubeFunction") {
         std::string file = 
           dynamic_cast<MapCubeFunction2*>(srcFuncs["SpatialDist"])->fitsFile();
         xmlBase::Dom::addAttribute(spatialElt, "file", file);
      } else if (type == "RadialProfile") {
         std::string file = 
            dynamic_cast<RadialProfile *>(srcFuncs["SpatialDist"])->templateFile();
         xmlBase::Dom::addAttribute(spatialElt, "file", file);
      }
      DiffuseSource * diffuseSource
         = dynamic_cast<DiffuseSource *>(srcFuncs["SpatialDist"]);
      if (diffuseSource !=0 && diffuseSource->mapBasedIntegral()) {
         xmlBase::Dom::addAttribute(spatialElt, "map_based_integral", "true");
      }
      srcFuncs["SpatialDist"]->appendParamDomElements(m_doc, spatialElt);
      optimizers::Dom::appendChild(srcElt, spatialElt);
   }
}

} // namespace Likelihood
