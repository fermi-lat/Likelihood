/**
 * @file SourceModelBuilder.cxx
 * @brief Implementation for class to provide methods to write xml
 * files for the Likelihood package source models.
 * @author J. Chiang
 *
 * $Header$
 */

#include <fstream>
#include <sstream>

#include "facilities/Util.h"

#include "Likelihood/SpatialMap.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceModelBuilder.h"

namespace {
   DomElement * createElement(DomDocument * doc, const std::string & name) {
      DomElement * elt = new DOM_Element();
      *elt = doc->createElement(name.c_str());
      return elt;
   }
}

namespace Likelihood {

SourceModelBuilder::SourceModelBuilder(const std::string &functionLibrary,
                                       const std::string &srcLibTitle) 
   : XmlBuilder() {
   m_srcLib = ::createElement(m_doc, "source_library");
   if (functionLibrary != "") {
      xml::Dom::addAttribute(*m_srcLib, "function_library", functionLibrary);
   }
   xml::Dom::addAttribute(*m_srcLib, "title", srcLibTitle);
}

SourceModelBuilder::~SourceModelBuilder() {
   delete m_srcLib;
}

void SourceModelBuilder::addSource(Source & src) {
   m_srcLib->appendChild(*likelihoodSource(src));
}

void SourceModelBuilder::write(std::string xmlFile) {
   facilities::Util::expandEnvVar(&xmlFile);
   std::ofstream outFile(xmlFile.c_str());
   outFile << "<?xml version='1.0' standalone='no'?>\n"
           << "<!DOCTYPE source_library SYSTEM "
           << "\"$(LIKELIHOODROOT)/xml/A1_Sources.dtd\" >\n";

   xml::Dom::prettyPrintElement(*m_srcLib, outFile, std::string(""));
}

DomElement * SourceModelBuilder::likelihoodSource(Source & src) {
   DomElement * srcElt = ::createElement(m_doc, "source");
   xml::Dom::addAttribute(*srcElt, "name", src.getName());
   srcElt->appendChild(*spectralPart(src));
   addSpatialPart(srcElt, src);
   return srcElt;
}

DomElement * SourceModelBuilder::spectralPart(Source & src) {
   Source::FuncMap & srcFuncs = src.getSrcFuncs();

   DomElement * specElt = ::createElement(m_doc, "spectrum");
   xml::Dom::addAttribute(*specElt, "type",
                          srcFuncs["Spectrum"]->genericName());
   srcFuncs["Spectrum"]->appendParamDomElements(*m_doc, *specElt);
   return specElt;
}

void SourceModelBuilder::addSpatialPart(DomElement * srcElt, Source & src) {
   Source::FuncMap & srcFuncs = src.getSrcFuncs();

   DomElement * spatialElt = ::createElement(m_doc, "spatialModel");
   if (srcFuncs.count("Position")) {
      xml::Dom::addAttribute(*srcElt, "type", "PointSource");
      xml::Dom::addAttribute(*spatialElt, "type", "SkyDirFunction");
      srcFuncs["Position"]->appendParamDomElements(*m_doc, *spatialElt);
      srcElt->appendChild(*spatialElt);
   } else if (srcFuncs.count("SpatialDist")) {
      xml::Dom::addAttribute(*srcElt, "type", "DiffuseSource");
      std::string type = srcFuncs["SpatialDist"]->genericName();
      xml::Dom::addAttribute(*spatialElt, "type", type);
      if (type == "SpatialMap") {
         std::string file = 
            dynamic_cast<SpatialMap *>(srcFuncs["SpatialDist"])->fitsFile();
         xml::Dom::addAttribute(*spatialElt, "file", file);
      }
      srcFuncs["SpatialDist"]->appendParamDomElements(*m_doc, *spatialElt);
      srcElt->appendChild(*spatialElt);
   }
}

} // namespace Likelihood
