/**
 * @file FluxModel.cxx
 * @brief Implementation for class to provide methods to write flux package
 * style xml files.
 * @author J. Chiang
 *
 * $Header$
 */

#include <string>
#include <utility>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "facilities/Util.h"

#include <xercesc/dom/DOM_Document.hpp>
#include <xercesc/dom/DOM_Element.hpp>
#include <xercesc/dom/DOM_NodeList.hpp>
#include <xercesc/dom/DOM_DOMException.hpp>

#include "xml/XmlParser.h"
#include "xml/Dom.h"

#include "optimizers/Function.h"

#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/TrapQuad.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/Source.h"
#include "Likelihood/FluxModel.h"

namespace {
   DomDocument * createDocument() {
      DomDocument * doc = new DOM_Document();
      *doc = DOM_Document::createDocument();
      return doc;
   }

   DomElement * createElement(DomDocument * doc, const std::string & name) {
      DomElement * elt = new DOM_Element();
      *elt = doc->createElement(name.c_str());
      return elt;
   }

   std::string basename(const std::string &path) {
      std::vector<std::string> names;
      facilities::Util::stringTokenize(path, "\\/", names);
      return *(names.end() - 1);
   }
}

namespace Likelihood {

FluxModel::FluxModel() {

   m_parser = new xml::XmlParser();

   m_doc = ::createDocument();

   m_srcLib = ::createElement(m_doc, "source_library");
   xml::Dom::addAttribute(*m_srcLib, "title", "Likelihood_model");

   m_allSrcsElt = ::createElement(m_doc, "source");

   makeEnergyGrid();
}

FluxModel::~FluxModel() {
   delete m_allSrcsElt;
   delete m_srcLib;
   delete m_doc;
   delete m_parser;
}

void FluxModel::addSource(Source & src) {
   DomElement * srcElt = fluxSource(src);

   if (srcElt) {
      m_srcLib->appendChild(*srcElt);
      
      DomElement * nestedSrc = ::createElement(m_doc, "nestedSource");
      std::string name = src.getName();
      addUnderscores(name);
      xml::Dom::addAttribute(*nestedSrc, "sourceRef", name.c_str());
      m_allSrcsElt->appendChild(*nestedSrc);
   }
}

void FluxModel::write(std::string xmlFile) {
   facilities::Util::expandEnvVar(&xmlFile);

   std::string name = std::string("all_in_") + ::basename(xmlFile.c_str());
   xml::Dom::addAttribute(*m_allSrcsElt, "name", name);
   m_srcLib->appendChild(*m_allSrcsElt);

   std::ofstream outFile(xmlFile.c_str());
   xml::Dom::prettyPrintElement(*m_srcLib, outFile, std::string(""));
}

DomElement * FluxModel::fluxSource(Source & src) {

   DomElement * srcElt = ::createElement(m_doc, "source");
   std::string name = src.getName();
   addUnderscores(name);
   xml::Dom::addAttribute(*srcElt, "name", name);

   std::string sourceType;
   getSourceType(src, sourceType);

   Source::FuncMap & srcFuncs = src.getSrcFuncs();
   TrapQuad fluxIntegral(srcFuncs["Spectrum"]);
   if (sourceType == "PointSource" || sourceType == "Isotropic") {
      xml::Dom::addAttribute(*srcElt, std::string("flux"),
                             fluxIntegral.integral(m_energies)/1e-4);
      DomElement * specElt = gammaSpectrum(*srcFuncs["Spectrum"]);
      if (sourceType == "PointSource") {
         specElt->appendChild(*srcDirection(*srcFuncs["Position"]));
      } else {
         specElt->appendChild(*solidAngle(-0.4, 1.0));
      }
      srcElt->appendChild(*specElt);
      return srcElt;
   } else if (sourceType == "GalDiffuse") {
      srcElt->appendChild(*galDiffuse(src));
      return srcElt;
   }
   return 0;
}

void FluxModel::getSourceType(Source &src, std::string & srcType) {
   Source::FuncMap & srcFuncs = src.getSrcFuncs();

   if (srcFuncs.count("Position")) {
      srcType = "PointSource";
   } else if (srcFuncs.count("SpatialDist") 
              && srcFuncs["SpatialDist"]->genericName() == "ConstantValue") {
      srcType = "Isotropic";
   } else if (srcFuncs.count("SpatialDist")
              && srcFuncs["SpatialDist"]->genericName() == "SpatialMap") {
      std::string fitsFile 
         = dynamic_cast<SpatialMap *>(srcFuncs["SpatialDist"])->fitsFile();
      std::string basefilename = ::basename(fitsFile.c_str());
      if (basefilename == "gas.cel") {
         srcType = "GalDiffuse";
      }      
   } else {
      throw Exception(std::string("Likelihood::FluxModel::getSourceType:\n")
                      + "unknown source type");
   }
   return;
}   

DomElement * FluxModel::gammaSpectrum(optimizers::Function & spectrum) {
   
   DomElement * specElt = ::createElement(m_doc, "spectrum");
   xml::Dom::addAttribute(*specElt, std::string("escale"), std::string("MeV"));

   DomElement * partElt = ::createElement(m_doc, "particle");
   xml::Dom::addAttribute(*partElt, std::string("name"), std::string("gamma"));

// Determine spectral type and set parameter values.
   DomElement * spectralTypeElt = ::createElement(m_doc, "power_law");
   
   xml::Dom::addAttribute(*spectralTypeElt, std::string("emin"), 
                          m_energies.front());
   xml::Dom::addAttribute(*spectralTypeElt, std::string("emax"), 
                          m_energies.back());

// It might be better here to have the Function objects set their own
// parameter value attributes, but this keeps the coupling between
// Functions and the xml interface looser.
   if (spectrum.genericName() == "PowerLaw") {
      xml::Dom::addAttribute(*spectralTypeElt, std::string("gamma"), 
                             -spectrum.getParamValue("Index"));
   } else if (spectrum.genericName() == "BrokenPowerLaw") {
      xml::Dom::addAttribute(*spectralTypeElt, std::string("gamma"), 
                             -spectrum.getParamValue("Index1"));
      xml::Dom::addAttribute(*spectralTypeElt, std::string("gamma2"), 
                             -spectrum.getParamValue("Index2"));
      xml::Dom::addAttribute(*spectralTypeElt, std::string("ebreak"), 
                             -spectrum.getParamValue("BreakValue"));
   }
   partElt->appendChild(*spectralTypeElt);
   specElt->appendChild(*partElt);
   return specElt;
}

DomElement * FluxModel::srcDirection(optimizers::Function & dir) {
   DomElement * dirElt = ::createElement(m_doc, "celestial_dir");
   xml::Dom::addAttribute(*dirElt, std::string("ra"), 
                          dynamic_cast<SkyDirFunction*>(&dir)->getDir().ra());
   xml::Dom::addAttribute(*dirElt, std::string("dec"), 
                          dynamic_cast<SkyDirFunction*>(&dir)->getDir().dec());
   return dirElt;
}

DomElement * FluxModel::solidAngle(double mincos, double maxcos) {
   DomElement * solidAngle = ::createElement(m_doc, "solid_angle");
   xml::Dom::addAttribute(*solidAngle, std::string("mincos"), mincos);
   xml::Dom::addAttribute(*solidAngle, std::string("maxcos"), maxcos);
   return solidAngle;
}

DomElement * FluxModel::galDiffuse(Source & src) {
   Source::FuncMap & srcFuncs = src.getSrcFuncs();

   DomElement * specElt = ::createElement(m_doc, "spectrum");
   xml::Dom::addAttribute(*specElt, std::string("escale"), std::string("MeV"));

   DomElement * specClassElt = ::createElement(m_doc, "SpectrumClass");
   xml::Dom::addAttribute(*specClassElt, std::string("name"), 
                          std::string("MapSpectrum"));

   if (srcFuncs["Spectrum"]->genericName() != "PowerLaw") {
      throw Exception(std::string("SourceModel::write_fluxXml:\n")
                      + "Galactic Diffuse spectral model is not a power-law.");
   } else {
// flux::MapSpectrum does not allow one to change the overall scaling
// of the FITS map used for modeling the Galactic diffuse emission. So
// we can only use the spectral index from our fit in this xml entry.
      std::ostringstream params;
      params << "100.," 
             << m_energies[0] << "," 
             << m_energies[m_energies.size()-1] << "," 
             << -srcFuncs["Spectrum"]->getParamValue("Index") << ","
             << "/sources/gas_gal.fits";
      xml::Dom::addAttribute(*specClassElt, std::string("params"), 
                             params.str());
   }
   specElt->appendChild(*specClassElt);
   DomElement * useSpecElt = ::createElement(m_doc, "use_spectrum");
   xml::Dom::addAttribute(*useSpecElt, std::string("frame"), 
                          std::string("galaxy"));
   specElt->appendChild(*useSpecElt);
   return specElt;
}

void FluxModel::makeEnergyGrid(unsigned int nee) {
   RoiCuts * roiCuts = RoiCuts::instance();
   std::pair<double, double> elims = roiCuts->getEnergyCuts();
   double estep = log(elims.second/elims.first)/(nee-1);
   m_energies.reserve(nee);
   for (unsigned int i = 0; i < nee; i++) {
      m_energies.push_back(elims.first*exp(estep*i));
   }
}

void FluxModel::addUnderscores(std::string &name) {
// Replace spaces with underscores.
   std::replace(name.begin(), name.end(), ' ', '_');

// Prepend underscore if name starts with an integer character.
   if (static_cast<int>(*name.begin()) >= '0'
       && static_cast<int>(*name.begin()) <= '9') {
      name = "_" + name;
   }
}

} // namespace Likelihood
