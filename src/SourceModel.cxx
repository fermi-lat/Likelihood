/** 
 * @file SourceModel.cxx
 * @brief SourceModel class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceModel.cxx,v 1.28 2003/11/07 02:27:10 jchiang Exp $
 */

#include <cmath>
#include <cassert>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include <xercesc/dom/DOM_Document.hpp>
#include <xercesc/dom/DOM_Element.hpp>
#include <xercesc/dom/DOM_NodeList.hpp>
#include <xercesc/dom/DOM_DOMException.hpp>

#include "facilities/Util.h"

#include "optimizers/Arg.h"
#include "optimizers/FunctionFactory.h"

#include "Likelihood/RoiCuts.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SourceFactory.h"
#include "Likelihood/TrapQuad.h"
#include "Likelihood/SourceModel.h"

namespace {
   void add_underscores(std::string &name) {
// Replace spaces with underscores.
      std::replace(name.begin(), name.end(), ' ', '_');
// Prepend underscore if name starts with an integer character.
      if (static_cast<int>(*name.begin()) >= '0'
          && static_cast<int>(*name.begin()) <= '9') {
         std::ostringstream new_name;
         new_name << "_" << name;
         name = new_name.str();
      }
   }
}

namespace Likelihood {

int SourceModel::s_refCount = 0;
std::vector<Source *> SourceModel::s_sources;

SourceModel::~SourceModel() {
   s_refCount--;
   if (s_refCount == 0) {
      std::vector<Source *>::iterator it = s_sources.begin();
      for (; it != s_sources.end(); it++)
         delete (*it);
   }
   s_sources.clear();
}

void SourceModel::setParam(const optimizers::Parameter &param, 
                           const std::string &funcName,
                           const std::string &srcName) 
   throw(optimizers::Exception) {
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      if (srcName == (*s_sources[i]).getName()) {
         Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
         if (srcFuncs.count(funcName)) {
            srcFuncs[funcName]->setParam(param.getName(), 
                                         param.getValue(),
                                         param.isFree());
// this seems inefficient, but necessary because of srcFuncs map(?)
            syncParams();
            return;
         }
      }
   }
   std::ostringstream errorMessage;
   errorMessage << "SourceModel::setParam:  Function " 
                << funcName << " for Source "
                << srcName << " was not found.\n";
   throw Exception(errorMessage.str());
}
 
std::vector<double>::const_iterator SourceModel::setParamValues_(
   std::vector<double>::const_iterator it) { 
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) 
         it = (*func_it).second->setParamValues_(it);
   }
   syncParams();
   return it;
}

std::vector<double>::const_iterator SourceModel::setFreeParamValues_(
   std::vector<double>::const_iterator it) { 
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) 
         it = (*func_it).second->setFreeParamValues_(it);
   }
   syncParams();
   return it;
}

optimizers::Parameter SourceModel::getParam(const std::string &paramName,
                                            const std::string &funcName,
                                            const std::string &srcName) const 
   throw(optimizers::Exception, optimizers::ParameterNotFound) {
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      if (srcName == (*s_sources[i]).getName()) {
         std::vector<optimizers::Parameter> params;
         Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
         if (srcFuncs.count(funcName)) {    //check for funcName
            srcFuncs[funcName]->getParams(params);
            for (unsigned int j = 0; j < params.size(); j++) {
               if (paramName == params[j].getName()) {
                  return params[j];
               }
            }
            throw optimizers::ParameterNotFound(paramName, funcName, 
                                    "SourceModel::getParam");
         }
         std::ostringstream errorMessage;
         errorMessage << "SourceModel::getParam:\n"
                      << "Function " << funcName 
                      << " was not found in Source " 
                      << srcName << "\n";
         throw Exception(errorMessage.str());
      }
   }
   std::ostringstream errorMessage;
   errorMessage << "SourceModel::getParam: "
                << "Source " << srcName 
                << " was not found.\n";
   throw Exception(errorMessage.str());
}

void SourceModel::setParamBounds(const std::string &paramName,
                                 const std::string &funcName,
                                 const std::string &srcName,
                                 double lower, double upper)
   throw(optimizers::ParameterNotFound, optimizers::OutOfBounds) {
   optimizers::Parameter my_param = getParam(paramName, funcName, srcName);
   my_param.setBounds(lower, upper);
   setParam(my_param, funcName, srcName);
   syncParams();
}

void SourceModel::setParamScale(const std::string &paramName,
                                const std::string &funcName,
                                const std::string &srcName,
                                double scale) 
   throw(optimizers::ParameterNotFound) {
   optimizers::Parameter my_param = getParam(paramName, funcName, srcName);
   my_param.setScale(scale);
   setParam(my_param, funcName, srcName);
   syncParams();
}

void SourceModel::setParamTrueValue(const std::string &paramName,
                                    const std::string &funcName,
                                    const std::string &srcName,
                                    double paramValue)
   throw(optimizers::ParameterNotFound, optimizers::OutOfBounds) {
   optimizers::Parameter my_param = getParam(paramName, funcName, srcName);
   my_param.setTrueValue(paramValue);
   setParam(my_param, funcName, srcName);
   syncParams();
}

void SourceModel::setParams_(std::vector<optimizers::Parameter> &params, 
                             bool setFree) 
   throw(optimizers::Exception, optimizers::ParameterNotFound) {
// ensure the number of Parameters matches
   unsigned int numParams;
   if (setFree) {
      numParams = getNumFreeParams();
   } else {
      numParams = getNumParams();
   }
   if (params.size() != numParams) {
      std::string errorMessage = std::string("SourceModel::setParams:\n") 
         + std::string("Inconsistent number of Parameters.");
      throw Exception(errorMessage);
   }
// assume ordering of Parameters in params matches that given by the
// ordering of the Sources and their Functions
   int k = 0;  // params' index
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
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
// loop over sources to ensure unique names
   for (unsigned int i = 0; i < s_sources.size(); i++) 
      assert((*src).getName() != (*s_sources[i]).getName());

// add a clone of this Source to the vector
   s_sources.push_back(src->clone());

// add the Parameters to the m_parameter vector 
// (would it be better just to use syncParams() here?)
   Source::FuncMap srcFuncs = (*src).getSrcFuncs();
   Source::FuncMap::iterator func_it = srcFuncs.begin();
   for (; func_it != srcFuncs.end(); func_it++) {
      std::vector<optimizers::Parameter> params;
      (*func_it).second->getParams(params);
      for (unsigned int ip = 0; ip < params.size(); ip++) 
         m_parameter.push_back(params[ip]);
   }      
}
 
void SourceModel::deleteSource(const std::string &srcName) 
   throw(optimizers::Exception) {
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      if (srcName == (*s_sources[i]).getName()) {
         s_sources.erase(s_sources.begin() + i, 
                         s_sources.begin() + i + 1);
         syncParams();
         return;
      }
   }
   std::ostringstream errorMessage;
   errorMessage << "SourceModel::deleteSource: " 
                << srcName << " was not found.\n";
   throw Exception(errorMessage.str());
}

void SourceModel::deleteAllSources() {
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   for (unsigned int i = 0; i < srcNames.size(); i++)
      deleteSource(srcNames[i]);
   m_parameter.clear();
}

Source * SourceModel::getSource(const std::string &srcName) {
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      if (srcName == s_sources[i]->getName())
         return s_sources[i];
   }
   return 0;
}

void SourceModel::getSrcNames(std::vector<std::string> &names) const {
// ensure names is empty
   if (!names.empty()) names.erase(names.begin(), names.end());

   for (unsigned int i = 0; i < s_sources.size(); i++) {
      names.push_back(s_sources[i]->getName());
   }
}

double SourceModel::evaluate_at(optimizers::Arg &x) const {
   double my_val = 0.;
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         my_val += (*func_it).second->value(x);
      }
   }
   return my_val;
}

// remake parameter vector from scratch 
void SourceModel::syncParams() {
   m_parameter.clear();

   for (unsigned int i = 0; i < s_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         std::vector<optimizers::Parameter> params;
         (*func_it).second->getParams(params);
         for (unsigned int ip = 0; ip < params.size(); ip++)
            m_parameter.push_back(params[ip]);
      }
   }
}

void SourceModel::fetchDerivs(optimizers::Arg &x, std::vector<double> &derivs, 
                              bool getFree) const {
   if (!derivs.empty()) derivs.clear();

   for (unsigned int i = 0; i < s_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
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
                          optimizers::FunctionFactory &funcFactory) {

// Expand any environment variables in the xmlFile name.
   facilities::Util::expandEnvVar(&xmlFile);

// Create a SourceFactory to read in the xml file.
   SourceFactory srcFactory;
   try {
      srcFactory.readXml(xmlFile, funcFactory);
   } catch (DOM_DOMException &eObj) {
      std::cout << "DOMException: " 
                << eObj.code << std::endl;
   }

// Loop over the sources that are now contained in srcFactory and add
// each one to the source model.
   std::vector<std::string> srcNames;
   srcFactory.fetchSrcNames(srcNames);

   std::vector<std::string>::iterator nameIt = srcNames.begin();
   for ( ; nameIt != srcNames.end(); nameIt++) {
      Source *src = srcFactory.create(*nameIt);
      std::cout << "adding source " << *nameIt << std::endl;
      addSource(src);
   }
   syncParams();
}

void SourceModel::writeXml(std::string xmlFile,
                           const std::string &functionLibrary) {
   xml::XmlParser *parser = new xml::XmlParser();

   DOM_Document doc = DOM_Document::createDocument();

   DOM_Element srcLib = doc.createElement("source_library");
   srcLib.setAttribute("function_library", functionLibrary.c_str());

// This attribute value need to be settable, either via data members
// or by hand.
   srcLib.setAttribute("title", "prototype sources");

// Loop over Sources.
   std::vector<Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); srcIt++) {
      DOM_Element srcElt = doc.createElement("source");
      std::string name = (*srcIt)->getName();
      srcElt.setAttribute("name", name.c_str());

// Add the xml data for the spatial and spectral Functions describing
// each source.
      Source::FuncMap srcFuncs = (*srcIt)->getSrcFuncs();
      if (srcFuncs.count("Spectrum")) {
         DOM_Element specElt = doc.createElement("spectrum");
         std::string name = srcFuncs["Spectrum"]->genericName();
         if (name != std::string("")) {
            specElt.setAttribute("type", name.c_str());
         } else {
            std::ostringstream errorMessage;
            errorMessage << "SourceModel::writeXml: "
                         << "genericName not found for the spectrum "
                         << "Function object of Source "
                         << (*srcIt)->getName() << "."
                         << std::endl;
            throw Exception(errorMessage.str());
         }
         srcFuncs["Spectrum"]->appendParamDomElements(doc, specElt);
         srcElt.appendChild(specElt);
      }         

      DOM_Element spatialElt = doc.createElement("spatialModel");
      if (srcFuncs.count("Position")) {
// This is a PointSource.
         srcElt.setAttribute("type", "PointSource");
         spatialElt.setAttribute("type", "SkyDirFunction");
         srcFuncs["Position"]->appendParamDomElements(doc, spatialElt);
         srcElt.appendChild(spatialElt);
      } else if (srcFuncs.count("SpatialDist")) {
// It's a DiffuseSource.
         srcElt.setAttribute("type", "DiffuseSource");
         std::string type = srcFuncs["SpatialDist"]->genericName();
         spatialElt.setAttribute("type", type.c_str());
         if (type == "SpatialMap") {
            std::string file = 
               dynamic_cast<SpatialMap *>(srcFuncs["SpatialDist"])->fitsFile();
            spatialElt.setAttribute("file", file.c_str());
         }
         srcFuncs["SpatialDist"]->appendParamDomElements(doc, spatialElt);
         srcElt.appendChild(spatialElt);
      }
      srcLib.appendChild(srcElt);
   }

// Expand any environment variables in the xmlFile name.
   facilities::Util::expandEnvVar(&xmlFile);

   std::ofstream outFile(xmlFile.c_str());
//    outFile << "<?xml version='1.0' standalone='no'?>\n"
//            << "<!DOCTYPE source_library SYSTEM \"A1_Sources.dtd\" >\n";
   xml::Dom::prettyPrintElement(srcLib, outFile, "");

   delete parser;
}

void SourceModel::write_fluxXml(std::string xmlFile) {
   xml::XmlParser *parser = new xml::XmlParser();

   DOM_Document doc = DOM_Document::createDocument();

   DOM_Element srcLib = doc.createElement("source_library");
   srcLib.setAttribute("title", "Likelihood_model");

// Create an energy vector for integrating the fluxes from Sources.
   RoiCuts * roiCuts = RoiCuts::instance();
   std::pair<double, double> elims = roiCuts->getEnergyCuts();
   int nee = 200;
   double estep = log(elims.second/elims.first)/(nee-1);
   std::vector<double> energies(nee);
   std::vector<double>::iterator enIt = energies.begin();
   for (int i=0 ; enIt != energies.end(); enIt++, i++) {
      *enIt = elims.first*exp(estep*i);
   }

// Create a nested source containing all sources.
   DOM_Element allSrcsElt = doc.createElement("source");
   std::ostringstream allSrcsName;
   allSrcsName << "all_in_" << xmlFile;
   allSrcsElt.setAttribute("name", allSrcsName.str().c_str());

// Loop over Sources.
   std::vector<Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); srcIt++) {
      DOM_Element srcElt = doc.createElement("source");
      std::string name = (*srcIt)->getName();
      ::add_underscores(name);
      srcElt.setAttribute("name", name.c_str());

// Get each Source's properties in terms of its functions.
      Source::FuncMap srcFuncs = (*srcIt)->getSrcFuncs();

// Consider only PointSources for now.
      if (srcFuncs.count("Position")) {

// Compute the flux integrated over the energy range given by the ROI.
         TrapQuad fluxIntegral(srcFuncs["Spectrum"]);
         std::ostringstream flux;
         flux << fluxIntegral.integral(energies)/1e-4;
         srcElt.setAttribute("flux", flux.str().c_str());

// Prepare the spectrum tag.
         DOM_Element specElt = doc.createElement("spectrum");
         specElt.setAttribute("escale", "MeV");

// particle tag.
         DOM_Element partElt = doc.createElement("particle");
         partElt.setAttribute("name", "gamma");

// Determine spectral type and set parameter values.
         DOM_Element spectralTypeElt = doc.createElement("power_law");
         std::ostringstream emin;
         emin << elims.first;
         spectralTypeElt.setAttribute("emin", emin.str().c_str());
         std::ostringstream emax;
         emax << elims.second;
         spectralTypeElt.setAttribute("emax", emax.str().c_str());
         if (srcFuncs["Spectrum"]->genericName() == "PowerLaw") {
            std::ostringstream gamma;
            gamma << -srcFuncs["Spectrum"]->getParamValue("Index");
            spectralTypeElt.setAttribute("gamma", gamma.str().c_str());
         } else if (srcFuncs["Spectrum"]->genericName() == "BrokenPowerLaw") {
            std::ostringstream gamma, gamma2, ebreak;
            gamma << -srcFuncs["Spectrum"]->getParamValue("Index1");
            spectralTypeElt.setAttribute("gamma", gamma.str().c_str());
            gamma2 << -srcFuncs["Spectrum"]->getParamValue("Index2");
            spectralTypeElt.setAttribute("gamma2", gamma2.str().c_str());
            ebreak << srcFuncs["Spectrum"]->getParamValue("BreakValue");
            spectralTypeElt.setAttribute("ebreak", ebreak.str().c_str());
         }
         partElt.appendChild(spectralTypeElt);
         specElt.appendChild(partElt);

// The source direction.
         DOM_Element dirElt = doc.createElement("celestial_dir");
         std::ostringstream ra, dec;
         ra << dynamic_cast<SkyDirFunction *>
            (srcFuncs["Position"])->getDir().ra();
         dirElt.setAttribute("ra", ra.str().c_str());
         dec << dynamic_cast<SkyDirFunction *>
            (srcFuncs["Position"])->getDir().dec();
         dirElt.setAttribute("dec", dec.str().c_str());
         specElt.appendChild(dirElt);

         srcElt.appendChild(specElt);
      }
      srcLib.appendChild(srcElt);

// Append the nested sources.
      DOM_Element nestedSrcElt = doc.createElement("nestedSource");
      nestedSrcElt.setAttribute("sourceRef", name.c_str());
      allSrcsElt.appendChild(nestedSrcElt);
   }

   srcLib.appendChild(allSrcsElt);

// Expand any environment variables in the xmlFile name.
   facilities::Util::expandEnvVar(&xmlFile);

   std::ofstream outFile(xmlFile.c_str());
//    outFile << "<?xml version='1.0' standalone='no'?>\n"
//            << "<!DOCTYPE source_library SYSTEM \"A1_Sources.dtd\" >\n";
   xml::Dom::prettyPrintElement(srcLib, outFile, "");

   delete parser;
}

} // namespace Likelihood
