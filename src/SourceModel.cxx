/** 
 * @file SourceModel.cxx
 * @brief SourceModel class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceModel.cxx,v 1.52 2004/09/13 15:30:39 jchiang Exp $
 */

#include <cassert>
#include <cmath>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "fitsio.h"

#include "xml/Dom.h"
#include "xml/XmlParser.h"

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "optimizers/Arg.h"
#include "optimizers/FunctionFactory.h"

#include "Likelihood/CountsMap.h"
#include "Likelihood/Exception.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SourceFactory.h"
#include "Likelihood/TrapQuad.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/FluxBuilder.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SourceModelBuilder.h"
#include "Likelihood/SourceModel.h"

namespace {
   void fitsReportError(FILE *stream, int status) {
      fits_report_error(stream, status);
      if (status != 0) {
         throw std::string("writeExposureFile: cfitsio error.");
      }
   }
} // unnamed namespace

namespace Likelihood {

int SourceModel::s_refCount = 0;
std::map<std::string, Source *> SourceModel::s_sources;

SourceModel::~SourceModel() {
   s_refCount--;
   if (s_refCount == 0) {
      std::map<std::string, Source *>::iterator it = s_sources.begin();
      for (; it != s_sources.end(); it++)
         delete (*it).second;
   }
   s_sources.clear();
}

void SourceModel::setParam(const optimizers::Parameter &param, 
                           const std::string &funcName,
                           const std::string &srcName) 
   throw(optimizers::Exception) {
   if (s_sources.count(srcName)) {
      Source::FuncMap srcFuncs = (*s_sources[srcName]).getSrcFuncs();
      if (srcFuncs.count(funcName)) {
         srcFuncs[funcName]->setParam(param.getName(), 
                                      param.getValue(),
                                      param.isFree());
         syncParams();
         return;
      }
   }
   std::string errorMessage = "SourceModel::setParam:\n Function " 
      + funcName + " for Source " + srcName + " was not found.";
   throw optimizers::Exception(errorMessage);
}
 
std::vector<double>::const_iterator SourceModel::setParamValues_(
   std::vector<double>::const_iterator it) { 
   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) 
         it = (*func_it).second->setParamValues_(it);
   }
   syncParams();
   return it;
}

std::vector<double>::const_iterator SourceModel::setFreeParamValues_(
   std::vector<double>::const_iterator it) { 
   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
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
   if (s_sources.count(srcName)) {
      std::vector<optimizers::Parameter> params;
      Source::FuncMap srcFuncs = s_sources[srcName]->getSrcFuncs();
      if (srcFuncs.count(funcName)) {    //check for funcName
         try {
            srcFuncs[funcName]->getParams(params);
         } catch (optimizers::Exception &eObj) {
            std::cerr << eObj.what() << std::endl;
            assert(false);
         }                 
         for (unsigned int j = 0; j < params.size(); j++) {
            if (paramName == params[j].getName()) {
               return params[j];
            }
         }
         throw optimizers::ParameterNotFound(paramName, funcName, 
                                             "SourceModel::getParam");
      }
      std::string errorMessage = "SourceModel::getParam:\n Function"
         + funcName + " was not found in Source " 
         + srcName + ".";
      throw optimizers::Exception(errorMessage);
   }
   std::string errorMessage = "SourceModel::getParam:\nSource "
      + srcName + " was not found.";
   throw optimizers::Exception(errorMessage);
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
// Ensure the number of Parameters matches.
   unsigned int numParams;
   if (setFree) {
      numParams = getNumFreeParams();
   } else {
      numParams = getNumParams();
   }
   if (params.size() != numParams) {
      std::string errorMessage = std::string("SourceModel::setParams:\n") 
         + "Inconsistent number of Parameters.";
      throw optimizers::Exception(errorMessage);
   }
// Assume ordering of Parameters in params matches that given by the
// ordering of the Sources and their Functions.
   int k = 0;  // params' index
   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
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
   if (!s_sources.count(src->getName())) {
      s_sources[src->getName()] = src->clone();
      syncParams();
   } else {
      throw Exception("Likelihood::SourceModel:\nSource named " 
                      + src->getName() + " alread exists.");
   }
}
 
Source * SourceModel::deleteSource(const std::string &srcName) 
   throw(optimizers::Exception) {
   std::map<std::string, Source *>::iterator it = s_sources.find(srcName);
   if (it != s_sources.end()) {
      Source * mySource = it->second;
      s_sources.erase(it);
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
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source * source = deleteSource(srcNames[i]);
      delete source;
   }
   m_parameter.clear();
}

Source * SourceModel::getSource(const std::string &srcName) {
   if (s_sources.count(srcName)) {
      return s_sources[srcName];
   }
   return 0;
}

void SourceModel::getSrcNames(std::vector<std::string> &names) const {
   names.clear();
   std::map<std::string, Source *>::iterator it = s_sources.begin();
   for ( ; it != s_sources.end(); ++it) {
      names.push_back(it->first);
   }
}

double SourceModel::value(optimizers::Arg &x) const {
   double my_val = 0.;
   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         my_val += (*func_it).second->value(x);
      }
   }
   return my_val;
}

void SourceModel::syncParams() { // remake parameter vector from scratch 
   m_parameter.clear();

   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         std::vector<optimizers::Parameter> params;
         (*func_it).second->getParams(params);
         for (unsigned int ip = 0; ip < params.size(); ip++)
            m_parameter.push_back(params[ip]);
      }
   }
}

void SourceModel::fetchDerivs(optimizers::Arg &x,
                              std::vector<double> &derivs, 
                              bool getFree) const {
   derivs.clear();

   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      Source::FuncMap srcFuncs = srcIt->second->getSrcFuncs();
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
                          optimizers::FunctionFactory &funcFactory,
                          bool requireExposure) {

// Expand any environment variables in the xmlFile name.
   facilities::Util::expandEnvVar(&xmlFile);

// Create a SourceFactory to read in the xml file.
   SourceFactory srcFactory;
   try {
      srcFactory.readXml(xmlFile, funcFactory, requireExposure);
   } catch (xml::DomException &eObj) {
      std::cout << "SourceModel::readXml:\n DomException: " 
                << eObj.what() << std::endl;
   }

// Loop over the sources that are now contained in srcFactory and add
// each one to the source model.
   std::vector<std::string> srcNames;
   srcFactory.fetchSrcNames(srcNames);

   std::vector<std::string>::iterator nameIt = srcNames.begin();
   for ( ; nameIt != srcNames.end(); nameIt++) {
      Source *src = srcFactory.create(*nameIt);
      if (m_verbose) std::cout << "adding source " << *nameIt << std::endl;
      addSource(src);
   }
   syncParams();
}

void SourceModel::reReadXml(std::string xmlFile) {

   facilities::Util::expandEnvVar(&xmlFile);

   xml::XmlParser *parser = new xml::XmlParser();

   DomDocument doc = parser->parse(xmlFile.c_str());

   if (doc == 0) {
      std::string errorMessage = "SourceFactory::readXml:\nInput xml file, "
         + xmlFile + ", not parsed successfully.";
      throw Exception(errorMessage);
   }

//*** direct Xerces API call...still available in Xerces 2.4.0 ***/
   DomElement source_library = doc.getDocumentElement();
   if (!xml::Dom::checkTagName(source_library, "source_library")) {
      throw Exception("SourceModel::reReadXml:\nsource_library not found");
   }

// Loop through source DOM_Elements and Source objects in parallel.
   std::vector<DomElement> srcs;
   xml::Dom::getChildrenByTagName(source_library, "source", srcs);

   for (unsigned int j = 0; j < srcs.size(); j++) {
      std::string srcName = xml::Dom::getAttribute(srcs[j], "name");
      Source * my_source(0);
      if (s_sources.count(srcName)) {
         my_source = s_sources[srcName];
      } else {
         continue;
      }
      std::string srcType = xml::Dom::getAttribute(srcs[j], "type");
// Get spectrum and spatialModel elements
      std::vector<DomElement> child;
      xml::Dom::getChildrenByTagName(srcs[j], "spectrum", child);
      DomElement spectrum = child[0];

      xml::Dom::getChildrenByTagName(srcs[j], "spatialModel", child);
      DomElement spatialModel = child[0];

      my_source->getSrcFuncs()["Spectrum"]->setParams(spectrum);
      if (srcType == "PointSource") {
// Extract (RA, Dec) from the parameter elements and use
// PointSource::setDir(...) method to ensure that exposure gets
// recalculated.
         double ra(0), dec(0);
         std::vector<DomElement> params;
         xml::Dom::getChildrenByTagName(spatialModel, "parameter", params);

         std::vector<DomElement>::const_iterator paramIt = params.begin();
         for ( ; paramIt != params.end(); paramIt++) {
            std::string name = xml::Dom::getAttribute(*paramIt, "name");
            if (name == "RA")
               ra = ::atof(xml::Dom::getAttribute(*paramIt, "value").c_str());
            if (name == "DEC")
               dec = ::atof(xml::Dom::getAttribute(*paramIt, "value").c_str());
         }
         astro::SkyDir newDir(ra, dec);
         double tol = 1e-4;
         if (newDir.difference(
                dynamic_cast<PointSource *>(my_source)->getDir() ) > tol) {
// Reset the direction, re-computing the PointSource exposure.
            my_source->setDir(newDir);
         }
      } else if (srcType == "DiffuseSource") {
         my_source->getSrcFuncs()["SpatialDist"]->setParams(spatialModel);
      } else {
         std::cerr << "SourceModel::reReadXml: "
                   << "Unknown Source type: " << srcType << std::endl;
         assert(false);
      }
   }
   syncParams();
   delete parser;
}

void SourceModel::writeXml(std::string xmlFile,
                           const std::string &functionLibrary,
                           const std::string &srcLibTitle) {

   SourceModelBuilder builder(functionLibrary, srcLibTitle);

   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); srcIt++) {
      builder.addSource(*(srcIt->second));
   }

   builder.write(xmlFile);
}

void SourceModel::write_fluxXml(std::string xmlFile) {

   FluxBuilder builder;

   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); srcIt++) {
      builder.addSource(*(srcIt->second));
   }

   builder.write(xmlFile);
}

bool SourceModel::hasSrcNamed(const std::string & srcName) const {
   std::vector<std::string> names;
   getSrcNames(names);
   std::vector<std::string>::iterator it 
      = std::find(names.begin(), names.end(), srcName);
   if (it == names.end()) {
      return false;
   }
   return true;
}

CountsMap * SourceModel::createCountsMap(const CountsMap & dataMap) const {

   if (ExposureCube::instance() == 0) {
      std::runtime_error("SourceModel::makeCountsMap:\n"
                         + std::string("Exposure cube not available."));
   }

   std::vector<double> map;

   std::vector<astro::SkyDir> pixelDirs;
   std::vector<double> pixelSolidAngles;
   getPixels(dataMap, pixelDirs, pixelSolidAngles);

   std::vector<double> energies;
   dataMap.getAxisVector(2, energies);

   computeModelMap(pixelDirs, pixelSolidAngles, energies, map);

   CountsMap * modelMap = new CountsMap(dataMap);
// @bug This will not work properly until
// tip::FitsExtensionManager::setImageDimensions is fixed.
   modelMap->setImage(map);
   return modelMap;
}

// void SourceModel::computeModelMap(const std::vector<astro::SkyDir> & pixelDirs,
//                                   const std::vector<double> & pixelSolidAngles,
//                                   const std::vector<double> & energies,
//                                   std::vector<double> & modelMap) {

//    modelMap.resize(pixelDirs.size()*(energies.size()-1));

//    std::map<std::string, Source *>::const_iterator src;
// // Loop over pixel directions.
//    for (unsigned int j = 0; j < pixelDirs.size(); j++) {
//       for (unsigned int k = 0; k < energies.size()-1; k++) {
//          int indx = k*pixelDirs.size() + j;
//          for (int evtType = 0; evtType < 2; evtType++) {
//             for (src = s_sources.begin(); src != s_sources.end(); ++src) {
//                Aeff aeff1(src->second, pixelDirs[j], energies[k], evtType);
//                double map_lower =
//                   ExposureCube::instance()->value(pixelDirs[j], aeff1);
//                Aeff aeff2(src->second, pixelDirs[j], energies[k+1], evtType);
//                double map_upper = 
//                   ExposureCube::instance()->value(pixelDirs[j], aeff2);
//                modelMap[indx] += (map_lower + map_upper)/2.*pixelSolidAngles[j]
//                   *(energies[k+1] - energies[k]);
//             } // src
//          } // evtType
//       } // k
//    } // j
// }

void SourceModel::computeModelMap(const std::vector<Pixel> & pixels,
                                  const std::vector<double> & energies,
                                  std::vector<double> & modelMap) {

   modelMap.reserve(pixels.size()*(energies.size()-1));
   for (unsigned int j = 0; j < pixels.size(); j++) {
      for (unsigned int k = 0; k < energies.size()-1; k++) {
         modelMap.push_back(pixels[j].modelCounts(energies[k], energies[k+1]));
      }
   }
}

void SourceModel::getPixels(const CountsMap & countsMap,
                            std::vector<Pixel> & pixels) {
   pixels.clear();
   std::vector<astro::SkyDir> pixelDirs;
   std::vector<double> solidAngles;
   getPixels(countsMap, pixelDirs, solidAngles);
   pixels.reserve(pixelDirs.size());
   for (unsigned int i = 0; i < pixelDirs.size(); i++) {
      pixels.push_back(Pixel(pixelDirs[i], solidAngles[i]));
   }
}

void SourceModel::getPixels(const CountsMap & countsMap, 
                            std::vector<astro::SkyDir> & pixelDirs,
                            std::vector<double> & pixelSolidAngles) {

   long nx = countsMap.imageDimension(0);
   long ny = countsMap.imageDimension(1);

   std::vector<double> longitudes;
   countsMap.getAxisVector(0, longitudes);
   std::vector<double> latitudes;
   countsMap.getAxisVector(1, latitudes);
   std::vector<double> energies;
   countsMap.getAxisVector(2, energies);

   pixelDirs.clear();
   pixelSolidAngles.clear();

   pixelDirs.reserve(nx*ny);
   pixelSolidAngles.reserve(nx*ny);
   std::vector<double>::const_iterator latIt = latitudes.begin();
   for ( ; latIt != latitudes.end() - 1; ++latIt) {
      double latitude = (*latIt + *(latIt+1))/2.;
      std::vector<double>::const_iterator lonIt = longitudes.begin();
      for ( ; lonIt != longitudes.end() - 1; ++lonIt) {
         double longitude = (*lonIt + *(lonIt+1))/2.;
         pixelDirs.push_back(astro::SkyDir(longitude, latitude, 
                                           countsMap.projection()));
         pixelSolidAngles.push_back(computeSolidAngle(lonIt, latIt, 
                                                      countsMap.projection()));
      }
   }
}

double SourceModel::computeSolidAngle(std::vector<double>::const_iterator lon,
                                      std::vector<double>::const_iterator lat,
                                      const astro::SkyProj & proj) {
   astro::SkyDir lower_left(*lon, *lat, proj);
   astro::SkyDir upper_right(*(lon+1), *(lat+1), proj);
   std::vector<double> theta(2);
   std::vector<double> phi(2);
   if (proj.isGalactic()) {
      phi[0] = lower_left.l()*M_PI/180.;
      theta[0] = lower_left.b()*M_PI/180.;
      phi[1] = upper_right.l()*M_PI/180.;
      theta[1] = upper_right.b()*M_PI/180.;
   } else {
      phi[0] = lower_left.ra()*M_PI/180.;
      theta[0] = lower_left.dec()*M_PI/180.;
      phi[1] = upper_right.ra()*M_PI/180.;
      theta[1] = upper_right.dec()*M_PI/180.;
   }
   return std::fabs((phi[1] - phi[0])*(sin(theta[1]) - sin(theta[0])));
}

// SourceModel::Aeff::Aeff(Source * src, const astro::SkyDir & appDir, 
//                         double energy, int type)
//    : m_src(src), m_appDir(appDir), m_energy(energy), m_type(type) {
//    PointSource * ptsrc = dynamic_cast<PointSource *>(src);
//    if (ptsrc == 0) {
//       m_separation = 90.;
//    } else {
//       m_separation = ptsrc->getDir().difference(appDir)*180./M_PI;
//    }
// }

// double SourceModel::Aeff::operator()(double costheta) const {
//    if (m_separation < 90.) {
//       double inclination = acos(costheta)*180./M_PI;
//       static double phi(0);
//       return m_src->fluxDensity(inclination, phi, m_energy, m_separation,
//                                 m_type);
//    }
//    return 0;
// }

} // namespace Likelihood
