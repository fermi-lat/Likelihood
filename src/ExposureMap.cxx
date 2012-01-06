/**
 * @file ExposureMap.cxx
 * @brief This class encapsulates exposure map information and makes
 * it available for use (primarily) by the DiffuseSource class.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/ExposureMap.cxx,v 1.45 2012/01/06 07:11:59 jchiang Exp $
 */

#include <cstdio>

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <utility>

#include "st_stream/StreamFormatter.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "facilities/Util.h"

#include "st_facilities/FitsImage.h"
#include "st_facilities/Util.h"

#include "Likelihood/ExposureMap.h"
#include "Likelihood/Observation.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/WcsMap2.h"

namespace Likelihood {

ExposureMap::~ExposureMap() {
   delete m_wcsmap;
}

void ExposureMap::readExposureFile(std::string exposureFile) {

   facilities::Util::expandEnvVar(&exposureFile);

   std::string extension;
   bool interpolate, enforceEnergyRange;
   m_wcsmap = new WcsMap2(exposureFile, extension="", interpolate=true,
                          enforceEnergyRange=true);

   st_facilities::FitsImage mapData(exposureFile);

   mapData.getCelestialArrays(m_ra, m_dec);

   readEnergyExtension(exposureFile, m_energies);

   mapData.getSolidAngles(m_solidAngles);

// Fill the vector of the planes of m_exposure for each true photon
// energy.
   m_exposure.clear();
   std::vector<double> exposure;
   mapData.getImageData(exposure);

   int indx = 0;
   int npixels = m_solidAngles.size();
   for (unsigned int k = 0; k < m_energies.size(); k++) {
      std::vector<double> expArray(npixels);
      for (int sa_indx = 0; sa_indx < npixels; sa_indx++) {
         expArray.at(sa_indx) = exposure.at(indx);
         indx++;
      }
      m_exposure.push_back(expArray);
   }
   m_haveExposureMap = true;
}

void ExposureMap::integrateSpatialDist(const std::vector<double> & energies,
                                       optimizers::Function * spatialDist,
                                       std::vector<double> & exposure) const {
   exposure.clear();
   exposure.reserve(energies.size());
   for (unsigned int k = 0; k < energies.size(); k++) {
      double srcExposure = 0;
// Find the index kk (of the energy array that describes the exposure
// map data) that will be used for interpolating at energies[k]; here
// we assume the ExposureMap energies are logarithmically spaced.
      std::vector<double>::const_iterator iterE;
      if (energies[k] < m_energies[0]) {
         iterE = m_energies.begin() + 1;
      } else if (energies[k] >= *(m_energies.end() - 1)) {
         iterE = m_energies.end() - 1;
      } else {
         iterE = std::upper_bound(m_energies.begin(), m_energies.end(), 
                                  energies[k]);
      }
      int kk = iterE - m_energies.begin();
      for (unsigned int j = 0; j < m_ra.size(); j++) {
         double expsr = log(energies[k]/(*(iterE-1)))/log(*iterE/(*(iterE-1)));
         expsr = expsr*(m_exposure[kk][j] - m_exposure[kk-1][j])
            + m_exposure[kk-1][j];
         astro::SkyDir skyDir(m_ra[j], m_dec[j]);
         SkyDirArg dir(skyDir, energies[k]);
         srcExposure += expsr*(*spatialDist)(dir)*m_solidAngles[j];
      }
      exposure.push_back(srcExposure);
   }
}

void ExposureMap::computeMap(std::string filename, 
                             const Observation & observation,
                             double sr_radius, int nlon, int nlat,
                             int nenergies, bool compute_submap,
                             int nlongmin, int nlongmax, 
                             int nlatmin, int nlatmax) {

   facilities::Util::expandEnvVar(&filename);

   const RoiCuts & roiCuts = observation.roiCuts();

   astro::SkyDir roiCenter = roiCuts.extractionRegion().center();

   std::pair<double, double> elims = roiCuts.getEnergyCuts();
   double estep = log(elims.second/elims.first)/(nenergies - 1);
   std::vector<double> energies;
   for (int k = 0; k < nenergies; k++) {
      energies.push_back(elims.first*exp(estep*k));
   }

   double crpix[] = {nlon/2 + 0.5, nlat/2 + 0.5, 1.};
   double crval[] = {roiCenter.ra(), roiCenter.dec(), energies.front()};
   double cdelt[] = {2.*sr_radius/(nlon-1), 2.*sr_radius/(nlat-1), estep};

   std::vector<long> naxes(3);
   naxes.at(0) = nlon;
   naxes.at(1) = nlat;
   naxes.at(2) = nenergies;

   astro::SkyProj proj("STG", crpix, crval, cdelt, 0, false);

   std::vector<double> exposure;
   exposure.reserve(nenergies);

   std::vector<float> expMap;
   expMap.resize(nenergies*nlat*nlon);

   st_stream::StreamFormatter formatter("ExposureMap", "computeMap", 2);

   formatter.info() << "Computing the ExposureMap";
   if (observation.expCube().haveFile()) {
      formatter.info() << " using " << observation.expCube().fileName()
                       << std::endl;
   } else {
      formatter.info() << " (no expCube file given) " << std::endl;
   }

   int ncount = 0;

   int imin(0);
   int imax(nlon);
   int jmin(0);
   int jmax(nlat);
   if (compute_submap) {
      imin = std::max(0, nlongmin);
      imax = std::min(nlon, nlongmax);
      jmin = std::max(0, nlatmin);
      jmax = std::min(nlat, nlatmax);
   }

   for (int j = jmin; j < jmax; j++) {
      for (int i = imin; i < imax; i++) {
         int step(((imax-imin)*(jmax-jmin))/20);
         if (step == 0) {
            step = 2;
         }
         if ((ncount % step) == 0) {
            formatter.warn() << ".";
         }

// NB: wcslib (via astro::SkyProj) starts indexing pixels at 1, not 0, 
// so apply correction here to avoid off-by-one error.
         std::pair<double, double> coords = proj.pix2sph(i+1, j+1);
         astro::SkyDir dir(coords.first, coords.second);

         bool verbose(false);
         if (observation.expCube().haveFile()) {
            PointSource::computeExposureWithHyperCube(dir, energies, 
                                                      observation, 
                                                      exposure, verbose);
         } else {
            PointSource::computeExposure(dir, energies, observation,
                                         exposure, verbose);
         }
         for (int k = 0; k < nenergies; k++) {
            int indx = (k*nlat + j)*nlon + i;
            expMap.at(indx) = exposure[k];
         }
         ncount++;
      }
   }
   formatter.warn() << "!" << std::endl;

   writeFitsFile(filename, naxes, crpix, crval, cdelt, energies, expMap);
}

void ExposureMap::writeFitsFile(const std::string & filename,
                                const std::vector<long> & naxes,
                                double * crpix, double * crval, double * cdelt,
                                const std::vector<double> & energies, 
                                const std::vector<float> & expMap) {
   
   if (st_facilities::Util::fileExists(filename)) {
      std::remove(filename.c_str());
   }
   std::string ext("PRIMARY");
   tip::IFileSvc::instance().appendImage(filename, ext, naxes);
   tip::Image * image = tip::IFileSvc::instance().editImage(filename, ext);

   tip::Header & header = image->getHeader();

   header["TELESCOP"].set("GLAST");
   header["INSTRUME"].set("LAT SIMULATION");
   header["DATE-OBS"].set("");
   header["DATE-END"].set("");

   header["CRVAL1"].set(crval[0]);
   header["CRPIX1"].set(crpix[0]);
   header["CDELT1"].set(cdelt[0]);
   header["CTYPE1"].set("RA---STG");

   header["CRVAL2"].set(crval[1]);
   header["CRPIX2"].set(crpix[1]);
   header["CDELT2"].set(cdelt[1]);
   header["CTYPE2"].set("DEC--STG");

   header["CRVAL3"].set(crval[2]);
   header["CRPIX3"].set(crpix[2]);
   header["CDELT3"].set(cdelt[2]);
   header["CTYPE3"].set("log_Energy");

   image->set(expMap);
   delete image;

   ext = "ENERGIES";
   tip::IFileSvc::instance().appendTable(filename, ext);
   tip::Table * table = tip::IFileSvc::instance().editTable(filename, ext);
   table->appendField("Energy", "1D");
// Repair field by removing TNULL keyword that is added by tip. The
// null value is usually ok for integers, but is inappropriate for
// floats and is not needed by either, so we remove it in every case.
   int fieldIndex = table->getFieldIndex("Energy") + 1;
   std::ostringstream nullkeyword;
   nullkeyword << "TNULL" << fieldIndex;
   try {
      table->getHeader().erase(nullkeyword.str());
   } catch (...) {
      // do nothing if tip fails us again
   }

   table->setNumRecords(energies.size());

   tip::Table::Iterator row = table->begin();
   tip::Table::Record & record = *row;

   std::vector<double>::const_iterator energy = energies.begin();
   for ( ; energy != energies.end(); ++energy, ++row) {
      record["Energy"].set(*energy);
   }

   delete table;
}

void ExposureMap::readEnergyExtension(const std::string & filename,
                                      std::vector<double> & energies) {
   std::auto_ptr<const tip::Table> 
      table(tip::IFileSvc::instance().readTable(filename, "ENERGIES"));
   energies.resize(table->getNumRecords());

   tip::Table::ConstIterator row = table->begin();
   tip::ConstTableRecord & record = *row;

   std::vector<double>::iterator energy = energies.begin();
   for ( ; energy != energies.end(); ++energy, ++row) {
      record["Energy"].get(*energy);
   }
}

double ExposureMap::operator()(const astro::SkyDir & dir,
                               double energy) const {
   return m_wcsmap->operator()(dir, energy);
}

double ExposureMap::operator()(const astro::SkyDir & dir,
                               int k) const {
   return m_wcsmap->operator()(dir, k);
}

bool ExposureMap::withinMapRadius(const astro::SkyDir & dir) const {
   return m_wcsmap->withinMapRadius(dir);
}

} // namespace Likelihood
