/**
 * @file BinnedExposure.cxx
 * @brief Integral of effective area over time for the entire sky at
 * various energies.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BinnedExposure.cxx,v 1.4 2004/11/28 06:58:21 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <memory>
#include <stdexcept>

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/Observation.h"
#include "Likelihood/ResponseFunctions.h"

#include "Verbosity.h"

namespace Likelihood {

#include "fitsio.h"

BinnedExposure::BinnedExposure() : m_observation(0) {}

BinnedExposure::BinnedExposure(const std::vector<double> & energies,
                               const Observation & observation) 
   : m_observation(&observation), m_energies(energies) {
   computeMap();
}

BinnedExposure::BinnedExposure(const std::string & filename) {
   std::auto_ptr<const tip::Image> 
      image(tip::IFileSvc::instance().readImage(filename, ""));
   std::vector<float> data;
   image->get(data);
   m_exposureMap.resize(data.size());
   std::copy(data.begin(), data.end(), m_exposureMap.begin());

   const tip::Header & header = image->getHeader();

   /// @bug tip doesn't allow direct access to NAXIS[1-3] keywords
   /// through the Header, so we are forced to use
   /// getImageDimensions() method.
   std::vector<tip::PixOrd_t> dims = image->getImageDimensions();

   unsigned int nra = dims[0];
   double ra0;
   double rastep;
   double refpix1;
   header["CRVAL1"].get(ra0);
   header["CDELT1"].get(rastep);
   header["CRPIX1"].get(refpix1);
   m_ras.resize(nra);
   for (unsigned int i = 0; i < nra; i++) {
      m_ras.at(i) = (i - refpix1)*rastep + ra0;
   }
   
   unsigned int ndec = dims[1];
   double dec0;
   double decstep;
   double refpix2;
   header["CRVAL2"].get(dec0);
   header["CDELT2"].get(decstep);
   header["CRPIX2"].get(refpix2);
   m_decs.resize(ndec);
   for (unsigned int i = 0; i < ndec; i++) {
      m_decs.at(i) = (i - refpix2)*decstep + dec0;
   }

   std::auto_ptr<const tip::Table>
      energies(tip::IFileSvc::instance().readTable(filename, "Energies"));

   m_energies.clear();
   tip::Table::ConstIterator it = energies->begin();
   tip::ConstTableRecord & row = *it;
   for ( ; it != energies->end(); ++it) {
      double value;
      row["Energy"].get(value);
      m_energies.push_back(value);
   }
}

double BinnedExposure::operator()(double energy, double ra, double dec) const {
   std::vector<double>::const_iterator ie = std::find(m_energies.begin(),
                                                      m_energies.end(),
                                                      energy);
   if (ie == m_energies.end()) {
      throw std::runtime_error("BinnedExposure::operator(): The energy " +
                               std::string("selected is not available."));
   }
   unsigned int k = ie - m_energies.begin();
   if (ra < 0) ra += 360.; 
   if (ra > 360.) ra = fmod(ra, 360.);
   unsigned int i = findIndex(m_ras.begin(), m_ras.end(), ra);
   unsigned int j = findIndex(m_decs.begin(), m_decs.end(), dec);
   unsigned int indx = (k*m_decs.size() + j)*m_ras.size() + i;
   return m_exposureMap.at(indx);
}

unsigned int 
BinnedExposure::findIndex(std::vector<double>::const_iterator begin,
                          std::vector<double>::const_iterator end,
                          double value) const {
   std::vector<double>::const_iterator it 
      = std::upper_bound(begin, end, value);
   if (it == end) {
      throw std::range_error("BinnedExposure::findIndex: out of range value.");
   }
   int indx = it - begin - 1;
   return indx;
}

void BinnedExposure::computeMap() {
   linearArray(0., 360., 361, m_ras);
   linearArray(-90., 90., 181, m_decs);

   m_exposureMap.resize(m_ras.size()*m_decs.size()*m_energies.size(), 0);
   int iter(0);
   if (print_output()) std::cerr << "Computing binned exposure map";
   for (unsigned int j = 0; j < m_decs.size(); j++) {
      for (unsigned int i = 0; i < m_ras.size(); i++) {
         if ( print_output() && 
              (iter % ((m_decs.size()*m_ras.size())/20)) == 0 ) {
            std::cerr << ".";
         }
         astro::SkyDir dir(m_ras[i], m_decs[j], astro::SkyDir::EQUATORIAL);
         for (unsigned int k = 0; k < m_energies.size(); k++) {
            unsigned int indx 
               = k*m_ras.size()*m_decs.size() + j*m_ras.size() + i;
            for (int evtType = 0; evtType < 2; evtType++) {
               Aeff aeff(m_energies[k], evtType, *m_observation);
               m_exposureMap.at(indx)
                  += m_observation->expCube().value(dir, aeff);
            }
         }
         iter++;
      }
   }
   if (print_output()) std::cerr << "!" << std::endl;
}

void BinnedExposure::linearArray(double xmin, double xmax, unsigned int npts,
                                 std::vector<double> &xx) const {
   double xstep = (xmax - xmin)/(npts - 1.);
   for (unsigned int i = 0; i < npts; i++) {
      xx.push_back(i*xstep + xmin);
   }
}

double BinnedExposure::Aeff::s_phi(0);

double BinnedExposure::Aeff::operator()(double cosTheta) const {
   double inclination = acos(cosTheta)*180./M_PI;
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
      = m_observation.respFuncs().begin();
   for ( ; respIt != m_observation.respFuncs().end(); ++respIt) {
      if (respIt->second->irfID() == m_evtType) {
         irfInterface::IAeff * aeff = respIt->second->aeff();
         double aeff_val = aeff->value(m_energy, inclination, s_phi);
         return aeff_val;
      }
   }
   return 0;
}

void BinnedExposure::writeOutput(const std::string & filename) const {
   fitsfile *fptr;
   int status = 0;
   
// Always overwrite an existing file.
   remove(filename.c_str());
   fits_create_file(&fptr, filename.c_str(), &status);
   fitsReportError(stderr, status);

   int bitpix = DOUBLE_IMG;
   long naxis = 3;
   long naxes[] = {m_ras.size(), m_decs.size(), m_energies.size()};
   fits_create_img(fptr, bitpix, naxis, naxes, &status);
   fitsReportError(stderr, status);

// Write the exposure map data.
   long group = 0;
   long dim1 = m_ras.size();
   long dim2 = m_decs.size();

   fits_write_3d_dbl(fptr, group, dim1, dim2, 
                     m_ras.size(), m_decs.size(), m_energies.size(),
                     const_cast<double *>(&m_exposureMap[0]), &status);
   fitsReportError(stderr, status);

// Write some keywords.
   double ra0 = m_ras[180];
   fits_update_key(fptr, TDOUBLE, "CRVAL1", &ra0, 
                   "RA of reference pixel", &status);
   fitsReportError(stderr, status);
   double dec0 = m_decs[90];
   fits_update_key(fptr, TDOUBLE, "CRVAL2", &dec0, 
                   "Dec of reference pixel", &status);
   fitsReportError(stderr, status);
   
   double rastep = m_ras[1] - m_ras[0];
   fits_update_key(fptr, TDOUBLE, "CDELT1", &rastep, 
                   "RA step at reference pixel", &status);
   fitsReportError(stderr, status);
   double decstep = m_decs[1] - m_decs[0];
   fits_update_key(fptr, TDOUBLE, "CDELT2", &decstep, 
                   "Dec step at reference pixel", &status);
   fitsReportError(stderr, status);
   
   float crpix1 = 180.;
   fits_update_key(fptr, TFLOAT, "CRPIX1", &crpix1, 
                   "reference pixel for RA coordinate", &status);
   fitsReportError(stderr, status);
   float crpix2 = 90.;
   fits_update_key(fptr, TFLOAT, "CRPIX2", &crpix2, 
                   "reference pixel for Dec coordinate", &status);
   fitsReportError(stderr, status);
   
   char * ctype1 = "RA---CAR";
   fits_update_key(fptr, TSTRING, "CTYPE1", ctype1, 
                   "right ascension", &status);
   fitsReportError(stderr, status);
   char * ctype2 = "DEC--CAR";
   fits_update_key(fptr, TSTRING, "CTYPE2", ctype2, 
                   "declination", &status);
   fitsReportError(stderr, status);

   double logEmin = log(m_energies[0]);
   fits_update_key(fptr, TDOUBLE, "CRVAL3", &logEmin,
                   "reference value for log_energy coordinate", &status);
   fitsReportError(stderr, status);

   int nee = m_energies.size();
   double estep = log(m_energies[nee-1]/m_energies[0])/(nee-1);
   fits_update_key(fptr, TDOUBLE, "CDELT3", &estep, 
                   "step in log_energy coordinate", &status);
   fitsReportError(stderr, status);

   float crpix3 = 1.;
   fits_update_key(fptr, TFLOAT, "CRPIX3", &crpix3,
                   "reference pixel for log_energy coordinate", &status);
   fitsReportError(stderr, status);

   char * ctype3 = "log_MeV";
   fits_update_key(fptr, TSTRING, "CTYPE3", ctype3,
                   "units for log_energy", &status);
   fitsReportError(stderr, status);

// Write the energy array as a binary table.  
/// @bug Only the Energy column is needed, but tip can't read binary
/// tables that have just one column, so we are forced to add a dummy
/// column.
   int nrows = m_energies.size();
   int tfields = 2;
   char * ttype[] = {"Energy", "Emax"};
   char * tform[] = {"1D", "1D"};
   char * tunit[] = {"MeV", "MeV"};
   char extname[] = "Energies";
   
   int firstrow  = 1;
   int firstelem = 1;
   
   fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform,
                   tunit, extname, &status);
   fitsReportError(stderr, status);
   
   fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, nrows, 
                  const_cast<double *>(&m_energies[0]), &status);
   fitsReportError(stderr, status);

   fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, nrows, 
                  const_cast<double *>(&m_energies[0]), &status);
   fitsReportError(stderr, status);
   
   fits_close_file(fptr, &status);
   fitsReportError(stderr, status);
   
   return;
}

void BinnedExposure::fitsReportError(FILE *stream, int status) const {
   fits_report_error(stream, status);
   if (status != 0) {
      throw std::runtime_error("BinnedExposure: cfitsio error.");
   }
}

} // namespace Likelihood
