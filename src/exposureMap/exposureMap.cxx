#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#include "fitsio.h"

#include "facilities/Util.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "astro/SkyDir.h"

#include "irfInterface/Irfs.h"

#include "map_tools/Exposure.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/ResponseFunctions.h"

using namespace Likelihood;

class Aeff : public map_tools::Exposure::Aeff {
public:
   Aeff(double energy, int evtType) : m_energy(energy), m_evtType(evtType),
      m_cutoff(0.25) {}
   virtual ~Aeff() {}
   virtual double operator()(double cosTheta) const {
//      return cosTheta < m_cutoff ? 0 : (cosTheta - m_cutoff)/(1. - m_cutoff);
      double inclination = acos(cosTheta);
      std::map<unsigned int, irfInterface::Irfs *>::iterator respIt 
         = ResponseFunctions::instance()->begin();
      for ( ; respIt != ResponseFunctions::instance()->end(); ++respIt) {
         if (respIt->second->irfID() == m_evtType) {  
            irfInterface::IAeff * aeff = respIt->second->aeff();
            double aeff_val = aeff->value(m_energy, inclination, s_phi);
            return aeff_val;
         }
      }
      return 0;
   }
private:
   double m_energy;
   int m_evtType;
   double m_cutoff;
   static double s_phi;
};

double Aeff::s_phi(0);

class exposureMap : public st_app::StApp {
public:
   exposureMap();
   virtual ~exposureMap() throw() {}
   virtual void run();
private:
   AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;

   std::vector<double> m_energies;

   void computeMap(std::string filename);

   void computeEnergies();

   void writeFitsFile(const std::string &filename, 
                      std::vector<double> &glon,
                      std::vector<double> &glat,
                      std::vector<double> &energies,
                      std::vector< std::vector<double> > &exposure,
                      double ra0, double dec0);

   void fitsReportError(FILE * stream, int status) const;

};

st_app::StAppFactory<exposureMap> myAppFactory;

exposureMap::exposureMap() 
   : st_app::StApp(), m_helper(0),
     m_pars(st_app::StApp::getParGroup("exposureMap")) {
   try {
      m_pars.Prompt();
      m_pars.Save();
      m_helper = new AppHelpers(m_pars);
//      ResponseFunctions::setEdispFlag(m_pars["use_energy_dispersion"]);
   } catch (std::exception & eObj) {
      std::cerr << eObj.what() << std::endl;
      std::exit(1);
   } catch (...) {
      std::cerr << "Caught unknown exception in exposureMap constructor." 
                << std::endl;
      std::exit(1);
   }
}

void exposureMap::run() {
   std::string expcube_file = m_pars["exposure_cube_file"];
   if (expcube_file == "none") {
      throw std::runtime_error("Please specify an exposure cube file.");
   }
   ExposureCube::readExposureCube(expcube_file);
   std::string output_file = m_pars["output_file_name"];
   computeEnergies();
   computeMap(output_file);
}

void exposureMap::computeEnergies() {
   double emin(20.);
   double emax(2e5);
   int nee(20);
   double estep = log(emax/emin)/(nee - 1.);
   m_energies.clear();
   m_energies.reserve(nee);
   for (int i = 0; i < nee; i++) {
      m_energies.push_back(emin*exp(estep*i));
   }
}

void exposureMap::computeMap(std::string filename) {
   facilities::Util::expandEnvVar(&filename);

   int nlon(360);
   double lonMin(-180.);
   double lonstep(1.);
   int nlat(180);
   double latMin(-90.);
   double latstep(1.);
   std::vector<double> lon;
   for (int i = 0; i < nlon; i++) {
      lon.push_back(lonstep*i + lonMin);
   }
   std::vector<double> lat;
   for (int j = 0; j < nlat; j++) {
      lat.push_back(latstep*j + latMin);
   }

   std::vector< std::vector<double> > exposureCube;
   int nenergies(m_energies.size());
   exposureCube.resize(nenergies);
   for (int k = 0; k < nenergies; k++) {
      exposureCube[k].resize(nlon*nlat, 0);
   }

   int indx = 0;
   for (int j = 0; j < nlat; j++) {
      for (int i = 0; i < nlon; i++) {
         if ( (indx % ((nlon*nlat)/20)) == 0 ) std::cerr << ".";
         for (int k = 0; k < nenergies; k++) {
            for (int evtType = 0; evtType < 2; evtType++) {
               Aeff aeff(m_energies[k], evtType);
               astro::SkyDir dir(lon[i], lat[j], astro::SkyDir::GALACTIC);
               exposureCube[k][indx] 
                  += ExposureCube::instance()->value(dir, aeff);
            }
         }
         indx++;
      }
   }
   std::cerr << "!" << std::endl;
   astro::SkyDir gc(0, 0, astro::SkyDir::GALACTIC);
   writeFitsFile(filename, lon, lat, m_energies, exposureCube, 
                 gc.l(), gc.b());
}

void exposureMap::writeFitsFile(const std::string &filename, 
                                std::vector<double> &glon,
                                std::vector<double> &glat,
                                std::vector<double> &energies,
                                std::vector< std::vector<double> > &exposure,
                                double ra0, double dec0) {
   fitsfile *fptr;
   int status = 0;
   
// Always overwrite an existing file.
   remove(filename.c_str());
   fits_create_file(&fptr, filename.c_str(), &status);
   fitsReportError(stderr, status);

   int bitpix = DOUBLE_IMG;
   long naxis = 3;
   long naxes[] = {glon.size(), glat.size(), energies.size()};
   fits_create_img(fptr, bitpix, naxis, naxes, &status);
   fitsReportError(stderr, status);

// Write the exposure map data.
   long group = 0;
   long dim1 = glon.size();
   long dim2 = glat.size();

// Repack exposure into a C array.
   double *exp_array = new double[glon.size()*glat.size()*energies.size()];
   int indx = 0;
   for (unsigned int k = 0; k < energies.size(); k++) {
      for (unsigned int i = 0; i < exposure[k].size(); i++) {
         exp_array[indx] = exposure[k][i];
         indx++;
      }
   }

   fits_write_3d_dbl(fptr, group, dim1, dim2, 
                     glon.size(), glat.size(), energies.size(),
                     exp_array, &status);
   delete[] exp_array;
   fitsReportError(stderr, status);

// Write some keywords.
   double l0 = glon[0];
   fits_update_key(fptr, TDOUBLE, "CRVAL1", &l0, 
                   "longitude of reference pixel", &status);
   fitsReportError(stderr, status);
   double b0 = glat[0];
   fits_update_key(fptr, TDOUBLE, "CRVAL2", &b0, 
                   "latitude of reference pixel", &status);
   fitsReportError(stderr, status);
   
   double lstep = glon[1] - glon[0];
   fits_update_key(fptr, TDOUBLE, "CDELT1", &lstep, 
                   "longitude step at ref. pixel", &status);
   fitsReportError(stderr, status);
   double bstep = glat[1] - glat[0];
   fits_update_key(fptr, TDOUBLE, "CDELT2", &bstep, 
                   "latitude step at ref. pixel", &status);
   fitsReportError(stderr, status);
   
   float crpix1 = 1.0;
   fits_update_key(fptr, TFLOAT, "CRPIX1", &crpix1, 
                   "reference pixel for longitude coordinate", &status);
   fitsReportError(stderr, status);
   float crpix2 = 1.0;
   fits_update_key(fptr, TFLOAT, "CRPIX2", &crpix2, 
                   "reference pixel for latitude coordinate", &status);
   fitsReportError(stderr, status);
   
   char *ctype1 = "GLON-CAR";
   fits_update_key(fptr, TSTRING, "CTYPE1", ctype1, 
                   "right ascension", &status);
   fitsReportError(stderr, status);
   char *ctype2 = "GLAT-CAR";
   fits_update_key(fptr, TSTRING, "CTYPE2", ctype2, 
                   "declination", &status);
   fitsReportError(stderr, status);

   double logEmin = log(energies[0]);
   fits_update_key(fptr, TDOUBLE, "CRVAL3", &logEmin,
                   "reference value for log_energy coordinate", &status);
   fitsReportError(stderr, status);

   int nee = energies.size();
   double estep = log(energies[nee-1]/energies[0])/(nee-1);
   fits_update_key(fptr, TDOUBLE, "CDELT3", &estep, 
                   "step in log_energy coordinate", &status);
   fitsReportError(stderr, status);

   float crpix3 = 1.;
   fits_update_key(fptr, TFLOAT, "CRPIX3", &crpix3,
                   "reference pixel for log_energy coordinate", &status);
   fitsReportError(stderr, status);

   char *ctype3 = "log_MeV";
   fits_update_key(fptr, TSTRING, "CTYPE3", ctype3,
                   "units for log_energy", &status);
   fitsReportError(stderr, status);

   fits_update_key(fptr, TDOUBLE, "ROI_RA", &ra0, "RA of ROI center", 
                   &status);
   fitsReportError(stderr, status);
   fits_update_key(fptr, TDOUBLE, "ROI_DEC", &dec0, "DEC of ROI center",
                   &status);
   fitsReportError(stderr, status);
   
// Write the energy array as a binary table.
   int nrows = energies.size();
   int tfields = 1;
   char *ttype[] = {"Energy"};
   char *tform[] = {"1D"};
   char *tunit[] = {"MeV"};
   char extname[] = "True Photon Energy Array";
   
   int firstrow  = 1;
   int firstelem = 1;
   
   fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform,
                   tunit, extname, &status);
   fitsReportError(stderr, status);
   
   fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, nrows, 
                  &energies[0], &status);
   fitsReportError(stderr, status);
   
   fits_close_file(fptr, &status);
   fitsReportError(stderr, status);
   
   return;
}

void exposureMap::fitsReportError(FILE *stream, int status) const {
   fits_report_error(stream, status);
   if (status != 0) {
      throw std::runtime_error("exposureMap: cfitsio error.");
   }
}
