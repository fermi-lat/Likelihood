/**
 * @file TsMap.cxx
 * @brief Prototype standalone application for producing
 * "test-statistic" maps.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/TsMap/TsMap.cxx,v 1.38 2007/07/03 22:48:19 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <sstream>

#include "fitsio.h"

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Util.h"

#include "optimizers/dArg.h"
#include "optimizers/Optimizer.h"
#include "optimizers/OptimizerFactory.h"
#include "optimizers/Exception.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/LogLike.h"

using namespace Likelihood;

/**
 * @class TsMap
 * @brief Class for encapsulating methods for creating a test-statistic
 * map.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/TsMap/TsMap.cxx,v 1.38 2007/07/03 22:48:19 jchiang Exp $
 */
class TsMap : public st_app::StApp {
public:
   TsMap();
   virtual ~TsMap() throw() {
      try {
         delete m_opt;
         delete m_helper;
      } catch (std::exception & eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
   }
   virtual void run();
   virtual void banner() const;
private:
   AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;
   LogLike * m_logLike;
   optimizers::Optimizer * m_opt;
   st_stream::StreamFormatter * m_formatter;

   std::vector<std::string> m_eventFiles;
   std::vector<double> m_lonValues;
   std::vector<double> m_latValues;
   std::vector< std::vector<double> > m_tsMap;
   std::string m_coordSys;
   void readSrcModel();
   void readEventData();
   void selectOptimizer();
   void setGrid();
   void computeMap();
   void writeFitsFile(const std::string &filename,
                      std::vector<double> &lon, 
                      std::vector<double> &lat,
                      std::vector< std::vector<double> > &map,
                      const std::string &coordSystem);
   void fitsReportError(FILE *stream, int status);
   void makeDoubleVector(double xmin, double xmax, int nx,
                         std::vector<double> &xVals);
   void setPointSourceSpectrum(PointSource &src);

   static std::string s_cvs_id;
};

st_app::StAppFactory<TsMap> myAppFactory("gttsmap");

TsMap::TsMap() : st_app::StApp(), m_helper(0), 
                 m_pars(st_app::StApp::getParGroup("gttsmap")),
                 m_logLike(0), m_opt(0),
                 m_formatter(new st_stream::StreamFormatter("gttsmap", "", 2)) 
{
   setVersion(s_cvs_id);
}

std::string TsMap::s_cvs_id("$Name:  $");

void TsMap::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void TsMap::run() {
   m_pars.Prompt();
   m_pars.Save();

   m_helper = new AppHelpers(&m_pars, "UNBINNED");
   m_helper->checkOutputFile();
   std::string coordSys = m_pars["coordsys"];
   m_coordSys = coordSys;
   std::string expCubeFile = m_pars["expcube"];
   if (expCubeFile != "" && expCubeFile != "none") {
      m_helper->observation().expCube().readExposureCube(expCubeFile);
   }
   st_facilities::Util::file_ok(m_pars["evfile"]);
   st_facilities::Util::resolve_fits_files(m_pars["evfile"], m_eventFiles);
   std::string ev_table = m_pars["evtable"];
   bool compareGtis(false);
   bool relyOnStreams(false);
   std::string respfunc = m_pars["irfs"];
   bool skipEventClassCuts(respfunc != "DSS");
   for (unsigned int i = 1; i < m_eventFiles.size(); i++) {
      AppHelpers::checkCuts(m_eventFiles[0], ev_table,
                            m_eventFiles[i], ev_table,
                            compareGtis, relyOnStreams, skipEventClassCuts);
   }
/// @bug Ascertain why this was called with just the first event file.
///   m_helper->setRoi(m_eventFiles[0]);
   m_helper->setRoi();
   m_helper->readScData();
   m_helper->readExposureMap();

   m_logLike = new LogLike(m_helper->observation());
   readEventData();
   readSrcModel();

   selectOptimizer();
   setGrid();

   computeMap();
   writeFitsFile(m_pars["outfile"], m_lonValues, m_latValues,
                 m_tsMap, m_coordSys);
}

void TsMap::readEventData() {
   std::vector<std::string>::const_iterator evIt = m_eventFiles.begin();
   for ( ; evIt != m_eventFiles.end(); evIt++) {
      st_facilities::Util::file_ok(*evIt);
      m_helper->observation().eventCont().getEvents(*evIt);
   }
}

void TsMap::readSrcModel() {
   std::string srcModelFile = m_pars["srcmdl"];
   if (srcModelFile != "" && srcModelFile != "none") {
      st_facilities::Util::file_ok(srcModelFile);
      m_logLike->readXml(srcModelFile, m_helper->funcFactory());
      m_logLike->computeEventResponses();
   }
}   

void TsMap::selectOptimizer() {
   std::string optimizer = m_pars["optimizer"];
   m_opt = optimizers::OptimizerFactory::instance().create(optimizer,
                                                           *m_logLike);
}

void TsMap::setGrid() {
   int nlon = m_pars["nx"];
   int nlat = m_pars["ny"];
   makeDoubleVector(m_pars["xref_min"], m_pars["xref_max"], nlon, m_lonValues);
   makeDoubleVector(m_pars["yref_min"], m_pars["yref_max"], nlat, m_latValues);
   m_tsMap.resize(nlat);
   for (int i = 0; i < nlat; i++) {
      m_tsMap.reserve(nlon);
   }
}

void TsMap::computeMap() {
   double ra(m_lonValues.front());
   double dec(m_latValues.front());
   std::string coordSys = m_pars["coordsys"];
   if (coordSys == "GAL") {
      astro::SkyDir my_dir(ra, dec, astro::SkyDir::GALACTIC);
      ra = my_dir.ra();
      dec = my_dir.dec();
   }
   PointSource testSrc(ra, dec, m_helper->observation());
   setPointSourceSpectrum(testSrc);
   testSrc.setName("testSource");

   int verbosity = m_pars["chatter"];
   verbosity -= 2;
   double tol = m_pars["ftol"];
   double logLike0;
   try {
      m_opt->find_min(verbosity, tol);
      logLike0 = m_logLike->value();
   } catch (...) {
      logLike0 = 0;
   }
   bool computeExposure(true);

   for (unsigned int jj = 0; jj < m_latValues.size(); jj++) {
      if ((jj % m_latValues.size()/20) == 0) {
         m_formatter->warn() << ".";
      }
      for (unsigned int ii = 0; ii < m_lonValues.size(); ii++) {
         if (m_coordSys == "CEL") {
            testSrc.setDir(m_lonValues[ii], m_latValues[jj], 
                           computeExposure, false);
         } else if (m_coordSys == "GAL") {
            testSrc.setGalDir(m_lonValues[ii], m_latValues[jj], 
                              computeExposure, false);
         } else {
            throw std::invalid_argument("Invalid coordinate system: "
                                        + m_coordSys);
         }
         m_logLike->addSource(&testSrc);
         try {
            m_opt->find_min(verbosity, tol);
            m_tsMap[jj].push_back(2.*(m_logLike->value() - logLike0));
         } catch (optimizers::Exception &eObj) {
            m_formatter->err() << eObj.what() << std::endl;
            // Default null value.
            m_tsMap[jj].push_back(0);
         }
         m_formatter->info(3) << m_lonValues[ii] << "  "
                              << m_latValues[jj] << "  "
                              << m_tsMap[jj][ii] << "  "
                              << m_helper->observation().eventCont().events().size() 
                              << std::endl;
         m_logLike->deleteSource(testSrc.getName());
      }
   }
   m_formatter->warn() << "!" << std::endl;
}

void TsMap::makeDoubleVector(double xmin, double xmax, int nx,
                             std::vector<double> &xVals) {
   xVals.clear();
   xVals.reserve(nx);
   double xstep = (xmax - xmin)/(nx-1);
   for (int i = 0; i < nx; i++) {
      xVals.push_back(xstep*i + xmin);
   }
}

void TsMap::setPointSourceSpectrum(PointSource &src) {
   optimizers::Function * pl = m_helper->funcFactory().create("PowerLaw");
   double parValues[] = {1., -2., 100.};
   std::vector<double> pars(parValues, parValues + 3);
   pl->setParamValues(pars);
   optimizers::Parameter indexParam = pl->getParam("Index");
   indexParam.setBounds(-3.5, -1.);
   pl->setParam(indexParam);
   optimizers::Parameter prefactorParam = pl->getParam("Prefactor");
   prefactorParam.setBounds(1e-10, 1e3);
   prefactorParam.setScale(1e-9);
   pl->setParam(prefactorParam);
   src.setSpectrum(pl);
}

void TsMap::writeFitsFile(const std::string &filename, 
                          std::vector<double> &lon,
                          std::vector<double> &lat,
                          std::vector< std::vector<double> > &map,
                          const std::string &coordSystem) {
   fitsfile *fptr;
   int status = 0;
   
// Always overwrite an existing file.
   remove(filename.c_str());
   fits_create_file(&fptr, filename.c_str(), &status);
   fitsReportError(stderr, status);

   int bitpix = DOUBLE_IMG;
   long naxis = 2;
   long naxes[] = {lon.size(), lat.size()};
   fits_create_img(fptr, bitpix, naxis, naxes, &status);
   fitsReportError(stderr, status);

// Repack exposure into a 1D vector
   std::vector<double> mapVector;
   mapVector.reserve(lon.size()*lat.size());
   for (unsigned int i = 0; i < lon.size(); i++) {
      for (unsigned int j = 0; j < lat.size(); j++) {
         mapVector.push_back(map[i][j]);
      }
   }

// Write the exposure map data.
   long group = 0;
   long firstelem = 1;
   long nelements = lon.size()*lat.size();
   fits_write_img_dbl(fptr, group, firstelem, nelements,
                      &mapVector[0], &status);
   fitsReportError(stderr, status);

// Write some keywords.
   double l0 = lon[0];
   fits_update_key(fptr, TDOUBLE, "CRVAL1", &l0, 
                   "longitude of reference pixel", &status);
   fitsReportError(stderr, status);
   double b0 = lat[0];
   fits_update_key(fptr, TDOUBLE, "CRVAL2", &b0, 
                   "latitude of reference pixel", &status);
   fitsReportError(stderr, status);
   
   double lstep = lon[1] - lon[0];
   fits_update_key(fptr, TDOUBLE, "CDELT1", &lstep, 
                   "longitude step at ref. pixel", &status);
   fitsReportError(stderr, status);
   double bstep = lat[1] - lat[0];
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
   
   if (coordSystem == "CEL") {
      char *ctype1 = "RA---CAR";
      fits_update_key(fptr, TSTRING, "CTYPE1", ctype1, 
                      "right ascension", &status);
      fitsReportError(stderr, status);
      char *ctype2 = "DEC--CAR";
      fits_update_key(fptr, TSTRING, "CTYPE2", ctype2, 
                      "declination", &status);
      fitsReportError(stderr, status);
   } else if (coordSystem == "GAL") {
      char *ctype1 = "GLON-CAR";
      fits_update_key(fptr, TSTRING, "CTYPE1", ctype1, 
                      "Galactic longitude", &status);
      fitsReportError(stderr, status);
      char *ctype2 = "GLAT-CAR";
      fits_update_key(fptr, TSTRING, "CTYPE2", ctype2, 
                      "Galactic latitude", &status);
      fitsReportError(stderr, status);
   }
   
   fits_close_file(fptr, &status);
   fitsReportError(stderr, status);
   
   return;
}

void TsMap::fitsReportError(FILE *stream, int status) {
   fits_report_error(stream, status);
   if (status != 0) {
      throw std::runtime_error("TsMap::writeFitsFile: cfitsio error.");
   }
}
