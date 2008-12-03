/**
 * @file LikeExposure.cxx
 * @brief Implementation of Exposure class for use by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/LikeExposure.cxx,v 1.32 2008/09/07 20:34:10 jchiang Exp $
 */

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "facilities/commonUtilities.h"
#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "fitsio.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "healpix/CosineBinner.h"
#include "healpix/HealpixArray.h"

#include "Likelihood/LikeExposure.h"
#include "Likelihood/RoiCuts.h"

namespace {
   bool compareFirst(const std::pair<double, double> & a, 
                     const std::pair<double, double> & b) {
      return a.first < b.first;
   }
   bool compareSecond(const std::pair<double, double> & a, 
                      const std::pair<double, double> & b) {
      return a.second < b.second;
   }
}

namespace Likelihood {

LikeExposure::
LikeExposure(double skybin, double costhetabin, 
             const std::vector< std::pair<double, double> > & timeCuts,
             const std::vector< std::pair<double, double> > & gtis,
             double zenmax)
   : map_tools::Exposure(skybin, costhetabin, std::cos(zenmax*M_PI/180.)), 
     m_costhetabin(costhetabin), m_timeCuts(timeCuts), m_gtis(gtis), m_numIntervals(0) {
   if (!gtis.empty()) {
      for (size_t i = 0; i < gtis.size(); i++) {
         if (i == 0 || gtis.at(i).first < m_tmin) {
            m_tmin = gtis.at(i).first;
         }
         if (i == 0 || gtis.at(i).second > m_tmax) {
            m_tmax = gtis.at(i).second;
         }
      }
   } else {
      throw std::runtime_error("LikeExposure::LikeExposure: GTIs are empty.\n"
                               "Cannot proceed with livetime calculation.");
   }
}

void LikeExposure::load(const tip::Table * scData, bool verbose) {
   st_stream::StreamFormatter formatter("LikeExposure", "load", 2);
   
   double ra, dec, ra_zenith, dec_zenith, start, stop, livetime;

   tip::Table::ConstIterator it(scData->end());
   tip::ConstTableRecord & row(*it);

// Count the rows within the user selected time interval (m_tmin, m_tmax)
// by counting inwards from the top and bottom of the scData.
   --it;
   tip::Index_t nrows(scData->getNumRecords());
   if (nrows == 0) {
      return;
   }
   for ( ; it != scData->begin(); --it, nrows--) {
      row["stop"].get(stop);
      if (stop < m_tmax) {
         break;
      }
   }

   it = scData->begin();
   for ( ; it != scData->end(); ++it, nrows--) {
      row["start"].get(start);
      if (start > m_tmin) {
         break;
      }
   }

// Reset to the FT2 interval to the one that preceeds the
// user-selected interval, if possible; and set the start time to that
// of the initial row.
   if (it != scData->begin()) {
      --it;
   }
   row["start"].get(start);

// Set the step size for the printing out the little progress dots.
   tip::Index_t istep(nrows/20);
   if (istep == 0) {
      istep = 1;
   }

   for (tip::Index_t irow = 0; it != scData->end() && start < m_tmax; ++it, ++irow) {
      if (verbose && (irow % istep) == 0 ) {
         formatter.warn() << "."; 
      }
      row["livetime"].get(livetime);
      row["start"].get(start);
      row["stop"].get(stop);
      double deltat = livetime;
      double fraction;
      if (acceptInterval(start, stop, m_timeCuts, m_gtis, fraction)) {
         row["ra_scz"].get(ra);
         row["dec_scz"].get(dec);
         row["ra_zenith"].get(ra_zenith);
         row["dec_zenith"].get(dec_zenith);
         fill(astro::SkyDir(ra, dec), astro::SkyDir(ra_zenith, dec_zenith), 
              deltat*fraction);
         m_numIntervals++;
      }
   }
   if (verbose) {
      formatter.warn() << "!" << std::endl;
   }
}

void LikeExposure::writeFile(const std::string & outfile) const {
   std::string dataPath = 
      facilities::commonUtilities::getDataPath("Likelihood");
   std::string templateFile = 
      facilities::commonUtilities::joinPath(dataPath, "LivetimeCubeTemplate");
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   fileSvc.createFile(outfile, templateFile);

   writeFilename(outfile);

   writeLivetimes(outfile);

   writeCosbins(outfile);
}

void LikeExposure::writeFilename(const std::string & outfile) const {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   tip::Image * phdu(fileSvc.editImage(outfile, ""));
   phdu->getHeader()["FILENAME"].set(facilities::Util::basename(outfile));
   delete phdu;
}

void LikeExposure::writeLivetimes(const std::string & outfile) const {
   setCosbinsFieldFormat(outfile);

   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   tip::Table * table = fileSvc.editTable(outfile, "EXPOSURE");
   table->setNumRecords(data().size());

   tip::Table::Iterator it(table->begin());
   tip::TableRecord & row(*it);

   healpix::HealpixArray<healpix::CosineBinner>::const_iterator 
      pixel(data().begin());

   for ( ; pixel != data().end(); ++pixel, ++it) {
      row["COSBINS"].set(*pixel);
      astro::SkyDir dir(data().dir(pixel));
      row["RA"].set(dir.ra());
      row["DEC"].set(dir.dec());
   }

   tip::Header & header(table->getHeader());
   header["PIXTYPE"].set("HEALPIX"); 
   header["ORDERING"].set("NESTED"); 
   header["COORDSYS"].set(data().healpix().galactic()? "GAL" : "EQU");
   header["NSIDE"].set(data().healpix().nside()); 
   header["FIRSTPIX"].set(0); 
   header["LASTPIX"].set(data().size()-1); 
   header["THETABIN"].set(healpix::CosineBinner::thetaBinning());
   header["NBRBINS"].set(healpix::CosineBinner::nbins());
   header["COSMIN"].set(healpix::CosineBinner::cosmin());

   delete table;
}

void LikeExposure::writeCosbins(const std::string & outfile) const {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   tip::Table * table = fileSvc.editTable(outfile, "CTHETABOUNDS");
   table->setNumRecords(healpix::CosineBinner::nbins());

   tip::Table::Iterator it(table->begin());
   tip::TableRecord & row(*it);
   
   std::vector<double> mubounds;
   computeCosbins(mubounds);

   for (size_t i(0); i < mubounds.size() -1; i++, ++it) {
      row["CTHETA_MIN"].set(mubounds.at(i+1));
      row["CTHETA_MAX"].set(mubounds.at(i));
   }
   delete table;
}
   
void LikeExposure::setCosbinsFieldFormat(const std::string & outfile) const {
   int status(0);

   fitsfile * fptr(0);
   std::string extfilename(outfile + "[EXPOSURE]");
   fits_open_file(&fptr, extfilename.c_str(), READWRITE, &status);
   fitsReportError(status, "LikeExposure::setCosbinsFieldFormat");
   
   int colnum(1); // by assumption
   fits_modify_vector_len(fptr, colnum, data().at(0).size(), &status);
   fitsReportError(status, "LikeExposure::setCosbinsFieldFormat");

   fits_close_file(fptr, &status);
   fitsReportError(status, "LikeExposure::setCosbinsFieldFormat");
}

bool LikeExposure::
acceptInterval(double start, double stop, 
               const std::vector< std::pair<double, double> > & timeCuts,
               const std::vector< std::pair<double, double> > & gtis,
               double & fraction) {
                                  
   std::pair<double, double> candidateInterval(start, stop);

   typedef std::vector< std::pair<double, double> > IntervalCont_t;
   IntervalCont_t::const_iterator it;

   for (it = timeCuts.begin(); it != timeCuts.end(); ++it) {
      if (!overlaps(*it, candidateInterval)) {
         fraction = 0;
         return false;
      }
   }
   
   double total(0);
   double maxTotal(candidateInterval.second - candidateInterval.first);
   
   IntervalCont_t::const_iterator gti_first = 
      std::upper_bound(gtis.begin(), gtis.end(), candidateInterval,
                       ::compareFirst);
   if (gti_first != gtis.begin()) {
      --gti_first;
   }

   IntervalCont_t::const_iterator gti_last = 
      std::upper_bound(gti_first, gtis.end(), candidateInterval,
                       ::compareSecond);

   if (gti_last != gtis.end()) {
      ++gti_last;
   }

   for (it = gti_first; it != gti_last; ++it) {
      double dt(overlap(*it, candidateInterval));
      total += dt;
   }
   if (total > maxTotal || gtis.size() == 0) {
      total = maxTotal;
   }
   fraction = total/(stop - start);
   return true;
}

bool LikeExposure::
overlaps(const std::pair<double, double> & interval1,
         std::pair<double, double> & interval2) {
   double start = std::max(interval1.first, interval2.first);
   double stop = std::min(interval1.second, interval2.second);
   if (start < stop) {
      interval2.first = start;
      interval2.second = stop;
      return true;
   }
   return false;
}

double LikeExposure::
overlap(const std::pair<double, double> & interval1,
        const std::pair<double, double> & interval2) {
   double start = std::max(interval1.first, interval2.first);
   double stop = std::min(interval1.second, interval2.second);
   if (start < stop) {
      return stop - start;
   }
   return 0;
}

void LikeExposure::
fitsReportError(int status, const std::string & routine) const {
   if (status == 0) {
      return;
   }
   fits_report_error(stderr, status);
   std::ostringstream message;
   message << routine << ": CFITSIO error " << status;
   throw std::runtime_error(message.str());
}

void LikeExposure::
computeCosbins(std::vector<double> & mubounds) const {
   bool sqrtbins(healpix::CosineBinner::thetaBinning() == "SQRT(1-COSTHETA)");
   double cosmin(healpix::CosineBinner::cosmin());
   size_t nbins(healpix::CosineBinner::nbins());
   mubounds.clear();
//   for (int i(nbins); i >= 0; i--) {
   for (size_t i(0); i < nbins+1; i++) {
      double factor(static_cast<double>(i)/nbins);
      if (sqrtbins) {
         factor *= factor;
      }
      mubounds.push_back(1. - factor*(1. - cosmin));
   }
}

} // namespace Likelihood
