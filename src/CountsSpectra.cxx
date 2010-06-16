/**
 * @file CountsSpectra.cxx
 * @brief Encapsulation of counts spectra for a Likelihood fit.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/CountsSpectra.cxx,v 1.7 2010/05/03 18:25:15 jchiang Exp $
 */

#include <cmath>
#include <cstdio>

#include <stdexcept>

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "Likelihood/CountsSpectra.h"
#include "Likelihood/Observation.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceMap.h"

namespace Likelihood {

void CountsSpectra::setEbounds(double emin, double emax, size_t nbounds) {
   if (m_binnedLike) {
      throw std::runtime_error("CountsSpectra: cannot set ebounds for "
                               "binned likelihood.");
   }
   double estep(std::log(emax/emin)/(nbounds - 1.));
   m_ebounds.clear();
   m_ebounds.reserve(nbounds);
   for (size_t k = 0; k < nbounds; k++) {
      m_ebounds.push_back(emin*std::exp(k*estep));
   }
}

void CountsSpectra::getSrcCounts(const std::string & srcName,
                                 std::vector<double> & srcCounts) const {
   const Source * src(const_cast<LogLike &>(m_logLike).getSource(srcName));
   if (!src) {
      throw std::runtime_error("CountsSpectra: Source " + srcName + 
                               " not found.");
   }
   srcCounts.clear();
   srcCounts.reserve(m_ebounds.size() - 1);
   if (m_binnedLike) {
/// @todo reimplement LogLike::Npred for BinnedLikelihood

/// @bug We cannot use the following const reference because of brain-dead
/// linkage problems on Windows.
///       const std::vector<double> & 
///          npreds(m_binnedLike->sourceMap(srcName).npreds());
/// so we are forced to use instead the getNpreds function that pollutes our
/// namespace.
      std::vector<double> npreds;
//      getNpreds(m_binnedLike->sourceMap(srcName), npreds);
      m_binnedLike->getNpreds(srcName, npreds);
      for (size_t k = 0; k < m_ebounds.size() - 1; k++) {
         srcCounts.push_back(src->pixelCounts(m_ebounds.at(k),
                                              m_ebounds.at(k+1),
                                              npreds.at(k), 
                                              npreds.at(k+1)));
      }
   } else {
      check_ebounds();
      for (size_t k = 0; k < m_ebounds.size() - 1; k++) {
         srcCounts.push_back(src->Npred(m_ebounds.at(k), m_ebounds.at(k+1)));
      }
   }
}

void CountsSpectra::getSrcFluxes(const std::string & srcName, 
                                 std::vector<double> & srcFluxes) const {
   check_ebounds();
   const Source * src(const_cast<LogLike &>(m_logLike).getSource(srcName));
   if (!src) {
      throw std::runtime_error("CountsSpectra: Source " + srcName + 
                               " not found.");
   }
   const PointSource * ptsrc = 
      dynamic_cast<PointSource *>(const_cast<Source *>(src));
   if (!ptsrc) {
      std::string message = 
         "CountsSpectra: Cannot obtain flux from a diffuse source, "
         + srcName;
      throw std::runtime_error(message);
   }
   srcFluxes.clear();
   srcFluxes.reserve(m_ebounds.size() - 1);
   for (size_t k = 0; k < m_ebounds.size() - 1; k++) {
      srcFluxes.push_back(ptsrc->flux(m_ebounds.at(k), m_ebounds.at(k+1)));
   }
}

void CountsSpectra::getTotalSrcCounts(std::vector<double> & srcCounts) const {
   srcCounts.clear();
   std::vector<std::string> sourceNames;
   m_logLike.getSrcNames(sourceNames);
   if (!sourceNames.empty()) {
      // Get counts due to first source.
      getSrcCounts(sourceNames.at(0), srcCounts);
      for (size_t i = 1; i < sourceNames.size(); i++) {
         std::vector<double> componentCounts(srcCounts.size(), 0.);
         // Get counts due to ith component.
         getSrcCounts(sourceNames.at(i), componentCounts);
         // Add this component to the total.
         for (size_t j = 0; j < componentCounts.size(); j++) {
            srcCounts[j] += componentCounts[j];
         }
      }
   }
}

void CountsSpectra::writeTable(const std::string & outfile) const {
   std::remove(outfile.c_str());

   std::vector<std::string> sourceNames;
   m_logLike.getSrcNames(sourceNames);

   writeCounts(outfile, sourceNames);
   writeFluxes(outfile, sourceNames);
   writeEbounds(outfile);
}

void CountsSpectra::
writeCounts(const std::string & outfile, 
            const std::vector<std::string> & sourceNames) const {
   std::string extname("COUNTS_SPECTRA");

   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   fileSvc.appendTable(outfile, extname);
   tip::Table * counts(fileSvc.editTable(outfile, extname));

   counts->appendField("ObsCounts", "1D");
   std::vector<double> obsCounts;
   getObsCounts(obsCounts);

   std::vector<std::vector<double> > srcCnts;
   std::vector<std::string>::const_iterator source = sourceNames.begin();
   for ( ; source != sourceNames.end(); ++source) {
      counts->appendField(*source, "1D");
      std::vector<double> my_counts;
      getSrcCounts(*source, my_counts);
      srcCnts.push_back(my_counts);
   }

   counts->setNumRecords(m_ebounds.size() - 1);
   tip::Table::Iterator row = counts->begin();
   tip::Table::Record & record = *row;
   for (size_t k = 0; k < m_ebounds.size() - 1; ++row, k++) {
      record["ObsCounts"].set(obsCounts.at(k));
      for (size_t i = 0; i < sourceNames.size(); i++) {
         record[sourceNames.at(i)].set(srcCnts.at(i).at(k));
      }
   }
   delete counts;
}

void CountsSpectra::
writeFluxes(const std::string & outfile, 
            const std::vector<std::string> & sourceNames) const {
   std::string extname("FLUXES");

   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   fileSvc.appendTable(outfile, extname);
   tip::Table * fluxes(fileSvc.editTable(outfile, extname));

   std::map<std::string, std::vector<double> > srcFluxes;
   std::vector<std::string>::const_iterator source = sourceNames.begin();
   for (source = sourceNames.begin(); source != sourceNames.end(); ++source) {
      if (m_logLike.source(*source).getType() == "Point") {
         std::vector<double> my_fluxes;
         getSrcFluxes(*source, my_fluxes);
         srcFluxes[*source] = my_fluxes;
         fluxes->appendField(*source, "1D");
      }
   }

   if (srcFluxes.size() > 0) {
      fluxes->setNumRecords(m_ebounds.size() - 1);
      tip::Table::Iterator row = fluxes->begin();
      tip::Table::Record & record = *row;
      for (size_t k = 0; k < m_ebounds.size() - 1; ++row, k++) {
         for (source = sourceNames.begin(); source != sourceNames.end();
              ++source) {
            if (m_logLike.source(*source).getType() == "Point") {
               record[*source].set(srcFluxes[*source].at(k));
            }
         }
      }
   }
   delete fluxes;
}

void CountsSpectra::
writeEbounds(const std::string & outfile) const {
   std::string extname("EBOUNDS");

   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   fileSvc.appendTable(outfile, extname);
   tip::Table * ebounds(fileSvc.editTable(outfile, extname));

   ebounds->appendField("E_MIN", "1D");
   ebounds->appendField("E_MAX", "1D");
   ebounds->setNumRecords(m_ebounds.size() - 1);

   tip::Table::Iterator row = ebounds->begin();
   tip::Table::Record & record = *row;
   for (size_t k = 0; k < m_ebounds.size() - 1; ++row, k++) {
      record["E_MIN"].set(m_ebounds.at(k));
      record["E_MAX"].set(m_ebounds.at(k+1));
   }
   delete ebounds;
}

void CountsSpectra::getObsCounts(std::vector<double> & counts) const {
   if (m_binnedLike) {
      counts = m_binnedLike->countsSpectrum();
   } else {
      check_ebounds();
      counts = m_observation.eventCont().nobs(m_ebounds);
   }
}

void CountsSpectra::check_ebounds() const {
   if (m_ebounds.size() == 0) {
      throw std::runtime_error("CountsSpectra: energy bounds are not set.");
   }
}

} // namespace Likelihood
