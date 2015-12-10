/**
 * @file CountsMapBase.cxx
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/src/CountsMapBase.cxx,v 1.2 2015/03/03 06:00:00 echarles Exp $
 */

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"
#include "tip/tip_types.h"

#include "astro/SkyDir.h"
#include "astro/SkyProj.h"

#include "st_facilities/Environment.h"
#include "st_facilities/FitsImage.h"
#include "st_facilities/Util.h"

#include "facilities/commonUtilities.h"

#include "evtbin/Gti.h"
#include "evtbin/LinearBinner.h"
#include "evtbin/LogBinner.h"
#include "evtbin/OrderedBinner.h"

#include "Likelihood/CountsMapBase.h"
#include "Likelihood/HistND.h"

namespace Likelihood {

CountsMapBase::CountsMapBase(const std::string & event_file,
			     const std::string & ev_table,
			     const std::string & proj, 
			     bool use_lb,
			     double emin, double emax, unsigned long nenergies) 
   : DataProduct(event_file, ev_table, evtbin::Gti(event_file)), m_hist(0), 
     m_proj_name(proj), 
     m_use_lb(use_lb), m_proj(0) {
}

CountsMapBase::CountsMapBase(const std::string & event_file, 
			     const std::string & ev_table, 
			     const std::string & proj, 
			     bool use_lb,
			     const std::vector<double> & energies)
   : DataProduct(event_file, ev_table, evtbin::Gti(event_file)), 
     m_hist(0), m_proj_name(proj), 
     m_use_lb(use_lb), m_proj(0) {
}


CountsMapBase::CountsMapBase(const std::string & event_file, 
			     const std::string & ev_table, 
			     const std::string & proj, 
			     bool use_lb,
			     const std::vector<double> & emins,
			     const std::vector<double> & emaxs) 
   : DataProduct(event_file, ev_table, evtbin::Gti(event_file)), 
     m_hist(0), m_proj_name(proj), 
     m_use_lb(use_lb), m_proj(0) {
}

CountsMapBase::CountsMapBase(const std::string & countsMapFile) 
   : DataProduct(countsMapFile, "", evtbin::Gti(countsMapFile)),
     m_use_lb(false) {
}

void CountsMapBase::readEbounds(const std::string & countsMapFile, 
				std::vector<evtbin::Binner *> & binners) {
   std::auto_ptr<const tip::Table> 
      ebounds(tip::IFileSvc::instance().readTable(countsMapFile, "EBOUNDS"));
   tip::Table::ConstIterator it = ebounds->begin();
   tip::Table::ConstRecord & row = *it;
   std::vector<double> energies(ebounds->getNumRecords() + 1);
   double emax;
   for (int i = 0 ; it != ebounds->end(); ++it, i++) {
      row["E_MIN"].get(energies.at(i));
      row["E_MAX"].get(emax);
   }
   energies.back() = emax;

   m_energies.clear();
   for (size_t k(0); k < energies.size(); k++) {
      m_energies.push_back(energies[k]/1e3);
   }

   std::vector<evtbin::Binner::Interval> energy_intervals;
// Convert to MeV
   for (unsigned int i = 0; i < energies.size()-1; i++) {
      energy_intervals.push_back(evtbin::Binner::Interval(energies[i]/1e3, 
                                                          energies[i+1]/1e3));
   }
   binners.push_back(new evtbin::OrderedBinner(energy_intervals,
                                               "photon energy"));
}


CountsMapBase::CountsMapBase(const CountsMapBase & rhs) : 
  DataProduct(rhs),
  m_hist(rhs.m_hist->clone()),
  m_energies(rhs.m_energies),
  m_proj_name(rhs.m_proj_name),
  m_use_lb(rhs.m_use_lb),
  m_proj(rhs.m_proj ? rhs.m_proj->clone() : 0),
  m_refDir(rhs.m_refDir),
  m_tstart(rhs.m_tstart),
  m_tstop(rhs.m_tstop){
}

CountsMapBase::CountsMapBase(const CountsMapBase & rhs, 
			     unsigned int idim, unsigned int firstBin, unsigned int lastBin) :
  DataProduct(rhs),
  m_hist(rhs.m_hist->sumRange(idim,firstBin,lastBin)),
  m_proj_name(rhs.m_proj_name),
  m_use_lb(rhs.m_use_lb),
  m_proj(rhs.m_proj ? rhs.m_proj->clone() : 0),
  m_refDir(rhs.m_refDir),
  m_tstart(rhs.m_tstart),
  m_tstop(rhs.m_tstop){
  m_energies.clear();
  m_energies.push_back(rhs.energies()[firstBin]);
  m_energies.push_back(rhs.energies()[lastBin]);
}

CountsMapBase::~CountsMapBase() throw() { 
   try {
      delete m_proj;
      delete m_hist;
   } catch (std::exception & eObj) {
      std::cerr << eObj.what() << std::endl;
   } catch (...) {
   }
}



void CountsMapBase::setImage(const std::vector<float> & image) {
   m_hist->setData(image);
}
void CountsMapBase::setImage(const std::vector<double> & image) {
   m_hist->setData(image);
}

long CountsMapBase::imageDimension(int i) const {
   const evtbin::Hist::BinnerCont_t & binners = m_hist->getBinners();
   if (i < 0 || i > 2) {
      throw std::invalid_argument("CountsMapBase::imageDimension:\n"
                                  "Invalid image dimension value.");
   }
   return binners[i]->getNumBins();
}

void CountsMapBase::getAxisVector(int i, std::vector<double> & axisVector) const {
   const evtbin::Hist::BinnerCont_t & binners = m_hist->getBinners();
   if (i < 0 || i >= binners.size()) {
      throw std::invalid_argument("CountsMapBase::getAxisVector:\n"
                                  "Invalid image dimension value.");
   }
   axisVector.clear();
   for (long j = 0; j < binners[i]->getNumBins(); j++) {
      axisVector.push_back(binners[i]->getInterval(j).begin());
   }
   long jj = binners[i]->getNumBins() - 1;
   axisVector.push_back(binners[i]->getInterval(jj).end());
}

void CountsMapBase::getEnergies(std::vector<double>& energies) const {
  energies.clear();
  energies.reserve(m_energies.size());
  for ( std::vector<double>::const_iterator itr = m_energies.begin(); 
	itr != m_energies.end(); itr ++ ) {
    energies.push_back(*itr);
  }
  // EAC this doesn't work
  // std::copy(m_energies.begin(),m_energies.end(),energies.begin());
}

const std::vector<Pixel> & CountsMapBase::pixels() const {
   // EAC, this is not optimal for healpix, since all the pixels have the same area, 
   // but it works for now.
   if (m_pixels.empty()) {
      std::vector<astro::SkyDir> pixelDirs;
      std::vector<double> solidAngles;
      getPixels(pixelDirs, solidAngles);
      m_pixels.reserve(pixelDirs.size());
      for (unsigned int i = 0; i < pixelDirs.size(); i++) {
	  m_pixels.push_back(Pixel(pixelDirs[i], solidAngles[i], m_proj));
      }
   }
   return m_pixels;
}

void CountsMapBase::setRefDir(double val1, double val2) {
   m_refDir = astro::SkyDir(val1,val2,
			    m_use_lb ? astro::SkyDir::GALACTIC : astro::SkyDir::EQUATORIAL);
}

void CountsMapBase::setDataDir() {
   m_data_dir = st_facilities::Environment::dataPath("Likelihood");
}

void CountsMapBase::deleteBinners(std::vector<evtbin::Binner *> & binners) const {
   for (std::vector<evtbin::Binner *>::reverse_iterator it = binners.rbegin();
        it != binners.rend(); ++it) {
      delete *it;
   }
}

void CountsMapBase::readTimeKeywords(const std::string& event_file) {
  // Read TSTART and TSTOP keywords from event file header.
   const tip::Table * events = 
      tip::IFileSvc::instance().readTable(event_file, "EVENTS");
   const tip::Header & header(events->getHeader());
   header["TSTART"].get(m_tstart);
   header["TSTOP"].get(m_tstop);
   delete events;
}

} // namespace Likelihood
