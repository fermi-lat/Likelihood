/**
 * @file CountsMapBase.cxx
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/CountsMapBase.cxx,v 1.2 2017/09/29 01:44:18 echarles Exp $
 */

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <cstdio>
#include <sstream>

#include "facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"
#include "tip/tip_types.h"
#include "tip/Header.h"

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
#include "Likelihood/FileUtils.h"
#include "Likelihood/AppHelpers.h"

namespace Likelihood {

  void CountsMapBase::fillEnergyBinWidths(std::vector<double>& energyBinWidths,
					  const std::vector<double>& energyBinEdges) {
    size_t n_energies = energyBinEdges.size();
    if ( n_energies == 0 ) return;
    size_t nebins = n_energies -1;
    energyBinWidths.resize(nebins);
    for ( size_t k(0); k < nebins; k++) {
      energyBinWidths[k] = energyBinEdges[k+1] - energyBinEdges[k];
    }
  }

  void CountsMapBase::fillEnergyBinGeomCenters(std::vector<double>& energyBinCenters,
					       const std::vector<double>& energyBinEdges) {
    size_t n_energies = energyBinEdges.size();
    if ( n_energies == 0 ) return;
    size_t nebins = n_energies -1;
    energyBinCenters.resize(nebins);
    for ( size_t k(0); k < nebins; k++) {
      energyBinCenters[k] = sqrt( energyBinEdges[k+1]*energyBinEdges[k]);
    }
  }
  
  void CountsMapBase::getBMinAndSumMinOverK2(float& bmin, float& sumk2,
					     std::vector<std::vector<float>::const_iterator >& itrs) {
    bmin = 1e99;
    sumk2 = 0;
    std::vector<std::vector<float>::const_iterator >::iterator itr2 = itrs.begin();
    for ( ; itr2 != itrs.end(); itr2++ ) {
      std::vector<float>::const_iterator itr = *itr2;
      if ( bmin > *itr ) {
	bmin = *itr;
      }
    }
    if ( bmin <= 0 ) {
      sumk2 = 0;
      return;
    }
    for ( itr2 = itrs.begin(); itr2 != itrs.end(); itr2++ ) {
      std::vector<float>::const_iterator itr = *itr2;
      float addend = bmin / *itr;
      addend *= addend;	
      sumk2 += addend;
    }
  }
  
  void CountsMapBase::getAlphaWts(float& alpha, const float& epsilon2, 
				  std::vector<std::vector<float>::const_iterator >& itrs) {
    float bmin(0.);
    float sumk2(0.);
    getBMinAndSumMinOverK2(bmin, sumk2, itrs);
    alpha = 1 + epsilon2*bmin;
    alpha /= (1 + epsilon2*bmin*sumk2);
  }
  
  void CountsMapBase::getAlphaVector(std::vector<float>& alphaVect, const float& epsilon2, 
				     const std::vector<const std::vector<float>* >& beffVects) {
    std::vector<std::vector<float>::const_iterator > itrs;
    for ( std::vector<const std::vector<float>* >::const_iterator itrVect = beffVects.begin(); 
	  itrVect != beffVects.end(); itrVect++ ) {
      itrs.push_back((*itrVect)->begin());
    }

    int ndata = alphaVect.size();
    int nprint = ndata/20;
    size_t idx(0);
    std::cout << "Computing alpha map: " << std::flush;
    for ( std::vector<float>::iterator itr = alphaVect.begin(); itr!= alphaVect.end(); itr++, idx++ ) {
      if ( idx % nprint == 0 ) {
	std::cout << '.' << std::flush;
      }
      getAlphaWts(*itr, epsilon2, itrs);
      // This loop steps all the iterator of the beffVects
      for ( std::vector<std::vector<float>::const_iterator >::iterator itr2 = itrs.begin();
	    itr2 != itrs.end(); itr2++ ) {
	(*itr2)++;
      }      
    }
    std::cout << '!' << std::endl;
  }

  void CountsMapBase::getWts(std::vector<float>& wts,
			     const float& epsilon2, 
			     const std::vector<float>& alphaVect,
			     const std::vector<float>& beffVect) {
    size_t n = beffVect.size();
    int nprint = n/20;
    size_t idx(0);
    std::vector<float>::iterator itrWts = wts.begin();
    std::vector<float>::const_iterator itrAlpha = alphaVect.begin();
    std::vector<float>::const_iterator itrBeff = beffVect.begin();
    std::cout << "Computing wts map: " << std::flush;
    for ( ; itrWts != wts.end(); itrWts++, itrBeff++, idx++ ) {
      if ( idx % nprint == 0 ) {
	std::cout << '.' << std::flush;
      }
      *itrWts = 1 / ( 1 + (epsilon2 * (*itrBeff) ) );
      if ( itrAlpha != alphaVect.end() ) {
	*itrWts *= (*itrAlpha);
	itrAlpha++;
      }
    }
    std::cout << '!' << std::endl;
  }

  CountsMapBase* CountsMapBase::makeAlphaMap(const float& epsilon2, const std::vector<CountsMapBase*>& input_maps) {
    if (  input_maps.size() == 0 ) { return 0; }
    CountsMapBase* outMap = input_maps[0]->clone();
    std::vector<const std::vector<float>* > input_vects;
    for ( std::vector<CountsMapBase*>::const_iterator itr_in = input_maps.begin(); itr_in != input_maps.end();
	  itr_in++ ) {
      const CountsMapBase* in_map = *itr_in;
      input_vects.push_back(&in_map->data());      
    }
    getAlphaVector(outMap->data_access(), epsilon2, input_vects);
    return outMap;
  }


  CountsMapBase* CountsMapBase::makeAlphaMap(const float& epsilon2, const std::vector<std::string>& input_map_files) {
    std::vector<Likelihood::CountsMapBase*> inputMaps;
    for ( std::vector<std::string>::const_iterator itr = input_map_files.begin();
	 itr != input_map_files.end(); itr++ ) {
      CountsMapBase* cmap = Likelihood::AppHelpers::readCountsMap(*itr); 
      inputMaps.push_back(cmap);
    }
    return CountsMapBase::makeAlphaMap(epsilon2, inputMaps);
  }
  
  CountsMapBase* CountsMapBase::makeWtsMap(const float& epsilon2, 
					   const CountsMapBase* alphaMap, 
					   const CountsMapBase& beffMap) {
    CountsMapBase* outMap = beffMap.clone();
    static const std::vector<float> nullVector;
    const std::vector<float>& alphaVector =  alphaMap != 0 ? alphaMap->data() : nullVector;
    getWts(outMap->data_access(), epsilon2, alphaVector, beffMap.data());
    return outMap;
  }

  void CountsMapBase::copyAndUpdateDssKeywords(const std::string& infile,
					       const std::string& outfile,
					       AppHelpers* helper,
					       const std::string& irfs){
    
    dataSubselector::Cuts my_cuts(infile, "PRIMARY", false, false, false);
    // Ensure that the irfs used are written to the DSS keywords.
    if ( helper != 0 ) {
      my_cuts.addVersionCut("IRF_VERSION", helper->irfsName());
    }
    
    tip::Image * my_image = tip::IFileSvc::instance().editImage(outfile, "");
    if (irfs != "CALDB" && helper != 0 ) {
      helper->setBitMaskCuts(my_cuts);
    }
    my_cuts.writeDssKeywords(my_image->getHeader());
    delete my_image;
  }
  
  void CountsMapBase::addBkgEffKeywords(const std::string& outfile,
					const std::string& inputmap,
					const float& efact){
    tip::Image* my_image = tip::IFileSvc::instance().editImage(outfile, "");
    tip::Header& my_header = my_image->getHeader();
    my_header["INPUTMAP"].set(inputmap);
    my_header["EFACT"].set(efact);
    delete my_image;
  }
  
  void CountsMapBase::addAlphaMapKeywords(const std::string& outfile,
					  double epsilon,
					  const std::vector<std::string>& inputFiles) {
    tip::Image* my_image = tip::IFileSvc::instance().editImage(outfile, "");
    tip::Header& my_header = my_image->getHeader();
    my_header["EPSILON"].set(epsilon);
    for ( size_t i(0); i < inputFiles.size(); i++ ) {
      std::ostringstream cardName;
      cardName<<"BKGMAP"<<i+1;
      my_header[cardName.str()].set(inputFiles[i]);
    }
    delete my_image;    
  }
  
  void CountsMapBase::addWtsMapKeywords(const std::string& outfile,
					double epsilon,
					const std::string& bkgmap,
					const std::string& alphamap) {
    tip::Image* my_image = tip::IFileSvc::instance().editImage(outfile, "");
    tip::Header& my_header = my_image->getHeader();
    my_header["EPSILON"].set(epsilon);
    my_header["BKGMAP"].set(bkgmap);
    my_header["ALPHAMAP"].set(alphamap);
    delete my_image;    
  }
    

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

   FileUtils::read_ebounds_to_vector(countsMapFile, m_energies);

   std::vector<evtbin::Binner::Interval> energy_intervals;

// Convert to MeV
   for (unsigned int i = 0; i < m_energies.size()-1; i++) {
      energy_intervals.push_back(evtbin::Binner::Interval(m_energies[i], 
                                                          m_energies[i+1]));
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

void CountsMapBase::writeEmptyOutput(const std::string & creator, const std::string & out_file) const {

   createFile(creator, out_file, 
              facilities::commonUtilities::joinPath(m_data_dir,"LatCountsMapTemplate"));

   const evtbin::Hist::BinnerCont_t & binners = m_hist->getBinners();
   if ( binners.size() > 1 ) {
     writeEbounds(out_file,binners[1]);
   }
   writeGti(out_file);
}

void CountsMapBase::writeEnergies(const std::string & creator, 
				  const std::string & out_file,
				  int kmin, int kmax) const {
   std::remove(out_file.c_str());     
   std::string ext("PRIMARY");
   tip::ImageBase::PixelCoordinate null;
   tip::IFileSvc::instance().appendImage(out_file, ext, null);

   ext = "ENERGIES";
   tip::IFileSvc::instance().appendTable(out_file, ext);
   tip::Table * table = tip::IFileSvc::instance().editTable(out_file, ext);
   table->appendField("Energy", "1D");

   kmax = kmax < 0 ? m_energies.size() : kmax;
   size_t nebins = kmax - kmin;

   table->setNumRecords(nebins);

   tip::Table::Iterator row = table->begin();
   tip::Table::Record & record = *row;

   for ( size_t ie(kmin); ie != kmax; ++ie, ++row ) {
     record["Energy"].set(m_energies[ie]);
   }
   delete table;
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

void CountsMapBase::getEnergyBinGeomCenters(std::vector<double>& energyBinCenters) const {
  fillEnergyBinGeomCenters(energyBinCenters, m_energies);
}

void CountsMapBase::getEnergyBinWidths(std::vector<double>& energyBinWidths) const {
  fillEnergyBinWidths(energyBinWidths, m_energies);
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
