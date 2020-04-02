/**
 * @file CountsMapHealpix.cxx
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/CountsMapHealpix.cxx,v 1.5 2017/08/18 22:46:12 echarles Exp $
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
#include "evtbin/HealpixBinner.h"

#include "healpix/HealpixRegion.h"

#include "Likelihood/CountsMapHealpix.h"
#include "Likelihood/HistND.h"
#include "Likelihood/HealpixProjMap.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/FileUtils.h"

namespace Likelihood {


  CountsMapHealpix::CountsMapHealpix(const std::string & countsMapFile) 
    : CountsMapBase(countsMapFile),
      m_healpixProj(0),
      m_solidAngle(0.),
      m_pixelSize(0.),
      m_nPixels(0) {
    readKeywords(countsMapFile);
    std::vector<evtbin::Binner *> binners;  
    binners.push_back(m_hpx_binner); // ! watch out, we own this one, remove it before we delete the binners
    m_nPixels = m_hpx_binner->getNumBins();
    readEbounds(countsMapFile, binners); // This make the energy binner
    readImageData(countsMapFile, binners);
    setDataDir();
    //! Set the front of the binners list to a null pointer, so we don't accidentally delete the HEALPix binner
    binners[0] = 0;
    deleteBinners(binners);
    latchCacheData();
  }

  CountsMapHealpix::CountsMapHealpix(const CountsMapHealpix & rhs) 
    : CountsMapBase(rhs),
      m_solidAngle(rhs.m_solidAngle),
      m_pixelSize(rhs.m_pixelSize),
      m_nPixels(rhs.m_nPixels){
    m_healpixProj = static_cast<astro::HealpixProj*>(m_proj);
    m_hpx_binner = const_cast<evtbin::HealpixBinner*>(static_cast<const evtbin::HealpixBinner*>(m_hist->getBinners()[0]));    
  }

  CountsMapHealpix::CountsMapHealpix(const CountsMapHealpix & rhs,
				     unsigned int firstBin, unsigned int lastBin): 
    CountsMapBase(rhs,1,firstBin,lastBin),
    m_solidAngle(rhs.m_solidAngle),
    m_pixelSize(rhs.m_pixelSize),
    m_nPixels(rhs.m_nPixels){    
    m_healpixProj = static_cast<astro::HealpixProj*>(m_proj);
    m_hpx_binner = const_cast<evtbin::HealpixBinner*>(static_cast<const evtbin::HealpixBinner*>(m_hist-> getBinners()[0]));
  }

  CountsMapHealpix::CountsMapHealpix(const HealpixProjMap& projMap, 
				     const CountsMapHealpix& counts_map):
    CountsMapBase(counts_map),
    m_solidAngle(counts_map.m_solidAngle),
    m_pixelSize(counts_map.m_pixelSize),
    m_nPixels(counts_map.m_nPixels){    
    m_healpixProj = static_cast<astro::HealpixProj*>(m_proj);
    m_hpx_binner = const_cast<evtbin::HealpixBinner*>(static_cast<const evtbin::HealpixBinner*>(m_hist-> getBinners()[0]));
    std::vector<double> energyBinWidths(energies().size() - 1);
    for ( int k(0); k < energies().size(); k++ ) {
      energyBinWidths[k] = energies()[k+1] - energies()[k];
    }
    HealpixProjMap::convertToIntegral(m_hist->data_access(), projMap.image(), energyBinWidths, 
				      nPixels(), m_solidAngle);				      
  }

  CountsMapHealpix::~CountsMapHealpix() throw() {;}

  ProjMap* CountsMapHealpix::makeProjMap(CountsMapBase::ConversionType cType) const {
    return new HealpixProjMap(*this, cType);
  }

  CountsMapBase* CountsMapHealpix::makeBkgEffMap(const MeanPsf & psf, const float& efact) const {

    CountsMapHealpix* outMap = new CountsMapHealpix(*this);     
    ProjMap* projMap = makeProjMap();
    ProjMap* convMap = projMap->convolveAll(psf);
    delete projMap;
    HealpixProjMap* convMap_hpx = convMap->cast_healpix();
    
    const std::vector<double>& ebins = energies();
    std::vector<double> energyBinWidths(num_ebins());
    std::vector<double> energyBinMeans(num_ebins());
    CountsMapBase::fillEnergyBinWidths(energyBinWidths, ebins);
    CountsMapBase::fillEnergyBinGeomCenters(energyBinMeans, ebins);

    std::vector<float>& outData = outMap->m_hist->data_access();

    int idx_fill(0);
    size_t kStep = nPixels();

    for ( size_t k(0); k < num_ebins(); k++ ) {
      double mean_energy = energyBinMeans[k];
      double psf_peak = psf.peakValue(mean_energy);

      // The factor we need to convert back to counts is
      // double factor1 = energyBinWidths[k]*solidAngle();
      // (This is b/c the proj map used solid angle in sr)

      // The factor we need account for the PSF in the correct units is 
      // double factor2 = 1 / psf_peak*solidAngle();
      // (This is b/c this class is using solidAngle() in sr)

      // Combining these we get
      // double factor = factor1*factor2 = energyBinWidths[k]/psf_peak
      double factor = energyBinWidths[k]/psf_peak;
      const Healpix_Map<float>& conv_image = convMap_hpx->image()[k];
      for ( size_t i(0); i < kStep; i++, idx_fill++  ) {
	float addend = conv_image[i] * factor;
	// Zero out the output data from this pixel / energy.
	outData[idx_fill] = 0.;
	// Add this quantity to each of the energy layers below this one
	// Note that fillIt is counting DOWN towards zero	
	double emin_integ = mean_energy/efact;
	int kinteg(k);
	for ( int fillIt(idx_fill); fillIt >= 0; fillIt -= kStep, kinteg-=1 ) {
	  if ( energyBinMeans[kinteg] < emin_integ) break;
	  outData[fillIt] += addend;
	}
      }
    }
    // clean up    
    delete convMap_hpx;
    return outMap;
  }


  void CountsMapHealpix::binInput(tip::Table::ConstIterator begin, 
				  tip::Table::ConstIterator end) {

    const evtbin::Hist::BinnerCont_t & binners = m_hist->getBinners();

    std::string field1 = m_proj->isGalactic() ? "L" : "RA";
    std::string field2 = m_proj->isGalactic() ? "B" : "DEC";
   
    // Fill histogram, converting each RA/DEC to Sky X/Y on the fly:
    std::vector<double> values(2,0);
    for (tip::Table::ConstIterator itor = begin; itor != end; ++itor) {
      double s1 = (*itor)[field1].get();
      double s2 = (*itor)[field2].get();
      double energy = (*itor)["ENERGY"].get();
      values[0] = m_hpx_binner->computeIndex(s1,s2);
      values[1] = binners[1]->computeIndex(energy);  
      m_hist->fillBin(values);
    }
  }

  bool CountsMapHealpix::withinBounds(const astro::SkyDir & dir, double energy,
				      long border_size) const {
    
    double s1 = m_proj->isGalactic() ? dir.l() : dir.ra();
    double s2 = m_proj->isGalactic() ? dir.b() : dir.dec();
    const evtbin::Hist::BinnerCont_t& binners = m_hist->getBinners();

    static std::vector<double> values(2,0);
    values[0] = m_hpx_binner->computeIndex(s1,s2);
    values[1] = binners[1]->computeIndex(energy);
    long indx = m_hist->binIndex(values, border_size);
    return indx >= 0;
  }

  void CountsMapHealpix::writeOutput(const std::string & creator, 
				     const std::string & out_file) const {
    static const std::string ext="SKYMAP";
    createFile(creator, out_file, 
	       facilities::commonUtilities::joinPath(m_data_dir,
						     "LatHealpixTemplate"));
    tip::Table *table = tip::IFileSvc::instance().editTable(out_file, ext);  
    tip::Header & header(table->getHeader());
    FileUtils::replace_image_from_hist_hpx(*table, "CHANNEL", *m_hist, *m_hpx_binner);
    setKeywords(header);
    delete table;
    const evtbin::Hist::BinnerCont_t & binners = m_hist->getBinners();
    if ( binners.size() > 1 ) {
      writeEbounds(out_file,binners[1]);
    }
    writeGti(out_file);
  }

  void CountsMapHealpix::writeAsWeightsMap(const std::string & creator, 
					   const std::string & out_file) const {
    static const std::string ext="SKYMAP";
    createFile(creator, out_file, 
	       facilities::commonUtilities::joinPath(m_data_dir,
						     "LatHealpixWeightsTemplate"));
    tip::Table *table = tip::IFileSvc::instance().editTable(out_file, ext);  
    tip::Header & header(table->getHeader());
    FileUtils::replace_image_from_hist_hpx(*table, "ENERGY", *m_hist, *m_hpx_binner);
    setKeywords(header);
    delete table;

    std::vector<double> energyBinCenters;
    getEnergyBinGeomCenters(energyBinCenters);
    
    tip::Extension* energiesHdu = FileUtils::replace_energies(out_file, "ENERGIES", energyBinCenters);
    delete energiesHdu;
    writeGti(out_file);
    delete table;
  }

  void CountsMapHealpix::setKeywords(tip::Header & header) const {
    //these should perhaps stay in the Primary header
    astro::JulianDate current_time = st_facilities::Util::currentTime();
    header["DATE"].set(current_time.getGregorianDate());
    header["DATE-OBS"].set("");
    header["DATE-END"].set("");
    m_hpx_binner->setKeywords(header); 
  }



void CountsMapHealpix::
getBoundaryPixelDirs(std::vector<astro::SkyDir> & pixelDirs) const {
  // EAC_FIX, HEALPIX impl of getBoundaryPixelDirs missing throws std::runtime_error
  throw std::runtime_error("CountsMapHealpix::getBoundaryPixelDirs"
			   "is not implemented");  
}

void CountsMapHealpix::getPixels(std::vector<astro::SkyDir> & pixelDirs,
				 std::vector<double> & pixelSolidAngles) const {


  pixelDirs.clear();
  pixelSolidAngles.clear();
  
  pixelDirs.reserve(m_nPixels);
  pixelSolidAngles.reserve(m_nPixels);

  for ( unsigned int iPix(0); iPix < m_nPixels; iPix++ ) {
    try {
      // EAC, use the HealpixBinner account for the 
      // possiblity that we have a partial-sky CountsMap, 
      // In that case, we want to remap the index
      int gloPix = m_hpx_binner->allSky() ? iPix : m_hpx_binner->pixelIndices()[iPix];
      astro::SkyDir my_dir(double(gloPix), 0, projection());
      pixelDirs.push_back(my_dir);
      pixelSolidAngles.push_back(m_solidAngle);
    } catch (std::exception & eObj) {
      if (st_facilities::Util::
	  expectedException(eObj, "SkyProj wcslib error")) {
	// Attempted to create a pixel outside of the map boundary for
	// the current projection, so fill solid angle with zero
	// and place default direction in pixel dir.
	pixelDirs.push_back(astro::SkyDir(0, 0));
	pixelSolidAngles.push_back(0);
      } else {
	throw;
      }
    }
  }
}


void CountsMapHealpix::readKeywords(const std::string & countsMapFile) {
  static const std::string extName("SKYMAP");
  m_healpixProj = new astro::HealpixProj(countsMapFile,extName);
  m_hpx_binner = new evtbin::HealpixBinner(countsMapFile,extName);
  m_proj = m_healpixProj;
  m_proj_name = m_proj->projType();
  m_use_lb = m_healpixProj->isGalactic();
  latchCacheData();
}


void CountsMapHealpix::readImageData(const std::string & countsMapFile,
				     std::vector<evtbin::Binner *> & binners) {
  m_hist = new HistND(binners);
  static const std::string extName("SKYMAP");
  std::unique_ptr<const tip::Table> 
    table(tip::IFileSvc::instance().readTable(countsMapFile,extName));

  // This is a bit tricky, basically all the data we care about
  // are in the columns called "CHANNELx"
  // Note also that tip work in lowercase
  std::vector<tip::FieldIndex_t> dataColumns;
  const tip::Table::FieldCont& colNames = table->getValidFields();
  for ( tip::Table::FieldCont::const_iterator itr = colNames.begin(); 
        itr != colNames.end(); itr++ ) {
    if ( itr->find("channel") == 0 ) { 
      dataColumns.push_back( table->getFieldIndex(*itr) );     
    } else {
      continue;
    }
  }

  int ncol = dataColumns.size();
  tip::Index_t nrow = table->getNumRecords();
  m_nPixels = nrow;

  std::vector<float> oneColumn(nrow,0);
  double readVal(0.);
  // This keeps track of the energy slice we are filling.  
  // Starts with slice 0
  std::vector<unsigned int> ivalues(2,0);    
  for ( std::vector<tip::FieldIndex_t>::const_iterator itrData = dataColumns.begin();
        itrData != dataColumns.end(); itrData++, ivalues[1]++ ) {    
    const tip::IColumn* col = table->getColumn(*itrData);
    for ( tip::Index_t irow(0); irow < nrow; irow++ ) {
      col->get(irow,readVal);
      oneColumn[irow] = readVal;
    }    
    m_hist->setSlice(0,ivalues,oneColumn);
  }
}


int CountsMapHealpix::globalToLocalIndex(int glo) const {
  if ( allSky() ) return glo;
  std::map<int,int>::const_iterator itrFind = m_hpx_binner->pixGlobalToLocalMap().find(glo);
  return itrFind == m_hpx_binner->pixGlobalToLocalMap().end() ? -1 : itrFind->second;
}
   
int CountsMapHealpix::localToGlobalIndex(int loc) const {
  if ( allSky() ) return loc;
  if ( loc < 0 || loc >= m_hpx_binner->pixelIndices().size() ) return -1;
  return m_hpx_binner->pixelIndices()[loc];
}      
  

const std::vector<int>& CountsMapHealpix::pixelIndices() const {
  return m_hpx_binner->pixelIndices();
}


void CountsMapHealpix::latchCacheData() {
  // Total number of pixels
  int nPix = m_healpixProj->healpix().Npix();
  // Solid angle in SR = 4pi/npix
  // B/c of inconsisent interfaces,  we want the solid angle in SR, but the pixel size in deg.    
  m_solidAngle = ASTRO_4PI / float(nPix);
  m_pixelSize = astro::radToDeg(sqrt(m_solidAngle));

  double c1(0.);
  double c2(0.);  
  healpix::HealpixRegion* hpx_reg = const_cast<healpix::HealpixRegion*>(m_hpx_binner->region());
  if ( hpx_reg != 0 ) {
    hpx_reg->getRegionSize(m_mapRadius);
    hpx_reg->getRefDir(c1,c2);
  } else {
    m_mapRadius = 180.;
  }
  setRefDir(c1,c2);
}
  
} // namespace Likelihood
