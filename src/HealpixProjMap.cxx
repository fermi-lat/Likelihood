/**
 * @file HealpixProjMap.cxx
 * @brief A map with reference point centered on the image and that
 * uses WCS projections for indexing its internal representation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/HealpixProjMap.cxx,v 1.1 2015/12/10 00:58:00 echarles Exp $
 */

#include <cmath>

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "st_facilities/Util.h"

#include "astro/HealpixProj.h"

#include "Likelihood/BinnedExposureBase.h"
#include "Likelihood/Convolve.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/HealpixProjMap.h"
#include "Likelihood/ConvolveHealpix.h"

namespace {
//   class Image : public std::vector< std::vector<double> > {
   class Image : public std::vector<float> {
   public:
      Image() {}
      void normalize() {
         double total(0);
         for (unsigned int i = 0; i < this->size(); i++) {
	   total += this->at(i);
         }
         for (unsigned int i = 0; i < this->size(); i++) {
	   this->at(i) /= total;
         }
      }
   };

   double my_round(double x) {
      int xint = static_cast<int>(x);
      if (x - xint >= 0.5) {
         return xint + 1.;
      }
      return xint;
   }
} // unnamed namespace

namespace Likelihood {

HealpixProjMap::HealpixProjMap(const std::string & filename,
			       const std::string & extension,
			       bool interpolate,
			       bool enforceEnergyRange) 
  : ProjMap(filename,interpolate,enforceEnergyRange),
    m_healpixProj( new astro::HealpixProj(filename, extension) ) {

  const tip::Table* table = 
      tip::IFileSvc::instance().readTable(filename, extension.empty() ?  "SKYMAP" : extension);
  const tip::Header & header = table->getHeader();
  
  double refdir1(0.);
  double refdir2(0.);
  try {
    header["REFDIR1"].get(refdir1);
    header["REFDIR2"].get(refdir2);
  } catch (...) {
    ;
  }

  astro::SkyDir refDir(refdir1,refdir2,
		       m_healpixProj->isGalactic() ? astro::SkyDir::GALACTIC : astro::SkyDir::EQUATORIAL );
  setProjInfo(refDir,*m_healpixProj);
  latchCacheData();

  int naxes(0);
  int numEnergies(1);
  int numPixels(0);  
  header["NAXIS"].get(naxes);
  header["NAXIS1"].get(numEnergies);
  if ( naxes == 2 ) {
    header["NAXIS2"].get(numPixels);
    // EAC_FIX. There is a bug that we are actually writing 8* the number of energy planes into the FITS keywork
    numEnergies /= 8;
    ExposureMap::readEnergyExtension(filename, energies_access());
  } else {
    // In this case axis1 actually has the number of pixels
    numPixels = numEnergies;
    numEnergies = 1;
    energies_access().push_back(100.);
  }

  // This is a bit tricky, basically all the data we care about
  // are in the columns called "ENERGYx"
  // Note also that tip work in lowercase
  std::vector<tip::FieldIndex_t> dataColumns;
  const tip::Table::FieldCont& colNames = table->getValidFields();
  for ( tip::Table::FieldCont::const_iterator itr = colNames.begin(); 
	itr != colNames.end(); itr++ ) {
    if ( itr->find("energy") == 0 ) { 
      dataColumns.push_back( table->getFieldIndex(*itr) );     
    } else {
      continue;
    }
  }

  int ncol = dataColumns.size();
  tip::Index_t nrow = table->getNumRecords();

  m_image.clear();
  m_image.resize(ncol);
  int idx(0);
  for ( std::vector<tip::FieldIndex_t>::const_iterator itrData = dataColumns.begin();
	itrData != dataColumns.end(); itrData++, idx++ ) {
    ImagePlane_t& image = m_image[idx];
    // EAC_FIX. Until we implement partial sky HEALPix maps we can use this to build the map
    image.SetNside(m_healpixProj->healpix().Nside(),m_healpixProj->healpix().Scheme());
    const tip::IColumn* col = table->getColumn(*itrData);
    for ( tip::Index_t irow(0); irow < nrow; irow++ ) {
      col->get(irow,image[irow]);
    }    
  }
  // EAC_FIX. Until we implement partial sky HEALPix maps the radius is 180 degrees (all-sky)
  setMapRadius(180.);
  computeMapIntegrals();
}

HealpixProjMap::HealpixProjMap(const DiffuseSource & diffuseSource, 
			       int nside, Healpix_Ordering_Scheme scheme, const nside_dummy dummy,
			       double energy, bool use_lb,
			       double radius, double ra, double dec,
			       bool interpolate, bool enforceEnergyRange) 
  : ProjMap("",interpolate,enforceEnergyRange),
    m_healpixProj( new astro::HealpixProj(nside,scheme,dummy,use_lb) ) {
  
  astro::SkyDir refDir(ra,dec,
		       use_lb ? astro::SkyDir::GALACTIC : astro::SkyDir::EQUATORIAL );
  setProjInfo(refDir,*m_healpixProj);
  latchCacheData();
  const int nPix = nPixels();
  
  ImagePlane_t image_plane(m_healpixProj->healpix().Nside(),m_healpixProj->healpix().Scheme(),dummy);
  // EAC_FIX, until we implement partial sky HEALPix maps the radius is 180 degrees (all-sky)
  setMapRadius(180.);

  // Fill the image_plane by looping over the pixels
  for ( int iLoc(0); iLoc < nPix; iLoc++ ) {
    int iglo = localToGlobal(iLoc);
    if (m_healpixProj->testpix2sph(iglo,0.) == 0) {
      std::pair<double, double> coord = m_healpixProj->pix2sph(iglo,0.);
      astro::SkyDir dir(coord.first, coord.second, 
			use_lb ? astro::SkyDir::GALACTIC : astro::SkyDir::EQUATORIAL );
      SkyDirArg my_dir(dir, energy);
      image_plane[iLoc] = diffuseSource.spatialDist(my_dir);
    } else {
      continue;
    }
  }
  check_negative_pixels(image_plane);
  m_image.clear();
  m_image.push_back(image_plane);
  energies_access().push_back(energy);  
  computeMapIntegrals();
}


HealpixProjMap::HealpixProjMap(const DiffuseSource & diffuseSource, 
			       int order, Healpix_Ordering_Scheme scheme, 
			       double energy, bool use_lb,
			       double radius, double ra, double dec,
			       bool interpolate, bool enforceEnergyRange) 
  : ProjMap("",interpolate,enforceEnergyRange),
    m_healpixProj( new astro::HealpixProj(order,scheme,use_lb) ) {
  
  astro::SkyDir refDir(ra,dec,
		       use_lb ? astro::SkyDir::GALACTIC : astro::SkyDir::EQUATORIAL );
  setProjInfo(refDir,*m_healpixProj);
  latchCacheData();
  const int nPix = nPixels();
  
  ImagePlane_t image_plane(m_healpixProj->healpix().Order(),m_healpixProj->healpix().Scheme());
  // EAC_FIX. Until we implement partial sky HEALPix maps the radius is 180 degrees (all-sky)
  setMapRadius(180.);

  // Fill the image_plane by looping over the pixels
  for ( int iLoc(0); iLoc < nPix; iLoc++ ) {
    int iglo = localToGlobal(iLoc);
    if (m_healpixProj->testpix2sph(iglo,0.) == 0) {
      std::pair<double, double> coord = m_healpixProj->pix2sph(iglo,0.);
      astro::SkyDir dir(coord.first, coord.second, 
			use_lb ? astro::SkyDir::GALACTIC : astro::SkyDir::EQUATORIAL );
      SkyDirArg my_dir(dir, energy);
      image_plane[iLoc] = diffuseSource.spatialDist(my_dir);
    } else {
      continue;
    }
  }
  check_negative_pixels(image_plane);
  m_image.clear();
  m_image.push_back(image_plane);
  energies_access().push_back(energy);  
  computeMapIntegrals();
}

HealpixProjMap::~HealpixProjMap(){}

HealpixProjMap::HealpixProjMap(const HealpixProjMap & rhs) 
  :  ProjMap(rhs),
     m_image(rhs.m_image), 
     m_solidAngle(rhs.m_solidAngle),
     m_pixelSize(rhs.m_pixelSize){
  m_healpixProj = (astro::HealpixProj*)getProj();
}

HealpixProjMap::HealpixProjMap(const HealpixProjMap & rhs, const double& energy, const Healpix_Map<float>& image) 
  :  ProjMap(rhs),
     m_solidAngle(rhs.m_solidAngle),
     m_pixelSize(rhs.m_pixelSize){
  m_healpixProj = (astro::HealpixProj*)getProj();
  m_image.push_back(image);
}

HealpixProjMap & HealpixProjMap::operator=(const HealpixProjMap & rhs) {
   if (this != &rhs) {
      ProjMap::operator=(rhs);
      m_image = rhs.m_image;
      m_solidAngle = rhs.m_solidAngle;
      m_pixelSize = rhs.m_pixelSize;
      m_healpixProj = (astro::HealpixProj*)getProj();
   }
   return *this;
}

double HealpixProjMap::operator()(const astro::SkyDir & dir, int k) const {
   check_energy_index(k);

   double theta = astro::degToRad ( astro::latToTheta_Deg(  m_healpixProj->isGalactic() ? dir.b() : dir.dec() ) );
   double phi = astro::degToRad( m_healpixProj->isGalactic() ? dir.l() : dir.ra() );
   const pointing ang(theta,phi);
   
   try {
     if ( getInterpolate() ) {
       return m_image[k].interpolated_value(ang);
     } 
     int loc = globalToLocal( m_image[k].ang2pix(ang) );
     return m_image[k][loc];
   } catch (...) {
     ;
   }
   return 0;
}
 

double HealpixProjMap::operator()(const astro::SkyDir & dir, double energy) const {
  // I guess this is to catch default values
  if (energy < 0) {
    energy = energies().front();
  }
  check_energy(energy);

  int k(0);
  if ( energies().size() > 1) {
    k = std::upper_bound(energies().begin(), energies().end(), energy) - energies().begin() - 1;
    /// Extrapolate beyond highest energy.  This will only occur if
    /// m_enforceEnergyRange == false.
    if (k > static_cast<int>(energies().size() - 2)) {
      k = energies().size() - 2;
      extrapolated_access() += 1;
    }
  }
  // This is a bit inefficient, it might be better to get the pixel index
  // just once for both energies.
  // Leave it as is for now.
  double y1 = operator()(dir, k);
  if (energy == energies()[k]) { 
    return y1;
  }
  double y2 = operator()(dir, k+1);
  // EAC, FIXME, HEALPix can very slightly overshoot interpolation
  // this is a problem if the map has zeros in it, as you 
  // will get negative numbers and crash in interpolatePowerLaw
  static const double almost_zero(-1e-16);
  if ( y1 < 0 && y1 > almost_zero ) y1 = 0.;
  if ( y2 < 0 && y2 > almost_zero ) y2 = 0.;  
  
  double value = interpolatePowerLaw(energy, energies()[k],
				     energies()[k+1], y1, y2);
  return value;
}

ProjMap* HealpixProjMap::convolve(double energy, const MeanPsf & psf,
				  const BinnedExposureBase & exposure,
				  bool performConvolution,
				  int k) const {
   // Convolve for a single image plane.
   check_energy_index(k);

   // Compute unconvolved counts map by multiplying intensity image by exposure.
   Healpix_Map<float> counts(m_image[k]);
   for ( int i(0); i < counts.Npix(); i++ ) {
     if (getProj()->testpix2sph(i, 0) == 0) {
       std::pair<double, double> coord = getProj()->pix2sph(i, 0);
       astro::SkyDir dir(coord.first, coord.second, 
			 getProj()->isGalactic() ? 
			 astro::SkyDir::GALACTIC : astro::SkyDir::EQUATORIAL );
       counts[i] *= exposure(energy, dir.ra(), dir.dec());
     }
   }

   if (!performConvolution) {
     HealpixProjMap* my_image = new HealpixProjMap(*this,energy,counts);
     return my_image;
   }

   double pixelSize_image = pixelSize();
   int nside_used = m_healpixProj->healpix().Nside();   
   double pixelSize_psfMin = ConvolveHealpix::psfMinPixSize(energy);
   double psfSize = 100*pixelSize_psfMin;
   bool upgrade(false);

   if ( pixelSize_image > 2*pixelSize_psfMin ) {
     nside_used *= 2;
     upgrade = true;
   } 

   Healpix_Map<float> counts_image(nside_used,RING,SET_NSIDE);
   if ( upgrade ) {
     counts_image.Import_upgrade(counts);
   } else {
     counts_image.Import_nograde(counts);
   }
   Healpix_Map<float> psf_image(nside_used,RING,SET_NSIDE);
   ConvolveHealpix::fillMapWithPSF_pole(psf,energy,psf_image);

   Healpix_Map<float> conv_image(nside_used,RING,SET_NSIDE);

   ConvolveHealpix::convolve(counts_image,psf_image,conv_image);

   Healpix_Map<float> out_image(m_healpixProj->healpix().Nside(),RING,SET_NSIDE);
   if ( upgrade ) {
     out_image.Import_degrade(conv_image);
   } else {
     out_image.Import_nograde(conv_image);
   }
   
   HealpixProjMap* my_image = new HealpixProjMap(*this,energy,out_image);
   return my_image;
}

double HealpixProjMap::pixelValue(double ilon, double ilat, int k) const {
  // let's say that his expects values in the local (i.e., compact) numbering 
  check_energy_index(k);
  return m_image[k][ilon];
}


bool HealpixProjMap::insideMap(const astro::SkyDir & dir) const {
  // EAC_FIX, for now we have the whole sky, so this is always true
  return true;
}

std::pair<astro::SkyDir, astro::SkyDir> 
HealpixProjMap::minMaxDistPixels(const astro::SkyDir & dir) const {
   astro::SkyDir closest(skyDir(1, 1));
   astro::SkyDir farthest(skyDir(1, 1));
   // EAC_FIX, HEALPIX impl of minMaxDistPixels missing throws std::runtime_error
   throw std::runtime_error("CountsMapHealpix::minMaxDistPixels"
			    "is not implemented");     
   return std::make_pair(closest, farthest);
}


void HealpixProjMap::computeMapIntegrals() {
  std::vector<float>& values = mapIntegrals_access();
  double& totalValue = mapIntegral_access();
  values.clear();
  for (int k(0); k < nenergies(); k++) {
    values.push_back(0);
    for (int i(0); i < nPixels(); i++ ) {
      values.back() += m_solidAngle*m_image[k][i];
    }
  }
  // For a 2D map, the map integral is just the angle-integrated map.
  if (nenergies() == 1) {
    totalValue = values.at(0);
    return;
  }
  // 3D maps: integrate the angle-integrate maps over energy.
  totalValue = 0;
  for (size_t k(0); k < nenergies()-1; k++) {
    totalValue += ( (values[k+1]*energies()[k+1] +
		     values[k]*energies()[k])/2.
		    *std::log(energies()[k+1]/energies()[k]) );
  }
}

ProjMap* HealpixProjMap::rebin(unsigned int factor, bool average) {
  // We don't actually enforce that the nside is a power of two, 
  // so we just set the new nside to be smaller by factor
  if ( factor <= 1 ) return this;
  int newNside = m_healpixProj->healpix().Nside() / factor;
  if ( newNside == 0 ) {
    if ( m_healpixProj->healpix().Nside() == 1 ) return this;
    newNside = 1;
  }
  for (size_t k(0); k < nenergies(); k++) {
    // Make temp storage of the degraded map
    ImagePlane_t newPlane(newNside,m_healpixProj->healpix().Scheme());
    newPlane.Import_degrade(m_image[k],false);
    // Copy that into the image owned by this class
    m_image[k].Set(const_cast<arr<float>&>(newPlane.Map()),newPlane.Scheme());
  }
  return this;
}

int HealpixProjMap::globalToLocal(int glo) const {
  // EAC_FIX, This is no-op until we implement partial sky HEALPix maps
  return glo;
}
   
int HealpixProjMap::localToGlobal(int loc) const {
  // EAC_FIX, This is no-op until we implement partial sky HEALPix maps
  return loc;
}

int HealpixProjMap::nPixels() const {
  // EAC_FIX, Until we implement partial sky HEALPix maps, this just gives the number of 
  // pixels in the full sky
  return m_healpixProj->healpix().Npix();
}


void HealpixProjMap::latchCacheData() {
  // Total number of pixels
  int nPix = m_healpixProj->healpix().Npix();
  // Solid angle in SR = 4pi/npix
  double solidAngle_sr = ASTRO_4PI / float(nPix);
  // Do we want it in deg^2 or sr?
  m_solidAngle = astro::radToDeg( astro::radToDeg( solidAngle_sr ) );
  m_pixelSize = sqrt(m_solidAngle);
}

void HealpixProjMap::check_negative_pixels(const ImagePlane_t & image) const {
   for (size_t i(0); i < image.Npix(); i++) {
     if (image[i] < 0) {
       throw std::runtime_error("Image pixel value less than zero.");
     }
   }
}

}

