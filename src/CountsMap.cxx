/**
 * @file CountsMap.cxx
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/CountsMap.cxx,v 1.44 2010/02/17 19:00:29 jchiang Exp $
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

#include "st_facilities/Env.h"
#include "st_facilities/FitsImage.h"
#include "st_facilities/Util.h"

#include "facilities/commonUtilities.h"

#include "evtbin/Gti.h"
#include "evtbin/LinearBinner.h"
#include "evtbin/LogBinner.h"
#include "evtbin/OrderedBinner.h"

#include "Likelihood/CountsMap.h"
#include "Likelihood/HistND.h"

namespace Likelihood {

CountsMap::CountsMap(const std::string & event_file,
                     const std::string & ev_table,
                     const std::string & sc_file, 
                     const std::string & sc_table,
                     double ref_ra, double ref_dec, const std::string & proj, 
                     unsigned long num_x_pix, unsigned long num_y_pix, 
                     double pix_scale, double axis_rot, bool use_lb,
                     const std::string & ra_field, 
                     const std::string & dec_field, 
                     double emin, double emax, unsigned long nenergies) 
   : DataProduct(event_file, ev_table, evtbin::Gti(event_file)), m_hist(0), 
     m_proj_name(proj), 
     m_crpix(), m_crval(), m_cdelt(), m_axis_rot(axis_rot), 
     m_use_lb(use_lb), m_proj(0) {

   std::vector<evtbin::Binner *> binners;

// The astro::SkyDir::project method will convert ra and dec into 
// bin indices, so we set up the LinearBinners generically by index.
   binners.push_back(new evtbin::LinearBinner(0.5, num_x_pix + 0.5, 1., 
                                              ra_field));
   binners.push_back(new evtbin::LinearBinner(0.5, num_y_pix + 0.5, 1., 
                                              dec_field));

// Energy bins are set up normally.
   binners.push_back(new evtbin::LogBinner(emin, emax, nenergies, 
                                           "photon energy"));

   init(binners, event_file, ev_table, sc_file, sc_table,
        num_x_pix, num_y_pix, ref_ra, ref_dec, pix_scale, pix_scale,
        emin, emax, nenergies, use_lb, proj);

   deleteBinners(binners);
}

CountsMap::CountsMap(const std::string & event_file, 
                     const std::string & ev_table, 
                     const std::string & sc_file, 
                     const std::string & sc_table, 
                     double ref_ra, double ref_dec, const std::string & proj, 
                     unsigned long num_x_pix, unsigned long num_y_pix, 
                     double pix_scale, double axis_rot, bool use_lb,
                     const std::string & ra_field, 
                     const std::string & dec_field, 
                     const std::vector<double> & energies)
   : DataProduct(event_file, ev_table, evtbin::Gti(event_file)), 
     m_hist(0), m_proj_name(proj), 
     m_crpix(), m_crval(), m_cdelt(), m_axis_rot(axis_rot), 
     m_use_lb(use_lb), m_proj(0) {

   std::vector<evtbin::Binner *> binners;

// The astro::SkyDir::project method will convert ra and dec into 
// bin indices, so we set up the LinearBinners generically by index.
   binners.push_back(new evtbin::LinearBinner(0.5, num_x_pix + 0.5, 1., 
                                              ra_field));
   binners.push_back(new evtbin::LinearBinner(0.5, num_y_pix + 0.5, 1., 
                                              dec_field));

// Use custom energy bins.
   std::vector<evtbin::Binner::Interval> energy_intervals;
   for (unsigned int i = 0; i < energies.size()-1; i++) {
      energy_intervals.push_back(evtbin::Binner::Interval(energies[i], 
                                                          energies[i+1]));
   }
   binners.push_back(new evtbin::OrderedBinner(energy_intervals,
                                               "photon energy"));
   init(binners, event_file, ev_table, sc_file, sc_table, 
        num_x_pix, num_y_pix, ref_ra, ref_dec, pix_scale, pix_scale,
        energies.front(), energies.back(), 
        energies.size(), use_lb, proj);

   deleteBinners(binners);
}

CountsMap::CountsMap(const std::string & event_file, 
                     const std::string & ev_table, 
                     const std::string & sc_file, 
                     const std::string & sc_table, 
                     double ref_ra, double ref_dec, const std::string & proj, 
                     unsigned long num_x_pix, unsigned long num_y_pix, 
                     double x_pix_scale, double y_pix_scale,
                     double axis_rot, bool use_lb,
                     const std::string & ra_field, 
                     const std::string & dec_field, 
                     const std::vector<double> & energies)
   : DataProduct(event_file, ev_table, evtbin::Gti(event_file)), 
     m_hist(0), m_proj_name(proj), 
     m_crpix(), m_crval(), m_cdelt(), m_axis_rot(axis_rot), 
     m_use_lb(use_lb), m_proj(0) {

   std::vector<evtbin::Binner *> binners;

// The astro::SkyDir::project method will convert ra and dec into 
// bin indices, so we set up the LinearBinners generically by index.
   binners.push_back(new evtbin::LinearBinner(0.5, num_x_pix + 0.5, 1., 
                                              ra_field));
   binners.push_back(new evtbin::LinearBinner(0.5, num_y_pix + 0.5, 1., 
                                              dec_field));

// Use custom energy bins.
   std::vector<evtbin::Binner::Interval> energy_intervals;
   for (unsigned int i = 0; i < energies.size()-1; i++) {
      energy_intervals.push_back(evtbin::Binner::Interval(energies[i], 
                                                          energies[i+1]));
   }
   binners.push_back(new evtbin::OrderedBinner(energy_intervals,
                                               "photon energy"));
   init(binners, event_file, ev_table, sc_file, sc_table, 
        num_x_pix, num_y_pix, ref_ra, ref_dec, x_pix_scale, y_pix_scale,
        energies.front(), energies.back(), 
        energies.size(), use_lb, proj);

   deleteBinners(binners);
}

CountsMap::CountsMap(const std::string & event_file, 
                     const std::string & ev_table, 
                     const std::string & sc_file, 
                     const std::string & sc_table, 
                     double ref_ra, double ref_dec, const std::string & proj, 
                     unsigned long num_x_pix, unsigned long num_y_pix, 
                     double x_pix_scale, double y_pix_scale,
                     double axis_rot, bool use_lb,
                     const std::string & ra_field, 
                     const std::string & dec_field, 
                     const std::vector<double> & emins,
                     const std::vector<double> & emaxs) 
   : DataProduct(event_file, ev_table, evtbin::Gti(event_file)), 
     m_hist(0), m_proj_name(proj), 
     m_crpix(), m_crval(), m_cdelt(), m_axis_rot(axis_rot), 
     m_use_lb(use_lb), m_proj(0) {

   std::vector<evtbin::Binner *> binners;

// The astro::SkyDir::project method will convert ra and dec into 
// bin indices, so we set up the LinearBinners generically by index.
   binners.push_back(new evtbin::LinearBinner(0.5, num_x_pix + 0.5, 1., 
                                              ra_field));
   binners.push_back(new evtbin::LinearBinner(0.5, num_y_pix + 0.5, 1., 
                                              dec_field));

// Use custom energy bins.
   std::vector<evtbin::Binner::Interval> energy_intervals;
   std::vector<double>::const_iterator emin(emins.begin());
   std::vector<double>::const_iterator emax(emaxs.begin());
   for ( ; emin != emins.end() && emax != emaxs.end(); ++emin, ++emax) {
      energy_intervals.push_back(evtbin::Binner::Interval(*emin, *emax));
   }
   binners.push_back(new evtbin::OrderedBinner(energy_intervals,
                                               "photon energy"));
   init(binners, event_file, ev_table, sc_file, sc_table, 
        num_x_pix, num_y_pix, ref_ra, ref_dec, x_pix_scale, y_pix_scale,
        emins.front(), emaxs.back(), emins.size()+1, use_lb, proj);

   deleteBinners(binners);
}

CountsMap::CountsMap(const std::string & countsMapFile) 
   : DataProduct(countsMapFile, "", evtbin::Gti(countsMapFile)),
     m_use_lb(false) {
   readKeywords(countsMapFile);
   std::vector<evtbin::Binner *> binners;
   binners.push_back(new evtbin::LinearBinner(0.5, m_naxes[0]+0.5, 1., "RA"));
   binners.push_back(new evtbin::LinearBinner(0.5, m_naxes[1]+0.5, 1., "DEC"));
   readEbounds(countsMapFile, binners);
   readImageData(countsMapFile, binners);
   setRefDir();
   setDataDir();
   deleteBinners(binners);
   checkMapConforms();
}

void CountsMap::readImageData(const std::string & countsMapFile,
                              std::vector<evtbin::Binner *> & binners) {
   m_hist = new HistND(binners);
   std::auto_ptr<const tip::Image> 
      image(tip::IFileSvc::instance().readImage(countsMapFile, ""));
   std::vector<float> image_data;
   image->get(image_data);
   m_hist->setData(image_data);
}

void CountsMap::readKeywords(const std::string & countsMapFile) {
   std::auto_ptr<const tip::Image> 
      image(tip::IFileSvc::instance().readImage(countsMapFile, ""));
   const tip::Header & header = image->getHeader();
   typedef std::vector<tip::PixOrd_t> DimCont_t;
   DimCont_t dims = image->getImageDimensions();
   DimCont_t::size_type num_dims = dims.size();
   if (3 != num_dims) {
      throw std::runtime_error("CountsMap::readKeywords: "
                               + std::string("input CountsMap is not 3D."));
   }
   for (unsigned int i = 0; i < 3; i++) {
      m_naxes[i] = dims[i];
   }
   header["CRPIX1"].get(m_crpix[0]);
   header["CRPIX2"].get(m_crpix[1]);
   header["CRPIX3"].get(m_crpix[2]);
   header["CRVAL1"].get(m_crval[0]);
   header["CRVAL2"].get(m_crval[1]);
   header["CRVAL3"].get(m_crval[2]);
   header["CDELT1"].get(m_cdelt[0]);
   header["CDELT2"].get(m_cdelt[1]);
   header["CDELT3"].get(m_cdelt[2]);
   header["CROTA2"].get(m_axis_rot);
   std::string lon_coord;
   std::string lat_coord;
   header["CTYPE1"].get(lon_coord);
   header["CTYPE2"].get(lat_coord);
   std::vector<std::string> tokens;
   facilities::Util::stringTokenize(lon_coord, "-", tokens);
   if (tokens[0] == "GLON") {
      m_use_lb = true;
   }
   m_proj_name = tokens.back();
   m_proj = new astro::SkyProj(m_proj_name, m_crpix, m_crval, m_cdelt, 
                               m_axis_rot, m_use_lb);
}

void CountsMap::readEbounds(const std::string & countsMapFile, 
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

   std::vector<evtbin::Binner::Interval> energy_intervals;
// Convert to MeV
   for (unsigned int i = 0; i < energies.size()-1; i++) {
      energy_intervals.push_back(evtbin::Binner::Interval(energies[i]/1e3, 
                                                          energies[i+1]/1e3));
   }
   binners.push_back(new evtbin::OrderedBinner(energy_intervals,
                                               "photon energy"));
}

void CountsMap::init(std::vector<evtbin::Binner *> & binners, 
                     const std::string & event_file, 
                     const std::string & ev_table,
                     const std::string & sc_file, const std::string & sc_table,
                     unsigned long num_x_pix, unsigned long num_y_pix,
                     double ref_ra, double ref_dec, 
                     double x_pix_scale, double y_pix_scale, 
                     double emin, double emax, 
                     unsigned long nenergies, bool use_lb, 
                     const std::string & proj) {

   m_hist = new HistND(binners);

   m_naxes[0] = num_x_pix;
   m_naxes[1] = num_y_pix;
   m_naxes[2] = nenergies;
   m_crpix[0] = (num_x_pix + 1.) / 2.;
   m_crpix[1] = (num_y_pix + 1.) / 2.;
   m_crpix[2] = 1;
   m_crval[0] = ref_ra;
   m_crval[1] = ref_dec;
   m_crval[2] = emin;
   m_cdelt[0] = -x_pix_scale;
   m_cdelt[1] = y_pix_scale;
// This may not be correct if custom energy bins are used.
   m_cdelt[2] = log(emax/emin)/(nenergies-1.);
   
   m_proj = new astro::SkyProj(proj, m_crpix, m_crval, m_cdelt, 
                               m_axis_rot, use_lb);

   setRefDir();

   harvestKeywords(event_file, ev_table);

   adjustTimeKeywords(sc_file, sc_table);

   setDataDir();

   checkMapConforms();
}

CountsMap::CountsMap(const CountsMap & rhs) : DataProduct(rhs) {
   m_hist = rhs.m_hist->clone();
   m_proj_name = rhs.m_proj_name;
   for (int i = 0; i < 3; i++) {
      m_naxes[i] = rhs.m_naxes[i];
      m_crpix[i] = rhs.m_crpix[i];
      m_crval[i] = rhs.m_crval[i];
      m_cdelt[i] = rhs.m_cdelt[i];
   }
   m_axis_rot = rhs.m_axis_rot;
   m_use_lb = rhs.m_use_lb;
//   m_proj = new astro::SkyProj(*(rhs.m_proj));
   m_proj = new astro::SkyProj(m_proj_name, m_crpix, m_crval, m_cdelt, 
                               m_axis_rot, m_use_lb);
   m_refDir = rhs.m_refDir;
   m_conforms = rhs.m_conforms;
}

CountsMap::~CountsMap() throw() { 
   try {
      delete m_proj;
      delete m_hist;
   } catch (std::exception & eObj) {
      std::cerr << eObj.what() << std::endl;
   } catch (...) {
   }
}

bool CountsMap::withinBounds(const astro::SkyDir & dir, double energy) const {
   std::pair<double, double> coord = dir.project(*m_proj);
   double my_values[] = {coord.first, coord.second, energy};
   std::vector<double> values(my_values, my_values + 3);
   long indx = m_hist->binIndex(values);
   return indx >= 0;
}

void CountsMap::binInput(tip::Table::ConstIterator begin, 
                         tip::Table::ConstIterator end) {

   const evtbin::Hist::BinnerCont_t & binners = m_hist->getBinners();

// From each sky binner, get the name of its field, interpreted as ra
// and dec.
   std::string ra_field = binners[0]->getName();
   std::string dec_field = binners[1]->getName();
   
// Fill histogram, converting each RA/DEC to Sky X/Y on the fly:
   for (tip::Table::ConstIterator itor = begin; itor != end; ++itor) {
      double ra = (*itor)[ra_field].get();
      double dec = (*itor)[dec_field].get();
      
      std::pair<double, double> coord 
         = astro::SkyDir(ra, dec).project(*m_proj);

      double energy = (*itor)["ENERGY"].get();

      double my_values[] = {coord.first, coord.second, energy};
      std::vector<double> values(my_values, my_values + 3);

      m_hist->fillBin(values);
   }
}

void CountsMap::writeOutput(const std::string & creator, 
                            const std::string & out_file) const {

   createFile(creator, out_file, 
              facilities::commonUtilities::joinPath(m_data_dir,
						    "LatCountsMapTemplate"));
   
   std::auto_ptr<tip::Image> 
       output_image(tip::IFileSvc::instance().editImage(out_file, ""));

   typedef std::vector<tip::PixOrd_t> DimCont_t;
   DimCont_t dims = output_image->getImageDimensions();

   DimCont_t::size_type num_dims = dims.size();
   if (3 != num_dims) {
      throw std::runtime_error("CountsMap::writeOutput "
                               "cannot write a count map "
                               "to an image which is not 3D");
   }
   
   const evtbin::Hist::BinnerCont_t & binners = m_hist->getBinners();

   tip::Header & header = output_image->getHeader();
   setKeywords(header);

// Resize image dimensions to conform to the binner dimensions.
   for (DimCont_t::size_type index = 0; index != num_dims; ++index) {
      dims[index] = binners.at(index)->getNumBins();
   }

   output_image->setImageDimensions(dims);

// Copy bins into image.
   output_image->set(m_hist->data());

   writeEbounds(out_file, binners[2]);
   writeGti(out_file);
}

void CountsMap::setImage(const std::vector<float> & image) {
   m_hist->setData(image);
}

long CountsMap::imageDimension(int i) const {
   const evtbin::Hist::BinnerCont_t & binners = m_hist->getBinners();
   if (i < 0 || i > 2) {
      throw std::invalid_argument(std::string("CountsMap::imageDimension:\n")
                                  + "Invalid image dimension value.");
   }
   return binners[i]->getNumBins();
}

void CountsMap::getAxisVector(int i, std::vector<double> & axisVector) const {
   const evtbin::Hist::BinnerCont_t & binners = m_hist->getBinners();
   if (i < 0 || i > 2) {
      throw std::invalid_argument("CountsMap::getAxisVector:\n"
                                  "Invalid image dimension value.");
   }
   axisVector.clear();
   for (long j = 0; j < binners[i]->getNumBins(); j++) {
      axisVector.push_back(binners[i]->getInterval(j).begin());
   }
   long jj = binners[i]->getNumBins() - 1;
   axisVector.push_back(binners[i]->getInterval(jj).end());
}

void CountsMap::setKeywords(tip::Header & header) const {
   const evtbin::Hist::BinnerCont_t & binners = m_hist->getBinners();
   header["CRPIX1"].set(m_crpix[0]);
   header["CRPIX2"].set(m_crpix[1]);
   header["CRPIX3"].set(m_crpix[2]);
   header["CRVAL1"].set(m_crval[0]);
   header["CRVAL2"].set(m_crval[1]);
   header["CRVAL3"].set(m_crval[2]);
   header["CDELT1"].set(m_cdelt[0]);
   header["CDELT2"].set(m_cdelt[1]);
   header["CDELT3"].set(m_cdelt[2]);
   header["CROTA2"].set(m_axis_rot);
   if (m_use_lb) {
      header["CTYPE1"].set("GLON-" + m_proj_name);
      header["CTYPE2"].set("GLAT-" + m_proj_name);
   } else {
      header["CTYPE1"].set(binners[0]->getName() + "---" + m_proj_name);
      header["CTYPE2"].set(binners[1]->getName() + "--" + m_proj_name);
   }
   header["CTYPE3"].set(binners[2]->getName());
}

const std::vector<Pixel> & CountsMap::pixels() const {
   if (m_pixels.empty()) {
      std::vector<astro::SkyDir> pixelDirs;
      std::vector<double> solidAngles;
      getPixels(pixelDirs, solidAngles);
      m_pixels.reserve(pixelDirs.size());
      for (unsigned int i = 0; i < pixelDirs.size(); i++) {
         m_pixels.push_back(Pixel(pixelDirs[i], solidAngles[i]));
      }
   }
   return m_pixels;
}

void CountsMap::
getBoundaryPixelDirs(std::vector<astro::SkyDir> & pixelDirs) const {
   pixelDirs.clear();
   long nx = imageDimension(0);
   long ny = imageDimension(1);

   std::vector<double> longitudes;
   getAxisVector(0, longitudes);
   std::vector<double> latitudes;
   getAxisVector(1, latitudes);

   pixelDirs.reserve(2*nx + 2*ny);
   unsigned int i(0);
   unsigned int j(0);
   try {
      for (j = 0; j < latitudes.size(); j++) {
         astro::SkyDir dir(longitudes.at(i), latitudes.at(j), projection());
         pixelDirs.push_back(dir);
      }
      i = longitudes.size() - 1;
      for (j = 0; j < latitudes.size(); j++) {
         astro::SkyDir dir(longitudes.at(i), latitudes.at(j), projection());
         pixelDirs.push_back(dir);
      }
      j = 0;
      for (i = 0; i < longitudes.size(); i++) {
         astro::SkyDir dir(longitudes.at(i), latitudes.at(j), projection());
         pixelDirs.push_back(dir);
      }
      j = latitudes.size() - 1;
      for (i = 0; i < longitudes.size(); i++) {
         astro::SkyDir dir(longitudes.at(i), latitudes.at(j), projection());
         pixelDirs.push_back(dir);
      }
   } catch (std::exception & eObj) {
      if (st_facilities::Util::expectedException(eObj, 
                                                 "SkyProj wcslib error")) {
         throw std::runtime_error("The counts map appears to contain "
                                  "projection boundaries, e.g., the outer "
                                  "edge of a Hammer-Aitoff map.\n Please "
                                  "analyze a smaller region of the sky that "
                                  "does not include the projection boundary.");
      } else {
         throw;
      }
   }
}

void CountsMap::getPixels(std::vector<astro::SkyDir> & pixelDirs,
                          std::vector<double> & pixelSolidAngles) const {

   long nx = imageDimension(0);
   long ny = imageDimension(1);

   std::vector<double> longitudes;
   getAxisVector(0, longitudes);
   std::vector<double> latitudes;
   getAxisVector(1, latitudes);

   pixelDirs.clear();
   pixelSolidAngles.clear();

   pixelDirs.reserve(nx*ny);
   pixelSolidAngles.reserve(nx*ny);
   std::vector<double>::const_iterator latIt = latitudes.begin();
   for ( ; latIt != latitudes.end() - 1; ++latIt) {
      double latitude = (*latIt + *(latIt+1))/2.;
      std::vector<double>::const_iterator lonIt = longitudes.begin();
      for ( ; lonIt != longitudes.end() - 1; ++lonIt) {
         double longitude = (*lonIt + *(lonIt+1))/2.;
         try {
            astro::SkyDir my_dir(longitude, latitude, projection());
            double solid_angle(computeSolidAngle(lonIt, latIt, projection()));
            pixelDirs.push_back(my_dir);
            pixelSolidAngles.push_back(solid_angle);
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
}

double CountsMap::computeSolidAngle(std::vector<double>::const_iterator lon,
                                    std::vector<double>::const_iterator lat,
                                    const astro::SkyProj & proj) const {

   astro::SkyDir A(*lon, *lat, proj);
   astro::SkyDir B(*lon, *(lat+1), proj);
   astro::SkyDir C(*(lon+1), *(lat+1), proj);
   astro::SkyDir D(*(lon+1), *lat, proj);

   double solidAngle(st_facilities::FitsImage::solidAngle(A, B, C, D));

//    if (solidAngle == 0) {
//       std::cout << "CountsMap::computeSolidAngle: "
//                 << "solid angle = 0 at " << *lon << ", " << *lat << std::endl;
//    }
   return solidAngle;
}

void CountsMap::setRefDir() {
   if (m_use_lb) {
      m_refDir = astro::SkyDir(m_crval[0], m_crval[1], 
                               astro::SkyDir::GALACTIC);
   } else {
      m_refDir = astro::SkyDir(m_crval[0], m_crval[1],
                               astro::SkyDir::EQUATORIAL);
   }
}

void CountsMap::setDataDir() {
   m_data_dir = facilities::commonUtilities::getDataPath("Likelihood");
}

void CountsMap::deleteBinners(std::vector<evtbin::Binner *> & binners) const {
   for (std::vector<evtbin::Binner *>::reverse_iterator it = binners.rbegin();
        it != binners.rend(); ++it) {
      delete *it;
   }
}

void CountsMap::checkMapConforms() {
   m_conforms = false;
   if (static_cast<int>(2.*m_crpix[0]) == static_cast<int>(m_naxes[0] + 1.) &&
       static_cast<int>(2.*m_crpix[1]) == static_cast<int>(m_naxes[1] + 1.) &&
       std::fabs(m_cdelt[0]) == std::fabs(m_cdelt[1])) {
      m_conforms = true;
   }
}

} // namespace Likelihood
