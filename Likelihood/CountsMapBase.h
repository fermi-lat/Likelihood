/**
 * @file CountsMapBase.h
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/Likelihood/CountsMapBase.h,v 1.2 2015/03/03 05:59:55 echarles Exp $
 */

#ifndef Likelihood_CountsMapBase_h
#define Likelihood_CountsMapBase_h

#include <string>

#include "astro/SkyDir.h"

#include "evtbin/DataProduct.h"

#include "Likelihood/HistND.h"
#include "Likelihood/Pixel.h"

namespace astro {
   class SkyProj;
}

namespace evtbin {
   class Binner;
}

namespace tip {
   class Header;
}

namespace Likelihood {

/**
 * @class CountsMapBase
 * 
 */
   
class CountsMapBase : public evtbin::DataProduct {

public:

   CountsMapBase(const std::string & event_file,
		 const std::string & ev_table,
		 const std::string & proj,
		 bool use_lb, double emin, double emax, unsigned long nenergies);

   CountsMapBase(const std::string & event_file,
		 const std::string & ev_table,
		 const std::string & proj,		 
		 bool use_lb, 
		 const std::vector<double> & energies);

   CountsMapBase(const std::string & event_file,
		 const std::string & ev_table,
		 const std::string & proj,
		 bool use_lb, 
		 const std::vector<double> & emins, 
		 const std::vector<double> & emaxs);

   CountsMapBase(const std::string & countsMapFile);

   CountsMapBase(const CountsMapBase & counts_map);

   CountsMapBase(const CountsMapBase & counts_map, 
		 unsigned int idim, unsigned int firstBin, unsigned int lastBin);

   virtual ~CountsMapBase() throw();

   virtual CountsMapBase* clone() const = 0;

   //virtual void binInput(tip::Table::ConstIterator begin, 
   //                      tip::Table::ConstIterator end) = 0;

   //virtual void writeOutput(const std::string & creator, 
   //                         const std::string & out_file) const = 0;

   virtual bool withinBounds(const astro::SkyDir & dir, double energy, long border_size=0) const = 0;
   virtual void setKeywords(tip::Header & header) const = 0;
   virtual void getBoundaryPixelDirs(std::vector<astro::SkyDir> & pixelDirs) const = 0;
   virtual double pixelSize() const = 0;

   void setImage(const std::vector<float> & image);
   void setImage(const std::vector<double> & image);
   
   long imageDimension(int idim) const;

   void getAxisVector(int idim, std::vector<double> & axisVector) const;

   void getEnergies(std::vector<double>& energies) const;

   inline const astro::ProjBase & projection() const {return *m_proj;}

   inline const std::vector<float> & data() const { return m_hist->data(); }

   inline const astro::SkyDir & refDir() const {return m_refDir;}

   inline const std::string & proj_name() const {return m_proj_name;}

   inline bool isGalactic() const {return m_use_lb;}

   inline const std::vector<double> & energies() const { return m_energies; }

   inline const std::string & filename() const { return m_event_file; }
  
   inline double tstart() const { return m_tstart; }

   inline double tstop() const { return m_tstop; }

   const std::vector<Pixel> & pixels() const;

protected:

   HistND * m_hist;
   std::string m_proj_name;
   
   bool m_use_lb;
   astro::ProjBase* m_proj;

   std::vector<double> m_energies;

   virtual double computeSolidAngle(std::vector<double>::const_iterator lon,
				    std::vector<double>::const_iterator lat,
				    const astro::ProjBase & proj) const = 0;

   virtual void getPixels(std::vector<astro::SkyDir> & pixelDirs,
			  std::vector<double> & pixelSolidAngles) const = 0;

   virtual void readKeywords(const std::string & countsMapFile) = 0;
 
   void readEbounds(const std::string & countsMapfile,
                    std::vector<evtbin::Binner *> & binners);

   virtual void readImageData(const std::string & countsMapfile,
			      std::vector<evtbin::Binner *> & binners) = 0;

   void setRefDir(double val1, double val2);  

   void setDataDir();

   void deleteBinners(std::vector<evtbin::Binner *> & binners) const;

   void readTimeKeywords(const std::string& event_file);
  

private:

   astro::SkyDir m_refDir;

   mutable std::vector<Pixel> m_pixels;

   double m_tstart;
   double m_tstop;

   CountsMapBase & operator=(const CountsMapBase &) {return *this;}


};

}

#endif // Likelihood_CountsMapBase_h
