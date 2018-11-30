/**
 * @file CountsMapBase.h
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/CountsMapBase.h,v 1.4 2017/09/29 01:44:15 echarles Exp $
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
  
  class ProjMap;
  class MeanPsf;
  class AppHelpers;

/**
 * @class CountsMapBase
 * 
 */
   
class CountsMapBase : public evtbin::DataProduct {

public:

  //! Typedef for conversion to ProjMap
  typedef enum { 
      //! Convert to an intensity map
      Intensity = 0,
      //! convert to a weights map
      Weights = 1 } ConversionType;

  // These need to be protected for swig
#ifndef SWIG

  static void fillEnergyBinWidths(std::vector<double>& energyBinWidths,
				  const std::vector<double>& energyBinEdges);

  static void fillEnergyBinGeomCenters(std::vector<double>& energyBinCenters,
				       const std::vector<double>& energyBinEdges);

  static void getBMinAndSumMinOverK2(float& bmin, float& sumk2,
				     std::vector<std::vector<float>::const_iterator >& itrs);

  static void getAlphaWts(float& alpha, const float& epsilon2, 
			  std::vector<std::vector<float>::const_iterator >& itrs);

  static void getAlphaVector(std::vector<float>& alphaVect, const float& epsilon2, 
			     const std::vector<const std::vector<float>* >& beffVects);

  static void getWts(std::vector<float>& wts,
		     const float& epsilon2, 
		     const std::vector<float>& alphaVect,
		     const std::vector<float>& beffVect);

  static CountsMapBase* makeAlphaMap(const float& epsilon2, const std::vector<CountsMapBase*>& input_maps);

#endif

  static CountsMapBase* makeAlphaMap(const float& epsilon2, const std::vector<std::string>& input_map_files);

  static CountsMapBase* makeWtsMap(const float& epsilon2, 
				   const CountsMapBase* alphaMap, 
				   const CountsMapBase& beffMap);
  
  static void copyAndUpdateDssKeywords(const std::string& infile,
				       const std::string& outfile,
				       AppHelpers* helper,
				       const std::string& irfs);
  
  static void addBkgEffKeywords(const std::string& outfile,
				const std::string& inputmap,
				const float& efact);

  static void addAlphaMapKeywords(const std::string& outfile,
				  double epsilon,
				  const std::vector<std::string>& inputFiles);

  static void addWtsMapKeywords(const std::string& outfile,
				double epsilon,
				const std::string& bkgmap,
				const std::string& alphamap);

  static void addPhasedExpMapKeywords(const std::string& outfile,
				      const std::string& phased_expmap);

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

  virtual ProjMap* makeProjMap(ConversionType cType = CountsMapBase::Intensity) const {
    // NB, this should be a pure virtual, but that messes up swig
    return 0;
  }

   virtual CountsMapBase* makeBkgEffMap(const MeanPsf & psf, const float& efact) const {
     // NB, this should be a pure virtual, but that messes up swig
     return 0;
   }

   //virtual void binInput(tip::Table::ConstIterator begin, 
   //                      tip::Table::ConstIterator end) = 0;

  virtual void writeOutput(const std::string & creator, 
			   const std::string & out_file) const {
    // NB, this should be a pure virtual, but that messes up swig
  }
  
  virtual void writeAsWeightsMap(const std::string & creator, 
				 const std::string & out_file) const {
    // NB, this should be a pure virtual, but that messes up swig
  }

   virtual void writeEmptyOutput(const std::string & creator, const std::string & out_file) const;

   virtual bool withinBounds(const astro::SkyDir & dir, double energy, long border_size=0) const = 0;
   virtual void setKeywords(tip::Header & header) const = 0;
   virtual void getBoundaryPixelDirs(std::vector<astro::SkyDir> & pixelDirs) const = 0;
   virtual double pixelSize() const = 0;
   virtual double mapRadius() const = 0;


   void writeEnergies(const std::string & creator, 
		      const std::string & out_file,
		      int kmin = 0, int kmax=-1) const;

   void setImage(const std::vector<float> & image);
   void setImage(const std::vector<double> & image);
   
   long imageDimension(int idim) const;

   void getAxisVector(int idim, std::vector<double> & axisVector) const;

   void getEnergies(std::vector<double>& energies) const;

   void getEnergyBinGeomCenters(std::vector<double>& eBinCenters) const;

   void getEnergyBinWidths(std::vector<double>& eBinCenters) const;   

   inline const astro::ProjBase & projection() const {return *m_proj;}

   inline const std::vector<float> & data() const { return m_hist->data(); }

   inline std::vector<float> & data_access() const { return m_hist->data_access(); }

   inline const astro::SkyDir & refDir() const {return m_refDir;}

   inline const std::string & proj_name() const {return m_proj_name;}

   inline bool isGalactic() const {return m_use_lb;}

   inline const std::vector<double> & energies() const { return m_energies; }

   inline const std::string & filename() const { return m_event_file; }
  
   inline double tstart() const { return m_tstart; }

   inline double tstop() const { return m_tstop; }

   const std::vector<Pixel> & pixels() const;

   inline size_t num_energies() const { return m_energies.size(); }

   inline size_t num_ebins() const { return m_energies.size() -1; }
   
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
