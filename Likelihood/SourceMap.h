/**
 * @file SourceMap.h
 * @brief Spatial distribution of a source folded through the
 *        instrument response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/SourceMap.h,v 1.61 2016/04/01 01:35:07 echarles Exp $
 */

#ifndef Likelihood_SourceMap_h
#define Likelihood_SourceMap_h

#include "st_facilities/libStApiExports.h"

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/Pixel.h"

#include "healpix_map.h"

namespace astro {
   class HealpixProj;
}

namespace st_stream {
   class StreamFormatter;
}

namespace Likelihood {
  
   // EAC, switch to using CountsMapBase and projection specific sub-classes
   class CountsMapBase;
   class CountsMap;
   class CountsMapHealpix;
   class DiffuseSource;
   class MeanPsf;
   class PointSource;
   class Source;
   class WcsMap2;

/*
 * @class SourceMap
 *
 */

#ifdef SWIG
class SourceMap {
#else
class  SCIENCETOOLS_API SourceMap {
#endif

public:

   static void fillHealpixFromWcsMap(const WcsMap2& inputMap,
				     const astro::HealpixProj& proj,
				     Healpix_Map<float>& hpm);

   static void fillHealpixFromWcsMap(const WcsMap2& inputMap,
                                     const astro::HealpixProj& proj,
				     const std::vector<int>& pixList,
                                     Healpix_Map<float>& hpm);

   SourceMap(Source * src, const CountsMapBase * dataMap,
             const Observation & observation,
             bool applyPsfCorrections=false,
             bool performConvolution=true,
             bool resample=true,
             double resamp_factor=2,
             double minbinsz=0.1,
             bool verbose=true,
	     const SourceMap* weights = 0);

   SourceMap(const ProjMap& weight_map,
	     const CountsMapBase * dataMap,
             const Observation & observation,
	     bool verbose=true);

   SourceMap(const std::string & sourceMapsFile,
             const std::string & srcName,
             const Observation & observation,
	     const SourceMap* weights = 0,
	     bool isWeights = false);

   ~SourceMap();

   const std::vector<float> & model() const {
      return m_model;
   }

   const std::vector<double> & npreds() const {
      return m_npreds;
   }

   const std::vector<std::pair<double,double> > & npred_weights() const {
      return m_npred_weights;
   }

//   void addMap(const std::vector<float> & other_model);
   
   void setImage(const std::vector<float>& model);

   double maxPsfRadius(PointSource * src) const;

   const std::string & srcType() const {
      return m_srcType;
   }

   const std::string & name() const {
      return m_name;
   }

   const SourceMap* weights() const { return m_weights; }

   void setWeights(const SourceMap* weights);

protected:

   void readImage(const std::string& sourceMapFile);

   void readTable_healpix(const std::string& sourceMapFile);

private:

   std::string m_name;

   /// @brief "Diffuse" or "Point"
   std::string m_srcType;

   const CountsMapBase * m_dataMap;
   
   const Observation & m_observation;

   st_stream::StreamFormatter * m_formatter;

   const SourceMap* m_weights;

   bool m_deleteDataMap;

   /// @brief m_models has the same size as the data in the dataMap plus
   /// one energy plane.
   ///
   /// @todo Keep track of event types included in a given SourceMap.
   std::vector<float> m_model;

   /// @brief Each entry is the angular integral over the energy plane.
   std::vector<double> m_npreds;

   /// @brief Weights for npreds.
   std::vector< std::pair<double,double> > m_npred_weights;

   /// @brief Scaling factor between the true and projected angular separation
   std::vector< std::vector< double > > m_pixelOffset;
   std::vector< double > m_pixelCoordX;
   std::vector< double > m_pixelCoordY;   

   std::string m_psfEstimatorMethod;

   double m_psfEstimatorFtol;
   double m_psfEstimatorPeakTh;

   bool haveMapCubeFunction(DiffuseSource * src) const;

   void getMapCorrections(PointSource * src, const MeanPsf & meanPsf,
                          const std::vector<Pixel> & pixels,
			  const std::pair<double,double>& srcCoord,
			  const std::vector< std::pair<double,double> > & pixCoords,
                          const std::vector<double> & energies,
                          std::vector<double> & mapCorrections) const;


   void computeNpredArray(bool isWeights=false);

   double computeResampFactor(const DiffuseSource & src,
                              const CountsMapBase & dataMap) const;
   
   void makeDiffuseMap(Source * src, 
                       const CountsMapBase * dataMap,
                       bool applyPsfCorrections,
                       bool performConvolution,
                       bool resample,
                       double resamp_factor,
                       double minbinsiz,
                       bool verbose);

   void makeDiffuseMap_wcs(Source * src,
			   const CountsMap * dataMap,
			   bool applyPsfCorrections,
			   bool performConvolution,
			   bool resample,
			   double resamp_factor,
			   double minbinsiz,
			   bool verbose);

   void makeDiffuseMap_healpix(Source * src,
			       const CountsMapHealpix * dataMap,
			       bool applyPsfCorrections,
			       bool performConvolution,
			       bool resample,
			       double resamp_factor,
			       double minbinsiz,
			       bool verbose);

   void makeDiffuseMap_native(Source * src,
                              const CountsMapHealpix * dataMap,
			      bool applyPsfCorrections,
			      bool performConvolution,
			      bool resample,
			      double resamp_factor,
			      double minbinsiz,
			      bool verbose);

   void makePointSourceMap(Source * src,
                           const CountsMapBase * dataMap,
                           bool applyPsfCorrections,
                           bool performConvolution,
                           bool verbose);

   void makePointSourceMap_wcs(Source * src,
			       const CountsMap * dataMap,
			       bool applyPsfCorrections,
			       bool performConvolution,
			       bool verbose);
   
   void makePointSourceMap_healpix(Source * src,
				   const CountsMapHealpix * dataMap,
				   bool applyPsfCorrections,
				   bool performConvolution,
				   bool verbose);

   void makeProjectedMap(const ProjMap& weight_map, bool& extrapolated);

   void applyPhasedExposureMap();

   void createOffsetMap(Source * src,
			const CountsMap * dataMap);

   double psfValueEstimate(const MeanPsf & meanPsf, double energy, 
                           const astro::SkyDir & srcDir, 
			   const Pixel & pixel,
			   const std::pair<double, double>& srcCoord,
			   const std::pair<double, double>& pixCoord) const;

   double integrate_psf_adaptive(const MeanPsf & meanPsf, double energy, 
				 const astro::SkyDir & srcDir, 
				 const Pixel & pixel,
				 const std::pair<double, double>& srcCoord,
				 const std::pair<double, double>& pixCoord) const;

   double integrate_psf(const MeanPsf & meanPsf, double energy,
                        const astro::SkyDir & srcDir, 
			const Pixel & pixel, size_t npts = 11) const;

   double integrate_psf_fast(const MeanPsf & meanPsf, double energy, 
			     const astro::SkyDir & srcDir, 
			     const Pixel & pixel, 
			     const std::pair<double, double>& srcCoord,
			     const std::pair<double, double>& pixCoord,
			     size_t npts) const;
};

} // namespace Likelihood

#endif // Likelihood_SourceMap_h
