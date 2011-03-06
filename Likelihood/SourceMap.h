/**
 * @file SourceMap.h
 * @brief Spatial distribution of a source folded through the
 *        instrument response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/SourceMap.h,v 1.49 2011/03/06 05:31:25 jchiang Exp $
 */

#ifndef Likelihood_SourceMap_h
#define Likelihood_SourceMap_h

#include "st_facilities/libStApiExports.h"

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/Pixel.h"

namespace st_stream {
   class StreamFormatter;
}

namespace Likelihood {

   class CountsMap;
   class DiffuseSource;
   class PointSource;
   class Source;

/*
 * @class SourceMap
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/SourceMap.h,v 1.49 2011/03/06 05:31:25 jchiang Exp $
 */

#ifdef SWIG
class SourceMap {
#else
class  SCIENCETOOLS_API SourceMap {
#endif

public:

   SourceMap(Source * src, const CountsMap * dataMap,
             const Observation & observation,
             bool applyPsfCorrections=false,
             bool performConvolution=true,
             bool resample=true,
             double resamp_factor=2,
             bool verbose=true);

   SourceMap(const std::string & sourceMapsFile, const std::string & srcName);

   ~SourceMap();

   const std::vector<float> & model() const {return m_model;}

   const std::vector<double> & npreds() const {return m_npreds;}

   void addMap(const std::vector<float> & other_model);
   
   static void setBinnedExposure(const std::string & filename);

   double maxPsfRadius(PointSource * src) const;

   const std::string & srcType() const {
      return m_srcType;
   }

   static void setBinnedExpMapName(const std::string & filename);

   static const std::string & binnedExpMap();

   static BinnedExposure & binnedExposure() {
      return *s_binnedExposure;
   }

private:

   std::string m_name;

   /// @brief "Diffuse" or "Point"
   std::string m_srcType;

   const CountsMap * m_dataMap;

   st_stream::StreamFormatter * m_formatter;

   bool m_deleteDataMap;

   /// @brief m_models has the same size as the data in the dataMap plus
   /// one energy plane.
   ///
   /// @todo Keep track of event types included in a given SourceMap.
   std::vector<float> m_model;

   /// @brief Each entry is the angular integral over the energy plane.
   std::vector<double> m_npreds;

/// @bug The binned exposure and mean psf handling is not thread safe.
   static std::string s_expMapFileName;

   static MeanPsf * s_meanPsf;
   static BinnedExposure * s_binnedExposure;
   static unsigned int s_refCount;

   static std::vector<double> s_phi;
   static std::vector<double> s_mu;
   static std::vector<double> s_theta;

   void prepareAngleArrays(int nmu=100, int nphi=50);

   bool haveMapCubeFunction(DiffuseSource * src) const;

   void getMapCorrections(PointSource * src, const MeanPsf & meanPsf,
                          const std::vector<Pixel> & pixels,
                          const std::vector<double> & energies,
                          std::vector<double> & mapCorrections) const;

   void computeExposureAndPsf(const Observation & observation);

   void computeNpredArray();

   double computeResampFactor(const DiffuseSource & src,
                              const CountsMap & dataMap) const;
   
   void makeDiffuseMap(Source * src, 
                       const CountsMap * dataMap,
                       const Observation & observation,
                       bool applyPsfCorrections,
                       bool performConvolution,
                       bool resample,
                       double resamp_factor,
                       bool verbose);

   void makePointSourceMap(Source * src, const CountsMap * dataMap,
                           const Observation & observation,
                           bool applyPsfCorrections, bool performConvolution,
                           bool verbose);

};

} // namespace Likelihood

#endif // Likelihood_SourceMap_h
