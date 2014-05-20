/**
 * @file SourceMap.h
 * @brief Spatial distribution of a source folded through the
 *        instrument response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceMap.h,v 1.53 2013/09/18 06:33:59 jchiang Exp $
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
   class MeanPsf;
   class PointSource;
   class Source;

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

   SourceMap(Source * src, const CountsMap * dataMap,
             const Observation & observation,
             bool applyPsfCorrections=false,
             bool performConvolution=true,
             bool resample=true,
             double resamp_factor=2,
             double minbinsz=0.1,
             bool verbose=true);

   SourceMap(const std::string & sourceMapsFile,
             const std::string & srcName,
             const Observation & observation);

   ~SourceMap();

   const std::vector<float> & model() const {
      return m_model;
   }

   const std::vector<double> & npreds() const {
      return m_npreds;
   }

//   void addMap(const std::vector<float> & other_model);
   
   double maxPsfRadius(PointSource * src) const;

   const std::string & srcType() const {
      return m_srcType;
   }

   const std::string & name() const {
      return m_name;
   }

private:

   std::string m_name;

   /// @brief "Diffuse" or "Point"
   std::string m_srcType;

   const CountsMap * m_dataMap;
   
   const Observation & m_observation;

   st_stream::StreamFormatter * m_formatter;

   bool m_deleteDataMap;

   /// @brief m_models has the same size as the data in the dataMap plus
   /// one energy plane.
   ///
   /// @todo Keep track of event types included in a given SourceMap.
   std::vector<float> m_model;

   /// @brief Each entry is the angular integral over the energy plane.
   std::vector<double> m_npreds;

   bool haveMapCubeFunction(DiffuseSource * src) const;

   void getMapCorrections(PointSource * src, const MeanPsf & meanPsf,
                          const std::vector<Pixel> & pixels,
                          const std::vector<double> & energies,
                          std::vector<double> & mapCorrections) const;


   void computeNpredArray();

   double computeResampFactor(const DiffuseSource & src,
                              const CountsMap & dataMap) const;
   
   void makeDiffuseMap(Source * src, 
                       const CountsMap * dataMap,
                       bool applyPsfCorrections,
                       bool performConvolution,
                       bool resample,
                       double resamp_factor,
                       double minbinsiz,
                       bool verbose);

   void makePointSourceMap(Source * src,
                           const CountsMap * dataMap,
                           bool applyPsfCorrections,
                           bool performConvolution,
                           bool verbose);

   void applyPhasedExposureMap();

   double psfValueEstimate(const MeanPsf & meanPsf, double energy,
                           double offset, double pixelSolidAngle) const;

};

} // namespace Likelihood

#endif // Likelihood_SourceMap_h
