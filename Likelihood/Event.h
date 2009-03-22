/** 
 * @file Event.h
 * @brief Event class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Event.h,v 1.43 2009/01/19 15:18:17 sfegan Exp $
 */

#ifndef Likelihood_Event_h
#define Likelihood_Event_h

#include <vector>
#include <string>
#include <map>

#include "astro/SkyDir.h"

#include "Likelihood/Exception.h"

namespace Likelihood {

   class DiffuseSource;
   class EquinoxRotation;
   class ResponseFunctions;
   class Source;

/** 
 * @class Event
 *
 * @brief A gamma-ray event --- apparent direction, energy, arrival time, etc.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Event.h,v 1.43 2009/01/19 15:18:17 sfegan Exp $
 */

class Event {
    
public:

   Event();

   Event(double ra, double dec, double energy, double time, 
         const astro::SkyDir & scZAxis, const astro::SkyDir & scXAxis, 
         double muZenith, bool useEdisp, const std::string & respName,
         int type=2);

   ~Event() {}

   const astro::SkyDir & getDir() const {return m_appDir;}
   const astro::SkyDir & getScDir() const {return m_scDir;}

   const astro::SkyDir & zAxis() const {return m_scDir;}
   const astro::SkyDir & xAxis() const {return m_scXDir;}

   /// some useful accessor functions
   double getEnergy() const {return m_energy;}
   double getArrTime() const {return m_arrTime;}
   double getMuZenith() const {return m_muZenith;}
   int getType() const {return m_type;}

   /// separation in units of radians
   double getSeparation(const astro::SkyDir &dir) const 
      {return m_appDir.SkyDir::difference(dir);}

   /// return the Event specific diffuse response function 
   /// for the named diffuse component
   double diffuseResponse(double energy, 
                          const std::string& srcName) const;
    
   void computeResponseGQ(std::vector<DiffuseSource *> & srcs, 
                          const ResponseFunctions & respFuncs,
                          bool useDummyValue=false);

   /// This method takes the spatial distribution of the emission for
   /// the DiffuseSource src and computes the event-specific response.
   /// See section 1 of 
   /// <a href="http://lheawww.gsfc.nasa.gov/~jchiang/SSC/like_3.ps>
   /// LikeMemo 3</a>.  The computed response is added to the
   /// m_respDiffuseSrcs map with the specified name.  sr_radius is the
   /// "source region" radius (in degrees) over which the spatial
   /// distribution of src will be integrated.
   void computeResponse(DiffuseSource &src, 
                        const ResponseFunctions & respFuncs, 
                        double sr_radius=30., double sr_radius2=80.) {
      std::vector<DiffuseSource *> srcs;
      srcs.push_back(&src);
      computeResponse(srcs, respFuncs, sr_radius, sr_radius2);
   }

   /// Compute the response integrals for a vector of DiffuseSources
   void computeResponse(std::vector<DiffuseSource *> &srcs, 
                        const ResponseFunctions & respFuncs, 
                        double sr_radius=30., double sr_radius2=80);

   /// Find the DiffuseSources that need to have diffuse response
   /// values computed.
   void getNewDiffuseSrcs(const std::vector<DiffuseSource *> & srcList,
                          std::vector<DiffuseSource *> & srcs) const;

   /// Write the diffuse responses for each source to a file.
   void writeDiffuseResponses(const std::string & filename);

   /// Set diffuse response for infinite energy resolution.
//    void setDiffuseResponse(const std::string& srcName, double value) {
//      const std::string & diffuseComponent = diffuseSrcName(srcName);
//       m_respDiffuseSrcs[diffuseComponent].clear();
//       m_respDiffuseSrcs[diffuseComponent].push_back(value);
//    }
   void setDiffuseResponse(const std::string & componentName, double value) {
      m_respDiffuseSrcs[componentName].clear();
      m_respDiffuseSrcs[componentName].push_back(value);
   }

   /// Set diffuse response for finite energy resolution.
   void setDiffuseResponse(const std::string& srcName, 
                           const std::vector<double> & gaussianParams);

   static void toLower(std::string & name);

   /// Direct access to vector of true energies.
   const std::vector<double> & trueEnergies() const {
      return m_trueEnergies;
   }

   /// Direct access to diffuse responses.
   const std::vector<double> & diffuseResponse(const std::string& srcName) const;

   void computeGaussianParams(const std::string & srcName, double & norm, 
                              double & mean, double & sigma) const;
   
   /// @return The srcName with the response function name + "::" prepended.
   ///         NB: Everything is converted to lower-case since that is
   ///         what tip returns when you ask for FITS table column names.
   /// @param srcName The source name.
   const std::string& diffuseSrcName(const std::string & srcName) const;

   /// Add or subtract contribution from a given source to m_modelSum.
   //typedef Source::CachedResponse CachedResponse; -- cannot do this!
   typedef std::pair<bool, double> CachedResponse;
   void updateModelSum(const Source & src, CachedResponse* cResp = 0);

   void resetModelSum();

   double modelSum() const {
      return m_modelSum;
   }

   void deleteSource(const std::string & srcName);

   void set_ctbclasslevel(int ctbclasslevel) {
      m_ctbclasslevel = ctbclasslevel;
   }

   int ctbclasslevel() const {
      return m_ctbclasslevel;
   }

private:

   /// apparent direction, energy, arrival time, and cosine(zenith angle)
   astro::SkyDir m_appDir;
   double m_energy;
   double m_arrTime;
   double m_muZenith;

   /// Event type (front vs back for now)
   int m_type;

   /// spacecraft info at event arrival time
   astro::SkyDir m_scDir;
   astro::SkyDir m_scXDir;

   bool m_useEdisp;
   std::string m_respName;

   double m_modelSum;

   std::map<std::string, double> m_fluxDensities;
   
   /// Vector of true energies.
   double m_estep;
   std::vector<double> m_trueEnergies;

   /// Response function data, unique to each event, and comprising an
   /// energy redistribution function for each diffuse source.
   typedef std::vector<double> diffuse_response;
   std::map<std::string, diffuse_response> m_respDiffuseSrcs;
   mutable std::map<std::string, std::string> m_diffSrcNames;

   /// Use this variable to keep track of Classification tree info.
   int m_ctbclasslevel;

   /// Compute Celestial direction from (phi, mu) in Equinox-centered
   /// coordinates.
   void getCelestialDir(double phi, double mu, EquinoxRotation & eqRot,
                        astro::SkyDir &dir);

   /// Angular arrays over the source region for the diffuse integrals.
   static std::vector<double> s_mu;
   static std::vector<double> s_phi;
   static std::vector<double> s_mu_2;
   static bool s_haveSourceRegionData;

   /// @todo Find a more rational way of defining the source region
   /// radius.  For now, use the default values passed from
   /// computeResponse(...).
   /// @param sr_radius Default to use for most events
   /// @param sr_radius2 Use this in cases where the apparent
   /// inclination lies outside the region of non-zero effective area
   /// (i.e., for theta > 78 degrees).
   void prepareSrData(double sr_radius, double sr_radius2);

   void fillMuArray(double sr_radius, int nmu, std::vector<double> & mu) const;

};

} // namespace Likelihood
#endif // Likelihood_Event_h
