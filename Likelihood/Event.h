/** 
 * @file Event.h
 * @brief Event class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Event.h,v 1.24 2004/06/06 22:43:40 jchiang Exp $
 */

#ifndef Likelihood_Event_h
#define Likelihood_Event_h

#include <vector>
#include <string>
#include <map>

#include "astro/SkyDir.h"
#include "Likelihood/FitsImage.h"
#include "Likelihood/Exception.h"

namespace Likelihood {

class DiffuseSource;

/** 
 * @class Event
 *
 * @brief A gamma-ray event --- apparent direction, energy, arrival time, etc.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Event.h,v 1.24 2004/06/06 22:43:40 jchiang Exp $
 */

class Event {
    
public:

   Event() {}
   Event(double ra, double dec, double energy, double time, 
         double sc_ra, double sc_dec, double muZenith, int type=2);

// The compiler-supplied copy constructor should suffice.
//   Event(const Event &);

   ~Event() {}

   astro::SkyDir getDir() const {return m_appDir;}
   astro::SkyDir getScDir() const {return m_scDir;}

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
                          std::string diffuseComponent) const
      throw(Exception);
    
   /// This method takes the spatial distribution of the emission for
   /// the DiffuseSource src and computes the event-specific response.
   /// See section 1 of 
   /// <a href="http://lheawww.gsfc.nasa.gov/~jchiang/SSC/like_3.ps>
   /// LikeMemo 3</a>.  The computed response is added to the
   /// m_respDiffuseSrcs map with the specified name.  sr_radius is the
   /// "source region" radius (in degrees) over which the spatial
   /// distribution of src will be integrated.
   void computeResponse(DiffuseSource &src, 
                        double sr_radius = 30.) {
      std::vector<DiffuseSource *> srcs;
      srcs.push_back(&src);
      computeResponse(srcs, sr_radius);
   }

   /// Compute the reponse integrals for a vector of DiffuseSources
   void computeResponse(std::vector<DiffuseSource *> &srcs, 
                        double sr_radius = 30.);

   /// Write the diffuse responses for each source to a file.
   void writeDiffuseResponses(const std::string & filename);

   /// Set diffuse response for infinite energy resolution.
   void setDiffuseResponse(const std::string & srcName, double value) {
      m_respDiffuseSrcs[srcName].clear();
      m_respDiffuseSrcs[srcName].push_back(value);
   }

   static void toLower(std::string & name);
   
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

   /// Vector of true energies.
   double m_estep;
   std::vector<double> m_trueEnergies;

   /// Response function data, unique to each event, and comprising an
   /// energy redistribution function for each diffuse source.
   typedef std::vector<double> diffuse_response;
   std::map<std::string, diffuse_response> m_respDiffuseSrcs;

   /// Compute Celestial direction from (phi, mu) in Equinox-centered
   /// coordinates.
   void getCelestialDir(double phi, double mu, FitsImage::EquinoxRotation
                        &eqRot, astro::SkyDir &dir);

   /// Angular arrays over the source region for the diffuse integrals.
   static std::vector<double> s_mu;
   static std::vector<double> s_phi;
   static FitsImage::EquinoxRotation s_eqRot;
   static bool s_haveSourceRegionData;

   /// @todo Find a rational way of defining the source region radius.
   /// For now, use the default values passed from
   /// computeResponse(...)
   void prepareSrData(double sr_region, int nmu=100, int nphi=50);

   void getNewDiffuseSrcs(const std::vector<DiffuseSource *> & srcList,
                          std::vector<DiffuseSource *> & srcs) const;

};

} // namespace Likelihood
#endif // Likelihood_Event_h
