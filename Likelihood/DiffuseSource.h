/** @file DiffuseSource.h
 * @brief DiffuseSource class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/DiffuseSource.h,v 1.2 2003/03/25 23:22:02 jchiang Exp $
 */

#ifndef DiffuseSource_h
#define DiffuseSource_h

#include "Likelihood/Source.h"
#include "Likelihood/Function.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/Event.h"
#include "Likelihood/SkyDirArg.h"

namespace Likelihood {

/** 
 * @class DiffuseSource
 *
 * @brief Non-point-like gamma-ray source.  The extragalactic and
 * Galactic diffuse emission components are the obvious examples, but
 * this class also includes discrete diffuse sources such as the LMC
 * or supernova remnants.
 *
 * This representation assumes that a single spectral model describes
 * the emission over the entire angular extent of the source.
 * Therefore, in order to have spectral variations across a source, it
 * must comprise sub-components that can be represented using this
 * class.
 *
 * Note that the member function double spatialDist(astro::SkyDir &)
 * returns the spatial distribution of the emission as a function of
 * direction.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/DiffuseSource.h,v 1.2 2003/03/25 23:22:02 jchiang Exp $ 
 *  
 */

class DiffuseSource : public Source {

public:

   //! A Function describing the spatial distribution of emission is 
   //! required for instantiation.
   DiffuseSource(Function *spatialDist);

   DiffuseSource(const DiffuseSource &rhs);

   virtual ~DiffuseSource() {delete m_spatialDist; delete m_spectrum;}

   //! Returns photons/cm^2-s-sr-MeV having been convolved through
   //! the LAT instrument response
   double fluxDensity(const Event &evt) const;

   //! Returns the derivative wrt to the named Parameter
   double fluxDensityDeriv(const Event &evt, std::string &paramName) const;

   //! Predicted number of photons given RoiCuts and ScData
   double Npred();

   //! Derivative of Npred wrt named Parameter
   double NpredDeriv(const std::string &paramName);

   //! Return the spatial distribution of the gamma-ray emission
   double spatialDist(astro::SkyDir &dir) {
      SkyDirArg SDarg(dir);
      return (*m_spatialDist)(SDarg);
   }

   //! Set the spectral model (should also check that the Parameter
   //! names do not conflict with "longitude" and "latitude" of m_dir)
   void setSpectrum(Function *spectrum) {
      m_spectrum = spectrum->clone();
      m_functions["Spectrum"] = m_spectrum;
   }

   virtual Source *clone() const {
      return new DiffuseSource(*this);
   }

protected:

   //! spatial model
   Function *m_spatialDist;

   //! spectral model
   Function *m_spectrum;

private:

   //! flag to indicate that static member data have been computed
   static bool s_haveStaticMembers;

   //! vector of energy values for Npred spectrum quadrature
   static std::vector<double> s_energies;

   //! method to create a logrithmically spaced grid given RoiCuts
   static void makeEnergyVector(int nee = 100);

   //! vector of angle integrated diffuse exposure as as a function of
   //! energy
   std::vector<double> m_exposure;

   //! disable these pure virtual functions inherited from Source
   void setDir(double, double, bool) {}
   void setDir(const astro::SkyDir &, bool) {}

};

} //namespace Likelihood

#endif // DiffuseSource_h
