/** @file DiffuseSource.h
 * @brief DiffuseSource class declaration
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef DiffuseSource_h
#define DiffuseSource_h

#include "Likelihood/Source.h"
#include "Likelihood/Function.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/Event.h"
#include "Likelihood/Arg.h"

namespace Likelihood {

/** 
 * @class DiffuseSource
 *
 * @brief Non-point-like gamma-ray source.  The extragalactic and
 * Galactic diffuse emission components are the obvious examples, but
 * this class also includes discrete diffuse sources such as the LMC
 * or supernova remnants.
 *
 * @author J. Chiang
 *    
 * $Header$ 
 * 
 */

class DiffuseSource : public Source {

public:

   //! A Function describing the spatial distribution of emission is 
   //! required for instantiation.
   DiffuseSource(Function *spatialDist);

   DiffuseSource(const DiffuseSource &rhs);

   virtual ~DiffuseSource() {delete m_spectrum;}

   //! Returns photons/cm^2-s-sr-MeV having been convolved through
   //! the LAT instrument response
   double fluxDensity(const Event &evt) const;

   //! Returns the derivative wrt to the named Parameter
   double fluxDensityDeriv(const Event &evt, std::string &paramName) const;

   //! Predicted number of photons given RoiCuts and ScData
   double Npred();

   //! Derivative of Npred wrt named Parameter
   double NpredDeriv(const std::string &paramName);

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

   //! Computes the exposure integrated over the solid angle subtended by
   //! the Source
   void computeExposure();

   //! spatial model
   Function *m_spatialDist;

   //! spectral model
   Function *m_spectrum;

private:

   //! flag to indicate that static member data has been computed
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
