/** 
 * @file DMFitFunction2.h
 * @brief Declaration for the DMFit Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/DMFitFunction2.h,v 1.1 2009/03/04 11:24:27 cohen Exp $
 */

#ifndef Likelihood_DMFitFunction2_h
#define Likelihood_DMFitFunction2_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

#include <string>

namespace Likelihood {

/** 
 * @class DMFitFunction2
 *
 * @brief A DMFit function for modeling Dark Matter spectra.
 *
 * @author J. Cohen-Tanugi, based on the DMFit package by S. Profumo and T. Jeltema
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/DMFitFunction2.h,v 1.1 2009/03/04 11:24:27 cohen Exp $
 */
    
class DMFitFunction2 : public optimizers::Function {

public:

  DMFitFunction2() : m_filename(""){
      init(1., 100., 1.0, 4, 1, 100.);
   }

   /// @param norm Normalization of the function
   /// @param mass The mass of the Dark Matter particle
   /// @param bratio The branching ratio between the 2 allowed final states
   /// @param channel0 : index of the first final state
   /// @param channel1 : index of the second final state   
   DMFitFunction2(double norm, double mass, double bratio, 
                 int channel0, int channel1, double lowerLimit) : 
    m_filename(""),m_ch0(channel0),m_ch1(channel1),m_emin(lowerLimit*0.001)
    {
      init(norm, mass, bratio, channel0, channel1, lowerLimit);
    }

   double value(optimizers::Arg&) const;

   double derivByParam(optimizers::Arg &x, const std::string &paramName) const;

   virtual Function *clone() const {
      return new DMFitFunction2(*this);
   }

   void readFunction(const std::string & filename);

   const std::string & filename() const {
      return m_filename;
   }

protected:

   double integral(optimizers::Arg &, optimizers::Arg &) const {
      return 0;
   }

private:
//    void updateCache(const double x, const double m, 
//                     const double b, const int ch0,
//                     const int ch1,  const double emin) const;
    
   void updateCache(double x, double m, 
                    double b, int ch0,
                    int ch1,  double emin) const;
    
   void init(double norm, double mass, double bratio, 
             int channel0, int channel1, double lowerLimit);

   std::string m_filename;

   int m_ch0, m_ch1;
   double m_emin;
   mutable double m_1mb;
   mutable double m_dmMass;
   mutable double m_int_dnde0;
   mutable double m_int_dnde1;
   mutable double m_int_dndm0;
   mutable double m_int_dndm1;
};

} // namespace Likelihood

#endif // Likelihood_DMFitFunction2_h
