#ifndef Source_h
#define Source_h

#include "Function.h"

#include "Constraint.h"

/** 
 * @class Source
 *
 * @brief Gamma-ray sources.
 *
 * @author J. Chiang, P. Nolan
 *    
 * $Header$
 */

class Source : public Function {


public:
    
   Source();
   ~Source();
#if 0
   //! sets a parameter to be free or fixed in the optimization process
   void setIsFixed(std::string paramName, bool fixed);
   bool getIsFixed(std::string paramName) const;

   //! access the free parameters as a group
   void setFreeParams(std::vector<double>);
   std::vector<double> getFreeParams() const;

   //! derivatives wrt to only free parameters
   virtual std::vector<double> getFreeDerivs() const;

   void setConstraintType(std::string paramName, Constraint::Type type);
   Constraint::Type getConstraintType(std::string paramName);
   void setConstraintMin(std::string paramName, double minValue);
   void setConstraintMax(std::string paramName, double maxValue);
   double getConstraintMin(std::string paramName);
   double getConstraintMax(std::string paramName);

private:
        
   //! flags to determine which parameters are fixed during optimization
   std::vector<bool> m_fixed;

   //! upper and lower bounds for valid parameter ranges
   std::vector<double> m_paramLims;

   //! vectors of pointers to functions that compose the spectral
   //! model for each source
   std::vector<double> *spectralFilter();       // Do we need a grammar for 
   std::vector<double> *spectralComponent();    // combining these things?
#endif
};

#endif // Source_h
