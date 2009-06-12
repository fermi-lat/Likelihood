/** 
 * @file FileFunction.h
 * @brief Declaration for the FileFunction class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FileFunction.h,v 1.1 2005/10/21 20:22:31 jchiang Exp $
 */

#ifndef Likelihood_FileFunction_h
#define Likelihood_FileFunction_h

#include <string>
#include <vector>

#include "optimizers/Arg.h"
#include "optimizers/Function.h"

namespace Likelihood {

/** 
 * @class FileFunction
 *
 * @brief This reads in a tabulated function from an ascii file containing
 * two columns of data.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FileFunction.h,v 1.1 2005/10/21 20:22:31 jchiang Exp $
 */
    
class FileFunction : public optimizers::Function {

public:

   FileFunction(double Normalization=1);
   
   double value(optimizers::Arg & x) const;

   double derivByParam(optimizers::Arg & x, 
                       const std::string & paramName) const;

   virtual Function * clone() const {
      return new FileFunction(*this);
   }

   void readFunction(const std::string & filename);

   const std::string & filename() const {
      return m_filename;
   }

protected:

// This is disabled.
   double integral(optimizers::Arg & xmin, optimizers::Arg & xmax) const;

private:

   double interpolateFlux(double logEnergy) const;

   std::vector<double> m_x;
   std::vector<double> m_y;

   std::string m_filename;

};

} // namespace Likelihood

#endif // Likelihood_FileFunction_h
