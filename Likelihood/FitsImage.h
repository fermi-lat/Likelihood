/** 
 * @file FitsImage.h
 * @brief Declaration of FitsImage class
 * @authors J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FitsImage.h,v 1.17 2004/09/28 04:32:25 jchiang Exp $
 *
 */

#ifndef Likelihood_FitsImage_h
#define Likelihood_FitsImage_h

#include <string>
#include <vector>

#include "astro/SkyDir.h"

#include "st_facilities/FitsImage.h"

#include "Likelihood/Exception.h"

namespace Likelihood {

/** 
 * @class FitsImage
 *
 * @brief A class for reading and storing FITS image data.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FitsImage.h,v 1.17 2004/09/28 04:32:25 jchiang Exp $
 *
 */

class FitsImage : public st_facilities::FitsImage {
    
public:

   FitsImage() {}
   FitsImage(const std::string &fitsfile);
   FitsImage(const FitsImage &rhs);
   virtual ~FitsImage() {
      delete m_eqRot;
   }

   virtual void getCelestialArrays(std::vector<double> &lonArray,
                                   std::vector<double> &latArray);
   
#ifndef SWIG
   FitsImage &operator=(const FitsImage &rhs);

/**
 * @class EquinoxRotation
 * @brief Nested class to perform the "Equinox Rotation" described in
 * <a href="http://lheawww.gsfc.nasa.gov/~jchiang/SSC/like_3.ps">
 * LikeMemo 3</a>.
 */
   class EquinoxRotation {
   public:
      EquinoxRotation() {}
      EquinoxRotation(double alpha0, double delta0);
      ~EquinoxRotation() {}
      void do_rotation(const astro::SkyDir &inDir, astro::SkyDir &outDir);
      EquinoxRotation *clone() const {
         return new EquinoxRotation(*this);
      }
   private:
      std::vector< std::vector<double> > rotMatrix;
   };                       
#endif

private:

   double m_roiRa, m_roiDec;

   EquinoxRotation * m_eqRot;

   bool haveRefCoord();

   void fitsReportError(FILE * stream, int status) const;

};

} // namespace Likelihood

#endif // Likelihood_FitsImage.h
