/** 
 * @file FitsImage.h
 * @brief Declaration of FitsImage class
 * @authors J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FitsImage.h,v 1.22 2005/02/17 23:22:30 jchiang Exp $
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

#include "fitsio.h"

/** 
 * @class FitsImage
 *
 * @brief A class for reading and storing FITS image data.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FitsImage.h,v 1.22 2005/02/17 23:22:30 jchiang Exp $
 *
 */

class FitsImage : public st_facilities::FitsImage {
    
public:

   FitsImage() : m_eqRot(0) {}
   FitsImage(const std::string &fitsfile);
   FitsImage(const FitsImage &rhs);
   virtual ~FitsImage() {
      delete m_eqRot;
   }

   virtual void getCelestialArrays(std::vector<double> &lonArray,
                                   std::vector<double> &latArray);

   void getPixelBounds(unsigned int naxis,
                       std::vector<double> & pixelBounds) const;

   const std::string & coordSys() const {
      return m_coordSys;
   }

   static void fitsReportError(int status, std::string routine="");

   static int findHdu(const std::string & fitsfile,
                      const std::string & extension);

#ifndef SWIG
   static void readColumn(fitsfile * fptr, const std::string & colname,
                          std::vector<double> & coldata);

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

   std::string m_coordSys;

   EquinoxRotation * m_eqRot;

   bool haveRefCoord();

   static std::string s_routineName;
};

} // namespace Likelihood

#endif // Likelihood_FitsImage.h
