/** @file FitsImage.h
 * @brief Declaration of FitsImage class
 * @authors J. Chiang
 *
 * $Header$
 *
 */

#ifndef FitsImage_h
#define FitsImage_h

#include <vector>
#include <string>
#include <valarray>
#include "astro/SkyDir.h"

namespace Likelihood {

/** 
 * @class FitsImage
 *
 * @brief A class for reading and storing FITS image data.
 *
 * @author J. Chiang
 *    
 * $Header$
 *
 */

class FitsImage {
    
public:

   FitsImage() {}
   FitsImage(std::string &fitsfile);
   virtual ~FitsImage() {
      if (m_haveRefCoord) delete m_eqRot;
   }

   //! A vector of the image axes dimensions
   void fetchAxisDims(std::vector<int> &axisDims);

   //! The names (CTYPEs) of the image axes
   void fetchAxisNames(std::vector<std::string> &axisNames);

   //! Fetch a vector filled with axis abscissa points for the naxis-th
   //! coordinate (We are perhaps dangerously ignoring header-specified
   //! projection effects, i.e., we assume a "plate-carree" projection, 
   //! as for EGRET maps, for all images).
   void fetchAxisVector(unsigned int naxis, std::vector<double> &axisVector);

   //! This method computes arrays of longitude and latitude obtained
   //! by traversing the image plane by column number then row.
   //! If m_refCoord == true, then these will be the coordinates in
   //! the unrotated coordinate system (see section 4 of 
   //! <a href="http://lheawww.gsfc.nasa.gov/~jchiang/SSC/like_3.ps">LikeMemo 3</a>.
   void fetchCelestialArrays(std::valarray<double> &lonArray,
                             std::valarray<double> &latArray);

   //! Fetch the pixel values.  They will be indexed by column, row,
   //! then plane, i.e., indx = i + j*NAXIS1 + k*NAXIS1*NAXIS2.  Note
   //! that each image plane is indexed starting at the lower left
   //! (South-East) corner.
   void fetchImageData(std::valarray<double> &imageData);

   //! This returns the pixel solid angles.  Use of this method assumes
   //! that m_axis[0] represents a longitudinal coordinate and that
   //! m_axis[1] represents a latitudinal coordinate.  The pixel values
   //! will be indexed by column then row, indx = i + j*NAXIS1.
   void fetchSolidAngles(std::valarray<double> &solidAngles);

/**
 * @class EquinoxRotation
 * @brief Nested class to perform the "Equinox Rotation" described in
 * <a href="http://lheawww.gsfc.nasa.gov/~jchiang/SSC/like_3.ps">LikeMemo 3</a>.
 */
   class EquinoxRotation {
   public:
      EquinoxRotation(double alpha0, double delta0);
      ~EquinoxRotation() {}
      void do_rotation(const astro::SkyDir &inDir, astro::SkyDir &outDir);
   private:
      std::vector< std::vector<double> > rotMatrix;
   };                       

protected:

/** 
 * @class AxisParams
 * @brief Nested n-tuple class to represent FITS image axis information
 */
   class AxisParams {
   public:
      AxisParams() {}
      ~AxisParams() {}
      int size;
      double refVal;
      double step;
      int refPixel;
      std::string axisType;
      std::string comment;
      bool logScale;

      //! Returns a vector of abscissa values based on the axis parameters.
      void computeAxisVector(std::vector<double> &axisVector);
   };

   //! Interface to cfitsio routines
   void read_fits_image(std::string &filename, std::vector<AxisParams> &axes,
                        std::valarray<double> &image);

   //! FITS file name
   std::string m_filename;

   //! Descriptions for each image axis
   std::vector<AxisParams> m_axes;

   //! Vectors of abscissa values for each axis in local (i.e., rotated)
   //! coordinate system
   std::vector< std::vector<double> > m_axisVectors;

   //! The FITS image data
   std::valarray<double> m_image;

   //! Hi-jack the LONPOLE and LATPOLE FITS keywords for use as the 
   //! true longitude and latitude of the reference pixel...will fix
   //! this later once we have LAT-specific FITS keywords.
   bool m_haveRefCoord;
   double m_lonpole, m_latpole;

   EquinoxRotation *m_eqRot;

};

} // namespace Likelihood

#endif // FitsImage.h
