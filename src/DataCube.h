#ifndef SciTools_DataCube_h
#define SciTools_DataCube_h

#include <vector>
#include <string>

/** 
 * @class DataCube
 *
 * @brief Base class for Science Tools DataCubes.  These are 
 * N-dimensional arrays containing such data as stacks of images 
 * (e.g., read from FITS files) or more general representations of the 
 * instrument response.
 *
 * @author J. Chiang
 *    
 * $Header$
 */

class DataCube {
    
public:
    
   DataCube();
   ~DataCube();

   double value(const double *location); /* interpolated value at 
					    desired location */

private:

   // dimension of data cube
   long m_naxes;

   // sizes of each dimension
   std::vector<long> m_naxis;

   // reference pixel
   std::vector<long> m_refpix;

   // step size along each axis at reference pixel
   std::vector<double> m_refstep;

   // coordinate values at reference pixel center
   std::vector<double> m_refval;

   // the data cube itself
   std::vector<double> m_data;

};

#endif // SciTools_DataCube_h
