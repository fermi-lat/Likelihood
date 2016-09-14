/**
 * @file FileUtils.h
 * @brief Functions to deal with PSF Integration and convolution
 * @author E. Charles, from code in SourceMap by J. Chiang and M. Wood.
 *
 *  This file contains a number of functions useful for PSF Integration and convolution.
 *
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FileUtils.h,v 1.1 2016/09/09 21:21:02 echarles Exp $
 */

#ifndef Likelihood_FileUtils_h
#define Likelihood_FileUtils_h

#include <vector>
#include <string>



namespace Likelihood {
  
  class CountsMapBase;
  class CountsMap;
  class CountsMapHealpix;


  namespace FileUtils {
   
    /* Read a FITS image from a file to a vector 

       filename  : The FITS file
       extension : The FITS HDU extension name
       vect      : The vector being filled

       return 0 for success, error code otherwise
     */
    int read_fits_image_to_float_vector(const std::string& filename, 
					const std::string& extension,
					std::vector<float>& vect);

    /* Read a healpix image stored as a FITS table from a file to a vector 

       filename  : The FITS file
       extension : The FITS HDU extension name
       vect      : The vector being filled

       return 0 for success, error code otherwise
     */
    int read_healpix_table_to_float_vector(const std::string& filename, 
					   const std::string& extension,
					   std::vector<float>& vect);

    /* Replace an image in a FITS file 

       filename   : The FITS file
       extension  : The FITS HDU extension name
       dataMap    : The CountsMap used to define the binning 
       imageData  : The data for the image we are writing
       is_src_map : If true, the image has one more energy plane than the CountsMap

       return 0 for success, error code otherwise
    */
    int replace_image_from_float_vector(const std::string& filename, 
					const std::string& extension,
					const CountsMapBase& dataMap,
					const std::vector<float>& imageData,
					bool is_src_map);


    /* Replace an image in a FITS file.  
       This version is for WCS-based images.

       filename   : The FITS file
       extension  : The FITS HDU extension name
       dataMap    : The CountsMap used to define the binning 
       imageData  : The data for the image we are writing
       is_src_map : If true, the image has one more energy plane than the CountsMap

       return 0 for success, error code otherwise
     */   
    int replace_image_from_float_vector_wcs(const std::string& filename, 
					    const std::string& extension,
					    const CountsMap& dataMap,
					    const std::vector<float>& imageData,
					    bool is_src_map);

    /* Replace an image in a FITS file.  
       This version is for HEALPix images.

       filename   : The FITS file
       extension  : The FITS HDU extension name
       dataMap    : The CountsMap used to define the binning 
       imageData  : The data for the image we are writing
       is_src_map : If true, the image has one more energy plane than the CountsMap

       return 0 for success, error code otherwise
     */   
    int replace_image_from_float_vector_healpix(const std::string& filename, 
						const std::string& extension,
						const CountsMapHealpix& dataMap,
						const std::vector<float>& imageData,
						bool is_src_map);
   
    /* Append an image to a FITS file 

       filename   : The FITS file
       extension  : The FITS HDU extension name
       dataMap    : The CountsMap used to define the binning 
       imageData  : The data for the image we are writing
       is_src_map : If true, the image has one more energy plane than the CountsMap

       return 0 for success, error code otherwise
    */
    int append_image_from_float_vector(const std::string& filename, 
				       const std::string& extension,
				       const CountsMapBase& dataMap,
				       const std::vector<float>& imageData,
				       bool is_src_map);
    
    /* Append an image to a FITS file 
       This version is for WCS-based images.

       filename   : The FITS file
       extension  : The FITS HDU extension name
       dataMap    : The CountsMap used to define the binning 
       imageData  : The data for the image we are writing
       is_src_map : If true, the image has one more energy plane than the CountsMap

       return 0 for success, error code otherwise
    */
    int append_image_from_float_vector_wcs(const std::string& filename, 
					   const std::string& extension,
					   const CountsMap& dataMap,
					   const std::vector<float>& imageData,
					   bool is_src_map);

    /* Append an image to a FITS file 
       This version is for HEALPix images.

       filename   : The FITS file
       extension  : The FITS HDU extension name
       dataMap    : The CountsMap used to define the binning 
       imageData  : The data for the image we are writing
       is_src_map : If true, the image has one more energy plane than the CountsMap

       return 0 for success, error code otherwise
    */
    int append_image_from_float_vector_healpix(const std::string& filename, 
					       const std::string& extension,
					       const CountsMapHealpix& dataMap,
					       const std::vector<float>& imageData,
					       bool is_src_map);
   
  } // namespace FileUtils

} // namespace Likelihood

#endif // Likelihood_FileUtils_h
