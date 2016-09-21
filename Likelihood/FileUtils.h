/**
 * @file FileUtils.h
 * @brief Functions to deal with PSF Integration and convolution
 * @author E. Charles, from code in SourceMap by J. Chiang and M. Wood.
 *
 *  This file contains a number of functions useful for PSF Integration and convolution.
 *
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FileUtils.h,v 1.1 2016/09/14 20:10:07 echarles Exp $
 */

#ifndef Likelihood_FileUtils_h
#define Likelihood_FileUtils_h

#include <vector>
#include <string>
#include <map>



namespace Likelihood {
  
  class CountsMapBase;
  class CountsMap;
  class CountsMapHealpix;


  namespace FileUtils {

    
    typedef enum { 
      //! Unkownn
      Unknown = -1,
      //! WCS-based
      WCS = 0,
      //! HEALPix, all-sky, implicit numbering
      HPX_AllSky = 1,
      //! HEALPix, all-sky, sparse mapping
      HPX_Sparse = 2,
      //! HEALPix, partial sky
      HPX_Partial = 3 } SrcMapType;

   
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


    /* Read a sparse healpix image stored as a FITS table from a file to a map 

       filename  : The FITS file
       extension : The FITS HDU extension name
       theMap    : The map being filled

       return 0 for success, error code otherwise
    */
    int read_healpix_table_to_float_map(const std::string& filename, 
					const std::string& extension,
					std::map<size_t,float>& theMap);

    /* Gets the type of HEALPix table containted in a FITS file
       
       filename  : The FITS file
       extension : The FITS HDU extension name
       
       return an enum of SrcMapType
     */
    SrcMapType get_src_map_type(const std::string& filename, 
				const std::string& extension);


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
   
    /* Replace an image in a FITS file.  
       This version is for HEALPix images stored as sparse maps

       filename   : The FITS file
       extension  : The FITS HDU extension name
       dataMap    : The CountsMap used to define the binning 
       imageData  : The data for the image we are writing

       return 0 for success, error code otherwise
     */   
    int replace_image_from_float_map_healpix(const std::string& filename, 
					     const std::string& extension,
					     const CountsMapHealpix& dataMap,
					     const std::map<size_t,float>& imageData,
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

    /* Append an image to a FITS file 
       This version is for HEALPix images.

       filename   : The FITS file
       extension  : The FITS HDU extension name
       dataMap    : The CountsMap used to define the binning 
       imageData  : The data for the image we are writing
       is_src_map : If true, the image has one more energy plane than the CountsMap

       return 0 for success, error code otherwise
    */
    int append_image_from_float_map_healpix(const std::string& filename, 
					    const std::string& extension,
					    const CountsMapHealpix& dataMap,
					    const std::map<size_t,float>& imageData,
					    bool is_src_map);
   
  } // namespace FileUtils

} // namespace Likelihood

#endif // Likelihood_FileUtils_h
