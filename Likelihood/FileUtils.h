/**
 * @file FileUtils.h
 * @brief Functions to deal with PSF Integration and convolution
 * @author E. Charles, from code in SourceMap by J. Chiang and M. Wood.
 *
 *  This file contains a number of functions useful for PSF Integration and convolution.
 *
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FileUtils.h,v 1.2 2016/09/21 22:42:39 echarles Exp $
 */

#ifndef Likelihood_FileUtils_h
#define Likelihood_FileUtils_h

#include <vector>
#include <string>
#include <map>
#include "Likelihood/SparseVector.h"

namespace tip {
  class Extension;
  class Table;
  class IColumn;
}

namespace optimizers {
  class Parameter;
  class Function;
}


namespace Likelihood {
  
  class CountsMapBase;
  class CountsMap;
  class CountsMapHealpix;
  class Source;
  class SourceModel;

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


    /* Read a sparse healpix image stored as a FITS table from a file to a sparse vector 

       filename  : The FITS file
       extension : The FITS HDU extension name
       vect      : The vector being filled

       return 0 for success, error code otherwise
    */
    int read_healpix_table_to_sparse_vector(const std::string& filename, 
					    const std::string& extension,
					    SparseVector<float>& vect);

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

       return ptr for success, NULL for failure
    */
    tip::Extension* replace_image_from_float_vector(const std::string& filename, 
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

       return ptr for success, NULL for failure
     */   
    tip::Extension* replace_image_from_float_vector_wcs(const std::string& filename, 
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

       return ptr for success, NULL for failure
     */   
    tip::Extension* replace_image_from_float_vector_healpix(const std::string& filename, 
							    const std::string& extension,
							    const CountsMapHealpix& dataMap,
							    const std::vector<float>& imageData,
							    bool is_src_map);
    
    /* Replace an image in a FITS file.  
       This version is for HEALPix images stored as sparse vector

       filename   : The FITS file
       extension  : The FITS HDU extension name
       dataMap    : The CountsMap used to define the binning 
       imageData  : The data for the image we are writing

       return ptr for success, NULL for failure
     */   
    tip::Extension* replace_image_from_sparse_vector_healpix(const std::string& filename, 
							     const std::string& extension,
							     const CountsMapHealpix& dataMap,
							     const SparseVector<float>& imageData,
							     bool is_src_map);
    

    /* Append an image to a FITS file 

       filename   : The FITS file
       extension  : The FITS HDU extension name
       dataMap    : The CountsMap used to define the binning 
       imageData  : The data for the image we are writing
       is_src_map : If true, the image has one more energy plane than the CountsMap

       return ptr for success, NULL for failure
    */
    tip::Extension* append_image_from_float_vector(const std::string& filename, 
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

       return ptr for success, NULL for failure
    */
    tip::Extension* append_image_from_float_vector_wcs(const std::string& filename, 
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

       return ptr for success, NULL for failure
    */
    tip::Extension* append_image_from_float_vector_healpix(const std::string& filename, 
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
    tip::Extension* append_image_from_sparse_vector_healpix(const std::string& filename, 
							    const std::string& extension,
							    const CountsMapHealpix& dataMap,
							    const SparseVector<float>& imageData,
							    bool is_src_map);
   


    /* Append a table to a FITs file.
       This version avoids LOTS of overhead that tip incurs for files 
       with many HDUs by not recomputing checksums and closing the file.

       This means that you need to make sure that you close the 
       file.
       
       This can be done, for example, by deleting one of the pointer return
       by the functions above. 
    */
    void append_table_only(const std::string& file_name,
			   const std::string& table_name);


    /* Write parameters to a tip::Table */
    tip::Extension* write_model_parameters_to_table(const std::string& file_name,
						    const std::string& table_name,
						    const std::vector<const Source*>& sources);
    
    /* Write a source's parameter to a tip::Table*/
    void write_source_parameters_to_table(const Source& source,
					  size_t& irec,
					  tip::IColumn& src_name_col,
					  tip::IColumn& func_name_col,
					  tip::IColumn& par_name_col,
					  tip::IColumn& par_value_col,
					  tip::IColumn& par_scale_col,
					  tip::IColumn& par_error_col,
					  tip::IColumn& par_free_col);
 
    /* Write a functions parameter to a tip::Table*/
    void write_function_parameters_to_table(const std::string& srcName,
					    const optimizers::Function& func,
					    size_t& irec,
					    tip::IColumn& src_name_col,
					    tip::IColumn& func_name_col,
					    tip::IColumn& par_name_col,
					    tip::IColumn& par_value_col,
					    tip::IColumn& par_scale_col,
					    tip::IColumn& par_error_col,
					    tip::IColumn& par_free_col);
   
    /* Write a single parameter to a tip::Table*/
    void write_parameter_to_table(const std::string& srcName,
				  const std::string& funcName,
				  const optimizers::Parameter& param,
				  size_t irec,
				  tip::IColumn& src_name_col,
				  tip::IColumn& func_name_col,
				  tip::IColumn& par_name_col,
				  tip::IColumn& par_value_col,
				  tip::IColumn& par_scale_col,
				  tip::IColumn& par_error_col,
				  tip::IColumn& par_free_col);


    /* Append a column to a FITs table */
    tip::IColumn& append_column(tip::Table& table,
				const std::string& colName,
				const std::string& colFormat);

    /* Read a column from FITs table */
    const tip::IColumn& get_column(const tip::Table& table,
				   const std::string& colName);
        
    /* Write all the parameter from at tip::Table to a SourceModel */
    void read_parameters_from_table(const std::string& file_name,
				    const std::string& table_name,
				    SourceModel& srcModel);
   


  } // namespace FileUtils

} // namespace Likelihood

#endif // Likelihood_FileUtils_h
