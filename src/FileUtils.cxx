/**
 * @file Convolve.cxx
 * @brief Functions to getting data to and from FITS files
 * @author E. Charles, from code in SourceMap by J. Chiang and M. Wood.
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PSFUtils.cxx,v 1.1 2016/09/09 21:21:03 echarles Exp $
 */

#include <memory>

#include "Likelihood/FileUtils.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"

#include "Likelihood/CountsMap.h"
#include "Likelihood/CountsMapHealpix.h"

namespace Likelihood {

  namespace FileUtils {

    int read_fits_image_to_float_vector(const std::string& filename, 
					const std::string& extension ,
					std::vector<float>& vect) {
      std::auto_ptr<const tip::Image>
	image(tip::IFileSvc::instance().readImage(filename,extension));
      vect.clear();
      image->get(vect);     
      return 0;
    }
    
    int read_healpix_table_to_float_vector(const std::string& filename, 
					   const std::string& extension,
					   std::vector<float>& vect) {
      std::auto_ptr<const tip::Table> 
	table(tip::IFileSvc::instance().readTable(filename,extension));
      
      // This is a bit tricky, basically all the data we care about
      // are in the columns called "CHANNELx"
      // Note also that tip work in lowercase
      std::vector<tip::FieldIndex_t> dataColumns;
      const tip::Table::FieldCont& colNames = table->getValidFields();
      for ( tip::Table::FieldCont::const_iterator itr = colNames.begin(); 
	    itr != colNames.end(); itr++ ) {
	if ( itr->find("channel") == 0 ) { 
	  dataColumns.push_back( table->getFieldIndex(*itr) );     
	} else {
	  continue;
	}
      }
      
      int ncol = dataColumns.size();
      tip::Index_t nrow = table->getNumRecords();
      vect.clear();
      vect.reserve(ncol*nrow);
      
      double readVal(0.);
      for ( std::vector<tip::FieldIndex_t>::const_iterator itrData = dataColumns.begin();
	    itrData != dataColumns.end(); itrData++ ) {    
	const tip::IColumn* col = table->getColumn(*itrData);
	for ( tip::Index_t irow(0); irow < nrow; irow++ ) {
	  col->get(irow,readVal);
	  vect.push_back(readVal);
	}
      }      
      return 0;
    }


    int replace_image_from_float_vector(const std::string& filename, 
					const std::string& extension,
					const CountsMapBase& dataMap,
					const std::vector<float>& imageData,
					bool is_src_map) {
      int status(0);
      switch ( dataMap.projection().method() ) {
      case astro::ProjBase::WCS:
	status = FileUtils::replace_image_from_float_vector_wcs(filename,extension,
								static_cast<const CountsMap&>(dataMap),
								imageData,is_src_map);
	return status;
      case astro::ProjBase::HEALPIX:
	status = FileUtils::replace_image_from_float_vector_healpix(filename,extension,
								    static_cast<const CountsMapHealpix&>(dataMap),
								    imageData,is_src_map);
	return status;
      default:
	break;
      }
      std::string errMsg("FileUtils did not recognize projection method used for CountsMap: ");
      errMsg += extension;
      throw std::runtime_error(errMsg);      
      return -1;
    }
    
    int replace_image_from_float_vector_wcs(const std::string& filename, 
					    const std::string& extension,
					    const CountsMap& /* dataMap */,
					    const std::vector<float>& imageData,
					    bool /* is_src_map */) {
      std::auto_ptr<tip::Image> 
	image(tip::IFileSvc::instance().editImage(filename, extension));  
      image->set(imageData);
      return 0;
    }

    int replace_image_from_float_vector_healpix(const std::string& filename, 
						const std::string& extension,
						const CountsMapHealpix& dataMap,
						const std::vector<float>& imageData,
						bool is_src_map) {
      std::auto_ptr<tip::Table> 
	table(tip::IFileSvc::instance().editTable(filename, extension));
      tip::Header& header = table->getHeader();
      dataMap.setKeywords(header);
      long nPix = dataMap.imageDimension(0);
      long nEBins = dataMap.energies().size();
      if ( ! is_src_map ) { 
	nEBins -= 1;
	header["NAXIS2"].set(nEBins);
      }
   
      if ( !dataMap.allSky() ) {
	tip::Header & header(table->getHeader());     
	header["INDXSCHM"].set("EXPLICIT");
	header["REFDIR1"].set(dataMap.isGalactic() ? dataMap.refDir().l() :  dataMap.refDir().ra() );
	header["REFDIR2"].set(dataMap.isGalactic() ? dataMap.refDir().b() :  dataMap.refDir().dec() );
	header["MAPSIZE"].set(dataMap.mapRadius());     
	std::string pixname("PIX");
	table->appendField(pixname, std::string("J"));
	tip::IColumn* col = table->getColumn(table->getFieldIndex(pixname));
	int writeValue(-1);
	for(int iloc(0); iloc < dataMap.nPixels(); iloc++ ) {
	  writeValue = dataMap.localToGlobalIndex(iloc);
	  col->set(iloc,writeValue);
	}
      }
  
      long idx(0);
      double writeValue(0);
      for (long e_index = 0; e_index != nEBins; e_index++ ) {
	std::ostringstream e_channel;
	e_channel<<"CHANNEL"<<e_index+1;
	// Check to see if the column already exist
	tip::FieldIndex_t col_idx = table->getFieldIndex(e_channel.str());
	if ( col_idx < 0 ) {
	  table->appendField(e_channel.str(), std::string("D"));
	  col_idx = table->getFieldIndex(e_channel.str());
	}
	tip::IColumn* col = table->getColumn(col_idx);
	for(long hpx_index = 0; hpx_index != nPix; ++hpx_index, idx++) {
	  writeValue = double(imageData[idx]);
	  col->set(hpx_index,writeValue);
	}
      }
      return 0;
    }
    
    int append_image_from_float_vector(const std::string& filename, 
				       const std::string& extension,
				       const CountsMapBase& dataMap,
				       const std::vector<float>& imageData,
				       bool is_src_map) {
      int status(0);
      switch ( dataMap.projection().method() ) {
      case astro::ProjBase::WCS:
	status = FileUtils::append_image_from_float_vector_wcs(filename,extension,
							       static_cast<const CountsMap&>(dataMap),
							       imageData,is_src_map);
	return status;
      case astro::ProjBase::HEALPIX:
	status = FileUtils::append_image_from_float_vector_healpix(filename,extension,
								   static_cast<const CountsMapHealpix&>(dataMap),
								   imageData,is_src_map);	
	return status;
      default:
	break;
      }
      std::string errMsg("FileUtils did not recognize projection method used for CountsMap: ");
      errMsg += extension;
      throw std::runtime_error(errMsg);
      return -1;
    }
    
    int append_image_from_float_vector_wcs(const std::string& filename, 
					   const std::string& extension,
					   const CountsMap& dataMap,
					   const std::vector<float>& imageData,
					   bool is_src_map) {
      std::vector<long> naxes;
      naxes.push_back(dataMap.imageDimension(0));
      naxes.push_back(dataMap.imageDimension(1));
      long nEBins = is_src_map ? dataMap.energies().size() : dataMap.energies().size() - 1;
      naxes.push_back(nEBins);
      
      tip::IFileSvc::instance().appendImage(filename, extension, naxes);
      int status = FileUtils::replace_image_from_float_vector_wcs(filename, extension,dataMap,imageData,is_src_map);
      return status;
    }

    int append_image_from_float_vector_healpix(const std::string& filename, 
					       const std::string& extension,
					       const CountsMapHealpix& dataMap,
					       const std::vector<float>& imageData,
					       bool is_src_map) {
      tip::IFileSvc::instance().appendTable(filename,extension);
      int status = FileUtils::replace_image_from_float_vector_healpix(filename,extension,
								      dataMap,imageData,is_src_map);
      return status;
    }

  } // namespace FileUtils
 
} // namespace Likelihood
