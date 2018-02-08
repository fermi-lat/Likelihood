/**
 * @file Convolve.cxx
 * @brief Functions to getting data to and from FITS files
 * @author E. Charles, from code in SourceMap by J. Chiang and M. Wood.
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/FileUtils.cxx,v 1.9 2017/06/26 22:02:29 echarles Exp $
 */

#include <memory>
#include <cstdio>

#include "Likelihood/FileUtils.h"

#include "st_facilities/Timer.h"
#include "st_facilities/Environment.h"
#include "facilities/commonUtilities.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/IColumn.h"
#include "tip/Table.h"
#include "tip/Image.h"
#include "fitsio.h"
#include "tip/FileSummary.h"
#include "tip/Image.h"
#include "tip/TipException.h"

#include "optimizers/Function.h"
#include "optimizers/Parameter.h"

#include "Likelihood/CountsMap.h"
#include "Likelihood/CountsMapHealpix.h"
#include "Likelihood/HistND.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceModel.h"

namespace Likelihood {

  namespace FileUtils {

    void makeTemplateFilepath(const std::string& basename, 
			      std::string& filepath){
      std::string data_dir = st_facilities::Environment::dataPath("Likelihood");
      filepath = facilities::commonUtilities::joinPath(data_dir,basename);      
    }

    bool fileHasExtension(const std::string& filename, 
			  const std::string& extension) {
      try {
	std::auto_ptr<const tip::Extension> ext(tip::IFileSvc::instance().readExtension(filename,extension));
      } catch (tip::TipException &) {
	return false;
      }
      return true;
    }


    void read_ebounds_to_vector(const std::string& filename,
				std::vector<double>& energies) {
      
      std::auto_ptr<const tip::Table> ebounds(tip::IFileSvc::instance().readTable(filename, "EBOUNDS"));
      tip::Table::ConstIterator it = ebounds->begin();
      tip::Table::ConstRecord & row = *it;
      energies.clear();
      energies.resize(ebounds->getNumRecords() + 1, 0.);
      double emax;
      for (int i = 0 ; it != ebounds->end(); ++it, i++) {
	row["E_MIN"].get(energies[i]);
	row["E_MAX"].get(emax);
	energies[i] /= 1000.;
      }
      energies.back() = emax / 1000.;
    }

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

    int read_healpix_table_to_sparse_vector(const std::string& filename, 
					    const std::string& extension,
					    size_t npix, 
					    SparseVector<float>& vect) {
      vect.clear_data();
      std::auto_ptr<const tip::Table> 
	table(tip::IFileSvc::instance().readTable(filename,extension));

      // There are two columns: key and value
      tip::FieldIndex_t pix_field = table->getFieldIndex("PIX");
      tip::FieldIndex_t chan_field = table->getFieldIndex("CHANNEL");
      tip::FieldIndex_t val_field = table->getFieldIndex("VALUE");
      const tip::IColumn* pix_col = table->getColumn(pix_field);
      const tip::IColumn* chan_col = table->getColumn(chan_field);
      const tip::IColumn* val_col = table->getColumn(val_field);

      

      tip::Index_t nrow = table->getNumRecords();

      std::vector<size_t> pix_vect(nrow, 0);
      std::vector<short> chan_vect(nrow, 0);
      std::vector<float> val_vect(nrow, 0.);
      std::vector<size_t> key_vect(nrow, 0.);

      for ( tip::Index_t irow(0); irow < nrow; irow++ ) {
	try {
	  pix_col->get(irow,pix_vect[irow]);
	  chan_col->get(irow,chan_vect[irow]);
	  val_col->get(irow,val_vect[irow]);
	} catch (...) {
	  std::cout << "Read failed at " << filename << ' ' << extension << ' ' << irow << std::endl;
	  return -1;
	}
      
	key_vect[irow] = chan_vect[irow] * npix + pix_vect[irow];
      }

      vect.fill_from_key_and_value(key_vect,val_vect);

      return 0;
    }


    SrcMapType get_src_map_type(const std::string& filename, 
				const std::string& extension) {

      try {
	std::auto_ptr<const tip::Image>
	image(tip::IFileSvc::instance().readImage(filename,extension));
	// It's a WCS Image
	return FileUtils::WCS;
      } catch (tip::TipException &) {
	;
      }

      std::auto_ptr<const tip::Table> 
	table(tip::IFileSvc::instance().readTable(filename,extension));

      const tip::Header& header = table->getHeader();     
      std::string indxschm;
      std::string pixtype; 

      try {
	header["PIXTYPE"].get(pixtype);
      } catch (tip::TipException &) {
	// NOT a HEALPIX map
	return FileUtils::Unknown;
      }	 

      if ( pixtype != "HEALPIX" ) {
	// NOT a HEALPIX map
	return FileUtils::Unknown;
      }

      try {
	header["INDXSCHM"].get(indxschm);
      } catch (tip::TipException &) {
	// NOT a HEALPIX map
	return FileUtils::Unknown;
      }	 

      if ( indxschm == "IMPLICIT" ) {
	// All-sky HEALPIX map	
	return FileUtils::HPX_AllSky;
      } else if ( indxschm == "SPARSE" ) {
	return FileUtils::HPX_Sparse;
      }	else if ( indxschm == "EXPLICIT" ) {
	return FileUtils::HPX_Partial;
     }     
      return FileUtils::Unknown;
    }

    tip::Extension* replace_energies(const std::string& filename, 
				     const std::string& extension,
				     const std::vector<double>& energies) {
      
      tip::Table* table = tip::IFileSvc::instance().editTable(filename, extension);
      
      // Resize table: number of records in output file must == the number of bins in the binner.
      table->setNumRecords(energies.size());
    
      // Need output table iterator.
      tip::FieldIndex_t col_idx = table->getFieldIndex("Energy");
      tip::Header& header = table->getHeader();
      tip::IColumn* col = table->getColumn(col_idx);

      for (size_t index(0); index < energies.size(); ++index) {
	float writeValue = energies[index];
	col->set(index,writeValue);
      }
      
      return table;
    }

  
    
    void replace_image_from_hist_wcs(tip::Image& image,
				     HistND& hist) {
      typedef std::vector<tip::PixOrd_t> DimCont_t;
      DimCont_t dims = image.getImageDimensions();

      DimCont_t::size_type num_dims = dims.size();
      if (3 != num_dims) {
	throw std::runtime_error("CountsMap::writeOutput "
				 "cannot write a count map "
				 "to an image which is not 3D");
      }
   
      const evtbin::Hist::BinnerCont_t & binners = hist.getBinners();
 
      // Resize image dimensions to conform to the binner dimensions.
      for (DimCont_t::size_type index = 0; index != num_dims; ++index) {
	dims[index] = binners.at(index)->getNumBins();
      }

      image.setImageDimensions(dims);

      // Copy bins into image.
      image.set(hist.data());
    }

    
    void replace_image_from_hist_hpx(tip::Table& table,
				     const std::string& colString,
				     HistND& hist,
				     evtbin::HealpixBinner& hpx_binner){
      
      const evtbin::Hist::BinnerCont_t & binners = hist.getBinners();
      long nPix = hpx_binner.getNumBins();
      table.setNumRecords(nPix);
      
      long nEBins = 1;
      if ( binners.size() > 1 ) {
	nEBins = binners[1]->getNumBins();
      } 
      
      // If the map is less than all-sky we need to write the pixel indices
      if ( ! hpx_binner.allSky() ) {
	std::string pixname("PIX");
	table.appendField(pixname, std::string("J"));
	tip::IColumn* col = table.getColumn(table.getFieldIndex(pixname));
	int writeValue(-1);
	for(long hpx_index = 0; hpx_index != hpx_binner.getNumBins(); ++hpx_index) {
	  writeValue = hpx_binner.pixelIndices()[hpx_index];
	  col->set(hpx_index,writeValue);
	}
      }
      
      std::vector<float> histData(nPix);
      double writeValue(0.);
      
      std::vector<unsigned int> ivalues(binners.size(),0);
      
      // std::cout << "Writing map: " << std::flush;
      for (long e_index = 0; e_index != nEBins; e_index++ ) {
	//std::cout << '.' << std::flush;
	std::ostringstream e_channel;
	e_channel << colString << e_index+1;
	//create new column
	table.appendField(e_channel.str(), std::string("D"));
	// get the column
	tip::IColumn* col = table.getColumn(table.getFieldIndex(e_channel.str()));
	// get the data slice from the underlying histogram
	if ( binners.size() > 1 ) {
	  hist.getSlice(0,ivalues,histData);
	} else {
	  histData = hist.data();
	}
	for(long hpx_index = 0; hpx_index != nPix; ++hpx_index) {
	  writeValue = double(histData[hpx_index]);
	  col->set(hpx_index,writeValue);
	}
	if ( ivalues.size() > 1 ) ivalues[1]++;
      }
      //std::cout << '!' << std::endl;    
    }

    tip::Extension* replace_image_from_float_vector(const std::string& filename, 
						    const std::string& extension,
						    const CountsMapBase& dataMap,
						    const std::vector<float>& imageData,
						    bool is_src_map, 
						    int kmin, int kmax) {
      tip::Extension* ptr(0);
      switch ( dataMap.projection().method() ) {
      case astro::ProjBase::WCS:
	ptr = FileUtils::replace_image_from_float_vector_wcs(filename,extension,
							     static_cast<const CountsMap&>(dataMap),
							     imageData,is_src_map,kmin,kmax);
	return ptr;
      case astro::ProjBase::HEALPIX:
	ptr = FileUtils::replace_image_from_float_vector_healpix(filename,extension,
								 static_cast<const CountsMapHealpix&>(dataMap),
								 imageData,is_src_map,kmin,kmax);
	return ptr;
      default:
	break;
      }
      std::string errMsg("FileUtils did not recognize projection method used for CountsMap: ");
      errMsg += extension;
      throw std::runtime_error(errMsg);      
      return ptr;
    }
    
    tip::Extension* replace_image_from_float_vector_wcs(const std::string& filename, 
							const std::string& extension,
							const CountsMap& dataMap,
							const std::vector<float>& imageData,
							bool is_src_map, 
							int kmin, int kmax) {
      tip::Image* image = tip::IFileSvc::instance().editImage(filename, extension);  
      long nEBins_data = is_src_map ? dataMap.num_energies() : dataMap.num_ebins();
      kmax = kmax < 0 ? nEBins_data : kmax;
      long nEBins = kmax - kmin;
      
      if ( nEBins == nEBins_data ) {
	// Include the whole map
	image->set(imageData);
      } else {
	// Only include some energy bins
	tip::ImageBase::PixelCoordinate imageDims(image->getImageDimensions());
	imageDims[2] = nEBins;
	image->setImageDimensions(imageDims);
	size_t npix = dataMap.pixels().size();
	std::vector<float> partial(nEBins*npix);
	std::copy(imageData.begin()+(kmin*npix),imageData.begin()+(kmax*npix),partial.begin());
      }
      return image;
    }

    tip::Extension* replace_image_from_float_vector_healpix(const std::string& filename, 
							    const std::string& extension,
							    const CountsMapHealpix& dataMap,
							    const std::vector<float>& imageData,
							    bool is_src_map, 
							    int kmin, int kmax) {
      tip::Table* table =tip::IFileSvc::instance().editTable(filename, extension);
      tip::Header& header = table->getHeader();
      dataMap.setKeywords(header);
      long nPix = dataMap.imageDimension(0);
      long nEBins_data = is_src_map ? dataMap.num_energies() : dataMap.num_ebins();
      kmax = kmax < 0 ? nEBins_data : kmax;
      long nEBins = kmax - kmin;
      
   
      if ( !dataMap.allSky() ) {
	tip::Header & header(table->getHeader());     
	header["INDXSCHM"].set("EXPLICIT");
	header["REFDIR1"].set(dataMap.isGalactic() ? dataMap.refDir().l() :  dataMap.refDir().ra() );
	header["REFDIR2"].set(dataMap.isGalactic() ? dataMap.refDir().b() :  dataMap.refDir().dec() );
	header["MAPSIZE"].set(dataMap.mapRadius());     
	std::string pixname("PIX");
	tip::FieldIndex_t col_idx(-1);
	try {
	  col_idx = table->getFieldIndex(pixname);
	} catch (tip::TipException &) {
	  table->appendField(pixname, std::string("J"));
	  col_idx = table->getFieldIndex(pixname);
	}	  
	tip::IColumn* col = table->getColumn(col_idx);
	int writeValue(-1);
	for(int iloc(0); iloc < dataMap.nPixels(); iloc++ ) {
	  writeValue = dataMap.localToGlobalIndex(iloc);
	  col->set(iloc,writeValue);
	}
      }
  
      long idx(0);
      double writeValue(0);
      for (long e_index = kmin; e_index != kmax; e_index++ ) {
	std::ostringstream e_channel;
	e_channel<<"CHANNEL"<< (e_index-kmin)+1;
	// Check to see if the column already exist
	tip::FieldIndex_t col_idx(-1);
	try {
	  col_idx = table->getFieldIndex(e_channel.str());
	} catch (tip::TipException &) {
	  table->appendField(e_channel.str(), std::string("D"));
	  col_idx = table->getFieldIndex(e_channel.str());
	}	
	tip::IColumn* col = table->getColumn(col_idx);
	for(long hpx_index = 0; hpx_index != nPix; ++hpx_index, idx++) {
	  writeValue = double(imageData[idx]);
	  col->set(hpx_index,writeValue);
	}
      }
      return table;
    }


    tip::Extension* replace_image_from_sparse_vector_healpix(const std::string& filename, 
							     const std::string& extension,
							     const CountsMapHealpix& dataMap,
							     const SparseVector<float>& imageData,
							     bool is_src_map, 
							     int kmin, int kmax) {
      
      tip::Table* table = tip::IFileSvc::instance().editTable(filename, extension);
      tip::Header& header = table->getHeader();
      dataMap.setKeywords(header);
      long nPix = dataMap.imageDimension(0);      
      long nEBins_data = is_src_map ? dataMap.num_energies() : dataMap.num_ebins();
      kmax = kmax < 0 ? nEBins_data : kmax;
      long nEBins = kmax - kmin;
   
      header["INDXSCHM"].set("SPARSE");
      header["REFDIR1"].set(dataMap.isGalactic() ? dataMap.refDir().l() :  dataMap.refDir().ra() );
      header["REFDIR2"].set(dataMap.isGalactic() ? dataMap.refDir().b() :  dataMap.refDir().dec() );
 
      int idx(0);
      std::vector<size_t> key_vect;
      std::vector<float> val_vect;

      
      if ( nEBins_data == nEBins ) {
	imageData.fill_key_and_value(key_vect,val_vect);
      } else {
	imageData.fill_key_and_value(key_vect,val_vect,kmin*nPix,kmax*nPix);
      }

      size_t nfilled = key_vect.size();

      tip::FieldIndex_t pix_idx(-1);
      try {
	pix_idx = table->getFieldIndex("PIX");
      } catch (tip::TipException &) {
	table->appendField("PIX", "J");
	pix_idx = table->getFieldIndex("PIX");
      }	   
      tip::IColumn* pix_col = table->getColumn(pix_idx);

      tip::FieldIndex_t chan_idx(-1);
      try {
	chan_idx = table->getFieldIndex("CHANNEL");
      } catch (tip::TipException &) {
	table->appendField("CHANNEL", "I");
	chan_idx = table->getFieldIndex("CHANNEL");
      }	   
      tip::IColumn* chan_col = table->getColumn(chan_idx);
  
      tip::FieldIndex_t val_idx(-1);
      try {
	val_idx = table->getFieldIndex("VALUE");
      } catch (tip::TipException &) {
	table->appendField("VALUE", "E");
	val_idx = table->getFieldIndex("VALUE");
      }	   
      tip::IColumn* val_col = table->getColumn(val_idx);
       
      std::vector<size_t> keys;
      std::vector<float> values;
          
      imageData.fill_key_and_value(keys, values);
      for ( size_t ikey(0); ikey < keys.size(); ikey++ ) {
	size_t key = keys[ikey];
	size_t ipix = key % nPix;
	short ie = key / nPix;
	pix_col->set(ikey, ipix);
	chan_col->set(ikey, ie);
	val_col->set(ikey, values[ikey]);
      }
      return table;
    }

    tip::Extension* append_image_from_float_vector(const std::string& filename, 
						   const std::string& extension,
						   const CountsMapBase& dataMap,
						   const std::vector<float>& imageData,
						   bool is_src_map, 
						   int kmin, int kmax) {
      tip::Extension* ptr(0);
      switch ( dataMap.projection().method() ) {
      case astro::ProjBase::WCS:
	ptr = FileUtils::append_image_from_float_vector_wcs(filename,extension,
							    static_cast<const CountsMap&>(dataMap),
							    imageData,is_src_map,kmin,kmax);
	return ptr;
      case astro::ProjBase::HEALPIX:
	ptr = FileUtils::append_image_from_float_vector_healpix(filename,extension,
								static_cast<const CountsMapHealpix&>(dataMap),
								imageData,is_src_map,kmin,kmax);	
	return ptr;
      default:
	break;
      }
      std::string errMsg("FileUtils did not recognize projection method used for CountsMap: ");
      errMsg += extension;
      throw std::runtime_error(errMsg);
      return ptr;
    }
    
    tip::Extension* append_image_from_float_vector_wcs(const std::string& filename, 
						       const std::string& extension,
						       const CountsMap& dataMap,
						       const std::vector<float>& imageData,
						       bool is_src_map, 
						       int kmin, int kmax) {
      std::vector<long> naxes;
      naxes.push_back(dataMap.imageDimension(0));
      naxes.push_back(dataMap.imageDimension(1));
      long nEBins = is_src_map ? dataMap.energies().size() : dataMap.energies().size() - 1;
      naxes.push_back(nEBins);
      
      tip::IFileSvc::instance().appendImage(filename, extension, naxes);
      return FileUtils::replace_image_from_float_vector_wcs(filename, extension,dataMap,imageData,is_src_map,kmin,kmax);
    }

    tip::Extension* append_image_from_float_vector_healpix(const std::string& filename, 
							   const std::string& extension,
							   const CountsMapHealpix& dataMap,
							   const std::vector<float>& imageData,
							   bool is_src_map, 
							   int kmin, int kmax) {
      tip::IFileSvc::instance().appendTable(filename,extension);
      return FileUtils::replace_image_from_float_vector_healpix(filename,extension,
								dataMap,imageData,is_src_map,kmin,kmax);
    }

    tip::Extension* append_image_from_sparse_vector_healpix(const std::string& filename, 
							    const std::string& extension,
							    const CountsMapHealpix& dataMap,
							    const SparseVector<float>& imageData,
							    bool is_src_map, 
							    int kmin, int kmax) {
      /* Add the table by hand to avoid lots of overhead from tip */
      // This causes memory corruption, so lets' not use it
      // append_table_only(filename,extension);
      tip::IFileSvc::instance().appendTable(filename,extension);
      return FileUtils::replace_image_from_sparse_vector_healpix(filename,extension,
								 dataMap,imageData,is_src_map,kmin,kmax);
    }


    void append_table_only(const std::string& file_name,
			   const std::string& table_name) {
      fitsfile * fp = 0;
      int status = 0;
      
      // Open the file create the file.
      fits_open_file(&fp, const_cast<char *>(file_name.c_str()), READWRITE, &status);
      if (0 != status) {
	throw tip::TipException(status, "File does not exist \"" + file_name + "\"");
      }    
      fits_create_tbl(fp, BINARY_TBL, 0, 0, 0, 0, 0, const_cast<char *>(table_name.c_str()), &status);
      if (0 != status) {
	throw tip::TipException(status, "Unable to create table named \"" + table_name + "\" in file \"" + file_name + "\"");
      }
    }

    /* Write parameters to a tip::Table */
    tip::Extension* write_model_parameters_to_table(const std::string& filename,
						    const std::string& extension,
						    const std::vector<const Source*>& sources) {

      tip::IFileSvc::instance().appendTable(filename,extension);
      tip::Table* table = tip::IFileSvc::instance().editTable(filename, extension);
      
      tip::IColumn& src_name_col = append_column(*table,"SourceName","32A");
      tip::IColumn& func_name_col = append_column(*table,"FunctionName","32A");
      tip::IColumn& par_name_col = append_column(*table,"ParamName","32A");
      tip::IColumn& par_value_col = append_column(*table,"ParamValue","D");
      tip::IColumn& par_scale_col = append_column(*table,"ParamScale","D");
      tip::IColumn& par_error_col = append_column(*table,"ParamError","D");
      tip::IColumn& par_free_col = append_column(*table,"ParamFree","L");

      size_t irec(0);
      for ( std::vector<const Source*>::const_iterator itr = sources.begin(); itr != sources.end(); itr++ ) {
	const Source* src = *itr;
	write_source_parameters_to_table(*src,
					 irec,src_name_col,
					 func_name_col,par_name_col,par_value_col,
					 par_scale_col,par_error_col,par_free_col);
      }
      return table;
    }

    void write_source_parameters_to_table(const Source& source,
					  size_t& irec,
					  tip::IColumn& src_name_col,
					  tip::IColumn& func_name_col,
					  tip::IColumn& par_name_col,
					  tip::IColumn& par_value_col,
					  tip::IColumn& par_scale_col,
					  tip::IColumn& par_error_col,
					  tip::IColumn& par_free_col ) {
      const std::string srcName = source.getName();
      const Source::FuncMap& funcMap =  source.getSrcFuncs();      
      for ( Source::FuncMap::const_iterator itr = funcMap.begin(); itr != funcMap.end(); itr++ ) {
	const optimizers::Function& func = *(itr->second);
	write_function_parameters_to_table(srcName,func,
					   irec,src_name_col,
					   func_name_col,par_name_col,par_value_col,
					   par_scale_col,par_error_col,par_free_col);
      }
    }

    void write_function_parameters_to_table(const std::string& srcName,
					    const optimizers::Function& func,
					    size_t& irec,
					    tip::IColumn& src_name_col,
					    tip::IColumn& func_name_col,
					    tip::IColumn& par_name_col,
					    tip::IColumn& par_value_col,
					    tip::IColumn& par_scale_col,
					    tip::IColumn& par_error_col,
					    tip::IColumn& par_free_col) {
      const std::string& funcName = func.getName();
      std::vector<std::string> parNames; 
      func.getParamNames(parNames);      
      for ( std::vector<std::string>::const_iterator itr = parNames.begin(); itr != parNames.end(); itr++ ) {
	const optimizers::Parameter param = func.getParam(*itr);
	write_parameter_to_table(srcName,funcName,param,
				 irec,src_name_col,
				 func_name_col,par_name_col,par_value_col,
				 par_scale_col,par_error_col,par_free_col);
	irec += 1;
      }
    }
 

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
				  tip::IColumn& par_free_col) {
      src_name_col.set(irec,srcName.c_str());
      func_name_col.set(irec,funcName.c_str());
      par_name_col.set(irec,param.getName().c_str());
      par_value_col.set(irec,param.getValue());
      par_scale_col.set(irec,param.getScale());
      par_error_col.set(irec,param.error());
      par_free_col.set(irec,param.isFree());
    }

    tip::IColumn& append_column(tip::Table& table,
				const std::string& colName,
				const std::string& colFormat) {
      table.appendField(colName,colFormat);
      tip::FieldIndex_t field = table.getFieldIndex(colName);
      tip::IColumn* col = table.getColumn(field);
      return *col;
    }

    /* Append a column to a FITs table */
    const tip::IColumn& get_column(const tip::Table& table,
				   const std::string& colName) {
      tip::FieldIndex_t field = table.getFieldIndex(colName);
      const tip::IColumn* col = table.getColumn(field);
      return *col;
    }
   
      /* Write a single parameter to a tip::Table*/
    void read_parameters_from_table(const std::string& file_name,
				    const std::string& table_name,
				    SourceModel& srcModel) {
      std::auto_ptr<const tip::Table> 
	table(tip::IFileSvc::instance().readTable(file_name,table_name));
      
      const tip::IColumn& src_name_col = get_column(*table,"SourceName");
      const tip::IColumn& func_name_col = get_column(*table,"FunctionName");
      const tip::IColumn& par_name_col = get_column(*table,"ParamName");
      const tip::IColumn& par_value_col = get_column(*table,"ParamValue");
      const tip::IColumn& par_scale_col = get_column(*table,"ParamScale");
      const tip::IColumn& par_error_col = get_column(*table,"ParamError");
      const tip::IColumn& par_free_col = get_column(*table,"ParamFree");

      tip::Index_t nrow = table->getNumRecords();

      std::string src_name;
      std::string func_name;
      std::string par_name;
      double par_value(0.);
      double par_scale(0.);
      double par_error(0.);
      bool par_free(0.);

      for ( tip::Index_t irow(0); irow < nrow; irow++ ) {
	src_name_col.get(irow,src_name);
	func_name_col.get(irow,func_name);
	par_name_col.get(irow,par_name);
	par_value_col.get(irow,par_value);
	par_scale_col.get(irow,par_scale);
	par_error_col.get(irow,par_error);
	par_free_col.get(irow,par_free);
	optimizers::Parameter param = srcModel.getParam(src_name,func_name,par_name);
	param.setValue(par_value);
	param.setScale(par_scale);
	param.setError(par_error);
	param.setFree(par_free);
	srcModel.setParam(param,func_name,src_name);
      }      
    }
  
  } // namespace FileUtils
 
} // namespace Likelihood
