/** 
 * @file Table.cxx
 * @brief Implementation of Table member functions
 * @authors T. Burnett, J. Chiang using code from Y. Ikebe
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Table.cxx,v 1.8 2003/05/29 21:52:50 jchiang Exp $
 */

#include "Likelihood/Table.h"
#include "LikelihoodException.h"
#include "fitsio.h"

#include <cstdio>
#include <iostream>
#include <sstream>

namespace Likelihood {

void Table::add_column(std::string colname)
{
    m_pars.push_back(Column(colname));
}

void Table::add_columns(std::string colnames)
{
    int offset = 0, length=colnames.size();
    do {
        int next = colnames.substr(offset).find_first_of(" ");
        if( next>0 ) next += offset; else next = length;
        add_column(colnames.substr(offset,next-offset));
        if( next==length) break;
        offset = next+1;
    }while(1);

}

void Table::add_columns(std::vector<std::string> columnNames) {
   for (unsigned int i = 0; i < columnNames.size(); i++) {
      add_column(columnNames[i]);
   }
}

void Table::read_FITS_table(std::string filename, int hdu) 
   throw(LikelihoodException) {
   int status = 0;
   fitsfile * fp = 0;
   
   fits_open_file(&fp, filename.c_str(), READONLY, &status);
   fits_report_error(stderr, status);
   if (status != 0) {
      throw LikelihoodException("Table::read_FITS_table: cfitsio error");
   }

   int hdutype = 0;
   fits_movabs_hdu(fp, hdu, &hdutype, &status);
   if( status == END_OF_FILE ) status = 0; // ok, I guess??
   fits_report_error(stderr, status);
   
   long lnrows=0;
   fits_get_num_rows(fp, &lnrows, &status);
   fits_report_error(stderr, status);
   
   int ncols=0;
   fits_get_num_cols(fp, &ncols, &status);
   fits_report_error(stderr, status);
   
   // match Column objects, with names already set, 
   // with entries in the FITS table
   for (int i = 0; i < npar(); i++) {
      fits_get_colnum(fp, CASEINSEN, par(i).colname,
                      & par(i).colnum, &status);
      fits_report_error(stderr, status);
      if (status == COL_NOT_FOUND) {
         std::ostringstream errorMessage;
         errorMessage << "Table::read_FITS_table:\n"
                      << "Reading file \""<< filename 
                      << "\"; Column with name " << par(i).colname 
                      << " not found" << std::endl;
         throw LikelihoodException(errorMessage.str());
      }
   }
   
   // now copy the data
   for (int i = 0; i < npar(); i++) {
      fits_get_coltype(fp, par(i).colnum, &par(i).typecode,
                       &par(i).dim, &par(i).width, &status);
      fits_report_error(stderr, status);
      
      par(i).dim *= lnrows;  // need this to get proper dimension of data
      
      par(i).val = new double[par(i).dim];
      
      int anynul=0, nulval=0;
      fits_read_col(fp, TDOUBLE, par(i).colnum, 1, 1, par(i).dim,
                    &nulval, par(i).val, &anynul, &status);
      fits_report_error(stderr, status);
   }

   fits_close_file(fp, &status);
   fits_report_error(stderr, status);
}

void Table::read_FITS_colnames(std::string &filename, int hdu, 
                               std::vector<std::string> &columnNames) 
   throw(LikelihoodException) {
   if (!columnNames.empty()) columnNames.clear();

   int status = 0;
   fitsfile * fp = 0;
   
   fits_open_file(&fp, filename.c_str(), READONLY, &status);
   fits_report_error(stderr, status);
   if( status != 0) {
      std::ostringstream errorMessage;
      errorMessage << "Table::read_FITS_colnames:\n "
                   << "Could not open FITS file " << filename << "\n";
      throw LikelihoodException(errorMessage.str());
   }
   
   int hdutype = 0;
   fits_movabs_hdu(fp, hdu, &hdutype, &status);
   if (status == END_OF_FILE) status = 0;
   fits_report_error(stderr, status);
   
// read in the number of columns 
   long ncols;
   char comment[72];
   fits_read_key_lng(fp, "TFIELDS", &ncols, comment, &status);
   fits_report_error(stderr, status);

// read in the column names
   for (int i = 0; i < ncols; i++) {
      char colname[10];
      sprintf(colname, "TTYPE%i ", i+1);
      char my_string[20];
      fits_read_key_str(fp, colname, my_string, comment, &status);
      fits_report_error(stderr, status);
      columnNames.push_back(std::string(my_string));
   }

   fits_close_file(fp, &status);
   fits_report_error(stderr, status);
}
   

} // namespace Likelihood

