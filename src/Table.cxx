/** @file Table.cxx
* @brief Implementation of Table member functions
*
* $Header:
*/

#include "../Likelihood/Table.h"
#include "fitsio.h"

#include <iostream>

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
void Table::read_FITS_table(std::string filename, int hdu) 
{
    int  status=0;
    fitsfile * fp=0;

    fits_open_file(&fp, filename.c_str(), READONLY, &status);
    if( status !=0) {
        std::cerr << "Could not open FITS file " << filename << std::endl;
        return;
    }
    fits_report_error(stderr, status);

    int hdutype=0;
    fits_movabs_hdu(fp, hdu, &hdutype, &status);
    if( status == END_OF_FILE ) status = 0; // ok, I guess??
    fits_report_error(stderr, status);

    // these are not used???  Oh, yes they are....
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
        if( status== COL_NOT_FOUND ) {
            std::cerr << "Reading file \""<< filename 
                << "\"; Column with name " << par(i).colname 
                << " not found" << std::endl;
            return;
        }
        fits_report_error(stderr, status);
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
}

} // namespace Likelihood

