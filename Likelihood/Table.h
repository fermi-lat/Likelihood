/** 
 * @file Table.h
 * @brief Declaration of Table class
 * @authors T. Burnett, J. Chiang using code from Y. Ikebe
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Table.h,v 1.9 2003/06/10 23:58:30 jchiang Exp $
 */

#ifndef Table_h
#define Table_h

#include <vector>
#include <string>
#include "Likelihood/LikelihoodException.h"

namespace Likelihood {

/** 
 * @class Table
 *
 * @brief Table class for storing FITS table data
 *
 * @author T. Burnett, J. Chiang using code from Y. Ikebe
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Table.h,v 1.9 2003/06/10 23:58:30 jchiang Exp $
 */

class Table {
    
public:
    
   Table() {}
   
   //! add a new column to the list
   void add_column(std::string colname);
   
   //! add several columns, a string delimited by single blanks
   void add_columns(std::string colnames);

   //! add several columns more sensibly using a vector of strings
   void add_columns(std::vector<std::string> columnNames);
   
   //! fill the table from the file
   //! expect that the column list is already set
   void read_FITS_table(std::string file, int hdu) throw(LikelihoodException);
   
   //! grab the column names from an existing FITS file
   void read_FITS_colnames(std::string &file, int hdu, 
                           std::vector<std::string> &columnNames)
      throw(LikelihoodException);

   /** @class Column
    * @brief Nested class represents a column entry
    */
   class Column {
   public:
      Column(std::string name):val(0){
         ::strncpy(colname, name.c_str(), sizeof(colname));
      }
      Column():val(0){}
      ~Column(){delete[] val;}
      
      char   colname[32]; // name of the Ntuple column
      int    colnum;      // column index
      int    typecode;    // data type code    
      long   dim;         //
      long   width;       //
      double *val;        // pointer to array of values
   };
   
   //! external read-write access to column contents
   Column& operator[](int i) { return m_pars[i];}
   
   //! read-only access to column contents
   const Column& operator[](int i)const { return m_pars[i];}
   
   int npar() const { return m_pars.size(); }
   
private:
   Column& par(int i) { return m_pars[i];}
   
   typedef std::vector<Column> ColumnVec;
   typedef ColumnVec::const_iterator const_iterator;
   ColumnVec m_pars;
};

} // namespace Likelihood

#endif // Table.h
