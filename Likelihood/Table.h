/** @file Table.h
 * @brief Declaration of Table class
 * $Header:
 */

#ifndef Table_h
#define Table_h

#include <vector>
#include <string>

namespace Likelihood {

/** 
 * @class Table
 *
 * @brief Table class for storing FITS table data
 *
 * @author
 *    
 * $Header:
 */

class Table {
    
public:
    
    Table(){};

    //! add a new column to the list
    void add_column(std::string colname);

    //! add several columns, a string delimited by single blanks
    void add_columns(std::string colnames);

    //! fill the table from the file
    //! expect that the column list is already
    void read_FITS_table(std::string file, int hdu);

    /** @class Column
    * @brief Nested class represents a column entry
    */
    class Column {
    public:
        Column(std::string name):val(0){
            ::strncpy(colname,name.c_str(), sizeof(colname));
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
