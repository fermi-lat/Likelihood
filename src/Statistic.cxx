/** @file Statistic.cxx
 * @brief Statistic class implementation
 *
 * $Header:
 */

#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include "../Likelihood/Statistic.h"

#include "../Likelihood/Table.h"

namespace Likelihood {

//! read in the event data
void Statistic::readEventData(const std::string &eventFile, 
                              const std::string &colnames, int hdu) {
   m_eventFile = eventFile;
   m_eventHdu = hdu;
   m_eventData.add_columns(colnames);
   m_eventData.read_FITS_table(m_eventFile, m_eventHdu);
}

//! return pointer to data columns
std::pair<long, double*> 
Statistic::m_getColumn(const Table &tableData, 
                       const std::string &colname) const{
   std::pair<long, double*> my_column(0, 0);
   std::string colnames;

// loop over column names, return the matching one
   for (int i = 0; i < tableData.npar(); i++) {
      if (colname == std::string(tableData[i].colname)) {
         my_column.first = tableData[i].dim;
         my_column.second = tableData[i].val;
         return my_column;
      }
      colnames += " "; colnames += tableData[i].colname;
   }
   std::cerr << "Column " << colname << " was not found in event data.\n";
   std::cerr << "Valid names are \n" << colnames << std::endl;
   std::cerr << std::endl;
   return my_column;
}

} // namespace Likelihood
