/** 
 * @file Statistic.cxx
 * @brief Statistic class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Statistic.cxx,v 1.9 2003/05/29 20:10:46 jchiang Exp $
 */

#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <sstream>

#include "Likelihood/Statistic.h"
#include "Likelihood/Table.h"

namespace Likelihood {

//! read in the event data
void Statistic::readEventData(const std::string &eventFile, 
                              const std::string &colnames, int hdu) {
   m_eventFile = eventFile;
   m_eventHdu = hdu;
//   m_eventData.add_columns(colnames);
   std::vector<std::string> columnNames;
   m_eventData.read_FITS_colnames(m_eventFile, m_eventHdu, columnNames);
   m_eventData.add_columns(columnNames);
   std::cerr << "Columns in " << m_eventFile 
             << ", HDU " << m_eventHdu 
             << ": \n";
   for (unsigned int i = 0; i < columnNames.size(); i++) {
      std::cerr << columnNames[i] << "  ";
   }
   std::cerr << std::endl;
   m_eventData.read_FITS_table(m_eventFile, m_eventHdu);
}

//! return pointer to data columns
std::pair<long, double*> Statistic::getColumn(const Table &tableData, 
                                              const std::string &colname) const
   throw(LikelihoodException) {
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
   std::ostringstream errorMessage;
   errorMessage << "Statistic::getColumn:\n"
                << "Column " << colname << " was not found in event data.\n"
                << "Valid names are \n" << colnames << "\n";
   throw LikelihoodException(errorMessage.str());
}

} // namespace Likelihood
