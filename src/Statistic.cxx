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
#include "PowerLaw.h"

namespace Likelihood {

//! objective function as a function of the free parameters
//! do EML (specialized to Gaussian functions for now)
double Statistic::value(const std::vector<double> &paramVec) {
   setParamValues(paramVec);
   
   double my_value = 0;
   
   for (int j = 0; j < m_eventData[0].dim; j++) {
      double src_sum = 0.;
      for (unsigned int i = 0; i < getNumSrcs(); i++)
// NB: Here evaluate_at(double) is inherited from SourceModel and
// evaluates as a function of the data variable.  This will need
// generalization in SourceModel, i.e., an Event class object should
// be passed instead.
         src_sum += evaluate_at(m_eventData[0].val[j]);
      my_value += log(src_sum);
   }
   for (unsigned int i = 0; i < getNumSrcs(); i++) {
// should typedef this map...
      std::map<std::string, Function *>::iterator 
         func_it = (*m_sources[i]).m_functions.begin();
      for (; func_it != (*m_sources[i]).m_functions.end(); func_it++)
         my_value -= (*func_it).second->integral(-1e3, 1e3);
   }
   
   return my_value;
}

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
