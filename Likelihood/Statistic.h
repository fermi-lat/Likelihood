/** @file Statistic.h
 * @brief Declaration of Statistic class
 * $Header:
 */

#ifndef Statistic_h
#define Statistic_h

#include "Table.h"
#include "SourceModel.h"

namespace Likelihood {

/** 
 * @class Statistic
 *
 * @brief Objective functions used for parameter estimation.  Augments
 * the SourceModel class by adding event and spacecraft data.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/Likelihood/Statistic.h,v 1.2 2003/02/23 22:28:59 jchiang Exp $
 */

class Statistic : public SourceModel {
    
public:

   virtual ~Statistic(){};

   //! return the objective function value taking the free parameters 
   //! as the function argument
   virtual double value(const std::vector<double> &paramVec);
   virtual double operator()(const std::vector<double> &paramVec) 
      {return value(paramVec);};

   void readEventData(const std::string &eventFile, 
                      const std::string &colnames, int hdu);

   std::pair<long, double*> getEventColumn(const std::string &colname) const
      {return m_getColumn(m_eventData, colname);};

private:

   //! generalized column access
   std::pair<long, double*> m_getColumn(const Table &tableData, 
                                        const std::string &colname) const;

   //! Event data; read from m_eventFile, stored in Table form
   std::string m_eventFile;
   int m_eventHdu;
   Table m_eventData;
};

} // namespace Likelihood

#endif // Statistic_h
