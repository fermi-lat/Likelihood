/** @file Statistic.h
 * @brief Declaration of Statistic class
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Statistic_h
#define Statistic_h

#include "Likelihood/SourceModel.h"
#include "Likelihood/Table.h"

namespace Likelihood {

/** 
 * @class Statistic
 *
 * @brief Absract base class for objective functions used for
 * parameter estimation.  Augments the SourceModel class by adding
 * facility for reading and storing data.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Statistic.h,v 1.8 2003/03/16 21:53:26 jchiang Exp $ */

class Statistic : public SourceModel {
    
public:

   Statistic(){}
   virtual ~Statistic(){}

   //! return the objective function value taking the free parameters 
   //! as the function argument
   virtual double value(const std::vector<double> &paramVec) = 0;
   virtual double operator()(const std::vector<double> &paramVec) 
      {return value(paramVec);}

   //! non-argument version of getFreeDerivs
   virtual void getFreeDerivs(std::vector<double> &) = 0;

   void readEventData(const std::string &eventFile, 
                      const std::string &colnames, int hdu);

   std::pair<long, double*> getEventColumn(const std::string &colname) const
      {return m_getColumn(m_eventData, colname);}

protected:

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
