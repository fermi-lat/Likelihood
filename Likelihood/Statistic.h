/** 
 * @file Statistic.h
 * @brief Declaration of Statistic class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Statistic.h,v 1.11 2003/03/25 23:22:02 jchiang Exp $
 */

#ifdef _MSC_VER
#pragma warning(disable:4290)
#endif

#ifndef Statistic_h
#define Statistic_h

#include "Likelihood/SourceModel.h"
#include "Likelihood/Table.h"

namespace Likelihood {

/** 
 * @class Statistic
 *
 * @brief Abstract base class for objective functions used for
 * parameter estimation.  Augments the SourceModel class by adding
 * facility for reading and storing data.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Statistic.h,v 1.11 2003/03/25 23:22:02 jchiang Exp $ 
 */

class Statistic : public SourceModel {
    
public:

   virtual ~Statistic() {}

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
      {return getColumn(m_eventData, colname);}

protected:

   //! generalized column access
   std::pair<long, double*> getColumn(const Table &tableData, 
                                      const std::string &colname) const
      throw(LikelihoodException);

   //! Event data; read from m_eventFile, stored in Table form
   std::string m_eventFile;
   int m_eventHdu;
   Table m_eventData;
};

} // namespace Likelihood

#endif // Statistic_h
