#ifndef EventNtuple_h
#define EventNtuple_h

#include "FitsNtuple.h"

/** 
 * @class EventNtuple
 *
 * @brief NTuple for LAT photon event data.  Augments FitsNtuple 
 * by allowing for an energy dispersion function that is unique
 * to each event to be included as a class member.
 *
 * @author J. Chiang
 *    
 * $Header$
 */

class EventNtuple : FitsNtuple {
    
public:
    
    EventNtuple();
    ~EventNtuple();

    double energyDispersion(double true_energy, double observed_energy);
    
private:

    //! number of energy values
    long m_neTrue;

    //! energy abscissa values
    std::vector<double> m_trueEnergy;

    //! number of energy values
    long m_neObs;

    //! energy ordinate values
    std::vector<double> m_obsEnergy;

    //! distribution function; this needs to be an m_neTrue by m_neObs matrix
    std::vector<double> m_energyDisp;
        
};

#endif // EventNtuple_h
