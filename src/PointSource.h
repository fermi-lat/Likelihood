#ifndef PointSource_h
#define PointSource_h

#include "Source.h"
#include "astro/Skydir.h"
/** 
 * @class PointSource
 *
 * @brief A gamma-ray point source, e.g., pulsars, AGNs, X-ray binaries
 *
 * @author J. Chiang
 *    
 * $Header$
 */


class PointSource : public Source {
    
public:
    
    PointSource();
    ~PointSource();

   void setDir();
   astro::SkyDir getDir()const;
    
private:

   //! location on the Celestial sphere (assume the existence of a SkyDir class)
    astro::SkyDir m_phat;
        
};

#endif // PointSource_h
