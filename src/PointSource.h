#ifndef PointSource_h
#define PointSource_h

#include "Source.h"
#include "astro/SkyDir.h"
/** 
 * @class PointSource
 *
 * @brief A gamma-ray point source, e.g., pulsars, AGNs, X-ray binaries
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/PointSource.h,v 1.1.1.1 2003/01/30 23:23:03 burnett Exp $
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
