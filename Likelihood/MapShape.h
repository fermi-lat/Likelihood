#ifndef Likelihood_MapShape_h
#define Likelihood_MapShape_h

#include <string>
#include <vector>

#include "astro/SkyDir.h"

namespace Likelihood {

class MapShape {

public:

//   MapShape(const std::string & fitsFile);

   MapShape(const std::vector<double> & ra, 
            const std::vector<double> & dec,
            const std::vector<double> & energy) 
      : m_x(ra), m_y(dec), m_z(energy), m_coordType(astro::SkyDir::EQUATORIAL)
      {}

   unsigned int nx() const {return m_x.size();}
   unsigned int ny() const {return m_y.size();}
   unsigned int nz() const {return m_z.size();}

   const std::vector<double> & x_vector() const {return m_x;}
   const std::vector<double> & y_vector() const {return m_y;}
   const std::vector<double> & z_vector() const {return m_z;}

   astro::SkyDir::CoordSystem coordType() const {return m_coordType;}

private:
   
   std::vector<double> m_x;
   std::vector<double> m_y;
   std::vector<double> m_z;
   astro::SkyDir::CoordSystem m_coordType;

};

} // namespace Likelihood

#endif // Likelihood_MapShape_h
