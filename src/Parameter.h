#ifndef Parameter_h
#define Parameter_h

#include <vector>
#include <string>

/** 
 * @class Parameter
 *
 * @brief Model parameters are identified by a name.  The same name
 * can be given to different parameters of different source model
 * components.  This has the effect of "tying" parameters together.
 * Parameters with null names, "", are anonymous and are treated as
 * "untied".  This should really be implemented as an STL multimap.
 *
 * @authors J. Chiang
 *    
 * $Header$ */

class Parameter {
    
public:
    
   Parameter();
   ~Parameter();

private:
   
   //! parameter name
   std::string m_name;

   //! its value
   double m_value;

   //! lower bound
   double m_minValue;

   //! upper bound
   double m_maxValue;

   //! flag to indicate free or fixed
   bool m_free;

};

#endif // Parameter_h
