#ifndef Constraint_h
#define Constraint_h


/** 
 * @class Constraint
 *
 * @brief Parameter constraints
 *
 * @author P. Nolan, J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/Constraint.h,v 1.1.1.1 2003/01/30 23:23:02 burnett Exp $
 */

class Constraint {
    
public:
   //? Constraint( enum Constraint::Type type = EQ , double min=0, double max=0) : m_type(type){};
   ~Constraint();

   enum Type { EQ, GT, LT, RANGE};

   Type getType() const;
   void setType(Type);
   void setMin(double);
   void setMax(double);
   double getMin() const;
   double getMax() const;
   bool satisfied(double x) const;
    
private:
   Type m_type;
   double m_min, m_max; 
        
};

#endif // Constraint_h
