#ifndef SciTools_Function_h
#define SciTools_Function_h

#include <vector>
#include <string>
#include <map>

/** 
 * @class Function
 *
 * @brief Base class for Science Tools Functions, i.e., things that are
 * not data containers per se, such as source models, fit statistics, etc.
 *
 * @author J. Chiang, P. Nolan
 *    
 * $Header:$
 */

class Function {

public:
    
    Function(){};
    virtual ~Function(){};

   void setParam(std::string paramName, double paramValue);
   double getParam(std::string paramName) const;
   
   //! function value at current parameters pure virtual, must be implemented by subclass
   virtual double value() const =0;
   
   //! allow any subclass to behave like a function object (with no arguments)
   double operator()()const { return value();}

   //! partial derivative wrt param (virtual, with default)
   virtual double derivByParam(std::string paramName) const {return 0;};

   //! derivatives as a group
   virtual std::vector<double> getAllDerivs() const {return std::vector<double> ();}

   //! provide a string identifier
   void setMyName(std::string);
   std::string getMyName() const;

   //! parameter access in groups
   int getNumParams() const {return m_param.size();}
   std::vector<std::string> getParamNames() const;
   std::vector<double> getAllParams() const;
   void setAllParams(std::vector<double> paramList);

private:

   std::vector<std::string> m_paramNames;
   std::vector<double> m_param;
    
   // might be easier to implement this way
   std::map<std::string, double> m_param_map;
        
};

#endif // SciTools_Function_h
