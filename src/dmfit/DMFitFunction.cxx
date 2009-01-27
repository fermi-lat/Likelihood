/** 
 * @file DMFitFunction.cxx
 * @brief Implementation for the DMFitFunction class
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/dmfit/DMFitFunction.cxx,v 1.4 2008/10/12 21:25:06 cohen Exp $
 */

#include <cmath>

#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"
#include "optimizers/f2c_types.h"
typedef double doublereal;

#include "Likelihood/DMFitFunction.h"

#include "facilities/Util.h"
#include "st_facilities/Util.h"

#ifdef __cplusplus
extern "C" {
#endif

//#include "f2c/f2c.h"
doublereal yieldget_ (integer *zi, integer *mxi, integer *ch);
doublereal llg_ (doublereal *x, doublereal *mx, doublereal *ml);
integer dmfit_load__ (char *filename, ftnlen filename_len);
doublereal dmfit_de__ (doublereal *mx, integer *ch, doublereal *ee);
doublereal dmfit_dm__ (doublereal *mx, integer *ch, doublereal *ee);

#ifdef __cplusplus
}
#endif

namespace Likelihood {

void DMFitFunction::init(double norm, double mass, double bratio,  
                         int channel0, int channel1) {
   setMaxNumParams(5);
   //should add a check that ch0/1 are integer between 1 and 9

   addParam("norm", norm, true);
   addParam("mass", mass, true);
   addParam("bratio", bratio, true);
   addParam("channel0", channel0, false);
   addParam("channel1", channel1, false);

   setParamAlwaysFixed("channel0");
   setParamAlwaysFixed("channel1");

   m_funcType = Addend;
   m_argType = "dArg";

   m_genericName = "DMFitFunction";
   m_normParName = "norm";

}

double DMFitFunction::value(optimizers::Arg &xarg) const {

  double escale=1.e-3;//GeV
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue()*escale;
   double n=m_parameter[0].getTrueValue();
   double m=m_parameter[1].getTrueValue();
   double b=m_parameter[2].getTrueValue();
//    long int ch0 = m_parameter[3].getTrueValue();
//    long int ch1 = m_parameter[4].getTrueValue();
   integer ch0 = m_parameter[3].getTrueValue();
   integer ch1 = m_parameter[4].getTrueValue();


   double my_value=n*(b*dmfit_de__(&m,&ch0,&x)
                      +(1.-b)*dmfit_de__(&m,&ch1,&x))/1000.;
   return my_value;
}

double DMFitFunction::derivByParam(optimizers::Arg & xarg,
                          const std::string & paramName) const {
   int iparam(-1);
   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
         iparam = i;
      }
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "DMFitFunction::derivByParam");
   }

   enum ParamTypes {norm,mass,bratio,channel0,channel1};
   double escale=1.e-3;//GeV
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue()*escale;
   double n=m_parameter[0].getTrueValue();
   double m=m_parameter[1].getTrueValue();
   double b=m_parameter[2].getTrueValue();
//    long int ch0 = m_parameter[3].getTrueValue();
//    long int ch1 = m_parameter[4].getTrueValue();
   integer ch0 = m_parameter[3].getTrueValue();
   integer ch1 = m_parameter[4].getTrueValue();

   double value(0);
   switch(iparam) {
   case norm:
     return m_parameter[0].getScale()*
       (b*dmfit_de__(&m,&ch0,&x)+(1.-b)*dmfit_de__(&m,&ch1,&x))/1000.;
     break;
   case mass:
     value= m_parameter[1].getScale()*
       n*(b*dmfit_dm__(&m,&ch0,&x)+(1.-b)*dmfit_dm__(&m,&ch1,&x))/1000.;
     //std::cout<<m<<" "<<ch0<<" "<<x<<" "<<dmfit_dm__(&m,&ch0,&x)<<std::endl;
     return value;
     break;
   case bratio:
     return m_parameter[2].getScale()*
       n*(dmfit_de__(&m,&ch0,&x) - dmfit_de__(&m,&ch1,&x))/1000.;
     break;
   }
   return 0;
}

  void DMFitFunction::readFunction(const std::string & filename) {
    m_filename = filename;
    facilities::Util::expandEnvVar(&m_filename);
    st_facilities::Util::file_ok(m_filename);
    //std::cout<<"Loading file "<<m_filename<<std::endl; 
    dmfit_load__(const_cast<char *>(m_filename.c_str()),m_filename.size());
  }

} // namespace Likelihood
