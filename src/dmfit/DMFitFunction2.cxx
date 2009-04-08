/** 
 * @file DMFitFunction2.cxx
 * @brief Implementation for the DMFitFunction2 class
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/dmfit/DMFitFunction2.cxx,v 1.3 2009/03/19 22:24:55 cohen Exp $
 */

#include <cmath>

#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"
#include "optimizers/f2c_types.h"
typedef double doublereal;

#include "Likelihood/DMFitFunction2.h"

#include "facilities/Util.h"
#include "st_facilities/Util.h"

#ifdef __cplusplus
extern "C" {
#endif

doublereal yieldget_ (integer *zi, integer *mxi, integer *ch);
doublereal llg_ (doublereal *x, doublereal *mx, doublereal *ml);
integer dmfit_load__ (char *filename, ftnlen filename_len);
doublereal dmfit_de__ (doublereal *mx, integer *ch, doublereal *ee);
doublereal dmfit_dm__ (doublereal *mx, integer *ch, doublereal *ee);
doublereal dmfit_deint__(doublereal *mx, integer *ch, doublereal *ee);
doublereal dmfit_dmint__(doublereal *mx, integer *ch, doublereal *ee);

#ifdef __cplusplus
}
#endif

namespace Likelihood {

void DMFitFunction2::init(double norm, double mass, double bratio,  
                         int channel0, int channel1,
                         double lowerLimit) {
   setMaxNumParams(6);
   //should add a check that ch0/1 are integer between 1 and 9

   addParam("norm", norm, true);
   addParam("mass", mass, true);
   addParam("bratio", bratio, true);
   addParam("channel0", channel0, false);
   addParam("channel1", channel1, false);
   addParam("LowerLimit", lowerLimit, false);

   setParamAlwaysFixed("channel0");
   setParamAlwaysFixed("channel1");
   setParamAlwaysFixed("LowerLimit");

   m_funcType = Addend;
   m_argType = "dArg";

   m_genericName = "DMFitFunction2";
   m_normParName = "norm";

}

double DMFitFunction2::value(optimizers::Arg &xarg) const {

  double escale=1.e-3;//GeV
  double x = dynamic_cast<optimizers::dArg &>(xarg).getValue()*escale;
  double n=m_parameter[0].getTrueValue();
  double m=m_parameter[1].getTrueValue();
  double b=m_parameter[2].getTrueValue();
  integer ch0 = static_cast<integer>(m_parameter[3].getTrueValue());
  integer ch1 = static_cast<integer>(m_parameter[4].getTrueValue());
  double emin = m_parameter[5].getTrueValue()*escale;

  updateCache(x, m, b, ch0, ch1, emin);

  double N = (b*dmfit_de__(&m,&ch0,&x)
              +m_1mb*dmfit_de__(&m,&ch1,&x))*escale;

  //the integrated spectrum
  double D = b * m_int_dnde0 + m_1mb * m_int_dnde1;

  double my_value=n*N/D;

  return my_value;
}

double DMFitFunction2::derivByParam(optimizers::Arg & xarg,
                          const std::string & paramName) const {
   int iparam(-1);
   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
         iparam = i;
      }
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "DMFitFunction2::derivByParam");
   }

   enum ParamTypes {norm,mass,bratio,channel0,channel1};
   double escale=1.e-3;//GeV
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue()*escale;
   double n=m_parameter[0].getTrueValue();
   double m=m_parameter[1].getTrueValue();
   double b=m_parameter[2].getTrueValue();
   integer ch0 = static_cast<integer>(m_parameter[3].getTrueValue());
   integer ch1 = static_cast<integer>(m_parameter[4].getTrueValue());
   double emin = m_parameter[5].getTrueValue()*escale;

   updateCache(x, m, b, ch0, ch1, emin);


   double N = (b*dmfit_de__(&m,&ch0,&x)
               + m_1mb*dmfit_de__(&m,&ch1,&x))*escale;
   double D = b * m_int_dnde0 + m_1mb * m_int_dnde1;
   
   double value(0);
   switch(iparam) {
   case norm:
     return m_parameter[0].getScale()*N/D;
     break;
   case mass:
     value= m_parameter[1].getScale()* n*
       (
        (b*dmfit_dm__(&m,&ch0,&x)+m_1mb*dmfit_dm__(&m,&ch1,&x))* escale
        - N/D*( b*m_int_dndm0 + m_1mb*m_int_dndm1 )
        )/D;
     return value;
     break;
   case bratio:
     return m_parameter[2].getScale()*
       n*(
          (dmfit_de__(&m,&ch0,&x) - dmfit_de__(&m,&ch1,&x))* escale
          - N/D*( m_int_dnde0 - m_int_dnde1 )
          )/D;
     break;
   }
   return 0;
}

  void DMFitFunction2::readFunction(const std::string & filename) {
// Save data member version, preserving any environment variables.
    m_filename = filename; 
// Use expanded local copy for reading the data.
    std::string expanded_filename = filename;
    facilities::Util::expandEnvVar(&expanded_filename);
    st_facilities::Util::file_ok(expanded_filename);
    //std::cout<<"Loading file "<<expanded_filename<<std::endl; 
    dmfit_load__(const_cast<char *>(expanded_filename.c_str()),
                 expanded_filename.size());

    //Check and enforce XML input non-degeneracy
    integer ch0 = static_cast<integer>(m_parameter[3].getTrueValue());
    integer ch1 = static_cast<integer>(m_parameter[4].getTrueValue());
    if((ch0==ch1) && getParam("bratio").isFree())
      {
        std::cout<<"bratio set fixed to 1, as channels 0 and 1 are identical"<<std::endl;
        setParam("bratio",1.);
        parameter("bratio").setFree(false);
      }
  }


  void DMFitFunction2::updateCache(double x,
                                  double m, 
                                  double b,
                                  int ch0,
                                  int ch1,
                                  double emin)     
    const {
     (void)(x);
    if(m != m_dmMass)
      {
        m_int_dnde0 = dmfit_deint__(&m,&ch0,&emin);
        m_int_dnde1 = dmfit_deint__(&m,&ch1,&emin);
        m_int_dndm0 = dmfit_dmint__(&m,&ch0,&emin);
        m_int_dndm1 = dmfit_dmint__(&m,&ch1,&emin);
        m_dmMass = m;
      }
    if(1.-b != m_1mb)
      {
        m_1mb = 1.-b;
      }
  }
  
} // namespace Likelihood
