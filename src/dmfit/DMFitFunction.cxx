/** 
 * @file DMFitFunction.cxx
 * @brief Implementation for the DMFitFunction class
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/dmfit/DMFitFunction.cxx,v 1.15 2013/02/07 19:04:24 cohen Exp $
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

doublereal yieldget_ (integer *zi, integer *mxi, integer *ch);
doublereal llg_ (doublereal *x, doublereal *mx, doublereal *ml);
integer dmfit_load__ (char *filename, ftnlen filename_len);
doublereal dmfit_de__ (doublereal *mx, integer *ch, doublereal *ee);
doublereal dmfit_dm__ (doublereal *mx, integer *ch, doublereal *ee);

#ifdef __cplusplus
}
#endif

namespace Likelihood {

static  enum  {norm,sigmav,mass,bratio,channel0,channel1} ParamType;

DMFitFunction::
DMFitFunction(double norm, double sigmav, double mass, double bratio,  
              int channel0, int channel1) 
   : optimizers::Function("DMFitFunction", 6, "sigmav", "dArg", Addend),
     m_filename(""), m_8pi(8.*M_PI) {
   //should add a check that ch0/1 are integer between 1 and 9
   addParam("norm", norm, false);
   addParam("sigmav", sigmav, true);
   addParam("mass", mass, true);
   addParam("bratio", bratio, true);
   addParam("channel0", channel0, false);
   addParam("channel1", channel1, false);

   //setParamAlwaysFixed("norm");
   setParamAlwaysFixed("channel0");
   setParamAlwaysFixed("channel1");
}

double DMFitFunction::value(optimizers::Arg &xarg) const {

  double escale=1.e-3;//GeV
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue()*escale;
   double n=m_parameter[norm].getTrueValue();
   double sv=m_parameter[sigmav].getTrueValue();
   double m=m_parameter[mass].getTrueValue();
   double b=m_parameter[bratio].getTrueValue();
   integer ch0 = static_cast<integer>(m_parameter[channel0].getTrueValue());
   integer ch1 = static_cast<integer>(m_parameter[channel1].getTrueValue());


   double my_value=n*sv/m/m*(b*dmfit_de__(&m,&ch0,&x)
                      +(1.-b)*dmfit_de__(&m,&ch1,&x))/1000.;
   return my_value/m_8pi;
}

double DMFitFunction::derivByParamImp(optimizers::Arg & xarg,
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

   double escale=1.e-3;//GeV
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue()*escale;
   double n=m_parameter[norm].getTrueValue();
   double sv=m_parameter[sigmav].getTrueValue();
   double m=m_parameter[mass].getTrueValue();
   double b=m_parameter[bratio].getTrueValue();
   integer ch0 = static_cast<integer>(m_parameter[channel0].getTrueValue());
   integer ch1 = static_cast<integer>(m_parameter[channel1].getTrueValue());
   
   double m2=m*m;
   double value(0);
   switch(iparam) {
   case norm:
     value = m_parameter[norm].getScale()*sv/m2*
       (b*dmfit_de__(&m,&ch0,&x)+(1.-b)*dmfit_de__(&m,&ch1,&x))/1000.;
     break;
   case sigmav:
     value = n*m_parameter[sigmav].getScale()/m2*
       (b*dmfit_de__(&m,&ch0,&x)+(1.-b)*dmfit_de__(&m,&ch1,&x))/1000.;
     break;
   case mass:
     value = m_parameter[1].getScale()* n*sv/m2 
       * ( b*dmfit_dm__(&m,&ch0,&x)+(1.-b)*dmfit_dm__(&m,&ch1,&x)
          -2./m*(b*dmfit_de__(&m,&ch0,&x)+(1.-b)*dmfit_de__(&m,&ch1,&x)) 
           )/1000.;
     break;
   case bratio:
     value = m_parameter[2].getScale() * n*sv/m2 *
       (dmfit_de__(&m,&ch0,&x) - dmfit_de__(&m,&ch1,&x))/1000.;
     break;
   }
   return value/m_8pi;
}

  void DMFitFunction::readFunction(const std::string & filename) {
    // Save data member version, preserving any environment variables.
    m_filename = filename; 
    if(m_filename.empty()){
      m_filename="$(LIKELIHOODDATAPATH)/gammamc_dif.dat";
    }
    // Use expanded local copy for reading the data.
    std::string expanded_filename = m_filename;
    facilities::Util::expandEnvVar(&expanded_filename);
    //std::cout<<filename<<" "<<expanded_filename<<std::endl;
    st_facilities::Util::file_ok(expanded_filename);
    dmfit_load__(const_cast<char *>(expanded_filename.c_str()),
                 expanded_filename.size());
  }

} // namespace Likelihood
