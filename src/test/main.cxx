// test program for Likelihood

#include <cstring>
#include <cmath>

//  include everything for the compiler to test
//  #include "Source.h"
//  #include "DiffuseParametric.h"
//  #include "DiffuseSource.h"
//  #include "PointSource.h"
//  #include "Response.h"
//  #include "SpectralFilter.h"
//  #include "Spectrum.h"
//  #include "Statistic.h"

#include "../Likelihood/Parameter.h"
#include "../Likelihood/Function.h"
#include "../Likelihood/SourceModel.h" 
#include "../Likelihood/Table.h"
#include "MyFun.h"
#include "PowerLaw.h"

using namespace Likelihood;   // for testing purposes only

void test_Parameter_class();
void test_Function_class();
void test_PowerLaw_class();
void test_SourceModel_class();
void test_Table_class();
void report_SrcModel_values(const SourceModel &SrcModel);

int main(int argc, char * argv[]){

   test_Parameter_class();
   test_Function_class();
   test_PowerLaw_class();
   test_SourceModel_class();
//   test_Table_class();

   return 0;
}

/*************************/
/* Parameter class tests */
/*************************/
void test_Parameter_class() {
   std::cout << "Parameter class tests: " << std::endl;

   std::vector<Parameter> my_params;
   Parameter my_param("William", 42.);
   my_params.push_back(my_param);
   
   my_param.setName("Tecumseh");
   my_param.setValue(13.7);
   my_param.setBounds(0., 20.);
   my_params.push_back(my_param);
   
   my_param.setName("Sherman");
   my_param.setValue(25.);
   std::pair<double, double> myBounds(-10, 30);
   my_param.setBounds(myBounds);
   my_params.push_back(my_param);

/* exercise the iterator....requires Parameter copy constructor! */

   std::vector<Parameter>::iterator iter = my_params.begin();
   for (; iter != my_params.end(); iter++) {
      std::cout << (*iter).getName() << ":  " 
                << (*iter).getValue() << "  "
                << (*iter).getBounds().first << "  "
                << (*iter).getBounds().second << "  "
                << (*iter).isFree()
                << std::endl;
   }
   std::cout << std::endl;
} // Parameter class tests

/************************/
/* Function class tests */
/************************/
void test_Function_class() {
   MyFun f; 
   std::cout << "Giving your function a name...." << std::endl;
   f.setMyName("Hal");
   
   std::cout << "Its name is " << f.getMyName()
             << ".  Hi, " << f.getMyName() << "!"
             << std::endl;

/* setting and getting parameter names and values */
   std::cout << "Naming and setting " << f.getMyName()
             << "'s parameters..." << std::endl;

   f.setParam(std::string("Ruthie"), 1.);
   f.setParam(std::string("Mary"), 2.);
   f.setParam(std::string("Jane"), 3e-5);
   f.setParam(std::string("Plain"), 3.1415);

   std::vector<string> my_paramNames = f.getParamNames();
      
   std::cout << "Here they are: " << std::endl;
   for (int i = 0; i < f.getNumParams(); i++) {
      std::cout << my_paramNames[i] << ":  " 
                << f.getParamValue(my_paramNames[i]) 
                << std::endl;
   }
   std::cout << "f(3) = " << f(3) << std::endl;

/* try to access a parameter not named in the function */
   std::cout << f.getParamValue("foo") << std::endl;   

/* reset all of the parameters in one shot */
   std::cout << "Resetting these guys in one shot..." << std::endl;
   std::vector<double> inputVec;
   for (int i=0; i < f.getNumParams(); i++) { 
      inputVec.push_back(double(i*i));
   }
   f.setParamValues(inputVec);

/* change the value of an existing parameter */
   f.setParam(string("Ruthie"), 10.);

/* attempt to change the value of a non-existing parameter */
   f.setParam(string("Oscar"), 5.);
      
   std::cout << "The current set of values: " << std::endl;
   std::vector<double> my_params = f.getParamValues();
   for (unsigned int i = 0; i < my_params.size(); i++) {
      std::cout << my_params[i] << " ";
   }
   std::cout << " f(2) = " << f(2) << std::endl;

/* get derivatives wrt parameters */
   std::cout << "getting derivatives one-by-one:" << std::endl;
   for (unsigned int i = 0; i < my_paramNames.size(); i++) {
      std::cout << my_paramNames[i] << ":  "
                << f.derivByParam(2, my_paramNames[i]) << std::endl;
   }

   std::cout << "all derivatives in one shot:" << std::endl;
   std::vector<double> my_derivs = f.getDerivs(2);
   for (unsigned int i = 0; i < my_paramNames.size(); i++) {
      std::cout << f.derivByParam(2, my_paramNames[i]) << "  ";
   }
   std::cout << std::endl;

/* test of pointers to Parameter */
   Parameter *ptrP = f.getParam(string("Mary"));
   if (ptrP != NULL) {
      std::cout << ptrP->getName() << ":  " 
                << ptrP->getValue() << std::endl;
   }
      
   ptrP = f.getParam(string("Joan"));
   if (ptrP != NULL) {
      std::cout << ptrP->getName() << ":  " 
                << ptrP->getValue() << "\n" << std::endl;
   }
   std::cout << std::endl;

/* test the Function copy constructor */

   MyFun f2 = f;

   for (double x = 0; x < 100.; x += 5.) {
      std::cout << x << "  " 
                << f(x) << "  " 
                << f2(x) << "  " 
                << std::endl;
   }

} // Function class (MyFun) tests

/************************/
/* PowerLaw class tests */
/************************/
void test_PowerLaw_class() {
   PowerLaw pl;

/* inspect the default parameters */
   std::vector<Parameter> my_params = pl.getParams();

   std::vector<Parameter>::iterator iter = my_params.begin();
   for (; iter != my_params.end(); iter++) {
      std::cout << (*iter).getName() << ":  " 
                << (*iter).getValue() << ", "
                << (*iter).isFree() << std::endl;
   }
   std::cout << std::endl;

/* get the free ones */
   my_params = pl.getFreeParams();

   iter = my_params.begin();
   for (; iter != my_params.end(); iter++) {
      std::cout << (*iter).getName() << ":  " 
                << (*iter).getValue() << ", "
                << (*iter).isFree() << std::endl;
   }
   std::cout << std::endl;

/* reset the parameters and compute some values */
   pl.setParam("Prefactor", 2.);
   pl.setParam("Index", -2.2);
   double x;
   for (x = 1.05; x < 1e3; x *= x) {
      std::cout << x << "   " << pl(x) << std::endl;
   }

/* get the derivatives and compare to numerical estimates */
   x = 10.;

   std::vector<double> pl_derivs = pl.getDerivs(x);

   std::vector<std::string> paramNames = pl.getParamNames();
// save current parameter values
   std::vector<double> params_save = pl.getParamValues();

   std::cout << "\nDerivatives: " << std::endl;
   for (unsigned int i = 0; i < pl_derivs.size(); i++) {
      std::cout << pl_derivs[i] << ":  ";

// compute the numerical derivative wrt this parameter
      double delta_param = fabs(params_save[i]/1e5);
      double num_deriv = pl(x);
      pl.setParam(paramNames[i], params_save[i]+delta_param);
      num_deriv -= pl(x);
      num_deriv /= delta_param;
      std::cout << -num_deriv << std::endl;

// reset the parameters for next time around
      pl.setParamValues(params_save);
   }

/* free derivatives */
   paramNames = pl.getFreeParamNames();
   pl_derivs = pl.getFreeDerivs(x);

   std::cout << "\nFree Derivatives: " << std::endl;
   for (unsigned int i = 0; i < pl_derivs.size(); i++) {
      std::cout << paramNames[i] << ":  "
                << pl_derivs[i] << ":  "
                << pl.derivByParam(x, paramNames[i])
                << std::endl;
   }
   std::cout << std::endl;

/* instantiate another power-law and compare its output to that of pl */
   PowerLaw pl2(1., -2.1, 1.);

   for (x = 1.05; x < 1e3; x *= x) {
      std::cout << x << "   " 
                << pl(x) << "  " 
                << pl2(x) << std::endl;
   }
   std::cout << std::endl;

/* check the copy constructor */
   PowerLaw pl3(pl2);

   my_params = pl2.getParams();
   std::cout << "Parameters for pl2:" << std::endl;

   iter = my_params.begin();
   for (; iter != my_params.end(); iter++) {
      std::cout << (*iter).getName() << ":  " 
                << (*iter).getValue() << ", "
                << (*iter).isFree() << std::endl;
   }
   std::cout << std::endl;

   my_params = pl3.getParams();
   std::cout << "Parameters for pl3:" << std::endl;

   iter = my_params.begin();
   for (; iter != my_params.end(); iter++) {
      std::cout << (*iter).getName() << ":  " 
                << (*iter).getValue() << ", "
                << (*iter).isFree() << std::endl;
   }
   std::cout << std::endl;


} // PowerLaw class tests

/***************************/
/* SourceModel class tests */
/***************************/
void test_SourceModel_class() {
      
   SourceModel SrcModel;

/* add some initial functions */
   PowerLaw pl1(1., -2., 1.);
   PowerLaw pl2(3., -2.5, 10.);
   PowerLaw pl3(5., -1.9, 5.);

   SrcModel.addSource(&pl1, string("Joe"));
   SrcModel.addSource(&pl2, string("Blow"));
   SrcModel.addSource(&pl3, string("John"));

// should see 1  -2  3  -2.5  5  -1.9
//            Joe  Blow  John
   report_SrcModel_values(SrcModel);

/* delete a function */
   SrcModel.deleteSource(string("Blow"));

// should see 1  -2  5  -1.9
//            Joe  John
   report_SrcModel_values(SrcModel);

/* add another function */
   PowerLaw pl4(4.2, -1.7, 3.);

/* make its scale factor free */
   pl4.setParam("Scale", pl4.getParamValue("Scale"), true);
   SrcModel.addSource(&pl4, string("Doe"));

// should see 1  -2  5  -1.9  4.2  -1.7  3
//            Joe  John  Doe
   report_SrcModel_values(SrcModel);

/* set Joe's spectral index to be fixed */
   Function *my_func = SrcModel.getFunc("Joe");
   (*my_func).setParam("Index", (*my_func).getParamValue("Index"), false);

   SrcModel.deleteSource("Joe");
   SrcModel.addSource(my_func, "Joe");

// should see 5  -1.9  4.2  -1.7  3  1
//            John  Doe  Joe
   report_SrcModel_values(SrcModel);

/* get the parameter names */

// should see Prefactor  Index  Prefactor  Index  Scale  Prefactor
   std::vector<std::string> paramNames = SrcModel.getParamNames();
   for (unsigned int i = 0; i < paramNames.size(); i++) {
      std::cout << paramNames[i] << "   ";
   }
   std::cout << std::endl;
   std::cout << std::endl;

/* reset Doe's prefactor */
   Parameter my_param("Prefactor", 13.7, true);
   SrcModel.setParam(my_param, "Doe");

// should see 5  -1.9  13.7  -1.7  3  1
//            John  Doe  Joe
   report_SrcModel_values(SrcModel);

/* derivative tests */

   double x(20.);
   std::vector<double> params_save = SrcModel.getParamValues();
   std::vector<double> sm_derivs = SrcModel.getDerivs(x);
   std::vector<double> params = params_save;

   std::cout << "SrcModel Derivatives: " << std::endl;
   for (unsigned int i = 0; i < sm_derivs.size(); i++) {
      std::cout << sm_derivs[i] << ":  ";

// compute the numerical derivative wrt this parameter
      double delta_param = fabs(params_save[i]/1e7);
      double num_deriv = SrcModel(x);
      params[i] += delta_param;
      SrcModel.setParamValues(params);
      num_deriv -= SrcModel(x);
      num_deriv /= delta_param;
      std::cout << -num_deriv << std::endl;

// reset the parameters for next time around
      SrcModel.setParamValues(params_save);
      params = params_save;
   }

} // Source Model Class tests

void report_SrcModel_values(const SourceModel &SrcModel) {

/* retrieve parameter values */
   std::vector<double> paramValues = SrcModel.getParamValues();
   for (unsigned int i = 0; i < paramValues.size(); i++) {
      std::cout << paramValues[i] << "  ";
   }
   std::cout << std::endl;
   
/* retrieve function names */
   std::vector<std::string> srcNames = SrcModel.getSrcNames();
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      std::cout << srcNames[i] << "  ";
   }
   std::cout << std::endl;
   std::cout << std::endl;
}

void test_Table_class() {

/* read in PSF parameters */

   std::string psf_file("psf_lat.fits");
   Table psf_data;

   psf_data.add_columns("ENERGY THETA SIG1_F SIG2_F W");
   psf_data.read_FITS_table(psf_file, 2);

   int nenergy = psf_data[0].dim;
   double *energy = psf_data[0].val;
   int ntheta = psf_data[1].dim;
   double *theta = psf_data[1].val;

   std::cout << "From " << psf_file << ": \n";

   std::cout << "energies: ";
   for (int i = 0; i < nenergy; i++) 
      std::cout << energy[i] << "  ";
   std::cout << std::endl;

   std::cout << "theta values: ";
   for (int i = 0; i < ntheta; i++) 
      std::cout << theta[i] << "  ";
   std::cout << std::endl;
   std::cout << std::endl;

/* read in AEFF parameters */

   std::string aeff_file("aeff_lat.fits");
   Table aeff_data;

   aeff_data.add_columns("ENERGY THETA AEFF_F");
   aeff_data.read_FITS_table(aeff_file, 2);

   nenergy = aeff_data[0].dim;
   energy = aeff_data[0].val;
   ntheta = aeff_data[1].dim;
   theta = aeff_data[1].val;
   
   std::cout << "From " << aeff_file << ": \n";
   std::cout << "energies: ";
   for (int i = 0; i < nenergy; i++) 
      std::cout << energy[i] << "  ";
   std::cout << std::endl;

   std::cout << "theta values: ";
   for (int i = 0; i < ntheta; i++) 
      std::cout << theta[i] << "  ";
   std::cout << std::endl;
   std::cout << std::endl;

/* read in EVENT file */

   std::string event_file("1_src_0000");
   Table event_data;

   event_data.add_columns("RA DEC energy time SC_x SC_y SC_z zenith_angle");
   event_data.read_FITS_table(event_file, 2);

   int nRA = event_data[0].dim;
   double* RA = event_data[0].val;
   int nDEC = event_data[1].dim;
   double* DEC = event_data[1].val;
   nenergy = event_data[2].dim;
   energy = event_data[2].val;
   int ntime = event_data[3].dim;
   double* time = event_data[3].val;

   int nmax = 5;

   std::cout << "Some of the data from " << event_file << ":\n";
   std::cout << nRA << "  " << nDEC << "  "
             << nenergy << "  " << ntime << std::endl;
   for (int i = 0; i < nmax; i++) {
      std::cout << RA[i] << "  "
                << DEC[i] << "  "
                << energy[i] << "  "
                << time[i] << std::endl;
   }
   std::cout << std::endl;

/* Spacecraft data file */

   std::string sc_file("1_src_sc_0000");
   Table sc_data;

   sc_data.add_columns("time SC_x SC_y SC_z SAA_flag");
   sc_data.read_FITS_table(sc_file, 2);

   ntime= sc_data[0].dim;
   time = sc_data[0].val;
   int nSC_x = sc_data[1].dim;
   double* SC_x = sc_data[1].val;
   int nSC_y = sc_data[2].dim;
   double* SC_y = sc_data[2].val;
   int nSC_z = sc_data[3].dim;
   double* SC_z = sc_data[3].val;

   std::cout << "Some of the data from " << sc_file << ":\n";
   std::cout << ntime << "  " << nSC_x << "  "
             << nSC_y << "  " << nSC_z << std::endl;
   for (int i = 0; i < nmax; i++) {
      std::cout << time[i] << "  "
                << SC_x[i] << "  "
                << SC_y[i] << "  "
                << SC_z[i] << std::endl;
   }
   std::cout << std::endl;

} // Table class tests
