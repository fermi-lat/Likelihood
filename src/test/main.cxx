// test program for Likelihood

#include <cstring>
#include <cmath>

//  include everything for the compiler to test
//  #include "DiffuseParametric.h"
//  #include "DiffuseSource.h"
//  #include "SpectralFilter.h"
//  #include "Spectrum.h"

#include "../Likelihood/Parameter.h"
#include "../Likelihood/Function.h"
#include "../Likelihood/SourceModel.h" 
#include "../Likelihood/Table.h"
#include "../Likelihood/Statistic.h"
#include "../Likelihood/Event.h"
#include "../Likelihood/Source.h"
#include "../Likelihood/PointSource.h"
#include "../Likelihood/Response.h"
#include "../Likelihood/Aeff.h"
#include "../Likelihood/Psf.h"
#include "MyFun.h"
#include "PowerLaw.h"
#include "Gaussian.h"

using namespace Likelihood;   // for testing purposes only

void test_Parameter_class();
void test_Function_class();
void test_PowerLaw_class();
void test_SourceModel_class();
void report_SrcModel_values(const SourceModel &SrcModel);
void test_Table_class();
void test_Statistic_class();
void test_Event_class();
void test_PointSource_class();
void test_Aeff_class();
void test_Psf_class();

std::string test_path;

int main(){
/* get root path to test data */   
   const char * root = ::getenv("LIKELIHOODROOT");
   if (!root) {  //use relative path from cmt directory
      test_path = "../src/test/";
   } else {
      test_path = std::string(root) + "/src/test/";
   }

   test_Parameter_class();
   test_Function_class();
   test_PowerLaw_class();
   test_SourceModel_class();
   test_Table_class();
   test_Statistic_class();
   test_Event_class();
   test_PointSource_class();
   test_Aeff_class();
   test_Psf_class();
   return 0;
}

/********************/
/* Psf class tests */
/********************/
void test_Psf_class() {
   Psf *psf = Psf::instance();

   std::string psf_file = test_path + "/CALDB/psf_lat.fits";
   int psf_hdu = 4;

   psf->readPsfData(psf_file, psf_hdu);

   int nenergy = 10;
   double emin = 0.03;
   double emax = 30.;
   double estep = log(emax/emin)/(nenergy - 1);

   std::cout << "\nPoint Spread Function data:" << std::endl;
   for (double inc = 0.; inc <= 70.; inc += 10.) {
      for (int i = 0; i < nenergy; i++) {
	 double energy = emin*exp(estep*i);
	 std::cout << energy << "  ";

	 std::vector<double> psf_params;
	 psf->fillPsfParams(energy, inc, psf_params);
	 std::cout << psf_params[0] << "  "
		   << psf_params[1] << "  "
		   << psf_params[2] << "\n";
      }
      std::cout << std::endl;
   }
} // Psf class tests

/********************/
/* Aeff class tests */
/********************/
void test_Aeff_class() {
   Aeff *aeff = Aeff::instance();

   std::string aeff_file = test_path + "/CALDB/aeff_lat.fits";
   int aeff_hdu = 4;

   aeff->readAeffData(aeff_file, aeff_hdu);

   int nenergy = 10;
   double emin = 0.03;
   double emax = 30.;
   double estep = log(emax/emin)/(nenergy - 1);

   std::cout << "\nEffective area data:" << std::endl;
   for (int i = 0; i < nenergy; i++) {
      double energy = emin*exp(estep*i);
      std::cout << energy << "  ";
      for (double inc = 0.; inc <= 70.; inc += 10.) {
	 std::cout << (*aeff)(energy, inc) << "  ";
      }
      std::cout << std::endl;
   }
   std::cout << std::endl;
} // Aeff class tests

/***************************/
/* PointSource class tests */
/***************************/
void test_PointSource_class() {

/* instantiate the Psf and read in the Psf and spacecraft data */

   Psf * psf = Psf::instance();

   std::string psf_file = test_path + "/CALDB/psf_lat.fits";
   int psf_hdu = 4;

   psf->readPsfData(psf_file, psf_hdu);

   std::string sc_file = test_path + "/Data/one_src_sc_0000";
   int sc_hdu = 2;

   psf->readScData(sc_file, sc_hdu);

/* instantiate the Aeff and read in the Aeff data */

   Aeff * aeff = Aeff::instance();

   std::string aeff_file = test_path + "/CALDB/aeff_lat.fits";
   int aeff_hdu = 4;

   aeff->readAeffData(aeff_file, aeff_hdu);

/* put this source at the center of the extraction region for
   Data/one_src_0000 (ra, dec) = (193.98, -5.82) */
   
   PointSource my_ptsrc(193.98, -5.82);

   my_ptsrc.setName("3C_279");

/* create a PowerLaw object for the source spectrum */

   PowerLaw source_pl(1., -2.1, 0.1);  // f(E) = 1.*(E/0.1 GeV)^(-2.1)
   my_ptsrc.setSpectrum(&source_pl);

/* compute the source flux at various energies and aspect angles */

   int nenergy = 10;
   double emin = 0.03;
   double emax = 30.;
   double estep = log(emax/emin)/(nenergy - 1);

   std::cout << "\nThe on-source counts spectrum of "
             << my_ptsrc.getName() 
             << " as a function of time: " 
             << std::endl;

   double maxtime = 95.*60.;
   double tstep = 1e2;
   double photon_flux;
   for (double time = 0; time < maxtime; time += tstep) {
      std::cout << "Time = " << time << std::endl;
      for (int i = 0; i < nenergy; i++) { // loop over energies
         double energy = emin*exp(estep*i);
         if (i == 0) 
            photon_flux = my_ptsrc.fluxDensity(energy, time, my_ptsrc.getDir());
         if (photon_flux > 0) {
            std::cout << energy << "  "
                      << my_ptsrc.fluxDensity(energy, time, my_ptsrc.getDir())
                      << std::endl;
         }
      }
      if (photon_flux > 0) std::cout << std::endl;
   }

} // PointSource class tests

/*********************/
/* Event class tests */
/*********************/
void test_Event_class() {

// instantiate a Statistic, read in the event data, then
// stuff them into an Event class vector

   Statistic logLike;

/* read in EVENT file */
   std::string event_file = test_path + "/Data/one_src_0000";

   logLike.readEventData(event_file,
                         "RA DEC energy time SC_x SC_y SC_z zenith_angle", 2);

   typedef std::pair<long, double*> tableColumn;
   tableColumn ra = logLike.getEventColumn("RA");
   tableColumn dec = logLike.getEventColumn("DEC");
   tableColumn energy = logLike.getEventColumn("energy");
   tableColumn time = logLike.getEventColumn("time");
   tableColumn sc_x = logLike.getEventColumn("SC_x");
   tableColumn sc_y = logLike.getEventColumn("SC_y");
   tableColumn sc_z = logLike.getEventColumn("SC_z");
   tableColumn zenangle = logLike.getEventColumn("zenith_angle");

   std::vector<Event> my_events;

   for (int i = 0; i < ra.first; i++) {
// compute sc_ra and sc_dec from direction cosines (it would be nice
// if SkyDir could do this for me....)
      double sc_ra = atan2(sc_y.second[i], sc_x.second[i])*180./M_PI;
      double sc_dec = asin(sc_z.second[i])*180./M_PI;

      Event thisEvent(ra.second[i], dec.second[i], energy.second[i],
                      time.second[i], sc_ra, sc_dec,
                      cos(zenangle.second[i]*M_PI/180.));
      my_events.push_back(thisEvent);
   }

   astro::SkyDir _3c279_dir(193.98, -5.82);

/* write out what's available for nmax of these guys */
// since I know these photons were simulated for 3C 279, write out the
// difference for later binning to get the psf

   int nmax = 20;
   for (int i = 0; i < nmax; i++) {
      astro::SkyDir thisEventDir = my_events[i].getDir();
      std::cout << thisEventDir.ra() << "  " 
                << thisEventDir.dec() << "  "
                << my_events[i].getEnergy() << "  "
                << my_events[i].getArrTime() << "  "
                << my_events[i].getSeparation(_3c279_dir)*180./M_PI
                << std::endl;
   }
} // Event class tests

/*************************/
/* Statistic class tests */
/*************************/
void test_Statistic_class() {
   Statistic logLike;

/* SourceModel function access */

   PowerLaw pl1(2., -2.5, 1.);
   logLike.addSource(&pl1, "power_law_1");

   PowerLaw pl2(10., -1.7, 1.);
   logLike.addSource(&pl2, "power_law_2");

   PowerLaw pl3(1., -2., 20.);
   logLike.addSource(&pl3, "power_law_3");

   report_SrcModel_values(logLike);

/* read in EVENT file */
   std::string event_file = test_path + "Data/normal_dist.fits";

   logLike.readEventData(event_file, "xi", 2);

   std::pair<long, double*> xi = logLike.getEventColumn("xi");

   std::cout << "Some of the data from the event file: \n";
   std::cout << xi.first << std::endl;
   int nmax = 10;
   for (int i = 0; i < nmax; i++) {
      std::cout << xi.second[i] 
		<< std::endl;
   }
   std::cout << std::endl;

/* try to access a column that doesn't exist */

   xi = logLike.getEventColumn("foo");

/* delete existing functions and replace with a Gaussian */

   std::vector<std::string> fNames;
   logLike.getSrcNames(fNames);
   for (unsigned int i = 0; i < fNames.size(); i++) 
      logLike.deleteSource(fNames[i]);
   
// these are the numbers used to generate the data in "Data/normal_dist.fits"
   Gaussian gauss(1e4, 20., 5.);

   logLike.addSource(&gauss, "Gauss_model");

   report_SrcModel_values(logLike);

/* evaluate the integral over x */

   std::cout << gauss.integral(-1e3, 1e3) << std::endl << std::endl;

/* Compute the log-likelihood of this model */

// first get the vector of parameters
   std::vector<double> params;
   logLike.getParamValues(params);

// now, vary over the first parameter and compute the log-likelihood
// at each step

   double xmin = 5e3;
   double xmax = 1.5e4;
   int nx = 20;
   double xstep = log(xmax/xmin)/(nx-1.);

   for (int i = 0; i < nx; i++) {
      params[0] = xmin*exp(xstep*i);
      std::cout << params[0] << "  " << logLike(params) << std::endl;
   }
} // Statistic class tests

/*********************/
/* Table class tests */
/*********************/
void test_Table_class() {

/* read in PSF parameters */
   std::string psf_file = test_path + "/CALDB/psf_lat.fits";

   Table psf_data;
      
   psf_data.add_columns("ENERGY THETA SIG1_F SIG2_F W");
   psf_data.read_FITS_table(psf_file, 2);

   int nenergy = psf_data[0].dim;
   double *energy = psf_data[0].val;
   int ntheta = psf_data[1].dim;
   double *theta = psf_data[1].val;

   std::cout << "\nFrom " << psf_file << ": \n";

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
   std::string aeff_file = test_path + "CALDB/aeff_lat.fits";

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
   std::string event_file = test_path + "/Data/one_src_0000";

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
   std::string sc_file = test_path + "/Data/one_src_sc_0000";

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
   std::vector<std::string> paramNames;
   SrcModel.getParamNames(paramNames);
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
   std::vector<double> params_save;
   SrcModel.getParamValues(params_save);
   std::vector<double> sm_derivs;
   SrcModel.getDerivs(x, sm_derivs);
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
   std::cout << std::endl;

} // Source Model Class tests

void report_SrcModel_values(const SourceModel &SrcModel) {

/* retrieve parameter values */
   std::vector<double> paramValues;
   SrcModel.getParamValues(paramValues);
   for (unsigned int i = 0; i < paramValues.size(); i++) {
      std::cout << paramValues[i] << "  ";
   }
   std::cout << std::endl;
   
/* retrieve function names */
   std::vector<std::string> srcNames;
   SrcModel.getSrcNames(srcNames);
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      std::cout << srcNames[i] << "  ";
   }
   std::cout << std::endl;
   std::cout << std::endl;
}

/************************/
/* PowerLaw class tests */
/************************/
void test_PowerLaw_class() {
   PowerLaw pl;

/* inspect the default parameters */
   std::vector<Parameter> my_params;
   pl.getParams(my_params);

   std::vector<Parameter>::iterator iter = my_params.begin();
   for (; iter != my_params.end(); iter++) {
      std::cout << (*iter).getName() << ":  " 
                << (*iter).getValue() << ", "
                << (*iter).isFree() << std::endl;
   }
   std::cout << std::endl;

/* get the free ones */
   pl.getFreeParams(my_params);

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

   std::vector<double> pl_derivs;
   pl.getDerivs(x, pl_derivs);

   std::vector<std::string> paramNames;
   pl.getParamNames(paramNames);
// save current parameter values
   std::vector<double> params_save;
   pl.getParamValues(params_save);

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
   pl.getFreeParamNames(paramNames);
   pl.getFreeDerivs(x, pl_derivs);

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

   pl2.getParams(my_params);
   std::cout << "Parameters for pl2:" << std::endl;

   iter = my_params.begin();
   for (; iter != my_params.end(); iter++) {
      std::cout << (*iter).getName() << ":  " 
                << (*iter).getValue() << ", "
                << (*iter).isFree() << std::endl;
   }
   std::cout << std::endl;

   pl3.getParams(my_params);
   std::cout << "Parameters for pl3:" << std::endl;

   iter = my_params.begin();
   for (; iter != my_params.end(); iter++) {
      std::cout << (*iter).getName() << ":  " 
                << (*iter).getValue() << ", "
                << (*iter).isFree() << std::endl;
   }
   std::cout << std::endl;


} // PowerLaw class tests

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

   std::vector<string> my_paramNames;
   f.getParamNames(my_paramNames);
      
   std::cout << "Here they are: " << std::endl;
   for (unsigned int i = 0; i < f.getNumParams(); i++) {
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
   for (unsigned int i=0; i < f.getNumParams(); i++) { 
      inputVec.push_back(double(i*i));
   }
   f.setParamValues(inputVec);

/* change the value of an existing parameter */
   f.setParam(string("Ruthie"), 10.);

/* attempt to change the value of a non-existing parameter */
   f.setParam(string("Oscar"), 5.);
      
   std::cout << "The current set of values: " << std::endl;
   std::vector<double> my_params;
   f.getParamValues(my_params);
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
   std::vector<double> my_derivs;
   f.getDerivs(2, my_derivs);
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
