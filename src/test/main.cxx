// test program for Likelihood

//#define HAVE_OPTIMIZERS

#include <cstring>
#include <cmath>

//  include everything for the compiler to test
//  #include "DiffuseParametric.h"
//  #include "DiffuseSource.h"

#include "Likelihood/Parameter.h"
#include "Likelihood/Function.h"
#include "Likelihood/SourceModel.h" 
#include "Likelihood/Table.h"
#include "Likelihood/Statistic.h"
#include "Likelihood/Event.h"
#include "Likelihood/Source.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/Response.h"
#include "Likelihood/Aeff.h"
#include "Likelihood/Psf.h"
#include "Likelihood/ScData.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/dArg.h"
#include "Likelihood/SumFunction.h"
#include "Likelihood/ProductFunction.h"
#include "lbfgs.h"
#include "OptPP.h"
#include "logLike_gauss.h"
#include "logLike_ptsrc.h"
#include "MyFun.h"
#include "PowerLaw.h"
#include "Gaussian.h"
#include "AbsEdge.h"
#include "Rosen.h"

using namespace Likelihood;   // for testing purposes only

void read_SC_Response_data();
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
void test_logLike_ptsrc();
void fit_3C279();
void fit_anti_center();
void test_CompositeFunction();
void test_OptPP();

std::string test_path;

int main(){
   read_SC_Response_data();
//     test_Parameter_class();
//     test_Function_class();
//     test_PowerLaw_class();
//     test_SourceModel_class();
//     test_Table_class();
//     test_Statistic_class();
//     test_Event_class();
//     test_PointSource_class();
//     test_Aeff_class();
//     test_Psf_class();
//     test_logLike_ptsrc();
   test_CompositeFunction();
//     test_OptPP();
//     fit_3C279();
//     fit_anti_center();
   return 0;
}

/*******************/
/* Test OptPP code */
/*******************/
void test_OptPP() {
// Test the OptPP code using my standard 2D Rosenbrock test with
// bounds constraints
   Rosen my_rosen;

   std::vector<Parameter> params;
   my_rosen.getParams(params);
   params[0].setValue(2.);
   params[0].setBounds(5./4., 10);
   params[1].setValue(2.);
   params[1].setBounds(-4., 10);
   my_rosen.setParams(params);

#ifdef HAVE_OPTIMIZERS

   int verbose = 3;

// try lbfgs_bcm method first   
//   lbfgs my_lbfgsObj(my_rosen);
//   my_lbfgsObj.find_min(verbose);

// now restart and try OptPP
   params[0].setValue(2.);
   params[1].setValue(2.);
   my_rosen.setParams(params);
   OptPP my_OptppObj(my_rosen);
   my_OptppObj.find_min(verbose);

#endif  //HAVE_OPTIMIZERS
   
   my_rosen.getParams(params);
   for (unsigned int i = 0; i < params.size(); i++) 
      std::cout << params[i].getName() << ": "
                << params[i].getValue() << std::endl;
}
// test_OptPP

void test_CompositeFunction() {

   Gaussian Carl(1., 1., 0.1);
   Gaussian Fredrich(0.5, 2., 0.3);

   SumFunction gauss_sum(Carl, Fredrich);

   Gaussian Johann(.75, 3., 0.5);

   SumFunction another_gauss_sum(Johann, gauss_sum);

   int nmax = 100;
   double xstep = 5./(nmax-1);
//     for (int i = 0; i < nmax; i++) {
//        double x = xstep*i;
//        dArg xarg(x);
//        std::cout << x << "  " << another_gauss_sum(xarg) << std::endl;
//     }

// add a PowerLaw and Gaussian together 

   PowerLaw pl_continuum(1., -2., 1.);
   Gaussian emission_line(10., 1., 0.1);

   SumFunction spectrum(pl_continuum, emission_line);

// multiply by an absorption edge

   AbsEdge edge(2., 3.);
   ProductFunction absorbed_spec(edge, spectrum);

   double xmin = 0.03;
   double xmax = 30.;
   xstep = log(xmax/xmin)/(nmax-1);
   for (int i = 0; i < nmax; i++) {
      double x = xmin*exp(xstep*i);
      dArg xarg(x);
      std::cout << x << "  " << absorbed_spec(xarg) << std::endl;
   }
   
// check the derivatives of AbsEdge
   std::vector<double> params;
   edge.getParamValues(params);
   dArg xarg(4.);
   double edge_value = edge(xarg);
   std::vector<double> derivs;
   edge.getDerivs(xarg, derivs);

   double eps = 1e-7;
   for (unsigned int i = 0; i < params.size(); i++) {
      std::vector<double> new_params = params;
      double dparam = eps*new_params[i];
      new_params[i] += dparam;
      edge.setParamValues(new_params);
      double new_edge_value = edge(xarg);
      std::cerr << derivs[i] << "  "
                << (new_edge_value - edge_value)/dparam
                << std::endl;
   }
}

void fit_anti_center() {

// set the ROI cuts
   RoiCuts::setCuts(83.57, 22.01, 20.);

/* read in the spacecraft data */
   std::string sc_file = test_path + "Data/anti_center_sc_0000";
   int sc_hdu = 2;
   ScData::readData(sc_file, sc_hdu);

// add the sources to the model

   logLike_ptsrc logLike;

// Crab Pulsar
   double ra = 83.57;
   double dec = 22.01;
   PointSource Crab(ra, dec);

//   double Prefactor = 6.46e-4;
//   double Prefactor = 6.46e1;
   double Prefactor = 27.;
   double Gamma = -2.19;
   double Escale = 100.;

   PowerLaw Crab_pl(Prefactor, Gamma, Escale);

//set limits on index
   Parameter *indexParam = Crab_pl.getParam("Index");
   indexParam->setBounds(-3.5, -1.);  
// set limits on normalization
   Parameter *prefactorParam = Crab_pl.getParam("Prefactor");
//   prefactorParam->setBounds(1e-8, 1e-2);
   prefactorParam->setBounds(1e-3, 1e3);
   prefactorParam->setScale(1e-9);

   Crab.setSpectrum(&Crab_pl);
   Crab.setName("Crab Pulsar");
   logLike.addSource(&Crab);

// Geminga Pulsar
   ra = 98.49;
   dec = 17.86;
   PointSource Geminga(ra, dec);

//   Prefactor = 4.866e-5;
   Prefactor = 23.29;
   Gamma = -1.66;
   Escale = 100.;

   PowerLaw Geminga_pl(Prefactor, Gamma, Escale);
//set limits on index
   indexParam = Geminga_pl.getParam("Index");
   indexParam->setBounds(-3.5, -1.);  
// set limits on normalization
   prefactorParam = Geminga_pl.getParam("Prefactor");
//   prefactorParam->setBounds(1e-8, 1e-2);  
   prefactorParam->setBounds(1e-3, 1e3);  
   prefactorParam->setScale(1e-9);

   Geminga.setSpectrum(&Geminga_pl);
   Geminga.setName("Geminga Pulsar");
   logLike.addSource(&Geminga);

// PKS 0528+134
   ra = 82.74;
   dec = 13.38;
   PointSource _0528(ra, dec);

//   Prefactor = 1.135e-3;
   Prefactor = 13.6;
   Gamma = -2.46;
   Escale = 100.;

   PowerLaw _0528_pl(Prefactor, Gamma, Escale);
//set limits on index
   indexParam = _0528_pl.getParam("Index");
   indexParam->setBounds(-3.5, -1.);  
// set limits on normalization
   prefactorParam = _0528_pl.getParam("Prefactor");
//   prefactorParam->setBounds(1e-8, 1e-2);
   prefactorParam->setBounds(1e-3, 1e3);
   prefactorParam->setScale(1e-9);

   _0528.setSpectrum(&_0528_pl);
   _0528.setName("PKS 0528+134");
   logLike.addSource(&_0528);

// read in the event data
   std::string event_file = test_path + "Data/anti_center_0000";
   logLike.getEvents(event_file, 2);

// some derivative tests
   std::vector<double> derivs;
   logLike.getFreeDerivs(derivs);
   
   std::vector<double> params;
   logLike.getFreeParamValues(params);
   
   double statval0 = logLike(params);
   
   for (unsigned int i = 0; i < params.size(); i++) {
      std::vector<double> new_params = params;
      double dparam = params[0]*1e-7;
      new_params[i] += dparam;
      double statval = logLike(new_params);
      std::cout << derivs[i] << "  " 
                << (statval - statval0)/dparam << std::endl;
   }

#ifdef HAVE_OPTIMIZERS

// do the fit
//   lbfgs myOptimizer(logLike);
   OptPP myOptimizer(logLike);

   int verbose = 3;
   myOptimizer.find_min(verbose);

#endif  //HAVE_OPTIMIZERS

   std::vector<Parameter> parameters;
   logLike.getParams(parameters);

   for (unsigned int i = 0; i < parameters.size(); i++)
      std::cout << parameters[i].getName() << ": "
                << parameters[i].getValue() << std::endl;
}

void fit_3C279() {

   double ra = 193.98;
   double dec = -5.82;
   
// set the ROI cuts
   RoiCuts::setCuts(ra, dec, 30.);
   
/* read in the spacecraft data */
   std::string sc_file = test_path + "Data/one_src_sc_0000";
   int sc_hdu = 2;
   ScData::readData(sc_file, sc_hdu);
   
   logLike_ptsrc logLike;

   PointSource _3c279(ra, dec);
   
//    double Prefactor = 5.92e-5;
//    double Gamma = -1.96;
//    double Escale = 1;

   double Prefactor = 4.;
   double Gamma = -2.1;
   double Escale = 100.;
   
   PowerLaw pl(Prefactor, Gamma, Escale);

//set limits on index
   Parameter *indexParam = pl.getParam("Index");
   indexParam->setBounds(-3.5, -1.);  

// set limits on normalization
   Parameter *prefactorParam = pl.getParam("Prefactor");
//   prefactorParam->setBounds(1e-8, 1e-2);  
   prefactorParam->setBounds(1e-3, 1e3);  
   prefactorParam->setScale(1e-9);

   _3c279.setSpectrum(&pl);

   _3c279.setName("3C 279");

   logLike.addSource(&_3c279);

// read in the data
   std::string event_file = test_path + "Data/one_src_0000";
   logLike.getEvents(event_file, 2);

 // some derivative tests
   std::vector<double> derivs;
   logLike.getFreeDerivs(derivs);
   
   std::vector<double> params;
   logLike.getFreeParamValues(params);
   
   double statval0 = logLike(params);
   
   for (unsigned int i = 0; i < params.size(); i++) {
      std::vector<double> new_params = params;
      double dparam = params[0]*1e-7;
      new_params[i] += dparam;
      double statval = logLike(new_params);
      std::cout << derivs[i] << "  " 
                << (statval - statval0)/dparam << std::endl;
   }

#ifdef HAVE_OPTIMIZERS

// do the fit using lbfgs_bcm
//   lbfgs myOptimizer(logLike);

// do the fit using OptPP
   OptPP myOptimizer(logLike);
   
   int verbose = 3;
   myOptimizer.find_min(verbose);

#endif  //HAVE_OPTIMIZERS
   
   std::vector<Parameter> parameters;
   logLike.getParams(parameters);
   
   for (unsigned int i = 0; i < parameters.size(); i++)
      std::cout << parameters[i].getName() << ": "
                << parameters[i].getValue() << std::endl;
}

/***********************/
/* logLike_ptsrc tests */
/***********************/
void test_logLike_ptsrc() {

   double ra = 193.98;
   double dec = -5.82;

// set the ROI cuts
   RoiCuts::setCuts(ra, dec, 30.);

/* read in the spacecraft data */
   std::string sc_file = test_path + "Data/one_src_sc_0000";
   int sc_hdu = 2;
   ScData::readData(sc_file, sc_hdu);

   logLike_ptsrc logLike;

   PointSource _3c279(ra, dec);

   double Prefactor = 7.;
   double Gamma = -1.96;
   double Escale = 100.;

   PowerLaw pl(Prefactor, Gamma, Escale);
   Parameter *prefactorParam = pl.getParam("Prefactor");
   prefactorParam->setScale(1e-9);

   _3c279.setSpectrum(&pl);

   _3c279.setName("3C 279");

   logLike.addSource(&_3c279);

   std::string event_file = test_path + "Data/one_src_0000";
   logLike.getEvents(event_file, 2);

   report_SrcModel_values(logLike);

   std::vector<Parameter> params;
   logLike.getFreeParams(params);

// vary over the Prefactor of the Spectrum and compute the
// log-likelihood at each step

   double xmin = 3.;
   double xmax = 10.;
   int nx = 20;
   double xstep = log(xmax/xmin)/(nx-1.);

   unsigned int j;
   for (int i = 0; i < nx; i++) {
      for (j = 0; j < params.size(); j++) {
         if (params[j].getName() == "Prefactor") {
            params[j].setValue(xmin*exp(xstep*i));
            break;
         }
      }
      std::vector<double> param_vals;
      for (unsigned int k = 0; k < params.size(); k++)
         param_vals.push_back(params[k].getValue());
      std::cout << params[j].getValue() << "  " 
                << logLike(param_vals) << std::endl;
   }
}
// logLike_ptsrc tests

/********************/
/* Psf class tests */
/********************/
void test_Psf_class() {
   Psf *psf = Psf::instance();

   int nenergy = 10;
   double emin = 30;
   double emax = 3e4;
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

   int nenergy = 10;
   double emin = 30;
   double emax = 3e4;
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

// set the ROI cuts
   RoiCuts::setCuts(193.98, -5.82, 30.);

/* read in the spacecraft data */
   std::string sc_file = test_path + "Data/one_src_sc_0000";
   int sc_hdu = 2;
   ScData::readData(sc_file, sc_hdu);

/* put this source at the center of the extraction region for
   Data/one_src_0000 (ra, dec) = (193.98, -5.82) */
   
   PointSource my_ptsrc(193.98, -5.82);

   my_ptsrc.setName("3C_279");

/* create a PowerLaw object for the source spectrum */

   PowerLaw source_pl(1., -2.1, 100.); // f(E) = (E/0.1 GeV)^(-2.1)
   my_ptsrc.setSpectrum(&source_pl);

/* compute the source flux at various energies and aspect angles */

   int nenergy = 10;
   double emin = 30.;
   double emax = 3e4;
   double estep = log(emax/emin)/(nenergy - 1);

   std::cout << "\nThe on-source counts spectrum of "
             << my_ptsrc.getName() 
             << " as a function of time: " 
             << std::endl;

   double photon_flux;

   double maxtime = 95.*60.;
   double tstep = 1e2;
   for (double time = 0; time < maxtime; time += tstep) {
      std::cout << "Time = " << time << std::endl;
      for (int i = 0; i < nenergy; i++) { // loop over energies
         double energy = emin*exp(estep*i);
         if (i == 0) photon_flux = 
                        my_ptsrc.fluxDensity(energy, time, my_ptsrc.getDir());
         if (photon_flux > 0) {
            std::cout << energy << "  ";
            std::cout << my_ptsrc.fluxDensity(energy, time, my_ptsrc.getDir())
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

// read in the event data, then stuff them into an Event class vector

   Table evt_table;

/* read in EVENT file */
   std::string event_file = test_path + "Data/one_src_0000";

   evt_table.add_columns("RA DEC energy time SC_x SC_y SC_z zenith_angle");
   evt_table.read_FITS_table(event_file, 2);

   typedef std::pair<long, double*> tableColumn;
   tableColumn ra(evt_table[0].dim, evt_table[0].val);
   tableColumn dec(evt_table[1].dim, evt_table[1].val);
   tableColumn energy(evt_table[2].dim, evt_table[2].val);
   tableColumn time(evt_table[3].dim, evt_table[3].val);
   tableColumn sc_x(evt_table[4].dim, evt_table[4].val);
   tableColumn sc_y(evt_table[5].dim, evt_table[5].val);
   tableColumn sc_z(evt_table[6].dim, evt_table[6].val);
   tableColumn zenangle(evt_table[7].dim, evt_table[7].val);

   std::vector<Event> my_events;

// get pointer to RoiCuts
   RoiCuts *roi_cuts = RoiCuts::instance();
   RoiCuts::setCuts(193.98, -5.82, 30);

   unsigned int nReject = 0;

   for (int i = 0; i < ra.first; i++) {
// compute sc_ra and sc_dec from direction cosines (it would be nice
// if SkyDir could do this for me....)
      double sc_ra = atan2(sc_y.second[i], sc_x.second[i])*180./M_PI;
      double sc_dec = asin(sc_z.second[i])*180./M_PI;

      Event thisEvent(ra.second[i], dec.second[i], energy.second[i],
                      time.second[i], sc_ra, sc_dec,
                      cos(zenangle.second[i]*M_PI/180.));
      if (roi_cuts->accept(thisEvent)) {
         my_events.push_back(thisEvent);
      } else {
         nReject++;
      }
   }

   std::cout << "Out of " << ra.first << " events, "
             << my_events.size() << " were accepted, and "
             << nReject << " were rejected.\n" << std::endl;

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
   std::cout << std::endl;
} // Event class tests

/*************************/
/* Statistic class tests */
/*************************/
void test_Statistic_class() {
   logLike_gauss logLike;

/* Try to analyze data from a 1D Gaussian distribution using the new
   Function and SourceModel interfaces (i.e., with iterators and
   Sources) */

/* for this test, the position is arbitrary --- use default (0, 0) --- */
   PointSource ptsrc;

/* but the "Spectrum" is not... */
   Gaussian gauss(1e4, 20., 5.);
   ptsrc.setSpectrum(&gauss);
   ptsrc.setName("Fredrich");

/* check the integral of this Gaussian */
//   std::cout << gauss.integral(Arg(-1e3), Arg(1e3)) << std::endl;

   logLike.addSource(&ptsrc);
   report_SrcModel_values(logLike);

/* read in "event" file */
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

/* evaluate the integral over x */

   dArg gmin(-1e3);
   dArg gmax(1e3);
   std::cout << gauss.integral(gmin, gmax)
             << std::endl << std::endl;

/* Compute the log-likelihood of this model */

// first get the vector of parameters
   std::vector<Parameter> params;
   logLike.getParams(params);

// vary over the first parameter and compute the log-likelihood at
// each step

   double xmin = 5e3;
   double xmax = 1.5e4;
   int nx = 20;
   double xstep = log(xmax/xmin)/(nx-1.);

   unsigned int j;
   for (int i = 0; i < nx; i++) {
      for (j = 0; j < params.size(); j++) {
         if (params[j].getName() == "Prefactor") {
            params[j].setValue(xmin*exp(xstep*i));
            break;
         }
      }
      std::vector<double> param_vals;
      for (unsigned int k = 0; k < params.size(); k++)
         param_vals.push_back(params[k].getValue());
      std::cout << params[j].getValue() << "  " 
                << logLike(param_vals) << std::endl;
   }
} // Statistic class tests

/*********************/
/* Table class tests */
/*********************/
void test_Table_class() {

/* read in PSF parameters */
   std::string psf_file = test_path + "CALDB/psf_lat.fits";

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
   std::string event_file = test_path + "Data/one_src_0000";

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
   std::string sc_file = test_path + "Data/one_src_sc_0000";

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
   
/* instantiate some point sources */

   PointSource _3c279;
   _3c279.setDir(193.98, -5.82);
   _3c279.setSpectrum(new PowerLaw(74.2, -1.96, 0.1));
   _3c279.setName("3C 279");

   PointSource _3c273;
   _3c273.setDir(187.25, 2.17);
   _3c273.setSpectrum(new PowerLaw(15.4, -2.58, 0.1));
   _3c273.setName("3C 273");

   PointSource Crab;
   Crab.setDir(83.57, 22.01);
   Crab.setSpectrum(new PowerLaw(226.2, -2.19, 0.1));
   Crab.setName("Crab Pulsar");

   PointSource Vela;
   Vela.setDir(128.73, -45.20);
   Vela.setSpectrum(new PowerLaw(834.3, -1.69, 0.1));
   Vela.setName("Vela Pulsar");

/* add these guys to the SouceModel */

   SrcModel.addSource(&_3c279);
   SrcModel.addSource(&_3c273);
   SrcModel.addSource(&Crab);
   SrcModel.addSource(&Vela);
   report_SrcModel_values(SrcModel);

/* delete a Source */
   SrcModel.deleteSource("Crab Pulsar");
   report_SrcModel_values(SrcModel);

/* add another function */
   PointSource Geminga;
   Geminga.setDir(98.49, 17.86);
   Geminga.setSpectrum(new PowerLaw(352.9, -1.66, 0.1));
   Geminga.setName("Geminga");
   SrcModel.addSource(&Geminga);

/* make its scale factor free (this could be made easier, e.g., by
   giving direct parameter access from PointSource) */

   Parameter param = *(SrcModel.getParam(std::string("Scale"), 
                                         std::string("Spectrum"), 
                                         std::string("Geminga")));
   param.setFree(true);
   SrcModel.setParam(param, std::string("Spectrum"), std::string("Geminga"));
   report_SrcModel_values(SrcModel);

/* derivative tests */

   dArg x(20.);
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

/* derivative tests for free parameters only */

   SrcModel.getFreeParamValues(params_save);
   SrcModel.getFreeDerivs(x, sm_derivs);
   params = params_save;

   std::cout << "SrcModel Derivatives for Free Parameters: " << std::endl;
   for (unsigned int i = 0; i < sm_derivs.size(); i++) {
      std::cout << sm_derivs[i] << ":  ";

// compute the numerical derivative wrt this parameter
      double delta_param = fabs(params_save[i]/1e7);
      double num_deriv = SrcModel(x);
      params[i] += delta_param;
      SrcModel.setFreeParamValues(params);
      num_deriv -= SrcModel(x);
      num_deriv /= delta_param;
      std::cout << -num_deriv << std::endl;

// reset the parameters for next time around
     SrcModel.setFreeParamValues(params_save);
     params = params_save;
  }
  std::cout << std::endl;

} // Source Model Class tests

void report_SrcModel_values(const SourceModel &SrcModel) {

/* retrieve Free parameter values */
  std::vector<double> paramValues;
  SrcModel.getFreeParamValues(paramValues);
  for (unsigned int i = 0; i < paramValues.size(); i++) {
     std::cout << paramValues[i] << "  ";
  }
  std::cout << std::endl;

/* retrieve all parameter values */
  SrcModel.getParamValues(paramValues);
  for (unsigned int i = 0; i < paramValues.size(); i++) {
     std::cout << paramValues[i] << "  ";
  }
  std::cout << std::endl;
   
/* retrieve Source names */
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
   for (double xx = 1.05; xx < 1e3; xx *= xx) {
      dArg xarg(xx);
      std::cout << xx << "   " << pl(xarg) << std::endl;
   }

/* get the derivatives and compare to numerical estimates */
   dArg x(10.);

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

   for (double x = 1.05; x < 1e3; x *= x) {
      dArg xarg(x);
      std::cout << x << "   " 
                << pl(xarg) << "  " 
                << pl2(xarg) << std::endl;
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

   std::vector<std::string> my_paramNames;
   f.getParamNames(my_paramNames);
      
   std::cout << "Here they are: " << std::endl;
   for (unsigned int i = 0; i < f.getNumParams(); i++) {
      std::cout << my_paramNames[i] << ":  " 
                << f.getParamValue(my_paramNames[i]) 
                << std::endl;
   }
   dArg x(3); 
   std::cout << "f(3) = " << f(x) << std::endl;

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
   f.setParam(std::string("Ruthie"), 10.);

/* attempt to change the value of a non-existing parameter */
   f.setParam(std::string("Oscar"), 5.);
      
   std::cout << "The current set of values: " << std::endl;
   std::vector<double> my_params;
   f.getParamValues(my_params);
   for (unsigned int i = 0; i < my_params.size(); i++) {
      std::cout << my_params[i] << " ";
   }
   x = dArg(2);
   std::cout << " f(2) = " << f(x) << std::endl;

/* get derivatives wrt parameters */
   std::cout << "getting derivatives one-by-one:" << std::endl;
   for (unsigned int i = 0; i < my_paramNames.size(); i++) {
      std::cout << my_paramNames[i] << ":  "
                << f.derivByParam(x, my_paramNames[i]) << std::endl;
   }

   std::cout << "all derivatives in one shot:" << std::endl;
   std::vector<double> my_derivs;
   f.getDerivs(x, my_derivs);
   for (unsigned int i = 0; i < my_paramNames.size(); i++) {
      std::cout << f.derivByParam(x, my_paramNames[i]) << "  ";
   }
   std::cout << std::endl;

/* test of pointers to Parameter */
   Parameter *ptrP = f.getParam(std::string("Mary"));
   if (ptrP != NULL) {
      std::cout << ptrP->getName() << ":  " 
                << ptrP->getValue() << std::endl;
   }
      
   ptrP = f.getParam(std::string("Joan"));
   if (ptrP != NULL) {
      std::cout << ptrP->getName() << ":  " 
                << ptrP->getValue() << "\n" << std::endl;
   }
   std::cout << std::endl;

/* test the Function copy constructor */

   MyFun f2 = f;

   for (double x = 0; x < 100.; x += 5.) {
      dArg xarg(x);
      std::cout << x << "  " 
                << f(xarg) << "  " 
                << f2(xarg) << "  " 
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


void read_SC_Response_data() {

/* get root path to test data */   
   const char * root = ::getenv("LIKELIHOODROOT");
   if (!root) {  //use relative path from cmt directory
      test_path = "../src/test/";
   } else {
      test_path = std::string(root) + "/src/test/";
   }

/* instantiate the Psf and read in its data */
   Psf * psf = Psf::instance();
   std::string psf_file = test_path + "CALDB/psf_lat.fits";
   psf->readPsfData(psf_file, Response::Combined);

/* instantiate the Aeff and read in its data */
   Aeff * aeff = Aeff::instance();
   std::string aeff_file = test_path + "CALDB/aeff_lat.fits";
   aeff->readAeffData(aeff_file, Response::Combined);
}
