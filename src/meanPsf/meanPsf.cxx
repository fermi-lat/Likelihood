#include <cstdlib>

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/MeanPsf.h"
//#include "Likelihood/ResponseFunctions.h"

using namespace Likelihood;

class meanPsf : public st_app::StApp {
public:
   meanPsf();
   virtual ~meanPsf() throw() {}
   virtual void run();
private:
   AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;

   MeanPsf * m_meanPsf;

   std::vector<double> m_energies;
   std::vector<double> m_thetas;

   void computeEnergies();
   void computeThetas();
   void writeFitsFile();
};

st_app::StAppFactory<meanPsf> myAppFactory;

meanPsf::meanPsf() 
   : st_app::StApp(), m_helper(0),
     m_pars(st_app::StApp::getParGroup("gtpsf")),
     m_meanPsf(0) {
   try {
      m_pars.Prompt();
      m_pars.Save();
      m_helper = new AppHelpers(m_pars);
   } catch (std::exception & eObj) {
      std::cerr << eObj.what() << std::endl;
      std::exit(1);
   } catch (...) {
      std::cerr << "Caught unknown exception in meanPsf constructor." 
                << std::endl;
      std::exit(1);
   }
}

void meanPsf::run() {
   std::string expcube_file = m_pars["exposure_cube_file"];
   if (expcube_file == "none") {
      throw std::runtime_error("Please specify an exposure cube file.");
   }
   ExposureCube::readExposureCube(expcube_file);
   
   computeEnergies();
   computeThetas();
   double ra = m_pars["ra"];
   double dec = m_pars["dec"];
   m_meanPsf = new MeanPsf(ra, dec, m_energies, m_helper->observation());
   writeFitsFile();
   delete m_meanPsf;
}

void meanPsf::computeEnergies() {
   double emin = m_pars["emin"];
   double emax = m_pars["emax"];
   long nee = m_pars["numenergies"];
   double estep = log(emax/emin)/(nee - 1.);
   m_energies.clear();
   m_energies.reserve(nee);
   for (int i = 0; i < nee; i++) {
      m_energies.push_back(emin*exp(estep*i));
   }
}

void meanPsf::computeThetas() {
   double theta_max = m_pars["thetamax"];
   int nth = m_pars["numthetas"];

   m_thetas.clear();
   m_thetas.reserve(nth);
   for (int i = 0; i < nth; i++) {
      m_thetas.push_back((float(i*i)/float(nth*nth))*theta_max);
   }
}

void meanPsf::writeFitsFile() {
   std::string output_file = m_pars["outfile"];
   std::string extension = m_pars["outtable"];
   bool clobber = m_pars["clobber"];
   if (st_facilities::Util::fileExists(output_file) && !clobber) {
      throw std::runtime_error(output_file
                               + " file exists and clobber is set to no.");
   }
   std::remove(output_file.c_str());
   tip::IFileSvc::instance().appendTable(output_file, extension);
   std::auto_ptr<tip::Table> 
      psf_table(tip::IFileSvc::instance().editTable(output_file, extension));
   psf_table->appendField("Energy", "1D");
   std::ostringstream format;
   format << m_thetas.size() << "D";
   psf_table->appendField("Psf", format.str());
   psf_table->setNumRecords(m_energies.size());
   std::vector<double> psf_values(m_thetas.size());
   tip::Table::Iterator row = psf_table->begin();
   tip::Table::Record & record = *row;
   for (std::vector<double>::const_iterator energy = m_energies.begin();
        energy != m_energies.end(); ++energy, ++row) {
      record["Energy"].set(*energy);
      tip::Table::Vector<double> Psf = record["Psf"];
      std::vector<double>::const_iterator theta = m_thetas.begin();
      for (unsigned int i = 0; theta != m_thetas.end(); ++theta, ++i) {
         Psf[i] = (*m_meanPsf)(*energy, *theta, 0);
      }
   }
   extension = "THETA";
   tip::IFileSvc::instance().appendTable(output_file, extension);
   std::auto_ptr<tip::Table> 
      theta_table(tip::IFileSvc::instance().editTable(output_file, extension));
   theta_table->appendField("Theta", "1D");
   theta_table->setNumRecords(m_thetas.size());
   row = theta_table->begin();
   tip::Table::Record & theta_record = *row;
   for (std::vector<double>::const_iterator theta = m_thetas.begin();
        theta != m_thetas.end(); ++theta, ++row) {
      theta_record["Theta"].set(*theta);
   }
}
