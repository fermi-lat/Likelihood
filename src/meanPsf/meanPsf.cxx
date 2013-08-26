/**
 * @file meanPsf.cxx
 * @brief Application to compute a spectrally weighted Psf averaged over
 * an observaition.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/meanPsf/meanPsf.cxx,v 1.15 2010/06/16 22:49:54 jchiang Exp $
 */

#include <cstdio>
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

#include "dataSubselector/Cuts.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/MeanPsf.h"

using namespace Likelihood;

class meanPsf : public st_app::StApp {
public:
   meanPsf();
   virtual ~meanPsf() throw() {}
   virtual void run();
   virtual void banner() const;
private:
   AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;

   MeanPsf * m_meanPsf;

   std::vector<double> m_energies;
   std::vector<double> m_thetas;

   void computeEnergies();
   void computeThetas();
   void writeFitsFile();

   static std::string s_cvs_id;
};

st_app::StAppFactory<meanPsf> myAppFactory("gtpsf");

std::string meanPsf::s_cvs_id("$Name: Likelihood-18-00-04 $");

meanPsf::meanPsf() 
   : st_app::StApp(), m_helper(0),
     m_pars(st_app::StApp::getParGroup("gtpsf")),
     m_meanPsf(0) {
   setVersion(s_cvs_id);
}

void meanPsf::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void meanPsf::run() {
   m_pars.Prompt();
   m_pars.Save();
   m_helper = new AppHelpers(&m_pars, "none");
   std::string expcube_file = m_pars["expcube"];
   if (expcube_file == "none") {
      throw std::runtime_error("Please specify an exposure cube file.");
   }
   ExposureCube & expCube = 
      const_cast<ExposureCube &>(m_helper->observation().expCube());
   expCube.readExposureCube(expcube_file);
   
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
   long nee = m_pars["nenergies"];
   double estep = log(emax/emin)/(nee - 1.);
   m_energies.clear();
   m_energies.reserve(nee);
   for (int i = 0; i < nee; i++) {
      m_energies.push_back(emin*exp(estep*i));
   }
}

void meanPsf::computeThetas() {
   double theta_max = m_pars["thetamax"];
   int nth = m_pars["ntheta"];

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
   psf_table->appendField("Exposure", "1D");
   std::ostringstream format;
   format << m_thetas.size() << "D";
   psf_table->appendField("Psf", format.str());
   psf_table->setNumRecords(m_energies.size());
   std::vector<double> psf_values(m_thetas.size());
   tip::Table::Iterator row = psf_table->begin();
   tip::Table::Record & record = *row;

   const std::vector<double> & exposure = m_meanPsf->exposure();

   std::vector<double>::const_iterator energy = m_energies.begin();
   for (int k = 0; energy != m_energies.end(); ++energy, ++row, ++k) {
      record["Energy"].set(*energy);
      record["Exposure"].set(exposure.at(k));
      tip::Table::Vector<double> Psf = record["Psf"];
      std::vector<double>::const_iterator theta = m_thetas.begin();
      for (unsigned int i = 0; theta != m_thetas.end(); ++theta, ++i) {
         Psf[i] = (*m_meanPsf)(*energy, *theta, 0);
      }
   }

   // Write irfs version using DSS keywords.
   dataSubselector::Cuts my_cuts;
   my_cuts.addVersionCut("IRF_VERSION", m_helper->irfsName());
   my_cuts.writeDssKeywords(psf_table->getHeader());

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
