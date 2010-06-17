/**
 * @file gtebl.cxx
 * @brief Return EBL optical depth as a function of energy for a given
 * redshift and EBL model.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/gtebl/gtebl.cxx,v 1.1 2009/08/30 00:21:22 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>

#include <fstream>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "eblAtten/EblAtten.h"

class Ebl : public st_app::StApp {

public:

   Ebl() : st_app::StApp(),
           m_pars(st_app::StApp::getParGroup("gtebl")) {
      try {
         setVersion(s_cvs_id);
      } catch (std::exception & eObj) {
         std::cerr << eObj.what() << std::endl;
         std::exit(1);
      } catch (...) {
         std::cerr << "Caught unknown exception in Ebl constructor." 
                   << std::endl;
         std::exit(1);
      }
   }

   virtual ~Ebl() throw() {
      try {
      } catch (std::exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }  
   }

   virtual void run();
   virtual void banner() const;

private:

   st_app::AppParGroup & m_pars;

   std::string par(const std::string & key) const;

   static std::string s_cvs_id;

};

std::string Ebl::s_cvs_id("$Name:  $");

st_app::StAppFactory<Ebl> myAppFactory("gtebl");

void Ebl::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

std::string Ebl::par(const std::string & key) const {
   std::string value = m_pars[key.c_str()];
   return value;
}

void Ebl::run() {
   m_pars.Prompt();
   m_pars.Save();

   st_stream::StreamFormatter formatter("gtebl", "", 2);

   int model_id = m_pars["model_id"];
   IRB::EblModel modelId = static_cast<IRB::EblModel>(model_id);

   double redshift = m_pars["redshift"];
   double emin = m_pars["emin"];
   double emax = m_pars["emax"];
   size_t nee = m_pars["nenergies"];

   double estep(std::log(emax/emin)/(nee - 1));

   std::string outfile = m_pars["outfile"];
   std::ofstream output(outfile.c_str());

   IRB::EblAtten tau(modelId);

   for (size_t k(0); k < nee; k++) {
      double energy(emin*std::exp(estep*k));
      output << energy << "  "
             << tau(energy, redshift) << std::endl;
   }
   output.close();
}
