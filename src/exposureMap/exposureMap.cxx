#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BinnedExposure.h"
#include "Likelihood/ExposureCube.h"
//#include "Likelihood/ResponseFunctions.h"

using namespace Likelihood;

class exposureMap : public st_app::StApp {
public:
   exposureMap();
   virtual ~exposureMap() throw() {}
   virtual void run();
private:
   AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;

   std::vector<double> m_energies;

   void computeEnergies();
};

st_app::StAppFactory<exposureMap> myAppFactory;

exposureMap::exposureMap() 
   : st_app::StApp(), m_helper(0),
     m_pars(st_app::StApp::getParGroup("exposureMap")) {
   try {
      m_pars.Prompt();
      m_pars.Save();
      m_helper = new AppHelpers(m_pars);
//      ResponseFunctions::setEdispFlag(m_pars["use_energy_dispersion"]);
   } catch (std::exception & eObj) {
      std::cerr << eObj.what() << std::endl;
      std::exit(1);
   } catch (...) {
      std::cerr << "Caught unknown exception in exposureMap constructor." 
                << std::endl;
      std::exit(1);
   }
}

void exposureMap::run() {
   std::string expcube_file = m_pars["exposure_cube_file"];
   if (expcube_file == "none") {
      throw std::runtime_error("Please specify an exposure cube file.");
   }
   ExposureCube::readExposureCube(expcube_file);
   std::string output_file = m_pars["output_file_name"];
   computeEnergies();
   BinnedExposure exposure(m_energies);
   exposure.writeOutput(output_file);
}

void exposureMap::computeEnergies() {
   double emin(20.);
   double emax(2e5);
   int nee(20);
   double estep = log(emax/emin)/(nee - 1.);
   m_energies.clear();
   m_energies.reserve(nee);
   for (int i = 0; i < nee; i++) {
      m_energies.push_back(emin*exp(estep*i));
   }
}
