/**
 * @file BinnedLikelihood.cxx
 *
 */

#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/SourceModel.h"

namespace Likelihood {

BinnedLikelihood::BinnedLikelihood(const CountsMap & dataMap)
   : m_dataMap(dataMap), m_modelMap(0) {
   getPixels(dataMap, m_pixelDirs, m_pixelSolidAngles);
   dataMap.getAxisVector(2, m_energies);
}

double BinnedLikelihood::value(optimizers::Arg &dummy) const {

   (void)(dummy);

   std::vector<double> model;
   computeModelMap(m_pixelDirs, m_pixelSolidAngles, m_energies, model);
   
   const std::vector<double> & data = m_dataMap.data();
   double my_value(0);

   for (unsigned int i = 0; i < data.size(); i++) {
      if (model[i] > 0) {
         my_value += data[i]*log(model[i]) - model[i];
      }
   }

   return my_value;
}



} // namespace Likelihood
