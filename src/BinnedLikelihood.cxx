/**
 * @file BinnedLikelihood.cxx
 *
 */

#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/SourceModel.h"

namespace Likelihood {

BinnedLikelihood::BinnedLikelihood(const CountsMap & dataMap)
   : m_dataMap(dataMap), m_modelIsCurrent(false) {
   getPixels(dataMap, m_pixels);
   dataMap.getAxisVector(2, m_energies);
}

double BinnedLikelihood::value(optimizers::Arg &dummy) const {

   (void)(dummy);

   computeModelMap(m_pixels, m_energies, m_model);
   m_modelIsCurrent = true;
   
   const std::vector<double> & data = m_dataMap.data();
   double my_value(0);

   for (unsigned int j = 0; j < data.size(); j++) {
      if (m_model[j] > 0) {
         my_value += data[j]*log(m_model[j]) - m_model[j];
      }
   }
   
   return my_value;
}

std::vector<double>::const_iterator 
BinnedLikelihood::setParamValues_(std::vector<double>::const_iterator it) {
   m_modelIsCurrent = false;
   return SourceModel::setParamValues_(it);
}

std::vector<double>::const_iterator 
BinnedLikelihood::setFreeParamValues_(std::vector<double>::const_iterator it) {
   m_modelIsCurrent = false;
   return SourceModel::setFreeParamValues_(it);
}
   
void BinnedLikelihood::getFreeDerivs(std::vector<double> & derivs) const {
   int nparams(getNumFreeParams());
   derivs.resize(nparams, 0);
   if (!m_modelIsCurrent) {
      computeModelMap(m_pixels, m_energies, m_model);
      m_modelIsCurrent = true;
   }
   int indx;
   const std::vector<double> & data = m_dataMap.data();
   for (unsigned int j = 0; j < m_pixels.size(); j++) {
      for (unsigned int k = 0; j < m_energies.size()-1; k++) {
         indx = k*m_pixels.size() + j;
         std::vector<double> my_derivs;
         m_pixels.getFreeDerivs(m_energies[k], m_energies[k+1], my_derivs);
         for (int i = 0; i < nparams; i++) {
            derivs[i] += (data[indx]/m_model[indx] - 1.)*my_derivs[i];
         }
      }               
   }
}

} // namespace Likelihood
