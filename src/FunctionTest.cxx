/**
 * @file FunctionTest.cxx
 * @brief Unit test implementation for Function class hierarchy
 * @author J. Chiang
 * $Header$
 */

#include "Likelihood/FunctionTest.h"

namespace Likelihood {

void FunctionTest::parameters(const std::vector<Parameter> &params) 
   throw(LikelihoodException) {

   std::vector<std::string> paramNames;
   m_func->getParamNames(paramNames);
   unsigned int numParams = m_func->getNumParams();
   assert(numParams == paramNames.size());

// set the Parameter values
   for (unsigned int i = 0; i < numParams && i < params.size(); i++) {
      m_func->setParam(paramNames[i], params[i].getValue());
      m_func->setParamBounds(paramNames[i], params[i].getBounds().first,
                             params[i].getBounds().second);
   }

// and test in groups and individually
   std::vector<double> paramValues;
   m_func->getParamValues(paramValues);
   for (unsigned int i = 0; i < paramValues.size(); i++) {
      assert(paramValues[i] == params[i].getValue());
      assert(m_func->getParamValue(paramNames[i]) == params[i].getValue());
   }

// set Parameters using Parameter class
   std::vector<Parameter> my_params;
   for (unsigned int i = 0; i < numParams && i < params.size(); i++) {
      my_params.push_back(params[i]);
// ensure names match up 
      my_params[i].setName(paramNames[i]);
   }

// set Parameters as a group
   m_func->setParams(my_params);

// test values again
   m_func->getParamValues(paramValues);
   for (unsigned int i = 0; i < paramValues.size(); i++) {
      assert(paramValues[i] == params[i].getValue());
      assert(m_func->getParamValue(paramNames[i]) == params[i].getValue());
   }

// reset Parameters individually with different values
   for (unsigned int i = 0; i < my_params.size(); i++) {
      my_params[i].setValue(params[i].getValue()/10.);
      m_func->setParam(my_params[i]);
   }

// test yet again
   m_func->getParamValues(paramValues);
   for (unsigned int i = 0; i < paramValues.size(); i++) {
      double value = params[i].getValue()/10.;
      assert(paramValues[i] == value);
      assert(m_func->getParamValue(paramNames[i]) == value);
   }

// restore the original Parameters
   m_func->setParams(m_originalParameters);
}

void FunctionTest::freeParameters(const std::vector<Parameter> &params) 
   throw(LikelihoodException) {

// Create a local copy of the test Parameters
   std::vector<Parameter> my_params = params;

// Make every other parameter free, the others fixed
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if ((i % 2) == 0) {
         my_params[i].setFree(true);
      } else {
         my_params[i].setFree(false);
      }
   }

// Sync up the parameter names
   std::vector<std::string> paramNames;
   m_func->getParamNames(paramNames);
   my_params.resize(paramNames.size());
   for (unsigned int i = 0; i < my_params.size(); i++) {
      my_params[i].setName(paramNames[i]);
   }

// set the all the Parameters
   m_func->setParams(my_params);

// retrieve the free ones
   std::vector<std::string> freeParamNames;
   m_func->getFreeParamNames(freeParamNames);
   std::vector<double> freeParamValues;
   m_func->getFreeParamValues(freeParamValues);

// and check against my_params
   int j = 0;
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if (my_params[i].isFree()) {
         assert(my_params[i].getName() == freeParamNames[j]);
         assert(my_params[i].getValue() == freeParamValues[j]);
         j++;
      }
   }

// reset the free Parameter values as a group
   for (unsigned int i = 0; i < freeParamValues.size(); i++) {
      freeParamValues[i] /= 5.;
   }
   m_func->setFreeParamValues(freeParamValues);

// check
   j = 0;
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if (my_params[i].isFree()) {
         assert(my_params[i].getValue()/5. 
                == m_func->getParamValue(freeParamNames[j]));
         j++;
      }
   }

// restore the original Parameters
   m_func->setParams(m_originalParameters);
}

void FunctionTest::funcEvaluations(const std::vector<Arg*> &arguments,
                                   const std::vector<double> &returnValues) 
   throw(LikelihoodException) {
   if (arguments.size() != returnValues.size()) {
      throw LikelihoodException(
         std::string("FunctionTest::funcEvaluations:\n") + 
         std::string("Sizes of arguments and (expected) returnValues") + 
         std::string(" do not match."));
   }

//    std::vector<double> params;
//    m_func->getParamValues(params);
//    for (unsigned int i = 0; i < params.size(); i++) 
//       std::cout << params[i] << std::endl;
   for (unsigned int i = 0; i < arguments.size(); i++) {
      double val = m_func->value(*arguments[i]);
      assert(val == returnValues[i]);
   }
}

void FunctionTest::derivatives(const std::vector<Arg*> &arguments,
                               double eps) throw(LikelihoodException) {
// loop over all Parameters
   for (unsigned int k = 0; k < arguments.size(); k++) {
      Arg *my_arg = arguments[k];

// Check the derivatives wrt all Parameters against numerical estimates 
      std::vector<double> derivs;
      m_func->getDerivs(*my_arg, derivs);
      std::vector<double> params;
      m_func->getParamValues(params);
      double f0 = m_func->value(*my_arg);
      for (unsigned int i = 0; i < params.size(); i++) {
         std::vector<double> new_params = params;
         double delta = new_params[i]*eps;
         new_params[i] += delta;
         m_func->setParamValues(new_params);
         double f1 = m_func->value(*my_arg);
         double my_deriv = (f1 - f0)/delta;
         assert(abs(my_deriv/derivs[i] - 1.) < eps*10.);
      }

// Check that the free derivatives are being accessed correctly.
      m_func->setParamValues(params);     // restore the Parameter values
      std::vector<Parameter> parameters;  // get all Parameter info
      m_func->getParams(parameters);
      std::vector<double> freeDerivs;     // and the free derivatives
      m_func->getFreeDerivs(*my_arg, freeDerivs);
      int j = 0;
      for (unsigned int i = 0; i < parameters.size(); i++) {
         if (parameters[i].isFree()) {
            assert(m_func->derivByParam(*my_arg, parameters[i].getName())
                   == freeDerivs[j]);
            j++;
         }
      }
   }

// restore the original Parameters
   m_func->setParams(m_originalParameters);
}

} //namespace Likelihood
