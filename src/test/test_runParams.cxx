/**
 * @file test_RunParams.cxx
 * @brief Standalone application for testing RunParams/HOOPS
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/likelihood.cxx,v 1.17 2003/12/06 03:52:41 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <iostream>
#include <cstring>
#include <cmath>
#include <cassert>

#include "Likelihood/RunParams.h"
#include "hoops/hoops_exception.h"

using namespace Likelihood;

int main(int iargc, char* argv[]) {

// Read in the command-line parameters using HOOPS
   std::string filename("test_runParams.par");
   delete argv[0];
   argv[0] = strdup(filename.c_str());

   RunParams params(iargc, argv);

// Read a file.
   std::string inputFile;
   params.getParam("input file", inputFile);
   std::cout << inputFile << std::endl;

   double value;
   params.getParam("double value", value);
   std::cout << value << std::endl;
}
