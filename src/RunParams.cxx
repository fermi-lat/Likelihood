/**
 * @file RunParams.cxx
 * @brief Implementation for the hoops wrapper class to retrieve
 * command-line parameters.
 * @author J. Chiang
 *
 * $Header$
 */

#include "Likelihood/RunParams.h"

namespace Likelihood {

RunParams::RunParams(int iargc, char* argv[]) {

   hoops::IParFile * pf = hoops::PILParFileFactory().NewIParFile(argv[0]);
   pf->Load();

   m_prompter = hoops::PILParPromptFactory().NewIParPrompt(iargc, argv);
   m_prompter->Prompt();

   pf->Group() = m_prompter->Group();
   pf->Save();

//   delete pf;
}

RunParams::~RunParams() {
   delete m_prompter;
//   delete pf;
}

} // namespace Likelihood
