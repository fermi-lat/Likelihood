/**
 * @file RunParams.cxx
 * @brief Implementation for the hoops wrapper class to retrieve
 * command-line parameters.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/RunParams.cxx,v 1.1 2003/11/05 03:41:24 jchiang Exp $
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

   delete pf;
}

RunParams::~RunParams() {
   delete m_prompter;
}

} // namespace Likelihood
