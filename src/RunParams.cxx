/**
 * @file RunParams.cxx
 * @brief Implementation for the hoops wrapper class to retrieve
 * command-line parameters.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/RunParams.cxx,v 1.3 2003/11/10 23:06:20 jchiang Exp $
 */

#include <fstream>

#include "facilities/Util.h"

#include "hoops/hoops_exception.h"

#include "Likelihood/RunParams.h"

namespace Likelihood {

RunParams::RunParams(int iargc, char* argv[]) {

   hoops::IParFile * pf = hoops::PILParFileFactory().NewIParFile(argv[0]);
   try {
      pf->Load();
   } catch (hoops::Hexception &eObj) {
      if (eObj.Code() == -3003) {
         std::cout << "Likelihood::RunParams: .par file " << argv[0]
                   << "is not found.  Check your PFILES directory." 
                   << std::endl;
         assert(eObj.Code() != -3003);
      }
   }

   m_prompter = hoops::PILParPromptFactory().NewIParPrompt(iargc, argv);
   m_prompter->Prompt();

   pf->Group() = m_prompter->Group();
   pf->Save();

   delete pf;
}

RunParams::~RunParams() {
   delete m_prompter;
}

void RunParams::resolve_fits_files(std::string filename, 
                                   std::vector<std::string> &files) {

   facilities::Util::expandEnvVar(&filename);
   files.clear();

// Read the first line of the file and see if the first 6 characters
// are "SIMPLE".  If so, then we assume it's a FITS file.
   std::ifstream file(filename.c_str());
   std::string firstLine;
   std::getline(file, firstLine, '\n');
   if (firstLine.find("SIMPLE") == 0) {
// This is a FITS file. Return that as the sole element in the files
// vector.
      files.push_back(filename);
      return;
   } else {
// filename contains a list of fits files.
      readLines(filename, files);
      return;
   }
}

void RunParams::readLines(std::string inputFile, 
                          std::vector<std::string> &lines) {

   facilities::Util::expandEnvVar(&inputFile);

   std::ifstream file(inputFile.c_str());
   lines.clear();
   std::string line;
   while (std::getline(file, line, '\n')) {
      lines.push_back(line);
   }
}

} // namespace Likelihood
