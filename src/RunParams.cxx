/**
 * @file RunParams.cxx
 * @brief Implementation for the hoops wrapper class to retrieve
 * command-line parameters.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/RunParams.cxx,v 1.7 2004/01/06 00:10:28 jchiang Exp $
 */

#include <fstream>

#include "facilities/Util.h"

#include "hoops/hoops_exception.h"

#include "Likelihood/RunParams.h"

namespace Likelihood {

RunParams::RunParams(int iargc, char* argv[]) {

   errorCodes();

   hoops::IParFile * pf = hoops::PILParFileFactory().NewIParFile(argv[0]);
   try {
      pf->Load();
      m_prompter = hoops::PILParPromptFactory().NewIParPrompt(iargc, argv);
      m_prompter->Prompt();
      
      pf->Group() = m_prompter->Group();
      pf->Save();

      delete pf;
   } catch (hoops::Hexception &eObj) {
      if (m_pil_errors.count(eObj.Code())) {
         std::cout << "PIL exception: " << eObj.Code() << "\n"
                   << m_pil_errors[eObj.Code()];
         if (eObj.Code() == -3003) {
// This is obviously a fatal exception.
            std::cout << argv[0] << ".par file is not found. "
                      << "Check your PFILES directory." 
                      << std::endl;
            assert(eObj.Code() != -3003);
         }
      }
      throw;
   }
}

RunParams::~RunParams() {
   delete m_prompter;
}

void RunParams::errorCodes() {
   m_pil_errors[-3000] = std::string("Null pointer passed as an argument, ")
      + "and function does not allow it.";
   m_pil_errors[-3001] = "Bad argument passed.";
   m_pil_errors[-3002] = "Not enough memory.";
   m_pil_errors[-3003] = "Cannot open/create file.";
   m_pil_errors[-3004] = "Read from file failed.";
   m_pil_errors[-3005] = "Write to file failed.";
   m_pil_errors[-3006] = "Unexpected end of string.";
   m_pil_errors[-3007] = "Invalid name of parameter.";
   m_pil_errors[-3008] = "Invalid type of parameter in parameter file.";
   m_pil_errors[-3009] = "Invalid mode of parameter in parameter file.";
   m_pil_errors[-3010] = "Invalid line in parameter file encountered.";
   m_pil_errors[-3011] = "Feature not implemented.";
   m_pil_errors[-3012] = "File does not exist.";
   m_pil_errors[-3013] = "File exists.";
   m_pil_errors[-3014] = "File is not readable.";
   m_pil_errors[-3015] = "File is not writable.";
   m_pil_errors[-3016] = "Blank line encountered.";
   m_pil_errors[-3017] = "Comment line encountered.";
   m_pil_errors[-3018] = "Invalid line encountered.";
   m_pil_errors[-3019] = "No such parameter.";
   m_pil_errors[-3020] = "PFILES environment variable too long.";
   m_pil_errors[-3021] = "PFILES environment variable is badly formatted.";
   m_pil_errors[-3022] = "Cannot (un)lock parameter file.";
   m_pil_errors[-3023] = "Bogus parameters found in command line.";
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
      if (line != "" && line != " ") { //skip (most) blank lines
         lines.push_back(line);
      }
   }
}

} // namespace Likelihood
