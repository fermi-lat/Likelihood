#include <cstdlib>

#include <iostream>

#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/TipException.h"

#include "st_facilities/Util.h"

#include "Likelihood/CountsMap.h"

int main(int iargc, char * argv[]) {
   if (iargc != 3) {
      std::cout << "usage: count_map.exe <event_file> <sc_file>" 
                << std::endl;
      std::exit(0);
   }
   std::string event_file(argv[1]);
   std::vector<std::string> eventFiles;
   st_facilities::Util::resolve_fits_files(event_file, eventFiles);
   std::string sc_file(argv[2]);
   std::vector<std::string> scDataFiles;
   st_facilities::Util::resolve_fits_files(sc_file, scDataFiles);

   std::vector<double> energies;
   energies.push_back(100.);
   energies.push_back(2e5);
   Likelihood::CountsMap cmap(eventFiles[0], scDataFiles[0], 0, 0, "CAR",
                              720, 360, 0.5, 0, true, "RA", "DEC", 
                              energies);
   
   for (unsigned int i = 0; i < eventFiles.size(); i++) {
      const tip::Table * events 
         = tip::IFileSvc::instance().readTable(eventFiles[i], "events");
      cmap.binInput(events->begin(), events->end());
      delete events;
   }
   cmap.writeOutput("counts_map", "counts_map.fits");
}
