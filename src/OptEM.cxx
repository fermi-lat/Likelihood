/**
 * @file OptEM.cxx
 * @brief Implementation for Expectation Maximization class.
 * @author P. L. Nolan
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/OptEM.cxx,v 1.6 2004/05/09 19:18:21 jchiang Exp $
 */

#include "Verbosity.h"
#include "Likelihood/OptEM.h"
#include "Likelihood/OneSourceFunc.h"
#include "optimizers/Lbfgs.h"
#include "optimizers/Drmngb.h"
#include "optimizers/Minuit.h"
#include "optimizers/Exception.h"
#include "optimizers/FunctionTest.h"
#include "optimizers/ParameterNotFound.h"
#include "optimizers/OutOfBounds.h"
#include <vector>

namespace Likelihood {

  void OptEM::findMin(const int verbose) {
    double oldLogL;
    double logL = 0.;

    //! Build array of weight factors
    std::vector<std::vector<double>* > Warray;
    for (unsigned int i = 0; i < getNumSrcs(); i++) {
      Warray.push_back(new std::vector<double>(m_events.size()));
    }

    unsigned int iteration = 0;
    int nPar;

    //! Main EM loop.  Repeat until log(L) converges.
    do {
      oldLogL = logL;
      nPar = 0;
      logL = 0.;

      //! The E step.  Find weight factors
      for (unsigned int j = 0; j < m_events.size(); j++) {
	double ztot = 0.;
        std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
        for (unsigned int i = 0 ; srcIt != s_sources.end(); ++srcIt, i++) {
	  double x = srcIt->second->fluxDensity(m_events[j]);
	  (*Warray[i])[j] = x;
	  ztot += x;
	}
	for (unsigned int i = 0; i < getNumSrcs(); i++) {
	  if (ztot > 0.) (*Warray[i])[j] /= ztot;
	}
      }

      //! The M step.  Optimize parameters of each source.
      std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
      for (unsigned int i = 0 ; srcIt != s_sources.end(); ++srcIt, ++i) {
	OneSourceFunc f(srcIt->second, m_events, Warray[i]); 
	f.setEpsF(1.e-20);
	f.setEpsW(1.e-2);
	nPar += f.getNumFreeParams();
	optimizers::Arg arg;
	//	optimizers::Lbfgs opt(f);
	optimizers::Drmngb opt(f);
	opt.setNeedCovariance(false);
	//optimizers::Minuit opt(f);
	try {
	  opt.find_min(verbose, .1 * chifunc(f.getNumFreeParams()),
		       optimizers::ABSOLUTE);}
	catch (optimizers::Exception e) {
           if (print_output()) {
              std::cerr << "Optimizer Exception " << e.code() << std::endl;
              std::cerr << "  " << e.what() << std::endl;
              std::cerr << " on iteration " << iteration+1
                        << " source " << i+1 << std::endl;
           }
	  //	  assert(0);
	}

	logL += f.value(arg);
	//	std::cout << "Function value " << f.value(arg) << endl;
      }
      iteration++;
      if (print_output() && verbose != 0)
	std::cout << "Iteration #" << iteration << ", logL = " << logL << 
	  ", old logL = " << oldLogL  << " params " << nPar << std::endl;
    } while (fabs(logL-oldLogL) > 0.1*chifunc(nPar) || oldLogL == 0.);

    //! Clean up before exit
    for (unsigned int i=0; i < getNumSrcs(); i++) {
      delete Warray[i];
    }
  }

  double chifunc(int ndof) {
    // chi-squared value whose cumulative probability is 0.683
    if (ndof <= 0) {
      throw optimizers::Exception("chifunc invalid input < 1");
      abort();
    }
    double val = 0.;
    if (ndof < 10) {
      val = -.2265+ndof*(1.2719+ndof*(-.0097));
    } else {
      double x = (ndof-255)/141.5951;
      val = 265.7118+x*(144.4737+x*(-.2998+x*(.2429+x*(-.1234))));
    }
    return val;
  }

} // namespace Likelihood
