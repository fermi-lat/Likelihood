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
	for (unsigned int i = 0; i < getNumSrcs(); i++) {
	  double x = s_sources[i]->fluxDensity(m_events[j]);
	  (*Warray[i])[j] = x;
	  ztot += x;
	}
	for (unsigned int i = 0; i < getNumSrcs(); i++) {
	  if (ztot > 0.) (*Warray[i])[j] /= ztot;
	}
      }

      //! The M step.  Optimize parameters of each source.
      for (unsigned int i = 0; i < getNumSrcs(); i++) {
	OneSourceFunc f(s_sources[i], m_events, Warray[i]); 
	nPar += f.getNumFreeParams();
	optimizers::Arg arg;
	//	optimizers::Lbfgs opt(f);
	optimizers::Drmngb opt(f);
	opt.setNeedCovariance(false);
	//optimizers::Minuit opt(f);
	try {
	  opt.find_min(verbose, .1 * f.getNumFreeParams() * 1.1 * .2,
		       optimizers::ABSOLUTE);}
	catch (optimizers::Exception e) {
	  std::cerr << "Optimizer Exception " << e.code() << std::endl;
	  std::cerr << "  " << e.what() << std::endl;
	  std::cerr << " on iteration " << iteration+1
		    << " source " << i+1 << std::endl;
	  //	  assert(0);
	}

	logL += f.value(arg);
	//	std::cout << "Function value " << f.value(arg) << endl;
      }
      iteration++;
      if (verbose != 0)
	cout << "Iteration #" << iteration << ", logL = " << logL << 
	  ", old logL = " << oldLogL  << " params " << nPar << endl;
    } while (abs(logL-oldLogL) > 0.11 * nPar || oldLogL == 0.);

    //! Clean up before exit
    for (unsigned int i=0; i < getNumSrcs(); i++) {
      delete Warray[i];
    }
  }

} // namespace Likelihood
