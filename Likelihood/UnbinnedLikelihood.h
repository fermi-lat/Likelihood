/**
 * @file UnbinnedLikelihood.h
 * @brief Unbinned likelihood analysis
 * @author J. Chiang <jchiang@slac.stanford.edu>
 * @author S. Fegan <sfegan@llr.in2p3.fr>
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Attic/UnbinnedLikelihood.h,v 1.1.2.3 2013/03/06 14:36:50 sfegan Exp $
 */

#ifndef Likelihood_UnbinnedLikelihood_h
#define Likelihood_UnbinnedLikelihood_h

#include "Likelihood/LikelihoodBase.h"
#include "Likelihood/Observation.h"
#include "Likelihood/Source.h"

namespace Likelihood {

class UnbinnedLikelihood : public LikelihoodBase {

public:

   UnbinnedLikelihood(const Observation & observation);

   virtual ~UnbinnedLikelihood();

   virtual double value() const;

   virtual void getFreeDerivs(std::vector<double> & derivs) const;

   /// Used by UnbinnedAnalysis
   void computeEventResponses(double sr_radius=30.0);
   void set_ebounds(double emin, double emax);
   void unset_ebounds();

protected:   

   EventResponseCache * computeSourceEventResponses(Source* src) const;

   virtual void addFreeSource(Source* src, const Source* callers_src);
   virtual void removeFreeSource(Source* src);

   virtual void addFixedSource(Source* src, const Source* callers_src);
   
   void modelDensity(KahanAccumulator & flux_acc, nsigned ievent) const;
   void modelLogDensityDerivs(std::vector<double> & log_flux_derivs,
			      unsigned ievent) const;

   bool                                    m_use_ebounds;
   double                                  m_emin;
   double                                  m_emax;

   std::vector<EventResponseCache *>       m_free_source_resp;
   Source *                                m_last_removed_source;
   EventResponses *                        m_last_removed_source_resp;

   unsigned                                m_nthread;
};

} // namespace Likelihood

#endif // Likelihood_UnbinnedLikelihood_h
