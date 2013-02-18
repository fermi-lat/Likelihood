/**
 * @file UnbinnedLikelihood.h
 * @brief Unbinned likelihood analysis
 * @author J. Chiang <jchiang@slac.stanford.edu>
 * @author S. Fegan <sfegan@llr.in2p3.fr>
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Attic/LikelihoodBase.h,v 1.1.2.2 2013/01/30 16:09:29 jchiang Exp $
 */

#ifndef Likelihood_UnbinnedLikelihood_h
#define Likelihood_UnbinnedLikelihood_h

#ifndef ST_LIKELIHOOD_NOTHREADS
#include <pthread.h>
#endif

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

   typedef std::vector<Source::Response> EventResponses;

   void computeSourceEventResponses(Source* src, EventResponses& resp) const;
   void partialSourceEventResponses(Source* src, EventResponses& resp,
				    unsigned ievent_begin, 
				    unsigned ievent_end) const;

   virtual void addFreeSource(Source* src, const Source* callers_src);
   virtual void removeFreeSource(Source* src);

   virtual void addFixedSource(Source* src, const Source* callers_src);
   
   void partialDataSum(KahanAccumulator& acc, 
		       unsigned ievent_begin, unsigned ievent_end) const;
   void modelDensity(KahanAccumulator & flux_acc, nsigned ievent) const;

   void partialDataDerivs(std::vector<KahanAccumulator> & acc, 
			  unsigned ievent_begin, unsigned ievent_end) const;

#ifndef ST_LIKELIHOOD_NOTHREADS
   struct PartialDataSumThreadArgs {
      pthread_t thread_id;
      UnbinnedLikelihood* like;
      KahanAccumulator acc;
      unsigned ievent_begin;
      unsigned ievent_end;
   };

   static void startPartialDataSumThread(void * pds_args);

   struct PartialDataDerivsThreadArgs {
      pthread_t thread_id;
      UnbinnedLikelihood* like;
      std::vector<KahanAccumulator> acc;
      unsigned ievent_begin;
      unsigned ievent_end;
   };

   static void startPartialDataDerivsThread(void * pdd_args);

   struct PartialSourceEventResponseThreadArgs {
      pthread_t thread_id;
      UnbinnedLikelihood* like;
      Source * src;
      EventResponses * resp;
      unsigned ievent_begin;
      unsigned ievent_end;
   };
   
   static void startPartialSourceEventResponseThread(void * pesr_args);
#endif


   bool                                    m_use_ebounds;
   double                                  m_emin;
   double                                  m_emax;

   std::vector<EventResponses *>           m_free_source_resp;
   Source *                                m_last_removed_source;
   EventResponses *                        m_last_removed_source_resp;

#ifndef ST_LIKELIHOOD_NOTHREADS
   unsigned                                m_nthread;
#endif
};

} // namespace Likelihood

#endif // Likelihood_UnbinnedLikelihood_h
