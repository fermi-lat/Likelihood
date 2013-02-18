/**
 * @file UnbinnedLikelihood.cxx
 * @brief Unbinned likelihood analysis
 * @author J. Chiang <jchiang@slac.stanford.edu>
 * @author S. Fegan <sfegan@llr.in2p3.fr>
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Attic/LikelihoodBase.h,v 1.1.2.2 2013/01/30 16:09:29 jchiang Exp $
 */

#include "UnbinnedLikelihood.h"

UnbinnedLikelihood::UnbinnedLikelihood(const Observation & observation):
   LikelihoodBase(observation) {

}

UnbinnedLikelihood::~UnbinnedLikelihood() {
  for(std::vector<EventResponses *>::iterator iresp =
	m_free_source_resp.begin(); iresp!=m_free_source_resp.end(); iresp++) {
     delete *iresp;
  }
  delete m_last_removed_source_resp;
}

double UnbinnedLikelihood::value() const {
   std::clock_t start = std::clock();
   if (m_use_ebounds) {
      std::pair<double, double> ebounds
         = m_observation.roiCuts().getEnergyCuts();
      if (m_emin >= ebounds.second || m_emax <= ebounds.first) {
         // Updated energy range selection excludes all of the original 
         // data so return zero for log-likelihood.
         return 0;
      }
   }

   KahanAccumulator acc;
   const std::vector<Event> & events = m_observation.eventCont().events();

   // Data Sum - can be threaded
#ifndef ST_LIKELIHOOD_NOTHREADS
   if(m_threads) {
      unsigned nevent_per_thread = (events.size()+m_threads-1) / m_threads;
      std::vector<PartialDataSumThreadArgs> args(nthread);
      for(unsigned ithread=0;ithread<nthread;ithread++) {
	 args[ithread].thread_id     = 0;
	 args[ithread].like          = this;
	 args[ithread].acc.reset();
	 args[ithread].ievent_begin  = 
	   std::min(ithread*nevent_per_thread, events.size());
	 args[ithread].ievent_end    = 
	   std::min((ithread+1)*nevent_per_thread, events.size());
	 pthread_create(&args[ithread].thread_id, NULL,
			&startPartialDataSumThread, &args[ithread]);
      }
      for(unsigned ithread=0;ithread<nthread;ithread++) {
	 pthread_join(args[ithread].thread_id, NULL);
	 acc.add(args[ithread].acc);
      }
   } else {
#endif
      partialDataSum(acc, 0, events.size());
#ifndef ST_LIKELIHOOD_NOTHREADS
   }
#endif
				
   // Model integral
   KahanAccumulator npred_acc(m_fixed_source_npred);
   for(std::vector<Source *>::iterator ifreesrc = m_free_sources.begin();
       ifreesrc != m_free_sources.end(); ifreesrc++) {
      npred_acc.add(ifreesrc->Npred());
   }
   acc.add(-npred_acc);

   // Priors (including the fixed sources)
   for (std::vector<optimizers::Parameter>::const_iterator 
	  par = m_parameter.begin(); par != m_parameter.end(); ++par) {
      acc.add(par->log_prior_value());
   }

   // Zero value
   acc.add(m_zero_value);

   double value = acc.sum();
   return value;
}

void UnbinnedLikelihood::partialDataSum(KahanAccumulator & acc,
					unsigned ievent_begin, 
					unsigned ievent_end)
{
   const std::vector<Event> & events = m_observation.eventCont().events();
   KahanAccumulator density;
   for (size_t ievent = ievent_begin; ievent < ievent_end; ievent++) {
      if (m_use_ebounds && 
          (events[ievent].getEnergy() < m_emin 
	   || events[ievent].getEnergy() > m_emax)) {
         continue;
      }
      modelDensity(density, ievent);
      acc.addLog(density); 
   }
}

void UnbinnedLikelihood::modelDensity(KahanAccumulator & flux_acc,
				      unsigned ievent) const {
   const Event & event = m_observation.eventCont().events()[ievent];
   flux_acc = m_fixed_source_flux_in_bin[ievent];
   unisgned nfree = m_free_sources.size();
   for(unsigned ifree = 0; ifree<nfree; ifree++) {
      Source * src = m_free_sources[ifree];
      double flux = 
	 src->fluxDensity(event, (*m_free_source_resp[ifree])[ievent]);
       ev_acc.add(flux)
   }
   return ev_acc;
}

void UnbinnedLikelihood::getFreeDerivs(std::vector<double> & derivs) const {
   std::clock_t start = std::clock();
   if (m_use_ebounds) {
      std::pair<double, double> ebounds
         = m_observation.roiCuts().getEnergyCuts();
      if (m_emin >= ebounds.second || m_emax <= ebounds.first) {
         // Updated energy range selection excludes all of the original 
         // data so return zero for log-likelihood.
         return 0;
      }
   }

   unsigned nfree_param = getNumFreeParams();
   std::vector<KahanAccumulator> acc(nfree_param);
   const std::vector<Event> & events = m_observation.eventCont().events();

   // Data Sum - can be threaded
#ifndef ST_LIKELIHOOD_NOTHREADS
   if(m_threads) {
      unsigned nevent_per_thread = (events.size()+m_threads-1) / m_threads;
      std::vector<PartialDataDerivsThreadArgs> args(nthread);
      for(unsigned ithread=0;ithread<nthread;ithread++) {
	 args[ithread].thread_id     = 0;
	 args[ithread].like          = this;
	 args[ithread].acc.resice(acc.size());
	 args[ithread].ievent_begin  = 
	   std::min(ithread*nevent_per_thread, events.size());
	 args[ithread].ievent_end    = 
	   std::min((ithread+1)*nevent_per_thread, events.size());
	 pthread_create(&args[ithread].thread_id, NULL,
			&startPartialDataDerivsThread, &args[ithread]);
      }
      for(unsigned ithread=0;ithread<nthread;ithread++) {
	 pthread_join(args[ithread].thread_id, NULL);
	 for(unsigned iacc=0;iacc<axx.size();iacc++) {
	    acc[iacc].add(args[ithread].acc[iacc]);
	 }
      }
   } else {
#endif
      partialDataDerivs(acc, 0, events.size());
#ifndef ST_LIKELIHOOD_NOTHREADS
   }
#endif

   // Model integral
   std::vector<double> derivs;
   unsigned iacc = 0;
   for(std::vector<Source *>::const_iterator isrc = m_free_sources.begin();
       isrc != m_free_sources.end(); isrc++) {
      std::vector<std::string> params;
      isrc->getFreeParamNames(params);
      for(std::vector<std::string>::const_iterator iparam = params.begin();
	  iparam != params.end(); iparam++, iacc++) {
	 acc[iacc].add(-isrc->NpredDeriv(*iparam));
      }
   }
   assert(iacc == acc.size());

   // Priors
   iacc = 0;
   for (std::vector<optimizers::Parameter>::const_iterator 
	  par = m_parameter.begin(); par != m_parameter.end(); ++par, ++iacc) {
      if(par->isFree()) {
	 acc[iacc].add(par->log_prior_deriv());
      }
   }
}

void UnbinnedLikelihood::partialDataDerivs(std::vector<KahanAccumulator> & acc,
					   unsigned ievent_begin, 
					   unsigned ievent_end)
{
   const std::vector<Event> & events = m_observation.eventCont().events();
   for (size_t ievent = ievent_begin; ievent < ievent_end; ievent++) {
      if (m_use_ebounds && 
          (events[ievent].getEnergy() < m_emin 
	   || events[ievent].getEnergy() > m_emax)) {
         continue;
      }
      KahanAccumulator flux_acc;
      modelDensity(flux_acc, ievent);
      double flux = flux_acc.sum();

      unsigned iacc = 0;
      for(std::vector<Source *>::const_iterator isrc = m_free_sources.begin();
	  isrc != m_free_sources.end(); isrc++) {
	 std::vector<std::string> params;
	 isrc->getFreeParamNames(params);
	 for(std::vector<std::string>::const_iterator iparam = params.begin();
	     iparam != params.end(); iparam++, iacc++) {
	    double deriv = fluxDensity(event, *iparam, 
				       (*m_free_source_resp[ifree])[ievent]);
	    acc[iacc].add(deriv/flux);
	 }
      }
      assert(iacc == acc.size());
   }
}

void UnbinnedLikelihood::computeEventResponses(double sr_radius=30.0) {

}

void UnbinnedLikelihood::set_ebounds(double emin, double emax) {
   m_use_ebounds = true;
   m_emin = emin;
   m_emax = emax;

   size_t nee(observation().roiCuts().energies().size());
   double estep(std::log(emax/emin)/(nee-1));
   std::vector<double> energies;
   for (size_t k(0); k < nee-1; k++) {
      energies.push_back(emin*std::exp(k*estep));
   }
   energies.push_back(emax);

   std::map<std::string, Source *>::iterator it(m_sources.begin());
   for ( ; it != m_sources.end(); ++it) {
      it->second->computeExposure(energies);
   }

   // Update the fixed source npred value
   m_fixed_source_npred.reset();
   for(std::map<Source *, std::vector<double> >::iterator isec =
	 m_fixed_source_param_values.begin(); 
       isrc != m_fixed_source_param_values.end(); isrc++) {
      double npred = isrc->first->src->Npred();
      m_fixed_source_npred.add(npred);
   }
}

void UnbinnedLikelihood::unset_ebounds() {
   m_use_ebounds = false;
   m_emin = 0;
   m_emax = 0;

   std::map<std::string, Source *>::iterator it(m_sources.begin());
   for ( ; it != m_sources.end(); ++it) {
      it->second->computeExposure();
   }

   // Update the fixed source npred value
   m_fixed_source_npred.reset();
   for(std::map<Source *, std::vector<double> >::iterator isec =
	 m_fixed_source_param_values.begin(); 
       isrc != m_fixed_source_param_values.end(); isrc++) {
      double npred = isrc->first->src->Npred();
      m_fixed_source_npred.add(npred);
   }
}

void UnbinnedLikelihood::
partialSourceEventResponses(Source* src, EventResponses& resp,
			    unsigned ievent_begin, unsigned ievent_end) {
   const std::vector<Event> & events = m_observation.eventCont().events();
   for (size_t ievent = ievent_begin; ievent < ievent_end; ievent++) {
      src->computeResponse(resp[ievent], event[ievent]);
   }
}

void UnbinnedLikelihood::
computeSourceEventResponses(Source* src, EventResponses& resp) {
   const std::vector<Event> & events = m_observation.eventCont().events();
   resp.resize(events.size());
   
#ifndef ST_LIKELIHOOD_NOTHREADS
   if(m_threads) {
      unsigned nevent_per_thread = (events.size()+m_threads-1) / m_threads;
      std::vector<PartialSourceEventResponseThreadArgs> args(nthread);
      for(unsigned ithread=0;ithread<nthread;ithread++) {
	 args[ithread].thread_id = 0;
	 args[ithread].like = this;
	 args[ithread].src = src;
	 args[ithread].resp = &resp;
	 args[ithread].ievent_begin = ithread*nevent_per_thread;
	 args[ithread].ievent_end = 
	   std::min((ithread+1)*nevent_per_thread, events.size());
	 pthread_create(&args[ithread].thread_id, NULL,
			&startPartialSourceEventResponseThread, &args[ithread]);
      }
      for(unsigned ithread=0;ithread<nthread;ithread++) {
	 pthread_join(args[ithread].thread_id, NULL);
      }
   } else {
#endif
      partialSourceEventResponses(src, resp, 0, events.size());
#ifndef ST_LIKELIHOOD_NOTHREADS
   }
#endif
}

void UnbinnedLikelihood::addFreeSource(Source* src, const Source* callers_src) {
   std::map<Source *, EventResponses *> resps;
   for(unsigned ifree = 0; ifree<m_free_sources.size(); ifree++) {
      resps[m_free_sources[ifree]] = m_free_source_resp[ifree];
   }
   LikelihoodBase::addFreeSource(src, callers_src);
   m_free_source_resp.resize(m_free_sources.size());
   for(unsigned ifree = 0; ifree<m_free_sources.size(); ifree++) {
      if(m_free_sources[ifree] == src) {
	 if(callers_src && m_last_removed_source == callers_src) {
	    m_free_source_resp[ifree] = m_last_removed_source_resp;
	    m_last_removed_source = 0;
	    m_last_removed_source_resp = 0;
	 } else {
	    EventResponses * resp = new EventResponses;
	    computeSourceEventResponses(src, *resp);
	    m_free_source_resp[ifree] = resp;
	 }
      } else {
	 m_free_source_resp[ifree] = resps.at(m_free_sources[ifree]);
      }
   }
}

void UnbinnedLikelihood::removeFreeSource(Source* src) {
   std::map<Source *, EventResponses *> resps;
   for(unsigned ifree = 0; ifree<m_free_sources.size(); ifree++) {
      resps[m_free_sources[ifree]] = m_free_source_resp[ifree];
   }
   LikelihoodBase::removeFreeSource(src);   
   m_free_source_resp.resize(m_free_sources.size());
   for(unsigned ifree = 0; ifree<m_free_sources.size(); ifree++) {
      m_free_source_resp[ifree] = resps.at(m_free_sources[ifree]);
   }
   if(m_last_removed_source_resp)delete m_last_removed_source_resp;
   m_last_removed_source = src;
   m_last_removed_source_resp = resps.at(src);
}

void UnbinnedLikelihood::
addFixedSource(Source* src, const Source* callers_src) {
   (void)(callers_src);
   const std::vector<Event> & events = m_observation.eventCont().events();
   EventResponses resp;
   computeSourceEventResponses(src, resp);

   m_fixed_source_flux_in_bin.resize(events.size());
   
   // Data sum
   for(unsigned ievent = 0; ievent<events.size(); ievent++) {
      double fluxDensity = src->fluxDensity(event[ievent], resp[ievent]);
      m_fixed_source_flux_in_bin[ievent].add(fluxDensity);
   }

   // Model integral - Npred
   m_fixed_source_npred.add(sec->Npred());
}

#ifndef ST_LIKELIHOOD_NOTHREADS
void UnbinnedLikelihood::startPartialDataSumThread(void * pds_args) {
   PartialDataSumTheadArgs * args(pds_args);
   args->like->partialDataSum(args->acc, args->ievent_begin, args->ievent_end);
}

void UnbinnedLikelihood::startPartialDataDerivsThread(void * pds_args) {
   PartialDataDerivsTheadArgs * args(pds_args);
   args->like->partialDataDerivs(args->acc, 
				 args->ievent_begin, args->ievent_end);
}

void UnbinnedLikelihood::
startPartialSourceEventResponseThread(void * pesr_args) {
   PartialSourceEventResponseThreadArgs * args(pesr_args);
   args->like->partialSourceEventResponses(args->src, *args->resp,
					   args->ievent_begin, 
					   args->ievent_end);
}
#endif

