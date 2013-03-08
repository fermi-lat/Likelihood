/**
 * @file UnbinnedLikelihood.cxx
 * @brief Unbinned likelihood analysis
 * @author J. Chiang <jchiang@slac.stanford.edu>
 * @author S. Fegan <sfegan@llr.in2p3.fr>
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/Attic/UnbinnedLikelihood.cxx,v 1.1.2.2 2013/03/06 14:36:52 sfegan Exp $
 */

#include "UnbinnedLikelihood.h"

UnbinnedLikelihood::UnbinnedLikelihood(const Observation & observation):
   LikelihoodBase(observation) {
#ifndef ST_LIKELIHOOD_NOTHREADS
   if (::getenv("ST_NTHREADS")) {
      m_nthread = atoi(::getenv("ST_NTHREADS"));
   }
#endif
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

   const std::vector<Event> & events = m_observation.eventCont().events();
   unsigned nevent = events.size();
   KahanAccumulator acc;

   // Data Sum - can be threaded

   // A quick note on the threading: the events are divided amongst
   // set of thereads with each thread accumulating the log densities
   // to a private per-thread accumulator (acc_p), which are then
   // added into the global accumulator in the critical section at the
   // end of the parallel task. The simplest scheme possible here
   // would be to fill an vector of event densities from the threads
   // and add them all up after. This would not require the use of
   // "acc_p" or the critical section (i.e. there would be only one
   // OMP directive) but it would require allocating a vector for the
   // event densities. Here that may not be a problem, but in
   // "getFreeDerivs" it might be, if there are a lot of free
   // parameters and a lot of events. So we use the slightly more
   // complex pattern outlined above to avoid having to allocate a
   // large array.

#pragma omp parallel num_threads(m_nthread)
   {
      KahanAccumulator acc_p;
      unsigned ievent;
#pragma omp for // schedule(dynamic,std::max(1,nevent/(50*m_nthread)))
      for(ievent=0; ievent<nevent; ievent++) {
	 const Event & event = events[ievent];	
	 if (!m_use_ebounds ||
	     (event.getEnergy() >= m_emin && event.getEnergy() <= m_emax)) {
	    KahanAccumulator density;
	    modelDensity(density, ievent);
	    acc_p.addLog(density);
	 }
      }
#pragma omp critical
      acc.add(acc_p);
   }

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

   const std::vector<Event> & events = m_observation.eventCont().events();
   unsigned nevent = events.size();
   unsigned nfree_param = getNumFreeParams();
   std::vector<KahanAccumulator> acc(nfree_param);

   // Data Sum - can be threaded
#pragma omp parallel num_threads(m_nthread)
   {
      std::vector<KahanAccumulator> acc_p(nfree_param);
      std::vector<double> log_density_derivs(nfree_param);
      unsigned ievent;
#pragma omp for // schedule(dynamic,std::max(1,nevent/(50*m_nthread)))
      for(ievent=0; ievent<nevent; ievent++) {
	 const Event & event = events[ievent];	
	 if (!m_use_ebounds ||
	     (event.getEnergy() >= m_emin && event.getEnergy() <= m_emax)) {
	    modelLogDensityDerivs(log_density_derivs, ievent);
	    for(unsigned iparam=0;iparam<nfree_param;iparam++) {
	       acc_p[iparam].add(log_density_derivs[iparam]);
	    }
	 }
      }
#pragma omp critical
      for(unsigned ifree_param=0;ifree_param<nfree_param;ifree_param++)
	acc[ifree_param].add(acc_p[ifree_param]);
   }

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

void UnbinnedLikelihood::modelDensity(KahanAccumulator & flux_acc,
				      unsigned ievent) const {
   const Event & event = m_observation.eventCont().events()[ievent];
   flux_acc = m_fixed_source_flux_in_bin[ievent];
   unisgned nfree = m_free_sources.size();
   for(unsigned ifree = 0; ifree<nfree; ifree++) {
      Source * src = m_free_sources[ifree];
      EventResponseCache * resp_cache = m_free_source_resp[ifree];
      EventResponse resp;
      if(resp_cache != 0) {
	 resp_cache->getResponse(resp, ievent);
      } else {
	 src->computeResponse(resp, event);
      }
      double flux = src->fluxDensity(event, resp);
      flux_acc.add(flux);
   }
}

void UnbinnedLikelihood::
modelLogDensityDerivs(std::vector<double> & log_flux_derivs,
		      unsigned ievent) const {
   const Event & event = m_observation.eventCont().events()[ievent];
   KahanAccumulator flux_acc = m_fixed_source_flux_in_bin[ievent];
   unisgned nfreesrc = m_free_sources.size();
   unsigned ifreeparam = 0;
   for(unsigned ifreesrc = 0; ifreesrc<nfreesrc; ifreesrc++) {
      Source * src = m_free_sources[ifreesrc];
      EventResponseCache * resp_cache = m_free_source_resp[ifreesrc];
      EventResponse resp;
      if(resp_cache != 0) {
	 resp_cache->getResponse(resp, ievent);
      } else {
	 src->computeResponse(resp, event);
      }
      double flux = src->fluxDensity(event, resp);
      flux_acc.add(flux);
     
      std::vector<std::string> params;
      src->getFreeParamNames(params);
      for(std::vector<std::string>::const_iterator iparam = 
	    params.begin(); iparam != params.end(); iparam++, ifreeparam++) {
	 double deriv = fluxDensity(event, *iparam, resp);
	 log_flux_derivs[ifreeparam] = deriv;
      }
   }
   assert(ifreeparam == log_flux_dervs.size());
   double flux = flux_acc.sum();
   for(std::vector<double>::iterator ilog_flux_deriv = log_flux_derivs.begin();
       ilog_flux_deriv != log_flux_derivs.begin(); ilog_flux_deriv++)
     *ilog_flux_deriv /= flux;
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
}

EventResponseCache * UnbinnedLikelihood::
computeSourceEventResponses(Source* src) const {
   const std::vector<Event> & events = m_observation.eventCont().events();
   unsigned nevent = events.size();
   EventResponseCache * resp_cache = 
     new EventResponseCache(src->edisp(), nevent);
   
   unsigned ievent;
#pragma omp parallel for private(ievent) num_threads(m_nthread) // schedule(dynamic,std::max(1,nevent/(50*m_nthread)))
   for(ievent=0; ievent<nevent; ievent++) {
      const Event & event = events[ievent];	
      EventResponse resp;
      src->computeResponse(resp, event);
#pragma omp critical
      resp_cache->setResponse(resp, ievent);
   }
   return resp_cache;
}

void UnbinnedLikelihood::addFreeSource(Source* src, const Source* callers_src) {
   std::map<Source *, EventResponseCache *> resps;
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
	    m_free_source_resp[ifree] = computeSourceEventResponses(src);
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
   unsigned nevent = events.size();

   // Data sum
   unsigned ievent;
#pragma omp parallel for private(ievent) num_threads(m_nthread) // schedule(dynamic,std::max(1,nevent/(50*m_nthread)))
   for(ievent=0; ievent<nevent; ievent++) {
      const Event & event = events[ievent];	
      EventResponse resp;
      src->computeResponse(resp, event);
      double fluxDensity = src->fluxDensity(event, resp);
      m_fixed_source_flux_in_bin[ievent].add(fluxDensity);
   }

   // Model integral - Npred
   m_fixed_source_npred.add(src->Npred());
}
