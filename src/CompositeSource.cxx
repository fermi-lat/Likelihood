/** 
 * @file CompositeSource.cxx
 * @brief A source made of a collection of sources.
 * @author E. Charles
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/CompositeSource.cxx,v 1.2 2017/09/29 01:55:18 echarles Exp $
 */

#include "Likelihood/CompositeSource.h"
#include "Likelihood/Accumulator.h"
#include "Likelihood/SourceMapCache.h"
#include "Likelihood/SourceMap.h"
#include "Likelihood/BinnedLikelihood.h"

namespace Likelihood {

  CompositeSource::CompositeSource(const Observation& observation)
    :Source(&observation),
     m_sourceModel(observation),
     m_srcMapCache(0),
     m_xmlFile(""){     
    m_srcType = Source::Composite;
    m_energies = m_observation->roiCuts().energies();
  }
    

  CompositeSource::CompositeSource(const  CompositeSource&rhs)
    :Source(rhs),
     m_sourceModel(rhs.m_sourceModel),
     m_srcMapCache(rhs.m_srcMapCache != 0 ? rhs.m_srcMapCache->clone() : 0){
    m_srcType = Source::Composite;
    m_functions["Spectrum"] = m_spectrum;
  }


  CompositeSource::CompositeSource(const Observation& observation,
				   const std::string& name,
				   const std::string& specFuncName)
    :Source(&observation),
     m_sourceModel(observation),
     m_srcMapCache(0){
    m_srcType = Source::Composite;
    m_energies = m_observation->roiCuts().energies();
    setName(name);
    setSpectrum(specFuncName);
  }


  CompositeSource::~CompositeSource() {
    delete m_srcMapCache;
  }

  auto CompositeSource::clone() const -> Source * {
    return new CompositeSource(*this);
  }


  void CompositeSource::addSource(Source *src, SourceMap* srcMap, 
				  bool fromClone) {
    m_sourceModel.addSource(src,fromClone);
    if ( srcMap == 0 ) return;
    if ( m_srcMapCache == 0 ) {
      throw std::runtime_error("Can't add a source model to a CompositeSource until buildSourceMapCache has been called");
      return;
    }
    Source* added = m_sourceModel.getSource(src->getName());
    srcMap->setSource(*added);
    m_srcMapCache->insertSourceMap(added->getName(),*srcMap);
  }
        

  auto CompositeSource::deleteSource(const std::string &srcName,
					SourceMap*& srcMap) -> Source* {
    srcMap = 0;
    if ( m_srcMapCache != 0 && m_srcMapCache->hasSourceMap(srcName) ) {
      srcMap = m_srcMapCache->removeSourceMap(srcName);
    }
    Source* src = m_sourceModel.deleteSource(srcName);
    return src;
  }

  void CompositeSource::deleteAllSources() {
    m_sourceModel.deleteAllSources();
    delete m_srcMapCache;
    m_srcMapCache = 0;
  }

  auto CompositeSource::fluxDensity(const Event & evt, 
				      CachedResponse* cResp) const -> double {
    optimizers::dArg energy_arg(evt.getEnergy() );
    double spec_factor = (*m_spectrum)(energy_arg);
    Kahan_Accumulator accum;
    for ( std::map<std::string, Source*>::const_iterator itr = m_sourceModel.sources().begin();
	  itr != m_sourceModel.sources().end(); itr++ ) {
      accum.add( itr->second->fluxDensity(evt,cResp) );
    }
    double retVal = accum.total();
    retVal *= spec_factor;
    return retVal;
  }
 
  auto CompositeSource::fluxDensity(double inclination, double phi, double energy, 
				      const astro::SkyDir & appDir, 
				      int evtType, double time, 
				      CachedResponse* cResp) const -> double {
    optimizers::dArg energy_arg(energy);
    double spec_factor = (*m_spectrum)(energy_arg);
    Kahan_Accumulator accum;
    for ( std::map<std::string, Source*>::const_iterator itr = m_sourceModel.sources().begin();
	  itr != m_sourceModel.sources().end(); itr++ ) {
      accum.add( itr->second->fluxDensity(inclination,phi,energy,appDir,evtType,time,cResp) );
    }
    double retVal = accum.total();
    retVal *= spec_factor;
    return retVal;
  }
 
  auto CompositeSource::fluxDensityDeriv(const Event & evt, 
					   const std::string & paramName,
					   CachedResponse* cResp) const -> double {
    optimizers::dArg energy_arg(evt.getEnergy() );
    double spec_factor = (*m_spectrum)(energy_arg);
    Kahan_Accumulator accum_pos;
    Kahan_Accumulator accum_neg;
    for ( std::map<std::string, Source*>::const_iterator itr = m_sourceModel.sources().begin();
	  itr != m_sourceModel.sources().end(); itr++ ) {
      double val = itr->second->fluxDensityDeriv(evt,paramName,cResp);
      Kahan_Accumulator& accum = val > 0 ? accum_pos : accum_neg;
      double addend = val > 0 ? val : -val;
      accum.add( addend );
    }
    double retVal = accum_pos.total() - accum_neg.total();
    retVal *= spec_factor;
    return retVal;
  }
  
  auto CompositeSource::fluxDensityDeriv(double inclination, double phi, 
					   double energy, const astro::SkyDir & appDir,
					   int evtType, double time, 
					   const std::string & paramName,
					   CachedResponse* cResp) const -> double {
    optimizers::dArg energy_arg(energy);
    double spec_factor = (*m_spectrum)(energy_arg);
    Kahan_Accumulator accum_pos;
    Kahan_Accumulator accum_neg;
    for ( std::map<std::string, Source*>::const_iterator itr = m_sourceModel.sources().begin();
	  itr != m_sourceModel.sources().end(); itr++ ) {
      double val = itr->second->fluxDensityDeriv(inclination,phi,energy,appDir,evtType,time,paramName,cResp);
      Kahan_Accumulator& accum = val > 0 ? accum_pos : accum_neg;
      double addend = val > 0 ? val : -val;
      accum.add( addend );
    }
    double retVal = accum_pos.total() - accum_neg.total();
    retVal *= spec_factor;
    return retVal;  
  }


  void CompositeSource::computeExposure(const std::vector<double> & energies,
					bool verbose) {
    throw std::runtime_error("CompositeSource::computeExposure not implemented");
  }
  
  auto CompositeSource::flux() const -> double {
    return computeEnergyIntegral(*m_spectrum, m_energies);
  }
    
  auto CompositeSource::fluxDeriv(const std::string & parName) const -> double {
    FluxDeriv my_functor(*m_spectrum, parName);
    return computeEnergyIntegral(my_functor, m_energies);
  }
  
  /// @return Photon flux integrated over the given energy range.
  /// Units are #/cm^2/s
  auto CompositeSource::flux(double emin, double emax, size_t npts) const -> double {
    return computeEnergyIntegral(*m_spectrum, emin, emax, npts);
  }
  
  /// @return Derivative of integrated photon flux wrt the named parameter
  /// over the given energy range.
  auto CompositeSource::fluxDeriv(const std::string & parName, double emin,
				    double emax, size_t npts) const -> double {
    FluxDeriv my_functor(*m_spectrum, parName);
    return computeEnergyIntegral(my_functor, emin, emax, npts);
  }
  
  /// @return Energy flux integrated over the ROI energy bounds. 
  /// Units are MeV/cm^2/s
  auto CompositeSource::energyFlux() const -> double {
    EnergyFlux my_functor(*m_spectrum);
    return computeEnergyIntegral(my_functor, m_energies);
  }
  
  /// @return Derivative of integrated energy flux wrt the named parameter
  auto CompositeSource::energyFluxDeriv(const std::string & parName) const -> double {
    EnergyFluxDeriv my_functor(*m_spectrum, parName);
    return computeEnergyIntegral(my_functor, m_energies); 
  }
  
  /// @return Energy flux integrated over the given energy range.
  /// Units are MeV/cm^2/s
  auto CompositeSource::energyFlux(double emin, double emax,
				     size_t npts) const -> double {
    EnergyFlux my_functor(*m_spectrum);
    return computeEnergyIntegral(my_functor, emin, emax, npts);
  }
  
  /// @return Derivative of integrated energy flux wrt the named parameter
  /// over the given energy range.
  auto CompositeSource::energyFluxDeriv(const std::string & parName, double emin,
					  double emax, size_t npts) const -> double {
    EnergyFluxDeriv my_functor(*m_spectrum, parName);
    return computeEnergyIntegral(my_functor, emin, emax, npts);
  }

  
  
  auto CompositeSource::steal_source(SourceModel& other,
					const std::string& srcName) -> Source* {
    
    Source* src = other.deleteSource(srcName);
    BinnedLikelihood* bl = dynamic_cast<BinnedLikelihood*>(&other);
    SourceMap* srcMap(0);
    if ( bl != 0 ) {
      srcMap = bl->removeSourceMap(srcName);
    }
    addSource(src,srcMap,false);
    return src;
  }
  
  auto CompositeSource::give_source(SourceModel& other,
				       const std::string& srcName) -> Source* {
    SourceMap* srcMap(0);
    if ( m_srcMapCache != 0 ) {
      srcMap = m_srcMapCache->removeSourceMap(srcName);
    }
    Source* src = m_sourceModel.give_source(other,srcName,srcMap);
    BinnedLikelihood* bl = dynamic_cast<BinnedLikelihood*>(&other);    
    if ( bl == 0 ) {      
      // Not a binned likelihood, so we should clean up the SourceMap
      delete srcMap;     
    }
    return src;
  }
 

  auto CompositeSource::buildSourceMapCache(const BinnedCountsCache& dataCache,
						       const std::string & srcMapsFile,
						       const Drm* drm) -> SourceMapCache* {
    if ( m_srcMapCache != 0 ) {
      std::cout << "WARNING, deleting existing SourceMapCache from CompositeSource::" << getName() << std::endl;
      delete m_srcMapCache;
    } 
    m_srcMapCache = new SourceMapCache(dataCache,*m_observation,srcMapsFile,m_config,drm);
    return m_srcMapCache;
  }


  void CompositeSource::fillSummedSourceMap(std::vector<float>& model, int kmin, int kmax) const {
    if ( m_srcMapCache == 0 ) {
      throw std::runtime_error("CompositeSource::buildSourceMapCache must be called before fillSummedSourceMap");
    }
    std::vector<std::string> srcNames;
    std::vector<const Source*> srcs;
    model.clear();
    model.resize(m_srcMapCache->dataCache().source_map_size(),0.);
    m_sourceModel.getSrcNames(srcNames);
    m_sourceModel.getSources(srcNames,srcs);
    m_srcMapCache->fillSummedSourceMap(srcs,model,kmin,kmax);
  }
  



} // namespace Likelihood
