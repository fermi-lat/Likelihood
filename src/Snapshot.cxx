/**
 * @file Snapshot.cxx
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/FitScanner.cxx,v 1.20 2016/06/25 00:10:33 echarles Exp $
 */


#include "Likelihood/Snapshot.h"

#include <cstdio>
#include <stdexcept>

#include "optimizers/Parameter.h"
#include "optimizers/Function.h"

#include "Likelihood/Source.h"
#include "Likelihood/SourceModel.h"


namespace Likelihood {

  void Snapshot_Source::delete_map(ParameterMap& theMap) {
    for ( ParameterMap::iterator itr = theMap.begin(); itr != theMap.end(); itr++ ) {
      delete itr->second;
    }
    theMap.clear();
  }

  bool Snapshot_Source::parameter_unchanged(const optimizers::Parameter& p1, 
					    const optimizers::Parameter& p2) {
    double v1 = p1.getTrueValue();
    double v2 = p2.getTrueValue();
    double abs_diff = std::abs(v1-v2);
    double sum = v1+v2;
    double tol = sum * 1e-5;
    return ( abs_diff > tol );

  }

  void Snapshot_Source::extract_norm(const Source& src,
				     optimizers::Parameter*& theParam,
				     bool clone) {
    const optimizers::Parameter& np = src.spectrum().normPar();
    theParam = clone ? new optimizers::Parameter(np) : const_cast<optimizers::Parameter*>(&np);    
  }


  void Snapshot_Source::extract_spectral(const Source& src,
					 std::string& model_name,
					 ParameterMap& theMap,
					 bool clone ) {
    theMap.clear();
    const optimizers::Function& func = src.spectrum();
    model_name = func.getName();
    std::vector<std::string> parNames;
    func.getParamNames(parNames);
    for ( std::vector<std::string>::const_iterator itr = parNames.begin(); itr != parNames.end(); itr++ ) {
      const optimizers::Parameter& par = func.getParam(*itr);
      optimizers::Parameter* ptr = clone ? new optimizers::Parameter(par) : const_cast<optimizers::Parameter*>(&par);
      theMap[*itr] = ptr;
    }
  }
  
  void Snapshot_Source::extract_spatial(const Source& src,
					std::string& model_name,
					ParameterMap& theMap,
					bool clone) {

    const std::string& spec_name = src.spectrum().getName();
    const std::map<std::string, optimizers::Function *>& srcFuncs = src.getSrcFuncs();
    for ( std::map<std::string, optimizers::Function *>::const_iterator itr = srcFuncs.begin();
	  itr != srcFuncs.end(); itr++ ) {
      if ( itr->first == spec_name ) continue;
      const optimizers::Function& func = *(itr->second); 
      std::vector<std::string> parNames;
      func.getParamNames(parNames);
      for ( std::vector<std::string>::const_iterator itr2 = parNames.begin(); itr2 != parNames.end(); itr2++ ) {
	const optimizers::Parameter& par = func.getParam(*itr2);
	optimizers::Parameter* ptr = clone ? new optimizers::Parameter(par) : const_cast<optimizers::Parameter*>(&par);
	theMap[*itr2] = ptr;
      }
    }
  }
  
  void Snapshot_Source::extract_ancillary(const Source& /*src*/,
					  StringMap& theMap) {
    theMap.clear();
  }   


  Snapshot_Status Snapshot_Source::compare_norm(const optimizers::Parameter& par1, const optimizers::Parameter& par2) {
    return parameter_unchanged(par1,par2) ? 
      Snapshot_Status(Snapshot_Status::Unchanged) : 
      Snapshot_Status(Snapshot_Status::Norm_Changed);
  }
   
  Snapshot_Status Snapshot_Source::compare_spectral(const std::string& name1, const std::string& name2,
						    const ParameterMap& map1, const ParameterMap& map2) {
    // Spectral model name has changed,
    // Don't bother comparing parameters
    if ( name1 != name2 ) return Snapshot_Status(Snapshot_Status::SpecModel_Changed);
    if ( map1.size() != map2.size() ) return Snapshot_Status(Snapshot_Status::SpecModel_Changed);
    ParameterMap::const_iterator itr1 = map1.begin();
    ParameterMap::const_iterator itr2 = map2.begin();
    Snapshot_Status status;
    for ( ; itr1 != map1.end(); itr1++, itr2++ ) {
      // Parameter name has changed
      // Don't bother comparing values
      // Question, should we throw an exception, this shouldn't really happen
      if ( itr1->first != itr2->first ) return Snapshot_Status(Snapshot_Status::SpecModel_Changed);
      const optimizers::Parameter* p1 = itr1->second;
      const optimizers::Parameter* p2 = itr2->second;
      if ( parameter_unchanged(*p1,*p2) ) continue;
      if ( p1->isFree() )  {
	status |= Snapshot_Status::SpecFree_Changed;
      } else {
	status |= Snapshot_Status::SpecFixed_Changed;
      }
    }
    return status;
  }
  
  Snapshot_Status Snapshot_Source::compare_spatial(const std::string& name1, const std::string& name2,
						   const ParameterMap& map1, const ParameterMap& map2) {
   // Spectral model name has changed,
    // Don't bother comparing parameters
    if ( name1 != name2 ) return Snapshot_Status(Snapshot_Status::SpatialModel_Changed);
    if ( map1.size() != map2.size() ) return Snapshot_Status(Snapshot_Status::SpatialModel_Changed);
    ParameterMap::const_iterator itr1 = map1.begin();
    ParameterMap::const_iterator itr2 = map2.begin();
    Snapshot_Status status;
    for ( ; itr1 != map1.end(); itr1++, itr2++ ) {
      // Parameter name has changed
      // Don't bother comparing values
      // Question, should we throw an exception, this shouldn't really happen
      if ( itr1->first != itr2->first ) return Snapshot_Status(Snapshot_Status::SpatialModel_Changed);
      const optimizers::Parameter* p1 = itr1->second;
      const optimizers::Parameter* p2 = itr2->second;
      if ( parameter_unchanged(*p1,*p2) ) continue;
      if ( p1->isFree() )  {
	status |= Snapshot_Status::SpatialParam_Changed;
      } else {
	status |= Snapshot_Status::SpatialParam_Changed;
      }
    }
    return status; 
  }
  
  Snapshot_Status Snapshot_Source::compare_ancillary(const StringMap& map1, const StringMap& map2) {
    if ( map1.size() != map2.size() ) return Snapshot_Status(Snapshot_Status::Ancillary_Changed);
    StringMap::const_iterator itr1 = map1.begin();
    StringMap::const_iterator itr2 = map2.begin();
    Snapshot_Status status;
    for ( ; itr1 != map1.end(); itr1++, itr2++ ) {
      if ( itr1->first != itr2->first ) status |= Snapshot_Status::Ancillary_Changed;
      if ( itr1->second != itr2->second ) status |= Snapshot_Status::Ancillary_Changed;
    }
    return status;
  }     
    
  Snapshot_Source::~Snapshot_Source() {
    if ( m_owned ) {
      delete m_norm_param;
      delete_map(m_spatial_params);
      delete_map(m_spectral_params);
    }
  }

  void Snapshot_Source::latch_source(const Source& src, bool owned) {
    if ( m_owned ) {
      delete m_norm_param;			
      m_norm_param = 0;
      delete_map(m_spatial_params);
      delete_map(m_spectral_params);
    }
    extract_norm(src,m_norm_param,owned);
    extract_spectral(src,m_spectral_model,m_spectral_params,owned);
    extract_spatial(src,m_spatial_model,m_spatial_params,owned);
    extract_ancillary(src,m_ancillary);
    m_owned = owned;
  }

  Snapshot_Status Snapshot_Source::compare(const Source& other) const {
    Snapshot_Source other_snapshot(other);
    Snapshot_Status status;
    status |= compare_norm(*m_norm_param,*(other_snapshot.m_norm_param));
    status |= compare_spectral(m_spectral_model,other_snapshot.m_spectral_model,
			       m_spectral_params,other_snapshot.m_spectral_params);
    status |= compare_spatial(m_spatial_model,other_snapshot.m_spatial_model,
			      m_spatial_params,other_snapshot.m_spatial_params);
    status |= compare_ancillary(m_ancillary,other_snapshot.m_ancillary);
    return status;
  }
  
  void Snapshot::latch_model(const SourceModel& model) {
    const std::map<std::string, Source *>& sources = model.sources();
    for ( std::map<std::string, Source *>::const_iterator itr = sources.begin(); itr != sources.end(); itr++ ) {
      m_sources.insert(std::make_pair<std::string,Snapshot_Source>(itr->first,Snapshot_Source(*(itr->second),true)));
    }
  }
  
  Snapshot_Status Snapshot::compare_source(const std::string& srcName, const Source& src) const {
    std::map<std::string,Snapshot_Source>::const_iterator itrFind = m_sources.find(srcName);
    if ( itrFind == m_sources.end() ) return Snapshot_Status(Snapshot_Status::Source_Added);
    return itrFind->second.compare(src);
  }

  Snapshot_Status Snapshot::compare_model(const SourceModel& model, std::vector<std::string>& changed_sources) const {
    changed_sources.clear();
    const std::map<std::string, Source *>& sources = model.sources();
    Snapshot_Status status;
    for ( std::map<std::string, Source *>::const_iterator itr = sources.begin(); itr != sources.end(); itr++ ) {
      Snapshot_Status src_comp = compare_source(itr->first,*(itr->second));
      if ( src_comp.changed() ) changed_sources.push_back(itr->first);
      status |= src_comp;
    }
    return status;
  }
 

} // namespace Likelihood
