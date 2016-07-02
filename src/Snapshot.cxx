/**
 * @file Snapshot.cxx
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Snapshot.cxx,v 1.7 2016/07/02 01:23:08 echarles Exp $
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

  void Snapshot_Source::deep_copy_param_map(ParameterMap& theMap, bool delete_old) {
    for ( ParameterMap::iterator itr = theMap.begin(); itr != theMap.end(); itr++ ) {
      optimizers::Parameter* old =  itr->second;
      itr->second = new optimizers::Parameter(*old);
      if ( delete_old ) delete old;
    }
  }


  bool Snapshot_Source::compare_priors(const optimizers::Parameter& p1, 
				       const optimizers::Parameter& p2){

    if ( p1.has_prior() != p2.has_prior() ) return true;
    if ( ! ( p1.has_prior() || p1.has_prior() ) ) return false;
    if ( p1.log_prior().getNumParams() != p2.log_prior().getNumParams() ) return true;
    std::vector<std::string> n1;
    std::vector<std::string> n2;
    std::vector<double> v1;
    std::vector<double> v2;
    p1.log_prior().getParamNames(n1);
    p2.log_prior().getParamNames(n2);
    p1.log_prior().getParamValues(v1);
    p2.log_prior().getParamValues(v2);

    for ( size_t i(0); i < n1.size(); i++ ) {
      if ( n1[i] != n2[i] ) return true;
      double abs_diff = std::abs(v1[i]-v2[i]);
      double abs_sum = std::abs(v1[i]) + std::abs(v2[i]);
      double tol = abs_sum * 1e-5;
      bool value_changed = ( abs_diff > tol );   
      if ( value_changed ) return true;
    }
    
    return false;
  }
  
  unsigned Snapshot_Source::parameter_status(const optimizers::Parameter& p1, 
					     const optimizers::Parameter& p2) {

    unsigned retVal(0);
    
    double v1 = p1.getTrueValue();
    double v2 = p2.getTrueValue();
    double abs_diff = std::abs(v1-v2);
    double abs_sum = std::abs(v1) + std::abs(v2);
    double tol = abs_sum * 1e-5;
    
    bool value_changed = ( abs_diff > tol );   
    bool status_changed = p1.isFree() != p2.isFree();
    bool prior_changed = compare_priors(p1,p2);

    if ( value_changed ) retVal |= Snapshot_Status::Value_Changed;
    if ( status_changed ) retVal |= Snapshot_Status::Status_Changed;
    if ( prior_changed ) retVal |= Snapshot_Status::Prior_Changed;
    return retVal;
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
    model_name = func.genericName();
    const std::string norm_par_name = func.normPar().getName();
    std::vector<std::string> parNames;
    func.getParamNames(parNames);
    for ( std::vector<std::string>::const_iterator itr = parNames.begin(); itr != parNames.end(); itr++ ) {
      const optimizers::Parameter& par = func.getParam(*itr);
      if ( par.getName() == norm_par_name ) continue;
      optimizers::Parameter* ptr = clone ? new optimizers::Parameter(par) : const_cast<optimizers::Parameter*>(&par);
      theMap[*itr] = ptr;
    }
  }
  
  void Snapshot_Source::extract_spatial(const Source& src,
					std::string& model_name,
					ParameterMap& theMap,
					bool clone) {

    static const std::string spec_name("Spectrum");
    const std::map<std::string, optimizers::Function *>& srcFuncs = src.getSrcFuncs();
    for ( std::map<std::string, optimizers::Function *>::const_iterator itr = srcFuncs.begin();
	  itr != srcFuncs.end(); itr++ ) {
      if ( itr->first == spec_name ) continue;
      model_name += itr->first;
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


  Snapshot_Status Snapshot_Source::compare_norm(const optimizers::Parameter& par1, 
						const optimizers::Parameter& par2,
						bool fixed) {
    unsigned norm_status = parameter_status(par1,par2);
    if ( norm_status ) {
      norm_status |= (fixed ? Snapshot_Status::Fixed_Changed :  Snapshot_Status::Free_Changed );
    }
    if ( ( norm_status & Snapshot_Status::Status_Changed ) != 0  ) {
      // Note that fixed is the _old_ status, so the effect (i.e., freeing or fixing) is the opposite
      norm_status |= (fixed ? Snapshot_Status::Source_Freed :  Snapshot_Status::Source_Fixed );
    }
    return Snapshot_Status(norm_status);
  }
   
  Snapshot_Status Snapshot_Source::compare_spectral(const std::string& name1, const std::string& name2,
						    const ParameterMap& map1, const ParameterMap& map2,
						    bool fixed) {
    // Spectral model name has changed,
    // Don't bother comparing parameters
    unsigned status(0);
    if ( (name1 != name2) || (map1.size() != map2.size()) ) {
      status |= Snapshot_Status::Spec_Model_Type_Changed;
      status |= ( fixed ? Snapshot_Status::Fixed_Spec_Model_Type_Changed : Snapshot_Status::Free_Spec_Model_Type_Changed) ;
      return Snapshot_Status(status);
    }
    ParameterMap::const_iterator itr1 = map1.begin();
    ParameterMap::const_iterator itr2 = map2.begin();
    for ( ; itr1 != map1.end(); itr1++, itr2++ ) {
      // Parameter name has changed
      // Don't bother comparing values
      // Question, should we throw an exception, this shouldn't really happen
      if ( itr1->first != itr2->first ) {
	status |= Snapshot_Status::Spec_Model_Type_Changed;
	status |= ( fixed ? Snapshot_Status::Fixed_Spec_Model_Type_Changed : Snapshot_Status::Free_Spec_Model_Type_Changed) ;
	return Snapshot_Status(status);
      }
      const optimizers::Parameter* p1 = itr1->second;
      const optimizers::Parameter* p2 = itr2->second;
      unsigned par_status = parameter_status(*p1,*p2);
      if ( par_status ){
	par_status |= (fixed ? Snapshot_Status::Fixed_Changed :  Snapshot_Status::Free_Changed );
      }
      status |= par_status;      
    }
    return Snapshot_Status(status << Snapshot_Status::Spec_Offset );
  }
  
  Snapshot_Status Snapshot_Source::compare_spatial(const std::string& name1, const std::string& name2,
						   const ParameterMap& map1, const ParameterMap& map2,
						   bool fixed) {
    // Spectral model name has changed,
    // Don't bother comparing parameters
    unsigned status(0);
    if ( ( name1 != name2 ) || (map1.size() != map2.size()) ) {
      status |= Snapshot_Status::Spatial_Model_Type_Changed;
      status |= ( fixed ? Snapshot_Status::Fixed_Spatial_Model_Type_Changed : Snapshot_Status::Free_Spatial_Model_Type_Changed) ;
      return Snapshot_Status(status);
    }
    ParameterMap::const_iterator itr1 = map1.begin();
    ParameterMap::const_iterator itr2 = map2.begin();
    for ( ; itr1 != map1.end(); itr1++, itr2++ ) {
      // Parameter name has changed
      // Don't bother comparing values
      // Question, should we throw an exception, this shouldn't really happen
      if ( itr1->first != itr2->first )  {
	status |= Snapshot_Status::Spatial_Model_Type_Changed;
	status |= ( fixed ? Snapshot_Status::Fixed_Spatial_Model_Type_Changed : Snapshot_Status::Free_Spatial_Model_Type_Changed) ;
	return Snapshot_Status(status);
      }
      const optimizers::Parameter* p1 = itr1->second;
      const optimizers::Parameter* p2 = itr2->second;
      unsigned par_status = parameter_status(*p1,*p2);
      if ( par_status ){
	par_status |= (fixed ? Snapshot_Status::Fixed_Changed :  Snapshot_Status::Free_Changed );
      }
      status |= par_status;
    }
    return Snapshot_Status(status << Snapshot_Status::Spatial_Offset );
  }
  
  Snapshot_Status Snapshot_Source::compare_ancillary(const StringMap& map1, const StringMap& map2,
						     bool fixed) {
    unsigned status(0);
    if ( map1.size() != map2.size() ) {
      status |= Snapshot_Status::Ancillary_Changed;
      status |= ( fixed ? Snapshot_Status::Fixed_Ancillary_Changed : Snapshot_Status::Free_Ancillary_Changed);
      return Snapshot_Status(status);
    }
    StringMap::const_iterator itr1 = map1.begin();
    StringMap::const_iterator itr2 = map2.begin();
    for ( ; itr1 != map1.end(); itr1++, itr2++ ) {
      if ( ( itr1->first != itr2->first ) || 
	   ( itr1->second != itr2->second ) ) {
	status |= Snapshot_Status::Ancillary_Changed;
	status |= ( fixed ? Snapshot_Status::Fixed_Ancillary_Changed : Snapshot_Status::Free_Ancillary_Changed);
      }
    }
    return status;
  }     

  void Snapshot_Source::print_param(const optimizers::Parameter& par, const std::string& indent) {
    std::cout << indent 
	      << par.getName() << " : "  
	      << par.getTrueValue() << ' ' 
	      << (par.isFree() ? "Free  " : "Fixed ") 
	      << (par.has_prior() ? "T" : "F")
	      << std::endl;		  
  }

  void Snapshot_Source::print_param_map(const ParameterMap& parMap, const std::string& indent) {
    std::string indent2(indent);
    indent2 += indent;
    for ( ParameterMap::const_iterator itr = parMap.begin(); itr != parMap.end(); itr++ ){
      print_param(*(itr->second),indent2);
    }
  }

  void Snapshot_Source::print_string_map(const StringMap& strMap, const std::string& indent) {
    for ( StringMap::const_iterator itr = strMap.begin(); itr != strMap.end(); itr++ ){
      std::cout << indent << ' ' << itr->first << " : " << itr->second << std::endl;
    }
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
    m_fixed = ! (m_norm_param->isFree());
    m_owned = owned;
  }


  void Snapshot_Source::deep_copy() {
    if ( m_owned ) {
      optimizers::Parameter* old_norm = m_norm_param;
      m_norm_param = new optimizers::Parameter(*old_norm);
      delete old_norm;
    }
    deep_copy_param_map(m_spatial_params,m_owned);
    deep_copy_param_map(m_spectral_params,m_owned);
    m_owned = true;
  }

  Snapshot_Status Snapshot_Source::compare(const Source& other) const {
    Snapshot_Source other_snapshot(other);
    Snapshot_Status status;
    status |= compare_norm(*m_norm_param,*(other_snapshot.m_norm_param),m_fixed);
    status |= compare_spectral(m_spectral_model,other_snapshot.m_spectral_model,
			       m_spectral_params,other_snapshot.m_spectral_params,m_fixed);
    status |= compare_spatial(m_spatial_model,other_snapshot.m_spatial_model,
			      m_spatial_params,other_snapshot.m_spatial_params,m_fixed);
    status |= compare_ancillary(m_ancillary,other_snapshot.m_ancillary,m_fixed);
    return status;
  }
  

  void Snapshot_Source::print() const {
    static const std::string indent("  ");
    std::cout << indent << (m_fixed ? "Fixed" : "Free") << std::endl;
    print_param(*m_norm_param,indent);
    std::cout << indent << "Spectral Model : " <<  m_spectral_model << std::endl; 
    print_param_map(m_spectral_params,indent);
    std::cout << indent << "Spatial Model : " <<  m_spatial_model << std::endl; 
    print_param_map(m_spatial_params,indent);
    print_string_map(m_ancillary,indent);    
  }


  bool Snapshot::source_free(const Source& src) {
    // FIXME, do we want to test all the parameters
    return src.spectrum().normPar().isFree();
  }

  void Snapshot::latch_model(const SourceModel& model, bool owned) {
    m_sources.clear();
    const std::map<std::string, Source *>& sources = model.sources();
    for ( std::map<std::string, Source *>::const_iterator itr = sources.begin(); itr != sources.end(); itr++ ) {
      m_sources[itr->first] = Snapshot_Source();
      m_sources[itr->first].latch_source(*(itr->second),owned);
    }
  }

  void Snapshot::print() const {
    for ( std::map<std::string,Snapshot_Source>::const_iterator itr = m_sources.begin(); itr != m_sources.end(); itr++ ) {
      std::cout << itr->first << std::endl;
      itr->second.print();
    }
  }

  void Snapshot::print_source(const std::string& srcName) const {
    std::map<std::string,Snapshot_Source>::const_iterator itr = m_sources.find(srcName);
    if ( itr == m_sources.end() ) {
      std::cout << "No source named " << srcName << std::endl;
    }
    std::cout << itr->first << std::endl;
    itr->second.print();
  }

  
  Snapshot_Status Snapshot::compare_source(const std::string& srcName, const Source& src) const {
    std::map<std::string,Snapshot_Source>::const_iterator itrFind = m_sources.find(srcName);
    if ( itrFind == m_sources.end() ) {
      bool src_is_free = Snapshot::source_free(src);
      return src_is_free ? 
	Snapshot_Status(Snapshot_Status::Free_Source_Added) :
	Snapshot_Status(Snapshot_Status::Free_Source_Added);
    }
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
