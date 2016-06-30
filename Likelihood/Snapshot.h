/**
 * @file Snapshot.h
 * @author Eric Charles
 *
 * A class to provide a Snapshot of a SourceModel and allow for comparisons
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Snapshot.h,v 1.1 2016/06/30 00:27:06 echarles Exp $
 */

#ifndef Likelihood_Scanshot_h
#define Likelihood_Scanshot_h

// stl includes
#include <string>
#include <vector>
#include <map>


// Fermi includes


namespace optimizers{
  class Parameter;
}


// This class lives in the Likelihood namespace
namespace Likelihood {

  /**
   * @class Snapshot
   * 
   */
   
  class SourceModel;
  class Source;

  struct Snapshot_Status {
  public:
    enum { 	
      Unchanged = 0,                //! All parameters match	
      Norm_Changed = 0x1,           //! The normalization has changed
      SpecFree_Changed = 0x2,       //! A free spectral parameter (besides the normalization) has changed
      SpecFixed_Changed = 0x4,      //! A fixed spectral parameter (besides the normalization) has changed       
      SpecParam_Shape_Mask = 0x6,   //! A spectral paramter (besides the normalization) has changed
      SpecParam_Any_Mask = 0x7,     //! A spectral parameter (including the normalization) has changed
      SpecModel_Changed = 0x8,      //! The spectral model type has changed
      Spectrum_Changed_Mask = 0xF,  //! The spectrum has changed
      SpatialParam_Changed = 0x10,  //! A spatial model parameter has changed
      SpatialModel_Changed = 0x20,  //! The spatial model type has changed
      Spatial_Changed_Mask = 0x30,  //! The spatial model has chagned
      Ancillary_Changed = 0x40,     //! Ancillary data has changed
      Source_Added = 0x100,         //! A source has been added
      Source_Removed = 0x200,       //! A source has been removed
      Source_Changed_Mask = 0x300   //! A source has been added or removed
    } Flags;
    
    Snapshot_Status(unsigned val=Unchanged) :m_val(val) {}
    ~Snapshot_Status(){ }
    
    inline bool unchanged() { return m_val == 0; }
    inline bool changed() { return m_val != 0; }
    inline bool norm_changed() { return (m_val & Norm_Changed) != 0; }
    inline bool spec_free_changed() { return (m_val & SpecFree_Changed) != 0; }
    inline bool spec_fixed_changed() { return (m_val & SpecFixed_Changed) != 0; }
    inline bool spec_shape_changed() { return (m_val & SpecParam_Shape_Mask) != 0; }
    inline bool spec_any_changed() { return (m_val & SpecParam_Any_Mask) != 0; }
    inline bool spec_model_changed() { return (m_val & SpecModel_Changed) != 0; }    
    inline bool spectrum_changed() { return (m_val & Spectrum_Changed_Mask) != 0; }
    inline bool spatial_param_changed() { return (m_val & SpatialParam_Changed) != 0; }
    inline bool spatial_model_changed() { return (m_val & SpatialModel_Changed) != 0; }
    inline bool spatial_changed() { return (m_val & Spatial_Changed_Mask) != 0; }   
    inline bool ancillary_changed() { return (m_val & Ancillary_Changed) != 0; }
    inline bool source_added() { return (m_val & Source_Added) != 0; }
    inline bool source_removed() { return (m_val & Source_Removed) != 0; }
    inline bool source_changed() { return (m_val & Source_Changed_Mask) != 0; }
     
    inline unsigned operator()() const { return m_val; }
    inline Snapshot_Status operator|(const Snapshot_Status& other) const {
      return Snapshot_Status( m_val | other.m_val );
    }
    inline Snapshot_Status operator|(unsigned val) const {
      return Snapshot_Status( m_val | val );
    }
    inline Snapshot_Status operator&(const Snapshot_Status& other) const {
      return Snapshot_Status( m_val & other.m_val );
    }
    inline Snapshot_Status operator&(unsigned val) const {
      return Snapshot_Status( m_val | val );
    }

    inline Snapshot_Status& operator|=(const Snapshot_Status& other)  {
      m_val |= other();
      return *this;
    }
    inline Snapshot_Status& operator|=(unsigned val)  {
      m_val |= val;
      return *this;
    }
    inline Snapshot_Status& operator&=(const Snapshot_Status& other)  {
      m_val &= other();
      return *this;
    }
    inline Snapshot_Status& operator&=(unsigned val)  {
      m_val &= val;
      return *this;
    }

    inline Snapshot_Status& operator=(const Snapshot_Status& other) {
      if ( this != &other ) { m_val = other.m_val; }
      return *this;
    }
    inline Snapshot_Status& operator=(unsigned val) {
      m_val = val;
      return *this;
    }

  private:
    unsigned m_val;
  };


  typedef std::map<std::string,optimizers::Parameter*> ParameterMap;
  typedef std::map<std::string,std::string> StringMap;

  class Snapshot_Source {
    /* A utility class to cache a Snapshot of a Source, and to allow for
       comparisions
    */
    
  public:
    
    static void delete_map(ParameterMap& theMap);

    static bool parameter_unchanged(const optimizers::Parameter& p1, 
				    const optimizers::Parameter& p2);

    static void extract_norm(const Source& src,
			     optimizers::Parameter*& theParam,
			     bool clone = false);

    static void extract_spectral(const Source& src,
				 std::string& model_name,
				 ParameterMap& theMap,
				 bool clone = false);

    static void extract_spatial(const Source& src,
				std::string& model_name,
				ParameterMap& theMap,
				bool clone = false);
  
    static void extract_ancillary(const Source& src,
				  StringMap& theMap);    
    

    static Snapshot_Status compare_norm(const optimizers::Parameter& par1, const optimizers::Parameter& par2);     

    static Snapshot_Status compare_spectral(const std::string& name1, const std::string& name2,
					    const ParameterMap& map1, const ParameterMap& map2); 
    
    static Snapshot_Status compare_spatial(const std::string& name1, const std::string& name2,
					   const ParameterMap& map1, const ParameterMap& map2); 

    static Snapshot_Status compare_ancillary(const StringMap& map1, const StringMap& map2);     
    

  public:

    Snapshot_Source():
      m_owned(false),
      m_norm_param(0){
    }

    Snapshot_Source(const Source& src, bool owned=false):
      m_owned(false),
      m_norm_param(0){      
      latch_source(src,owned);
    }
    
    virtual ~Snapshot_Source();

    void latch_source(const Source& src, bool owned=false);

    Snapshot_Status compare(const Source& other) const;

  private:
    
    bool m_owned;
    std::string m_spatial_model;
    std::string m_spectral_model;

    optimizers::Parameter* m_norm_param;
    ParameterMap m_spatial_params;
    ParameterMap m_spectral_params;
    StringMap m_ancillary;

  };


  class Snapshot {

  public:

    /* Build from a SourceModel 
       
     */
    Snapshot(const SourceModel& model) {
      latch_model(model);
    }

    /* D'tor, deletes stuff */
    ~Snapshot(){;}
    
    void latch_model(const SourceModel& model);

    Snapshot_Status compare_source(const std::string& srcName, const Source& src) const;

    Snapshot_Status compare_model(const SourceModel& model, std::vector<std::string>& changed_sources) const;

  private:

    std::map<std::string,Snapshot_Source>  m_sources;

  };

}

#endif // Likelihood_Snapshot_h
