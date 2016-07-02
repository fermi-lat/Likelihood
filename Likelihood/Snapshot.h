/**
 * @file Snapshot.h
 * @author Eric Charles
 *
 * A class to provide a Snapshot of a SourceModel and allow for comparisons
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Snapshot.h,v 1.3 2016/07/01 23:37:00 echarles Exp $
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
    /* Encapsulate the status changes of a model as a bit mask
       and provide bit manipulation tools */

  public:
    enum { 	
      Unchanged = 0,                    //! No change
      Norm_Offset = 0,                  //! Offset for normalization parameter status
      Spec_Offset = 5,                  //! Offset for spectral parameter status
      Spatial_Offset = 10,              //! Offset for spatial parameter status
      Spec_Model_Offset = 15,           //! Offset for spectral model status
      Spatial_Model_Offset = 18,        //! Offset for spatial model status
      Ancillary_Offset = 21,            //! Offset for ancillary status
      Fixed_Source_Added_Offset = 24,   //! Offset for marking that a fixed source was added
      Free_Source_Added_Offset = 25,    //! Offset for marking that a free source was added
      Fixed_Source_Removed_Offset = 26, //! Offset for marking that a fixed source was removed
      Free_Source_Removed_Offset = 27,  //! Offset for marking that a free source was removed
      Source_Freed_Offset = 28,         //! Offset for marking that a source was free
      Source_Fixed_Offset = 29,         //! Offset for marking that a source was fixed

      Value_Changed = 0x1,              //! Parameter value changed
      Status_Changed = 0x2,             //! Parameter status chagned
      Prior_Changed = 0x4,              //! Parameter prior changed
      Fixed_Changed = 0x8,              //! A fixed parameter was changed
      Free_Changed = 0x10,              //! A free parameter was changed
      Par_Changed_Mask = 0x7,           //! Something changed
      
      Model_Type_Changed = 0x1,         //! A model changed type
      Free_Model_Type_Changed = 0x2,    //! A free model changed type
      Fixed_Model_Type_Changed = 0x4,   //! A fixed model changed type

      Norm_Value_Changed = Value_Changed << Norm_Offset,
      Norm_Status_Changed = Status_Changed << Norm_Offset,
      Norm_Prior_Changed = Prior_Changed << Norm_Offset,
      Norm_Fixed_Changed = Fixed_Changed << Norm_Offset,
      Norm_Free_Changed = Free_Changed << Norm_Offset,

      Spec_Value_Changed = Value_Changed << Spec_Offset,
      Spec_Status_Changed = Status_Changed << Spec_Offset,
      Spec_Prior_Changed = Prior_Changed << Spec_Offset,
      Spec_Fixed_Changed = Fixed_Changed << Spec_Offset,
      Spec_Free_Changed = Free_Changed << Spec_Offset,

      Spatial_Value_Changed = Value_Changed << Spatial_Offset,
      Spatial_Status_Changed = Status_Changed << Spatial_Offset,
      Spatial_Prior_Changed = Prior_Changed << Spatial_Offset,
      Spatial_Fixed_Changed = Fixed_Changed << Spatial_Offset,
      Spatial_Free_Changed = Free_Changed << Spatial_Offset,

      Spec_Model_Type_Changed = Model_Type_Changed << Spec_Model_Offset,
      Free_Spec_Model_Type_Changed = Free_Model_Type_Changed << Spec_Model_Offset,
      Fixed_Spec_Model_Type_Changed = Fixed_Model_Type_Changed << Spec_Model_Offset,

      Spatial_Model_Type_Changed = Model_Type_Changed << Spatial_Model_Offset,
      Free_Spatial_Model_Type_Changed = Free_Model_Type_Changed << Spatial_Model_Offset,
      Fixed_Spatial_Model_Type_Changed = Fixed_Model_Type_Changed << Spatial_Model_Offset,

      Ancillary_Changed = Model_Type_Changed << Ancillary_Offset,
      Free_Ancillary_Changed = Free_Model_Type_Changed << Ancillary_Offset,
      Fixed_Ancillary_Changed = Fixed_Model_Type_Changed << Ancillary_Offset,
   
      Fixed_Source_Added = Value_Changed << Fixed_Source_Added_Offset,
      Free_Source_Added = Value_Changed << Free_Source_Added_Offset,
      Fixed_Source_Removed = Value_Changed << Fixed_Source_Removed_Offset,
      Free_Source_Removed = Value_Changed << Free_Source_Removed_Offset,  
      
      Source_Freed = Value_Changed << Source_Freed_Offset,
      Source_Fixed = Value_Changed << Source_Fixed_Offset,  

      Fixed_Source_Changed = Fixed_Source_Added | Fixed_Source_Removed,
      Free_Source_Changed = Free_Source_Added | Free_Source_Removed,

      Souce_Added_Mask = Fixed_Source_Added | Free_Source_Added,       //! Mask for seeing if any source have been added
      Souce_Removed_Mask = Fixed_Source_Removed | Free_Source_Removed, //! Mask for seeing if any source have been removed    

      Model_Free_Mask = Spec_Free_Changed | Spatial_Free_Changed,     
      Model_Fixed_Mask = Spec_Fixed_Changed | Spatial_Fixed_Changed,     

      Norm_Mask = Par_Changed_Mask << Norm_Offset,        //! Mask for normalization parameter status      
      Spec_Mask = Par_Changed_Mask << Spec_Offset,        //! Mask for spectral parameter status
      Spatial_Mask =  Par_Changed_Mask << Spatial_Offset, //! Mask for spatial parameter status
      Par_Changed =                                       //! Mask for any parameter changed
      Norm_Mask | Spec_Mask | Spatial_Mask,
      Free_Par_Changed =                                  //! Mask for any fixed parameter changed 
      Norm_Free_Changed | Spec_Free_Changed | Spatial_Free_Changed,
      Fixed_Par_Changed =                                 //! Mask for any fixed parameter changed 
      Norm_Fixed_Changed | Spec_Fixed_Changed | Spatial_Fixed_Changed,
      Model_Type_Mask =                                   //! Mask for changes to model type
      Spec_Model_Type_Changed | Spatial_Model_Type_Changed | Ancillary_Changed,
      Free_Model_Type_Mask =                              //! Mask for any free model changed type 
      Free_Spec_Model_Type_Changed | Free_Spatial_Model_Type_Changed | Free_Ancillary_Changed,
      Fixed_Model_Type_Mask =                             //! Mask for any fixed model changed type
      Fixed_Spec_Model_Type_Changed | Fixed_Spatial_Model_Type_Changed | Fixed_Ancillary_Changed,

      Free_Changed_Mask  =                                //! Any free source changed
      Free_Par_Changed | Free_Model_Type_Changed | Free_Source_Changed, 
      Fixed_Changed_Mask  =                               //! Any fixed source changed
      Fixed_Par_Changed | Fixed_Model_Type_Mask | Fixed_Source_Changed,
    } Flags;
    
    Snapshot_Status(unsigned val=Unchanged) :m_val(val) {}
    ~Snapshot_Status(){ }
    
    inline bool unchanged() const { return m_val == 0; }
    inline bool changed() const { return m_val != 0; }

    inline bool norm_changed() const { return (m_val & Norm_Mask) != 0; }
    inline bool norm_value_changed() const { return (m_val & Norm_Value_Changed) != 0; }
    inline bool norm_status_changed() const { return (m_val & Norm_Status_Changed) != 0; }
    inline bool norm_prior_changed() const { return (m_val & Norm_Prior_Changed) != 0; }
    inline bool norm_fixed_changed() const { return (m_val & Norm_Fixed_Changed) != 0; }
    inline bool norm_free_changed() const { return (m_val & Norm_Free_Changed) != 0; }

    inline bool spec_changed() const { return (m_val & Spec_Mask) != 0; }
    inline bool spec_free_changed() const { return (m_val & Spec_Free_Changed) != 0; }
    inline bool spec_fixed_changed() const { return (m_val & Spec_Fixed_Changed) != 0; }    

    inline bool spatial_changed() const { return (m_val & Spatial_Mask) != 0; }
    inline bool spatial_free_changed() const { return (m_val & Spatial_Free_Changed) != 0; }
    inline bool spatial_fixed_changed() const { return (m_val & Spatial_Fixed_Changed) != 0; }

    inline bool spec_model_type_changed() const { return (m_val & Spec_Model_Type_Changed) != 0; }
    inline bool spatial_model_type_changed() const { return (m_val & Spatial_Model_Type_Changed) != 0; }

    inline bool model_type_changed() const { return (m_val & Model_Type_Mask) != 0; }
    inline bool ancillary_changed() const { return (m_val & Ancillary_Changed) != 0; }
    inline bool model_free_changed() const { return (m_val & Model_Free_Mask) != 0; }    
    inline bool model_fixed_changed() const { return (m_val & Model_Fixed_Mask) != 0; }    
    inline bool free_source_changed() const { return (m_val & Free_Source_Changed) != 0; }
    inline bool free_source_added() const { return (m_val & Free_Source_Removed) != 0; }
    inline bool free_source_removed() const { return (m_val & Free_Source_Added) != 0; }
    inline bool fixed_source_changed() const { return (m_val & Fixed_Source_Changed) != 0; }
    inline bool fixed_source_added() const { return (m_val & Fixed_Source_Removed) != 0; }
    inline bool fixed_source_removed() const { return (m_val & Fixed_Source_Added) != 0; }     
    inline bool fixed_changed() const { return (m_val & Fixed_Changed_Mask) != 0; }
    
    inline bool source_freed() const { return (m_val & Source_Freed) != 0; }     
    inline bool source_fixed() const { return (m_val & Source_Fixed) != 0; }

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

    /* delete all the owned parameters in a ParameterMap */
    static void delete_map(ParameterMap& theMap);

    /* deep copy all the owned parameters in a ParameterMap */
    static void deep_copy_param_map(ParameterMap& theMap, bool delete_old = false);
    

    /* Compare the Priors on two Parameters.
       Returns true if they have changed */
    static bool compare_priors(const optimizers::Parameter& p1, 
			       const optimizers::Parameter& p2);

    /* Compare two Parameters.
       Returns a 4 bit-mask
       bit 0 : value has changed
       bit 1 : status (fixed/free) has changed
       bit 2 : prior has changed
       bit 3 : parameter was previously fixed
    */
    static unsigned parameter_status(const optimizers::Parameter& p1, 
				     const optimizers::Parameter& p2);

    /* Extract the normalization parameter from a Source */
    static void extract_norm(const Source& src,
			     optimizers::Parameter*& theParam,
			     bool clone = false);

    /* Fill a ParameterMap with the spectral parameters of a source 
       This does not include the normalization parameter */
    static void extract_spectral(const Source& src,
				 std::string& model_name,
				 ParameterMap& theMap,
				 bool clone = false);

    /* Fill a ParameterMap with the spatial parameters of a source */
    static void extract_spatial(const Source& src,
				std::string& model_name,
				ParameterMap& theMap,
				bool clone = false);
    
    /* Fill a StringMap with the ancillary data of a source */
    static void extract_ancillary(const Source& src,
				  StringMap& theMap);    
    
    /* Compare the normalization parameters */    
    static Snapshot_Status compare_norm(const optimizers::Parameter& par1, 
					const optimizers::Parameter& par2,
					bool fixed);
    
    /* Compare the spectral parameters */       
    static Snapshot_Status compare_spectral(const std::string& name1, const std::string& name2,
					    const ParameterMap& map1, const ParameterMap& map2,
					    bool fixed); 
    
    /* Compare the spatial parameters */       
    static Snapshot_Status compare_spatial(const std::string& name1, const std::string& name2,
					   const ParameterMap& map1, const ParameterMap& map2,
					   bool fixed); 
    
    /* Compare the ancillary data */
    static Snapshot_Status compare_ancillary(const StringMap& map1, const StringMap& map2,
					     bool fixed);         

    /* print a parameter */
    static void print_param(const optimizers::Parameter& par, const std::string& indent);

    /* print a parameter map */
    static void print_param_map(const ParameterMap& parMap, const std::string& indent);

    /* print a string map */
    static void print_string_map(const StringMap& strMap, const std::string& indent);
    

  public:

    /* c'tor, does nothing */
    Snapshot_Source()
      :m_owned(false),
       m_fixed(false),
       m_norm_param(0){
    }

    /* c'tor, build from a source
       
       if owned is true this will clone the Parameters out of the source
    */
    Snapshot_Source(const Source& src, bool owned=false)
      :m_owned(false),
       m_fixed(false),
       m_norm_param(0){      
      latch_source(src,owned);
    }
    
    /* D'tor, 
       
       if owned is true this will delete the owned Paramters 
    */
    virtual ~Snapshot_Source();


    /* Latch parameters from a source
       
       if owned is true this will clone the Parameters out of the source
    */
    void latch_source(const Source& src, bool owned=false);

    /* Deep copy all of the Parameters */
    void deep_copy();

    /* Compare this snapshot to the current status of a source */
    Snapshot_Status compare(const Source& other) const;

    /* print this Snapshot_Source */
    void print() const;

  private:
    
    bool m_owned;
    bool m_fixed;

    std::string m_spatial_model;
    std::string m_spectral_model;

    optimizers::Parameter* m_norm_param;
    ParameterMap m_spatial_params;
    ParameterMap m_spectral_params;
    StringMap m_ancillary;

  };


  class Snapshot {
 
    /* A utility class to cache a Snapshot of a SourceModel, and to allow for
       comparisions
    */

  public:

    /* Check to see if a source is 'free'
    
       For now this only check the normalization 
    */
    static bool source_free(const Source& src);

  public:

    /* Build from a SourceModel 

       If owned is true this will clone the Parameters       
     */
    Snapshot(const SourceModel& model, bool owned = false) {
      latch_model(model,owned);
    }

    /* D'tor, deletes stuff */
    ~Snapshot(){;}
    
    
    /* Latch parameters from a SourceModel
       
       If owned is true this will clone the Parameters
    */
    void latch_model(const SourceModel& model, bool owned = false);

    /* Print the entire model */
    void print() const;

    /* Print a single source */
    void print_source(const std::string& srcName) const;

    /* Compare a single source with this Snapshot */
    Snapshot_Status compare_source(const std::string& srcName, const Source& src) const;

    /* Compare an entire model with this Snapshot */    
    Snapshot_Status compare_model(const SourceModel& model, std::vector<std::string>& changed_sources) const;

  private:

    std::map<std::string,Snapshot_Source>  m_sources;

  };

}

#endif // Likelihood_Snapshot_h
