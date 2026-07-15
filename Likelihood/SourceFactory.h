/** 
 * @file SourceFactory.h
 * @brief Declaration of SourceFactory class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceFactory.h,v 1.28 2015/02/28 03:27:30 jchiang Exp $
 */

#ifndef Likelihood_SourceFactory_h
#define Likelihood_SourceFactory_h

#include <map>
#include <memory>
#include <string>
#include <vector>

#ifndef SWIG
// RapidXML-based XML framework (replaces Xerces-C)
#include "xmlBase/rapidxml.hpp"
#include <string_view>
#endif

#include "Likelihood/Exception.h"
#include "Likelihood/Observation.h"
#include "Likelihood/Source.h"

namespace st_stream {
   class StreamFormatter;
}

namespace optimizers {
   class Function;
   class FunctionFactory;
}

namespace Likelihood {

/** 
 * @class SourceFactory
 *
 * @brief This class implements the Prototype pattern to return
 * clones of various gamma-ray Sources.
 *
 * The design of this class is based on the Factory template class
 * of Hippodraw.
 *
 */
class SourceFactory {

public:
   /**
    * @brief Construct a SourceFactory
    * @param observation Reference to the observation data
    * @param verbose Enable verbose output
    */
   SourceFactory(const Observation& observation, bool verbose = false);

   ~SourceFactory();

#ifndef SWIG
   // Non-copyable - hidden from SWIG as it doesn't handle = delete well
   SourceFactory(const SourceFactory&) = delete;
   SourceFactory& operator=(const SourceFactory&) = delete;
#endif

   /**
    * @brief Create a source by name (returns a clone)
    * @param name The source name
    * @return Pointer to cloned source (caller takes ownership)
    */
   Source* create(const std::string& name);

   /**
    * @brief Release ownership of a source
    * @param name The source name
    * @return Pointer to the source (caller takes ownership)
    */
   Source* releaseSource(const std::string& name);

   /**
    * @brief Add a source prototype
    * @param name The name to register the source under
    * @param src Pointer to the source
    * @param fromClone If true, stores a clone; if false, takes ownership
    * @note Clients should almost always have fromClone = true; otherwise,
    *       the destructor will delete their Source, rather than a clone.
    */
   void addSource(const std::string& name, Source* src, bool fromClone = true);

   /**
    * @brief Replace an existing source prototype
    * @param src Pointer to the replacement source
    * @param fromClone If true, stores a clone; if false, takes ownership
    */
   void replaceSource(Source* src, bool fromClone = true);

   /**
    * @brief Read source definitions from an XML file
    * @param xmlFile Path to the XML file
    * @param funcFactory Function factory for creating spectral functions
    * @param requireExposure Whether exposure calculation is required
    * @param addPointSources Whether to add point sources
    * @param loadMaps Whether to load spatial maps
    */
   void readXml(const std::string& xmlFile,
                optimizers::FunctionFactory& funcFactory,
                bool requireExposure = true,
                bool addPointSources = true,
                bool loadMaps = true);

   /**
    * @brief Get all registered source names
    * @param srcNames Output vector to receive source names
    */
   void fetchSrcNames(std::vector<std::string>& srcNames) const;

   /**
    * @brief Get all registered source names (C++17 style)
    * @return Vector of source names
    */
   std::vector<std::string> fetchSrcNames() const;

#ifndef SWIG
   // ==================== RapidXML-specific interface ====================
   // Hidden from SWIG - these use RapidXML types

   /**
    * @brief Read source definitions from an XML node
    */
   void readXml(rapidxml::xml_node<>* source_library,
                const std::string& xmlFile,
                optimizers::FunctionFactory& funcFactory,
                bool requireExposure = true,
                bool addPointSources = true,
                bool loadMaps = true);

private:
   // XML document and buffer for parsing
   mutable rapidxml::xml_document<> m_xmlDoc;
   mutable std::vector<char> m_xmlBuffer;

   /**
    * @brief Create sources from XML source_library element
    */
   void makeSources(const std::string& xmlFile,
                    const rapidxml::xml_node<>* source_library,
                    std::vector<Source*>& sources,
                    optimizers::FunctionFactory& funcFactory,
                    bool requireExposure = true,
                    bool addPointSources = true,
                    bool loadMaps = true);

   /**
    * @brief Create a point source from XML elements
    */
   Source* makePointSource(const rapidxml::xml_node<>* spectrum,
                           const rapidxml::xml_node<>* spatialModel,
                           optimizers::FunctionFactory& funcFactory);

   /**
    * @brief Create a diffuse source from XML elements
    */
   Source* makeDiffuseSource(const rapidxml::xml_node<>* spectrum,
                             const rapidxml::xml_node<>* spatialModel,
                             optimizers::FunctionFactory& funcFactory,
                             bool loadMap = true);

   /**
    * @brief Create a composite source from XML elements
    */
   Source* makeCompositeSource(const std::string& xmlFile,
                               const rapidxml::xml_node<>* spectrum,
                               rapidxml::xml_node<>* source_library,
                               optimizers::FunctionFactory& funcFactory,
                               bool requireExposure = true,
                               bool addPointSources = true,
                               bool loadMap = true);

   /**
    * @brief Set the spectrum for a source from XML
    */
   void setSpectrum(Source* src, 
                    const rapidxml::xml_node<>* spectrum,
                    optimizers::FunctionFactory& funcFactory);

   /**
    * @brief Add parameters to MultipleBrokenPowerLaw function
    */
   void addParamsToMultipleBPL(optimizers::Function* spec,
                               const std::vector<rapidxml::xml_node<>*>& params,
                               const Source* src) const;

   /**
    * @brief Add parameters to PiecewisePowerLaw function
    */
   void addParamsToPiecewisePL(optimizers::Function* spec,
                               const std::vector<rapidxml::xml_node<>*>& params,
                               const Source* src) const;

   /**
    * @brief Get attribute value from XML node
    */
   static std::string getAttributeValue(
       const rapidxml::xml_node<>* node,
       const char* attrName,
       const std::string& defaultValue);

   /**
    * @brief Get attribute value as double from XML node
    */
   static double getAttributeValueAsDouble(
       const rapidxml::xml_node<>* node,
       const char* attrName,
       double defaultValue);

   /**
    * @brief Collect child elements with a specific name
    */
   static std::vector<rapidxml::xml_node<>*> collectChildren(
       const rapidxml::xml_node<>* parent,
       const char* childName);

   /**
    * @brief Check if a source is within acceptable ROI distance
    */
   void checkRoiDist(double ra, double dec) const;

#else
   // When SWIG is processing, provide minimal private section
private:
   void checkRoiDist(double ra, double dec) const;
#endif // SWIG

   // Members visible to both SWIG and C++ compiler
   // (but SWIG doesn't need to wrap private members)
private:
   bool m_verbose;
   std::map<std::string, Source*> m_prototypes;
   bool m_requireExposure;
   const Observation& m_observation;
   st_stream::StreamFormatter* m_formatter;
   std::string m_currentSrcName;
};

} // namespace Likelihood

#endif // Likelihood_SourceFactory_h
