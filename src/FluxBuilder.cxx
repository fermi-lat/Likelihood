/**
 * @file FluxBuilder.h
 * @brief Builder class for creating flux-style xml files from 
 * Likelihood::Sources.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FluxBuilder.h,v 1.4 2005/02/27 06:42:24 jchiang Exp $
 */

#ifndef Likelihood_FluxBuilder_h
#define Likelihood_FluxBuilder_h

#include <string>
#include <string_view>
#include <vector>
#include <sstream>
#include <iomanip>

#include "Likelihood/XmlBuilder.h"

// RapidXML-based XML framework (replaces Xerces-C)
#include "xmlBase/rapidxml.hpp"

namespace optimizers {
   class Function;
}

namespace Likelihood {

// Forward declarations
class Source;
class SourceModel;

/**
 * @class FluxBuilder
 * @brief This class provides methods for writing the source
 * information from Source objects as xml output appropriate for the
 * flux package.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FluxBuilder.h,v 1.4 2005/02/27 06:42:24 jchiang Exp $
 */
class FluxBuilder : public XmlBuilder {

public:
   /**
    * @brief Construct a FluxBuilder with energy range
    * @param emin Minimum energy in MeV
    * @param emax Maximum energy in MeV
    */
   FluxBuilder(double emin, double emax);

   virtual ~FluxBuilder();

   // Non-copyable
   FluxBuilder(const FluxBuilder&) = delete;
   FluxBuilder& operator=(const FluxBuilder&) = delete;

   /**
    * @brief Add all sources from a SourceModel
    * @param srcModel The source model containing sources to add
    */
   virtual void addSourceModel(SourceModel& srcModel);

   /**
    * @brief Add a single source
    * @param src The source to add
    */
   virtual void addSource(Source& src);
   
   /**
    * @brief Write the XML to a file
    * @param xmlFile Path to the output file
    */
   virtual void write(std::string xmlFile);
      
private:
   // XML nodes for the source library structure
   rapidxml::xml_node<>* m_srcLib;
   rapidxml::xml_node<>* m_allSrcsElt;

   // Energy grid
   std::vector<double> m_energies;

   /**
    * @brief Determine the source type string
    * @param src The source to classify
    * @param srcType Output string for the source type
    */
   void getSourceType(Source& src, std::string& srcType);

   /**
    * @brief Create a flux source element
    * @param src The source
    * @return Pointer to the created XML node
    */
   [[nodiscard]] rapidxml::xml_node<>* fluxSource(Source& src);

   /**
    * @brief Create a gamma spectrum element
    * @param spectrum The spectrum function
    * @return Pointer to the created XML node
    */
   [[nodiscard]] rapidxml::xml_node<>* gammaSpectrum(optimizers::Function& spectrum);

   /**
    * @brief Create a source direction element
    * @param dir The direction function
    * @return Pointer to the created XML node
    */
   [[nodiscard]] rapidxml::xml_node<>* srcDirection(optimizers::Function& dir);

   /**
    * @brief Create a solid angle element
    * @param mincos Minimum cosine value
    * @param maxcos Maximum cosine value
    * @return Pointer to the created XML node
    */
   [[nodiscard]] rapidxml::xml_node<>* solidAngle(double mincos, double maxcos);

   /**
    * @brief Create a galactic diffuse element
    * @param src The source
    * @return Pointer to the created XML node
    */
   [[nodiscard]] rapidxml::xml_node<>* galDiffuse(Source& src);

   /**
    * @brief Create a map cube source element
    * @param src The source
    * @return Pointer to the created XML node
    */
   [[nodiscard]] rapidxml::xml_node<>* mapCubeSource(Source& src);

   /**
    * @brief Build the energy grid
    * @param emin Minimum energy
    * @param emax Maximum energy
    * @param nee Number of energy bins (default 200)
    */
   void makeEnergyGrid(double emin, double emax, unsigned int nee = 200);

   /**
    * @brief Replace spaces with underscores in a name
    * @param name The name to modify (in-place)
    */
   void addUnderscores(std::string& name);

   // ==================== RapidXML Helper Methods ====================

   /**
    * @brief Create an XML element
    * @param name Element name
    * @return Pointer to the created node
    */
   [[nodiscard]] rapidxml::xml_node<>* createElement(const char* name);

   /**
    * @brief Add a string attribute to a node
    * @param node The node to add the attribute to
    * @param name Attribute name
    * @param value Attribute value
    */
   void addAttribute(rapidxml::xml_node<>* node, 
                     const char* name, 
                     const std::string& value);

   /**
    * @brief Add a double attribute to a node
    * @param node The node to add the attribute to
    * @param name Attribute name
    * @param value Attribute value
    * @param precision Decimal precision (default 10)
    */
   void addAttribute(rapidxml::xml_node<>* node, 
                     const char* name, 
                     double value,
                     int precision = 10);

   /**
    * @brief Append a child node to a parent
    * @param parent Parent node
    * @param child Child node to append
    */
   void appendChild(rapidxml::xml_node<>* parent, 
                    rapidxml::xml_node<>* child);

   /**
    * @brief Pretty print an element to an output stream
    * @param node The node to print
    * @param out The output stream
    * @param indent Current indentation string
    */
   void prettyPrintElement(rapidxml::xml_node<>* node,
                           std::ostream& out,
                           const std::string& indent = "");
};

} // namespace Likelihood

#endif // Likelihood_FluxBuilder_h
