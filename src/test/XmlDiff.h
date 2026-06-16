/**
 * @file XmlDiff.h
 * @brief Class to compare two XML files.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/XmlDiff.h,v 1.2 2004/11/11 04:32:24 jchiang Exp $
 */

#ifndef LikelihoodTests_XmlDiff_h
#define LikelihoodTests_XmlDiff_h

#include <map>
#include <memory>
#include <string>
#include <vector>

// RapidXML-based framework includes
#include "xmlBase/rapidxml.hpp"

/**
 * @class XmlDiff
 *
 * @brief A class for comparing two XML files.  Since the order of
 * appearance of tags in an XML file is arbitrary, this class
 * reserializes the selected tags, which must be immediate children of
 * the root element, writing the output to temporary files for a
 * subsequent line-by-line comparison.  It uses a map, keyed by a
 * common attribute of the selected tags, to determine a unique
 * ordering for the serialization.  To ensure a thorough comparison,
 * the selected attribute should have a value unique to each child
 * element.
 *
 * @bug This class tacitly assumes that the ordering of descendents
 * for each child element of the root element is the same for both
 * files.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/XmlDiff.h,v 1.2 2004/11/11 04:32:24 jchiang Exp $
 */

class XmlDiff {

public:
   // Type aliases for RapidXML types
   using XmlNode = rapidxml::xml_node<char>;
   using XmlDocument = rapidxml::xml_document<char>;
   using XmlAttribute = rapidxml::xml_attribute<char>;

   /// @param file1 The first XML file to be compared.
   /// @param file2 The second XML file.
   /// @param tagName A tag specifying an immediate child element of 
   ///        the document root element.
   /// @param attribute An identifying attribute of the specified
   ///        child elements. It should have a unique value for each 
   ///        element.
   XmlDiff(std::string file1, std::string file2, 
           const std::string& tagName, const std::string& attribute);

   ~XmlDiff();

   /// Compare the two XML files
   /// @return true if files are equivalent, false otherwise
   [[nodiscard]] bool compare();
   [[nodiscard]] bool sameKeys();
  
private:
   std::string m_tagName;
   std::string m_attribute;

   std::string m_file1;
   std::string m_file2;

   // Document storage (RapidXML requires buffers to persist during document lifetime)
   std::unique_ptr<XmlDocument> m_doc1;
   std::unique_ptr<XmlDocument> m_doc2;
   std::vector<char> m_buffer1;
   std::vector<char> m_buffer2;

   using DomMap = std::map<std::string, XmlNode*>;
   DomMap m_domMap1;
   DomMap m_domMap2;

   /// Load and parse an XML file
   /// @param filename Path to the XML file
   /// @param doc Output document (will be populated)
   /// @param buffer Output buffer (must persist during document lifetime)
   /// @return true on success, false on failure
   bool loadDocument(const std::string& filename,
                     std::unique_ptr<XmlDocument>& doc,
                     std::vector<char>& buffer);

   /// Create a map of elements keyed by the specified attribute
   /// @param rootElt Root element to search
   /// @param domMap Output map to populate
   void createDomElementMap(XmlNode* rootElt, DomMap& domMap);

   /// Write elements from a map to a file in sorted order
   /// @param filename Output file path
   /// @param domMap Map of elements to write
   void writeReserializedFile(const std::string& filename, 
                              const DomMap& domMap);

   /// Read all lines from a file
   /// @param inputFile Input file path
   /// @param lines Output vector of lines
   void readLines(const std::string& inputFile, 
                  std::vector<std::string>& lines);

   // Helper methods for XML operations
   static std::string getAttribute(XmlNode* node, const char* name);
   static std::string getTagName(XmlNode* node);
   static void prettyPrint(XmlNode* node, std::ostream& out, 
                           const std::string& indent = "");
};

#endif // LikelihoodTests_XmlDiff_h
