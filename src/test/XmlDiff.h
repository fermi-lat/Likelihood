/**
 * @file XmlDiff.h
 * @brief Class to compare two XML files.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef LikelihoodTests_XmlDiff_h
#define LikelihoodTests_XmlDiff_h

#include <map>
#include <string>
#include <vector>

#include "xml/Dom.h"

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
 * $Header$
 */

class XmlDiff {

public:

   /// @param file1 The first XML file to be compared.
   /// @param file2 The second XML file.
   /// @param tagName A tag specifying an immediate child element of 
   ///        the document root element.
   /// @param attribute An identifying attribute of the specified
   ///        child elements. It should have a unique value for each 
   ///        element.
   XmlDiff(std::string file1, std::string file2, 
           const std::string & tagName, const std::string & attribute);

   ~XmlDiff();

   bool compare();

private:

   std::string m_tagName;
   std::string m_attribute;

   std::string m_file1;
   std::string m_file2;

   typedef std::map<std::string, DomElement> DomMap;
   DomMap m_domMap1;
   DomMap m_domMap2;

   void createDomElementMap(const DomElement & rootElt, DomMap & domMap);

   void writeReserializedFile(const std::string & filename, 
                              const DomMap & domMap);

   void readLines(const std::string &inputFile, 
                  std::vector<std::string> lines);
};

#endif // LikelihoodTests_XmlDiff_h
