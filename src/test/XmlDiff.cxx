/**
 * @file XmlDiff.cxx
 * @brief Class to compare two XML files.
 * @author J. Chiang
 *
 * $Header$
 */

#include <fstream>

#include "facilities/Util.h"

#include "xml/XmlParser.h"

#include "XmlDiff.h"

XmlDiff::XmlDiff(std::string file1, std::string file2,
                 const std::string & tagName, const std::string & attribute) 
   : m_tagName(tagName), m_attribute(attribute) {
   facilities::Util::expandEnvVar(&file1);
   facilities::Util::expandEnvVar(&file2);

   xml::XmlParser * parser = new xml::XmlParser();

   DomDocument doc1 = parser->parse(file1.c_str());
   createDomElementMap(doc1.getDocumentElement(), m_domMap1);

   DomDocument doc2 = parser->parse(file2.c_str());
   createDomElementMap(doc2.getDocumentElement(), m_domMap2);

   m_file1 = file1 + "_reserialized";
   m_file2 = file2 + "_reserialized";

   writeReserializedFile(m_file1, m_domMap1);
   writeReserializedFile(m_file2, m_domMap2);

   delete parser;
}

XmlDiff::~XmlDiff() {
   std::remove(m_file1.c_str());
   std::remove(m_file2.c_str());
}

bool XmlDiff::compare() {
   std::vector<std::string> file1_lines;
   readLines(m_file1, file1_lines);
   std::vector<std::string> file2_lines;
   readLines(m_file2, file2_lines);
   if (file1_lines.size() != file2_lines.size()) return false;
   for (unsigned int i = 0; i < file1_lines.size(); i++) {
      if (file1_lines[i] != file2_lines[i]) return false;
   }
   return true;
}

void XmlDiff::createDomElementMap(const DomElement & rootElt, 
                                  DomMap & domMap) {
   std::vector<DomElement> elts;
   xml::Dom::getChildrenByTagName(rootElt, m_tagName, elts);
   domMap.clear();
   for (unsigned int i = 0; i < elts.size(); i++) {
      std::string name = xml::Dom::getAttribute(elts[i], m_attribute);
      domMap[name] = elts[i];
   }
}

void XmlDiff::writeReserializedFile(const std::string & filename, 
                                    const DomMap & domMap) {
   std::ofstream outfile(filename.c_str());
   DomMap::const_iterator it;
   for (it = domMap.begin(); it != domMap.end(); it++) {
      xml::Dom::prettyPrintElement(it->second, outfile, std::string(""));
   }
   outfile.close();
}

void XmlDiff::readLines(const std::string & inputFile, 
                        std::vector<std::string> lines) {
   std::ifstream file(inputFile.c_str());
   lines.clear();
   std::string line;
   while (std::getline(file, line, '\n')) {
      lines.push_back(line);
   }
}
