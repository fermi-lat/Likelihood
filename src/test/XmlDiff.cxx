/**
 * @file XmlDiff.cxx
 * @brief Class to compare two XML files.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/test/XmlDiff.cxx,v 1.5 2011/02/02 01:21:49 jchiang Exp $
 */

#include "XmlDiff.h"

#include <cstdio>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "facilities/Util.h"

// RapidXML includes
#include "xmlBase/rapidxml.hpp"

XmlDiff::XmlDiff(std::string file1, std::string file2,
                 const std::string& tagName, const std::string& attribute) 
   : m_tagName(tagName)
   , m_attribute(attribute)
{
   // Expand environment variables in file paths
   facilities::Util::expandEnvVar(&file1);
   facilities::Util::expandEnvVar(&file2);

   // Load and parse both documents
   if (!loadDocument(file1, m_doc1, m_buffer1)) {
      throw std::runtime_error("Failed to parse XML file: " + file1);
   }
   
   if (!loadDocument(file2, m_doc2, m_buffer2)) {
      throw std::runtime_error("Failed to parse XML file: " + file2);
   }

   // Get root elements
   XmlNode* root1 = m_doc1->first_node();
   XmlNode* root2 = m_doc2->first_node();
   
   if (!root1 || !root2) {
      throw std::runtime_error("One or both XML files have no root element");
   }

   // Create maps of elements keyed by the specified attribute
   createDomElementMap(root1, m_domMap1);
   createDomElementMap(root2, m_domMap2);

   // Create temporary filenames for reserialized content
   m_file1 = file1 + "_reserialized";
   m_file2 = file2 + "_reserialized";

   // Write reserialized files for comparison
   writeReserializedFile(m_file1, m_domMap1);
   writeReserializedFile(m_file2, m_domMap2);
}

XmlDiff::~XmlDiff() {
   // Remove temporary files
   std::remove(m_file1.c_str());
   std::remove(m_file2.c_str());
}

bool XmlDiff::loadDocument(const std::string& filename,
                           std::unique_ptr<XmlDocument>& doc,
                           std::vector<char>& buffer) {
   // Open file
   std::ifstream file(filename, std::ios::binary | std::ios::ate);
   if (!file.is_open()) {
      return false;
   }
   
   // Read file content into buffer
   auto size = file.tellg();
   file.seekg(0, std::ios::beg);
   
   buffer.resize(static_cast<size_t>(size) + 1);
   if (!file.read(buffer.data(), size)) {
      return false;
   }
   buffer[static_cast<size_t>(size)] = '\0';
   
   // Parse the document
   doc = std::make_unique<XmlDocument>();
   try {
      doc->parse<rapidxml::parse_default>(buffer.data());
      return true;
   } catch (const rapidxml::parse_error& e) {
      doc.reset();
      return false;
   }
}

std::string XmlDiff::getAttribute(XmlNode* node, const char* name) {
   if (!node) return "";
   XmlAttribute* attr = node->first_attribute(name);
   return attr ? std::string(attr->value(), attr->value_size()) : "";
}

std::string XmlDiff::getTagName(XmlNode* node) {
   if (!node || !node->name()) return "";
   return std::string(node->name(), node->name_size());
}

void XmlDiff::createDomElementMap(XmlNode* rootElt, DomMap& domMap) {
   domMap.clear();
   
   if (!rootElt) return;
   
   // Iterate through child elements matching the tag name
   for (XmlNode* child = rootElt->first_node(m_tagName.c_str()); 
        child; 
        child = child->next_sibling(m_tagName.c_str())) {
      
      // Skip non-element nodes
      if (child->type() != rapidxml::node_element) continue;
      
      // Get the attribute value to use as the map key
      std::string attrValue = getAttribute(child, m_attribute.c_str());
      if (!attrValue.empty()) {
         domMap[attrValue] = child;
      }
   }
}

void XmlDiff::prettyPrint(XmlNode* node, std::ostream& out, 
                          const std::string& indent) {
   if (!node) return;
   
   switch (node->type()) {
      case rapidxml::node_element: {
         // Opening tag
         out << indent << "<" << std::string(node->name(), node->name_size());
         
         // Attributes (sorted for consistent output)
         std::vector<std::pair<std::string, std::string>> attrs;
         for (XmlAttribute* attr = node->first_attribute(); 
              attr; 
              attr = attr->next_attribute()) {
            attrs.emplace_back(
               std::string(attr->name(), attr->name_size()),
               std::string(attr->value(), attr->value_size())
            );
         }
         std::sort(attrs.begin(), attrs.end());
         
         for (const auto& [name, value] : attrs) {
            out << " " << name << "=\"" << value << "\"";
         }
         
         // Check for children
         XmlNode* firstChild = node->first_node();
         if (!firstChild) {
            // Self-closing tag
            out << " />\n";
         } else {
            out << ">\n";
            
            // Process children
            for (XmlNode* child = firstChild; child; child = child->next_sibling()) {
               prettyPrint(child, out, indent + "  ");
            }
            
            // Closing tag
            out << indent << "</" << std::string(node->name(), node->name_size()) << ">\n";
         }
         break;
      }
      
      case rapidxml::node_data:
      case rapidxml::node_cdata: {
         // Text content (trim whitespace)
         std::string value(node->value(), node->value_size());
         // Trim leading/trailing whitespace
         size_t start = value.find_first_not_of(" \t\n\r");
         size_t end = value.find_last_not_of(" \t\n\r");
         if (start != std::string::npos && end != std::string::npos) {
            out << indent << value.substr(start, end - start + 1) << "\n";
         }
         break;
      }
      
      case rapidxml::node_comment: {
         out << indent << "<!--" << std::string(node->value(), node->value_size()) << "-->\n";
         break;
      }
      
      default:
         // Ignore other node types
         break;
   }
}

void XmlDiff::writeReserializedFile(const std::string& filename, 
                                    const DomMap& domMap) {
   std::ofstream outfile(filename);
   if (!outfile.is_open()) {
      throw std::runtime_error("Cannot open file for writing: " + filename);
   }
   
   // Write elements in sorted order (map is already sorted by key)
   for (const auto& [key, node] : domMap) {
      if (node) {
         prettyPrint(node, outfile, "");
      }
   }
   
   outfile.close();
}

void XmlDiff::readLines(const std::string& inputFile, 
                        std::vector<std::string>& lines) {
   lines.clear();
   
   std::ifstream file(inputFile);
   if (!file.is_open()) {
      throw std::runtime_error("Cannot open file for reading: " + inputFile);
   }
   
   std::string line;
   while (std::getline(file, line)) {
      lines.push_back(line);
   }
}

bool XmlDiff::compare() {
   // Read lines from both reserialized files
   std::vector<std::string> file1_lines;
   readLines(m_file1, file1_lines);
   
   std::vector<std::string> file2_lines;
   readLines(m_file2, file2_lines);
   
   // Compare number of lines
   if (file1_lines.size() != file2_lines.size()) {
      return false;
   }
   
   // Compare line by line
   for (size_t i = 0; i < file1_lines.size(); ++i) {
      if (file1_lines[i] != file2_lines[i]) {
         return false;
      }
   }
   
   return true;
}

bool XmlDiff::sameKeys() {
   if (m_domMap1.size() != m_domMap2.size()) {
      return false;
   }
   
   for (const auto& [key, node] : m_domMap1) {
      if (m_domMap2.find(key) == m_domMap2.end()) {
         return false;
      }
   }
   return true;
}

