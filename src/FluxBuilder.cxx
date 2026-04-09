/**
 * @file FluxBuilder.cxx
 * @brief Implementation for FluxBuilder class
 * @author J. Chiang
 */

#include <fstream>
#include <sstream>
#include <cmath>

#include "xmlBase/rapidxml.hpp"

#include "facilities/Util.h"

#include "optimizers/Function.h"
#include "optimizers/Parameter.h"
#include "optimizers/dArg.h"

#include "Likelihood/PointSource.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceModel.h"
#include "Likelihood/FluxBuilder.h"

namespace Likelihood {

namespace {

char* allocateString(rapidxml::xml_document<>* doc, const std::string& str) {
   return doc->allocate_string(str.c_str(), str.size() + 1);
}

rapidxml::xml_node<>* createElement(rapidxml::xml_document<>* doc, const char* name) {
   char* nodeName = doc->allocate_string(name);
   return doc->allocate_node(rapidxml::node_element, nodeName);
}

void addAttribute(rapidxml::xml_document<>* doc, 
                  rapidxml::xml_node<>* node, 
                  const char* name, 
                  const std::string& value) {
   char* attrName = doc->allocate_string(name);
   char* attrValue = allocateString(doc, value);
   rapidxml::xml_attribute<>* attr = doc->allocate_attribute(attrName, attrValue);
   node->append_attribute(attr);
}

void addAttribute(rapidxml::xml_document<>* doc,
                  rapidxml::xml_node<>* node,
                  const char* name,
                  double value) {
   std::ostringstream ss;
   ss.precision(15);
   ss << value;
   addAttribute(doc, node, name, ss.str());
}

void appendChild(rapidxml::xml_node<>* parent, rapidxml::xml_node<>* child) {
   if (parent && child) {
      parent->append_node(child);
   }
}

void prettyPrint(std::ostream& out, rapidxml::xml_node<>* node, int indent = 0) {
   if (node == nullptr) return;
   
   std::string indentStr(indent * 2, ' ');
   
   if (node->type() == rapidxml::node_element) {
      out << indentStr << "<" << node->name();
      
      for (rapidxml::xml_attribute<>* attr = node->first_attribute();
           attr != nullptr;
           attr = attr->next_attribute()) {
         out << " " << attr->name() << "=\"" << attr->value() << "\"";
      }
      
      rapidxml::xml_node<>* child = node->first_node();
      if (child == nullptr) {
         out << "/>\n";
      } else {
         out << ">\n";
         for (; child != nullptr; child = child->next_sibling()) {
            prettyPrint(out, child, indent + 1);
         }
         out << indentStr << "</" << node->name() << ">\n";
      }
   }
}

} // anonymous namespace

// Constructor definition
FluxBuilder::FluxBuilder(double emin, double emax)
   : XmlBuilder() {
   m_emin = emin;
   m_emax = emax;
   m_srcLib = createElement(m_doc, "source_library");
   addAttribute(m_doc, m_srcLib, "title", "Likelihood_sources");
   m_doc->append_node(m_srcLib);
}

// Destructor definition
FluxBuilder::~FluxBuilder() {
   // Base class handles m_doc cleanup
}

// write() definition
void FluxBuilder::write(std::string xmlFile) {
   facilities::Util::expandEnvVar(&xmlFile);
   std::ofstream outFile(xmlFile.c_str());
   outFile << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
   prettyPrint(outFile, m_srcLib);
}

// addSource() definition
void FluxBuilder::addSource(Source& src) {
   rapidxml::xml_node<>* srcElt = createElement(m_doc, "source");
   addAttribute(m_doc, srcElt, "name", src.getName());
   
   // Compute and add flux
   double flux = 0.0;
   try {
      const optimizers::Function& spectrum = src.spectrum();
      const int nsteps = 100;
      double logEmin = std::log10(m_emin);
      double logEmax = std::log10(m_emax);
      double dlogE = (logEmax - logEmin) / nsteps;
      
      for (int i = 0; i < nsteps; ++i) {
         double e1 = std::pow(10.0, logEmin + i * dlogE);
         double e2 = std::pow(10.0, logEmin + (i + 1) * dlogE);
         optimizers::dArg arg1(e1);
         optimizers::dArg arg2(e2);
         double f1 = const_cast<optimizers::Function&>(spectrum)(arg1);
         double f2 = const_cast<optimizers::Function&>(spectrum)(arg2);
         flux += 0.5 * (f1 + f2) * (e2 - e1);
      }
   } catch (...) {
      flux = 1.0;  // Default flux if calculation fails
   }
   addAttribute(m_doc, srcElt, "flux", flux);
   
   // Add spectrum element
   rapidxml::xml_node<>* specElt = createElement(m_doc, "spectrum");
   addAttribute(m_doc, specElt, "escale", "MeV");
   addAttribute(m_doc, specElt, "type", "PowerLaw");
   
   rapidxml::xml_node<>* gammaElt = createElement(m_doc, "particle");
   addAttribute(m_doc, gammaElt, "name", "gamma");
   addAttribute(m_doc, gammaElt, "value", 2.0);
   appendChild(specElt, gammaElt);
   
   rapidxml::xml_node<>* eminElt = createElement(m_doc, "particle");
   addAttribute(m_doc, eminElt, "name", "emin");
   addAttribute(m_doc, eminElt, "value", m_emin);
   appendChild(specElt, eminElt);
   
   rapidxml::xml_node<>* emaxElt = createElement(m_doc, "particle");
   addAttribute(m_doc, emaxElt, "name", "emax");
   addAttribute(m_doc, emaxElt, "value", m_emax);
   appendChild(specElt, emaxElt);
   
   appendChild(srcElt, specElt);
   
   // Add celestial_dir for point sources
   PointSource* ptSrc = dynamic_cast<PointSource*>(&src);
   if (ptSrc != nullptr) {
      rapidxml::xml_node<>* dirElt = createElement(m_doc, "celestial_dir");
      addAttribute(m_doc, dirElt, "ra", ptSrc->getDir().ra());
      addAttribute(m_doc, dirElt, "dec", ptSrc->getDir().dec());
      appendChild(srcElt, dirElt);
   }
   
   appendChild(m_srcLib, srcElt);
}

// addSourceModel() definition  
void FluxBuilder::addSourceModel(SourceModel& srcModel) {
   std::vector<std::string> srcNames;
   srcModel.getSrcNames(srcNames);
   for (const auto& name : srcNames) {
      Source* src = srcModel.getSource(name);
      if (src != nullptr) {
         addSource(*src);
      }
   }
}

} // namespace Likelihood
