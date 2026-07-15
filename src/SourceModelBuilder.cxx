/**
 * @file SourceModelBuilder.cpp
 * @brief Implementation for class to provide methods to write xml
 * files for the Likelihood package source models.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceModelBuilder.cxx,v 1.21 2015/03/03 21:41:31 jchiang Exp $
 */

#include <fstream>
#include <sstream>

#include "xmlBase/rapidxml.hpp"
#include "xmlBase/xml_printer.hpp"

#include "facilities/Util.h"

#include "Likelihood/CompositeSource.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/DMFitFunction.h"
#include "Likelihood/FileFunction.h"
#include "Likelihood/MapCubeFunction2.h"
#include "Likelihood/RadialProfile.h"
#include "Likelihood/ScaleFactor.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceModel.h"
#include "Likelihood/SourceModelBuilder.h"

namespace Likelihood {

namespace {

// Helper function to allocate a string in the document's memory pool
char* allocateString(rapidxml::xml_document<>* doc, const std::string& str) {
   return doc->allocate_string(str.c_str(), str.size() + 1);
}

// Helper function to create an element node
rapidxml::xml_node<>* createElement(rapidxml::xml_document<>* doc, const char* name) {
   char* nodeName = doc->allocate_string(name);
   return doc->allocate_node(rapidxml::node_element, nodeName);
}

// Helper function to add an attribute to a node
void addAttribute(rapidxml::xml_document<>* doc, 
                  rapidxml::xml_node<>* node, 
                  const char* name, 
                  const std::string& value) {
   char* attrName = doc->allocate_string(name);
   char* attrValue = allocateString(doc, value);
   rapidxml::xml_attribute<>* attr = doc->allocate_attribute(attrName, attrValue);
   node->append_attribute(attr);
}

// Helper function to add an attribute with double value
void addAttribute(rapidxml::xml_document<>* doc,
                  rapidxml::xml_node<>* node,
                  const char* name,
                  double value) {
   std::ostringstream ss;
   ss.precision(15);
   ss << value;
   addAttribute(doc, node, name, ss.str());
}

// Helper function to add an attribute with int value
void addAttribute(rapidxml::xml_document<>* doc,
                  rapidxml::xml_node<>* node,
                  const char* name,
                  int value) {
   addAttribute(doc, node, name, std::to_string(value));
}

// Helper function to add an attribute with bool value
void addAttribute(rapidxml::xml_document<>* doc,
                  rapidxml::xml_node<>* node,
                  const char* name,
                  bool value) {
   addAttribute(doc, node, name, value ? std::string("true") : std::string("false"));
}

// Helper function to add an attribute with char array value
void addAttribute(rapidxml::xml_document<>* doc,
                  rapidxml::xml_node<>* node,
                  const char* name,
                  const char* value) {
   addAttribute(doc, node, name, std::string(value));
}
  
 
// Helper function to append a child node
void appendChild(rapidxml::xml_node<>* parent, rapidxml::xml_node<>* child) {
   if (parent && child) {
      parent->append_node(child);
   }
}

// Helper to write parameter nodes for a Function
void appendParameterNodes(rapidxml::xml_document<>* doc,
                          rapidxml::xml_node<>* parent,
                          optimizers::Function& func) {
   std::vector<optimizers::Parameter> params;
   func.getParams(params);
   
   for (const auto& param : params) {
      rapidxml::xml_node<>* paramNode = createElement(doc, "parameter");
      addAttribute(doc, paramNode, "name", param.getName());
      addAttribute(doc, paramNode, "value", param.getValue());
      addAttribute(doc, paramNode, "scale", param.getScale());
      addAttribute(doc, paramNode, "min", param.getBounds().first);
      addAttribute(doc, paramNode, "max", param.getBounds().second);
      addAttribute(doc, paramNode, "free", param.isFree() ? 1 : 0);
      appendChild(parent, paramNode);
   }
}

// Pretty print helper
void prettyPrint(std::ostream& out, rapidxml::xml_node<>* node, int indent = 0) {
   if (node == nullptr) return;
   
   std::string indentStr(indent * 2, ' ');
   
   if (node->type() == rapidxml::node_element) {
      out << indentStr << "<" << node->name();
      
      // Print attributes
      for (rapidxml::xml_attribute<>* attr = node->first_attribute();
           attr != nullptr;
           attr = attr->next_attribute()) {
         out << " " << attr->name() << "=\"" << attr->value() << "\"";
      }
      
      // Check for children
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

SourceModelBuilder::SourceModelBuilder(const std::string &functionLibrary,
                                       const std::string &srcLibTitle) 
   : XmlBuilder() {
   m_srcLib = createElement(m_doc, "source_library");
   if (!functionLibrary.empty()) {
      addAttribute(m_doc, m_srcLib, "function_library", functionLibrary);
   }
   addAttribute(m_doc, m_srcLib, "title", srcLibTitle);
   m_doc->append_node(m_srcLib);
}

SourceModelBuilder::~SourceModelBuilder() {}

void SourceModelBuilder::addSourceModel(const SourceModel& srcModel) {   
   append_source_model(m_srcLib, srcModel);
}

void SourceModelBuilder::addSource(const Source & src) {
   append_source(m_srcLib, src);
}

void SourceModelBuilder::write(std::string xmlFile) {
   facilities::Util::expandEnvVar(&xmlFile);
   std::ofstream outFile(xmlFile.c_str());
   outFile << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
   prettyPrint(outFile, m_srcLib);
}

rapidxml::xml_node<>* SourceModelBuilder::likelihoodSource(const Source & src) {
   rapidxml::xml_node<>* srcElt = createElement(m_doc, "source");
   addAttribute(m_doc, srcElt, "name", src.getName());
   appendChild(srcElt, spectralPart(src));
   
   switch (src.srcType()) {
   case Source::Composite:
      addComposite(srcElt, src);
      break;
   case Source::Point:
      addSpatialPart(srcElt, src);
      break;
   case Source::Diffuse:
      addSpatialPart(srcElt, src);
      break;
   default:
      throw std::runtime_error("SourceModelBuilder::likelihoodSource: unknown source type.");     
   }
   
   return srcElt;
}

rapidxml::xml_node<>* SourceModelBuilder::spectralPart(const Source & src) {
   const Source::FuncMap & srcFuncs = src.getSrcFuncs();
   
   Source::FuncMap::const_iterator find_spectrum = srcFuncs.find("Spectrum");
   if (find_spectrum == srcFuncs.end()) {
      throw std::runtime_error("SourceModelBuilder::spectralPart: Spectrum function not found");
   }
   
   const optimizers::Function* spectrum = find_spectrum->second;
   
   rapidxml::xml_node<>* specElt = createElement(m_doc, "spectrum");
   addAttribute(m_doc, specElt, "type", spectrum->genericName());
   
   // Handle DMFitFunction special case
   const DMFitFunction* dmfit = dynamic_cast<const DMFitFunction*>(spectrum);
   if (dmfit != nullptr) {
      addAttribute(m_doc, specElt, "file", dmfit->filename());
   }
   
   // Handle FileFunction special case
   const FileFunction* fileFunc = dynamic_cast<const FileFunction*>(spectrum);
   if (fileFunc != nullptr) {
      addAttribute(m_doc, specElt, "file", fileFunc->filename());
   }
   
   // Check if spectrum is wrapped in ScaleFactor
   const ScaleFactor* scaleFactor = dynamic_cast<const ScaleFactor*>(spectrum);
   if (scaleFactor != nullptr) {
      // Get the wrapped spectrum - spectrum() returns a pointer
      const optimizers::Function* wrappedSpec = scaleFactor->spectrum();
      if (wrappedSpec != nullptr) {
         const FileFunction* wrappedFileFunc = dynamic_cast<const FileFunction*>(wrappedSpec);
         if (wrappedFileFunc != nullptr) {
            addAttribute(m_doc, specElt, "file", wrappedFileFunc->filename());
         }
      }
   }
   
   // If the source spectrum has a scaling function, add the filename
   // as the scaling_file attribute.
   const optimizers::Function* scalingFunc = spectrum->scalingFunction();
   if (scalingFunc != nullptr) {
      const FileFunction* scalingFileFunc = dynamic_cast<const FileFunction*>(scalingFunc);
      if (scalingFileFunc != nullptr) {
         addAttribute(m_doc, specElt, "scaling_file", scalingFileFunc->filename());
      }
   }

   optimizers::Function* nc_spectrum = const_cast<optimizers::Function*>(spectrum);
   appendParameterNodes(m_doc, specElt, *nc_spectrum);
   return specElt;
}

void SourceModelBuilder::addSpatialPart(rapidxml::xml_node<>* srcElt, const Source & src) {
   const Source::FuncMap & srcFuncs = src.getSrcFuncs();

   rapidxml::xml_node<>* spatialElt = createElement(m_doc, "spatialModel");
   
   Source::FuncMap::const_iterator find_pos = srcFuncs.find("Position");
   Source::FuncMap::const_iterator find_spatial = srcFuncs.find("SpatialDist");

   if (find_pos != srcFuncs.end()) {
      addAttribute(m_doc, srcElt, "type", "PointSource");
      addAttribute(m_doc, spatialElt, "type", "SkyDirFunction");
      optimizers::Function* nc_pos = const_cast<optimizers::Function*>(find_pos->second);
      appendParameterNodes(m_doc, spatialElt, *nc_pos);
      appendChild(srcElt, spatialElt);
   } else if (find_spatial != srcFuncs.end()) {
      addAttribute(m_doc, srcElt, "type", "DiffuseSource");
      const optimizers::Function* spatial = find_spatial->second;
      std::string type = spatial->genericName();
      addAttribute(m_doc, spatialElt, "type", type);
      
      if (type == "SpatialMap") {
         std::string file = dynamic_cast<const SpatialMap*>(spatial)->fitsFile();
         addAttribute(m_doc, spatialElt, "file", file);
      } else if (type == "MapCubeFunction") {
         std::string file = dynamic_cast<const MapCubeFunction2*>(spatial)->fitsFile();
         addAttribute(m_doc, spatialElt, "file", file);
      } else if (type == "RadialProfile") {
         std::string file = dynamic_cast<const RadialProfile*>(spatial)->templateFile();
         addAttribute(m_doc, spatialElt, "file", file);
      }
      
      const DiffuseSource* diffuseSource = dynamic_cast<const DiffuseSource*>(&src);
      if (diffuseSource != nullptr && diffuseSource->mapBasedIntegral()) {
         addAttribute(m_doc, spatialElt, "map_based_integral", "true");
      }
      
      optimizers::Function* nc_spatial = const_cast<optimizers::Function*>(spatial);
      appendParameterNodes(m_doc, spatialElt, *nc_spatial);
      appendChild(srcElt, spatialElt);
   }
}

void SourceModelBuilder::addComposite(rapidxml::xml_node<>* srcElt, const Source & src) {
   addAttribute(m_doc, srcElt, "type", "CompositeSource");
   
   const CompositeSource* comp = dynamic_cast<const CompositeSource*>(&src);
   if (comp == nullptr) {
      throw std::runtime_error("SourceModelBuilder::addComposite: source not a composite");
   }
   
   const std::string& compXmlFile = comp->xmlFile();
   rapidxml::xml_node<>* compElt = createElement(m_doc, "source_library");
   
   if (compXmlFile.empty()) {
      append_source_model(compElt, comp->sourceModel());    
   } else {
      addAttribute(m_doc, compElt, "xmlFile", compXmlFile);    
   }
   appendChild(srcElt, compElt);
}

void SourceModelBuilder::append_source(rapidxml::xml_node<>* parent, const Source & src) {
   appendChild(parent, likelihoodSource(src));
} 

void SourceModelBuilder::append_source_model(rapidxml::xml_node<>* parent, const SourceModel& srcModel) {
   for (const auto& [name, source] : srcModel.sources()) {
      append_source(parent, *source);
   }
}

} // namespace Likelihood
