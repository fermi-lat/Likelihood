// -*- mode: c++ -*-
%module Likelihood
%{
#include "optimizers/Lbfgs.h"
#include "optimizers/Drmngb.h"
#include "optimizers/Minuit.h"
#include "optimizers/Optimizer.h"
#include "optimizers/Parameter.h"
#include "optimizers/ParameterNotFound.h"
#include "optimizers/Function.h"
#include "optimizers/FunctionFactory.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/Source.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/EventArg.h"
#include "Likelihood/Event.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/FitsImage.h"
#include "Likelihood/logSrcModel.h"
#include "Likelihood/Npred.h"
#include "Likelihood/OneSourceFunc.h"
#include "Likelihood/OptEM.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SourceFactory.h"
#include "Likelihood/SourceModel.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/SrcArg.h"
#include "Likelihood/TrapQuad.h"
#include "Likelihood/Exception.h"
#include "Likelihood/LogLike.h"
#include <vector>
#include <string>
#include <exception>
using optimizers::Parameter;
using optimizers::ParameterNotFound;
using optimizers::Function;
using optimizers::Exception;
%}
%include stl.i
//%include /home/jchiang/ST/optimizers/v1r3/optimizers/Statistic.h
%include /home/jchiang/ST/optimizers/v1r3/optimizers/Function.h
%include /home/jchiang/ST/optimizers/v1r3/optimizers/FunctionFactory.h
%include /home/jchiang/ST/optimizers/v1r3/optimizers/Optimizer.h
%include ../Likelihood/Exception.h
%include ../Likelihood/ResponseFunctions.h
%include ../Likelihood/Event.h
%include ../Likelihood/Source.h
%include ../Likelihood/SourceModel.h
%include ../Likelihood/DiffuseSource.h
%include ../Likelihood/EventArg.h
%include ../Likelihood/ExposureMap.h
%include ../Likelihood/FitsImage.h
%include ../Likelihood/LogLike.h
%include ../Likelihood/logSrcModel.h
%include ../Likelihood/Npred.h
%include ../Likelihood/OneSourceFunc.h
%include ../Likelihood/OptEM.h
%include ../Likelihood/PointSource.h
%include ../Likelihood/RoiCuts.h
%include ../Likelihood/ScData.h
%include ../Likelihood/SkyDirArg.h
%include ../Likelihood/SkyDirFunction.h
%include ../Likelihood/SourceFactory.h
%include ../Likelihood/SpatialMap.h
%include ../Likelihood/SrcArg.h
%include ../Likelihood/TrapQuad.h
%template(DoubleVector) std::vector<double>;
%template(DoubleVectorVector) std::vector< std::vector<double> >;
%template(StringVector) std::vector<std::string>;
%extend Likelihood::SourceFactory {
   optimizers::FunctionFactory * funcFactory() {
      optimizers::FunctionFactory * myFuncFactory 
         = new optimizers::FunctionFactory;
      myFuncFactory->addFunc("SkyDirFunction", 
                             new Likelihood::SkyDirFunction(), false);
      myFuncFactory->addFunc("SpatialMap", new Likelihood::SpatialMap(), 
                             false);
      return myFuncFactory;
   }
}
%extend Likelihood::LogLike {
   void print_source_params() {
      std::vector<std::string> srcNames;
      self->getSrcNames(srcNames);
      std::vector<Parameter> parameters;
      for (unsigned int i = 0; i < srcNames.size(); i++) {
         Likelihood::Source *src = self->getSource(srcNames[i]);
         Likelihood::Source::FuncMap srcFuncs = src->getSrcFuncs();
         srcFuncs["Spectrum"]->getParams(parameters);
         std::cout << "\n" << srcNames[i] << ":\n";
         for (unsigned int i = 0; i < parameters.size(); i++)
            std::cout << parameters[i].getName() << ": "
                      << parameters[i].getValue() << std::endl;
         std::cout << "Npred: "
                   << src->Npred() << std::endl;
      }
   }
   void src_param_table() {
      std::vector<std::string> srcNames;
      self->getSrcNames(srcNames);
      std::vector<Parameter> parameters;
      for (unsigned int i = 0; i < srcNames.size(); i++) {
         Likelihood::Source *src = self->getSource(srcNames[i]);
         Likelihood::Source::FuncMap srcFuncs = src->getSrcFuncs();
         srcFuncs["Spectrum"]->getParams(parameters);
         for (unsigned int i = 0; i < parameters.size(); i++)
            std::cout << parameters[i].getValue() << "  ";
         std::cout << src->Npred() << "  ";
         std::cout << srcNames[i] << std::endl;
      }
   }
   int getNumFreeParams() {
      return self->getNumFreeParams();
   }
   void getFreeParamValues(std::vector<double> & params) {
      self->getFreeParamValues(params);
   }
   optimizers::Optimizer * Minuit() {
      return new optimizers::Minuit(*self);
   }
   optimizers::Optimizer * Lbfgs() {
      return new optimizers::Lbfgs(*self);
   }
   optimizers::Optimizer * Drmngb() {
      return new optimizers::Drmngb(*self);
   }
}
%extend Likelihood::Event {
   double ra() {
      return self->getDir().ra();
   }
   double dec() {
      return self->getDir().dec();
   }
   double scRa() {
      return self->getScDir().ra();
   }
   double scDec() {
      return self->getScDir().dec();
   }
}
%extend Likelihood::FitsImage {
   void fetchCelestialArrays(std::vector<double> &longitude,
                             std::vector<double> &latitude) {
      std::valarray<double> lon;
      std::valarray<double> lat;
      self->fetchCelestialArrays(lon, lat);
      longitude.reserve(lon.size());
      latitude.reserve(lat.size());
      for (unsigned int i = 0; i < lon.size(); i++) {
         longitude.push_back(lon[i]);
         latitude.push_back(lat[i]);
      }
   }
   void fetchImageData(std::vector<double> &image) {
      std::valarray<double> imageArray;
      self->fetchImageData(imageArray);
      image.reserve(imageArray.size());
      for (unsigned int i = 0; i < imageArray.size(); i++)
         image.push_back(imageArray[i]);
   }
}
