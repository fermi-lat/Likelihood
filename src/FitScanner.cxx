/**
 * @file FitScanner.cxx
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/src/FitScanner.cxx,v 1.2 2015/12/02 00:53:06 echarles Exp $
 */


#include "Likelihood/FitScanner.h"

#include <cstdio>
#include <stdexcept>
#include <cmath>
#include <vector>

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"
#include "tip/tip_types.h"

#include "astro/SkyDir.h"
#include "astro/SkyProj.h"
//#include "astro/HealpixProj.h"

#include "optimizers/Optimizer.h"

#include "evtbin/LinearBinner.h"
#include "evtbin/LogBinner.h"
#include "evtbin/HealpixBinner.h"
#include "evtbin/OrderedBinner.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/HistND.h"
#include "Likelihood/Source.h"
#include "Likelihood/ScanUtils.h"
#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/SummedLikelihood.h"
#include "Likelihood/FitUtils.h"
#include "Likelihood/SourceMap.h"

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

namespace Likelihood {

  TestSourceModelCache::TestSourceModelCache(const BinnedLikelihood& logLike,
					     const Source& source)
    :m_refModel(logLike.countsMap().data().size(),0),
     m_proj(logLike.countsMap().projection()),
     m_refDir(logLike.countsMap().refDir()) {    
    // latch the reference direction and the size of the model axes
    m_refPixel = m_refDir.project( m_proj );
    m_nx = logLike.countsMap().imageDimension(0);
    m_ny = logLike.countsMap().imageDimension(1);    
    m_ne = logLike.countsMap().imageDimension(2);
    // Extract the reference image
    FitUtils::extractModelFromSource(source,logLike,m_refModel,true);
  }
    
  int TestSourceModelCache::translateMap(const astro::SkyDir& newRef,
					 std::vector<float>& out_model) const {
    
    /*
    if ( m_proj.method() == astro::ProjBase::HEALPIX ) {
      // FIXME, get these values
      double d_theta(0.);
      double d_phi(0.);
      return translateMap_Healpix(d_theta,d_phi,out_model); 
    }
    */
    std::pair<double,double> newPix = newRef.project( m_proj );
    double dx = newPix.first - m_refPixel.first;
    double dy = newPix.second - m_refPixel.second;
    return translateMap_Wcs(dx,dy,out_model);
  }

  void TestSourceModelCache::writeTestSourceToFitsImage(const std::string& fits_file,
							const std::string& ext_name) const {
    std::vector<long> naxes;    
    naxes.push_back(m_nx);
    naxes.push_back(m_ny);
    naxes.push_back(m_ne);
    // Add an image to the file
    tip::IFileSvc::instance().appendImage(fits_file,ext_name,naxes);
    tip::Image* image(tip::IFileSvc::instance().editImage(fits_file,ext_name));
    // Set the image dimension and keywords
    image->setImageDimensions(naxes);
    tip::Header & header(image->getHeader());
    // FIXME, astro::ProjBase::setKeywords() should be a const function
    astro::ProjBase& nc_proj = const_cast<astro::ProjBase&>(m_proj);
    nc_proj.setKeywords(header);
    // These are just placeholders, since this is just a test, that is ok
    header.setKeyword("CRPIX3",1);
    header.setKeyword("CRVAL3",100.);
    header.setKeyword("CDELT3",100.);    
    header.setKeyword("CTYPE3","photon energy");
    // actually set the data
    image->set(m_currentModel);
    delete image;    
  }   
  
  int TestSourceModelCache::translateMap_Wcs(double dx, double dy, std::vector<float>& out_model) const {

    // First fill the output map with 1e-9 everywhere
    // We don't use 0 to avoid failed matrix inversions 
    // when the number of observed counts is very small
    // FIX (we should actually use the map symmetry to deal with this problem)
    FitUtils::setVectorValue(1e-9,out_model.begin(),out_model.end());    
    
    // Convert the deltas to integers
    // I can't find a built-in rounding function, so do this instead (it is really ugly)
    int delta_x = dx >= 0 ? ( std::fmod(dx,1.0) > 0.5 ? std::ceil(dx) : std::floor(dx) ) :
      ( std::fmod(dx,1.0) < -0.5 ? std::floor(dx) : std::ceil(dx) );
    
    int delta_y = dy >= 0 ? ( std::fmod(dy,1.0) > 0.5 ? std::ceil(dy) : std::floor(dy) ) :
      ( std::fmod(dy,1.0) < -0.5 ? std::floor(dy) : std::ceil(dy) );
    
    if ( false ) {
      std::cout << "Translate WCS: " << delta_x << ' ' << delta_y << std::endl;
    }
    
    // Compute the start and stop indices for the input vector
    // and the start indices for the output vector
    int start_x_in = delta_x < 0 ? -delta_x : 0;
    int stop_x_in = delta_x < 0 ? m_nx : m_nx - delta_x;
    int start_x_out = start_x_in + delta_x;
    
    int start_y_in = delta_y < 0 ? -delta_y : 0;
    int stop_y_in = delta_y < 0 ?  m_ny : m_ny - delta_y;
    int start_y_out = start_y_in + delta_y;
    
    // This is the stride from one energy layer to the next
    int npix = m_nx*m_ny;
    
    size_t map_offset(0);
    for ( size_t ie(0); ie < m_ne; ie++, map_offset+=npix ) {
      // These are the indices of the first pixel in the first row in the
      // input and output maps
      size_t row_start_in =  map_offset + (start_y_in * m_nx);
      size_t row_start_out = map_offset + (start_y_out * m_nx);
      // Here we loop over rows
      for ( size_t iy(start_y_in); iy < stop_y_in; iy++, row_start_in += m_nx, row_start_out += m_nx) {
	// This are the indices of the start and stop pixels in the
	// input and output maps
	size_t in_start = row_start_in + start_x_in;
	size_t in_stop = row_start_in + stop_x_in;
	size_t out_start = row_start_out + start_x_out; 
	// copy from input to output
	std::copy( m_refModel.begin() + in_start,
		   m_refModel.begin() + in_stop,
		   out_model.begin() + out_start);
      }
    }
    
    // copy to the local cache
    m_currentModel.resize(out_model.size(),0.);
    std::copy(out_model.begin(),out_model.end(),m_currentModel.begin());
    
    // This is just a check to make sure that we have preserved the 
    // normalization
    if ( true ) {
      float total_in(0.);
      float total_out(0.);
      FitUtils::sumVector(m_refModel.begin(),m_refModel.end(),total_in);
      FitUtils::sumVector(out_model.begin(),out_model.end(),total_out);
      if ( false ) {
	std::cout << "Translate " << dx << ' ' << delta_x << ' '
		  << dy << ' ' << delta_y << ' '
		  << total_in << ' ' << total_out << std::endl;
      }
    }
    return 0;
  }
  

  void FitScanMVPrior::negativeLogLikelihood(const CLHEP::HepVector& params, double& logLike) const {
    CLHEP::HepVector delta = params - m_centralVals;
    // This returns delta.T * ( m_hessian * delta ), which is exactly what we want
    logLike = m_hessian.similarity(delta);
  }

  void FitScanMVPrior::gradient(const CLHEP::HepVector& params,  CLHEP::HepVector& grad) const {
    CLHEP::HepVector delta = params - m_centralVals;
    grad = m_hessian*delta;
  }
 
  
  void FitScanMVPrior::update(const CLHEP::HepVector& centralVals,
			      const CLHEP::HepSymMatrix& covariance,
			      const std::vector<bool>& constrainPars,
			      bool includeTestSource) {
    m_constrainPars = constrainPars;
    m_includeTestSource = includeTestSource;
    if ( includeTestSource ) {
      m_centralVals = CLHEP::HepVector(centralVals.num_row() + 1,0);
      m_centralVals.sub(1,centralVals);
    } else {
      m_centralVals = centralVals;      
    } 
    m_covariance = covariance;
    latchReducedMatrix();
  }

  int FitScanMVPrior::latchReducedMatrix() {
    // FIXME, invert then reduce, or reduce then invert?
    std::vector<int> idx_red;
    
    int idx(0);
    for ( std::vector<bool>::const_iterator itr = m_constrainPars.begin(); 
	  itr != m_constrainPars.end(); itr++, idx++ ) {
      idx_red.push_back(idx);
    }
    CLHEP::HepSymMatrix tempCov(idx_red.size(),0);
    int idx1(0);
    for ( std::vector<int>::const_iterator itr1 = idx_red.begin(); itr1 != idx_red.end(); itr1++, idx1++ ) {
      int idx2(idx1);
      for ( std::vector<int>::const_iterator itr2 = itr1; itr2 != idx_red.end(); itr2++, idx2++ ) {
	tempCov[idx1][idx2] = m_covariance[*itr1][*itr2];      
      }
    }
    int status(0);
    tempCov.invert(status);
    if ( status ) {
      throw std::runtime_error("Could not invert covarience matrix given for FitScanMVPrior");
      return status;
    }
    
    if ( m_includeTestSource ) {
      m_hessian = CLHEP::HepSymMatrix(idx_red.size()+1,0);
      m_hessian.sub(1,tempCov);
    } else {
      m_hessian = tempCov;    
    }

    if ( false ) {
      FitUtils::printSymMatrix("Prior: ",m_hessian);
    }

    return 0;
  }

  FitScanModelWrapper_Binned::FitScanModelWrapper_Binned(BinnedLikelihood& binnedLike)
    :FitScanModelWrapper(),
     m_binnedLike(binnedLike){
    setDims(binnedLike.countsMap().pixels().size(),binnedLike.countsMap().energies().size()-1);
  }

  const std::vector<float>& FitScanModelWrapper_Binned::data() const {
    return m_binnedLike.countsMap().data();
  }
  

  void FitScanModelWrapper_Binned::extractModelFromSource(Source& aSrc,
							  std::vector<float>& model,
							  bool rescaleToNormOne) const {
    FitUtils::extractModelFromSource(aSrc,m_binnedLike,model,rescaleToNormOne);
  }

  void FitScanModelWrapper_Binned::extractModels(const std::string& test_name,
						 std::vector<std::vector<float> >& templates,		       
						 std::vector<float>& fixed,
						 std::vector<float>& test_source_model,
						 std::vector<float>& refPars) const {
    FitUtils::extractModels(m_binnedLike,test_name,templates,fixed,test_source_model,refPars);
  }
  
  double FitScanModelWrapper_Binned::value() const {
    return m_binnedLike.value();
  }

  void FitScanModelWrapper_Binned::addSource(Source* aSrc) {
    m_binnedLike.addSource(aSrc);
  }

  void FitScanModelWrapper_Binned::syncParams() {
    m_binnedLike.syncParams();
  }
    
  void FitScanModelWrapper_Binned::removeSource(const std::string& sourceName) {
     Source* delSrc = m_binnedLike.deleteSource(sourceName);
     if ( ! delSrc->fixedSpectrum() ) {
       m_binnedLike.eraseSourceMap(sourceName);
     }
     delete delSrc;
  } 

  int FitScanModelWrapper_Binned::writeFits_EnergyBins(const std::string& fitsFile, const evtbin::Binner* binner) const {
    m_binnedLike.countsMap().writeEbounds(fitsFile,binner);
    return 0;
  }
  
  int FitScanModelWrapper_Binned::writeFits_GTIs(const std::string& fitsFile) const {
    m_binnedLike.countsMap().writeGti(fitsFile);
  }

  const std::vector<double>& FitScanModelWrapper_Binned::energies() const {
    return m_binnedLike.energies();
  }

  int FitScanModelWrapper_Binned::shiftTestSource(const std::vector<TestSourceModelCache*>& modelCaches,
						  const astro::SkyDir& newDir,
						  std::vector<float>& targetModel) const {
    if ( modelCaches.size() != 1 ) {
      throw std::runtime_error("FitScanModelWrapper_Binned should only have a single modelCache");
      return -1;
    }
    modelCaches[0]->translateMap(newDir,targetModel);
  }

  void FitScanModelWrapper_Binned::set_klims(size_t kmin, size_t kmax) {
    m_binnedLike.set_klims(kmin,kmax);
  }
  

  double FitScanModelWrapper_Summed::findMinAndMatches(const std::vector< const std::vector<double>* >& energyBins,
						       std::vector<size_t>& localIdx,
						       std::vector<int>& matches,
						       float tol) {
    // First find the vector with the lowest current value;
    int best_idx(-1);
    double lowestValue(1.0e32);
    for ( int idx(0); idx < localIdx.size(); idx++ ) {
      if ( localIdx[idx] >= energyBins[idx]->size() ) continue;
      double testVal = (*energyBins[idx])[localIdx[idx]];
      if ( testVal < lowestValue ) {
	lowestValue = testVal;
	best_idx = idx;
      }						       
    }

    // We have reached the end of all the vectors
    if ( best_idx < 0 ) return -1;
    
    // Loop over the vectors and tag & increment all the ones that also have the lowest current value
    for ( int idx2(0); idx2 < localIdx.size(); idx2++ ) {
      matches[idx2] = -1;
      if ( localIdx[idx2] >= energyBins[idx2]->size() ) continue;
      double testVal = (*energyBins[idx2])[localIdx[idx2]];
      if ( fabs( lowestValue - testVal ) < tol ) {
	matches[idx2] = localIdx[idx2];
	localIdx[idx2] += 1;
      }
    }
    return lowestValue;
  }

  void FitScanModelWrapper_Summed::mergeVectors(const std::vector< const std::vector<float>* >& toMerge,
						const std::vector<size_t>& npixelsByComp,
						const std::vector< std::vector<int> >& energyBinLocal,
						std::vector<float>& mergedData) {

    std::vector<float>::iterator writePtr = mergedData.begin();    
    // Loop over energy bins    
    for ( size_t i(0); i < energyBinLocal.size()-1; i++ ) {
      // Loop over components
      for ( size_t j(0); j < npixelsByComp.size(); j++ ) {
	int ebinLocal = energyBinLocal[i][j];
	size_t nPixLocal = npixelsByComp[j];
	if ( ebinLocal < 0 ) {
	  for ( std::vector<float>::iterator zeroPtr = writePtr; zeroPtr != writePtr + nPixLocal; zeroPtr++ ) {
	    *zeroPtr = 0;
	  }
	} else {
	  std::vector<float>::const_iterator readStartPtr = toMerge[j]->begin() + ebinLocal*nPixLocal;
	  std::vector<float>::const_iterator readEndPtr = readStartPtr + nPixLocal;
	  std::copy(readStartPtr,readEndPtr,writePtr);
	}
	writePtr += nPixLocal;
      }
    }
  }
     

  FitScanModelWrapper_Summed::FitScanModelWrapper_Summed(SummedLikelihood& summedLike)
    :m_summedLike(summedLike),
     m_master(0),
     m_localData(),
     m_nPixelsByComp(summedLike.numComponents(),0),
     m_nEBinsByComp(summedLike.numComponents(),0),
     m_sizeByComp(summedLike.numComponents(),0),
     m_energiesMerged(),
     m_energyBinLocal() {

    // First Loop, get the pixel and the energy bins and counts maps for each component
    std::vector< const std::vector<double>* > energyBins; 
    std::vector< const std::vector<float>* > toMerge;
    size_t nPixTot(0);

    for ( size_t i(0); i < summedLike.numComponents(); i++ ) {
      LogLike* comp = summedLike.getComponent(i);
      BinnedLikelihood* binnedLike = dynamic_cast<BinnedLikelihood*>(comp);

      if ( binnedLike == 0 ) {
	throw std::runtime_error("FitScanModelWrapper can only be built from a sum of BinnedLikelihoods");
	continue;
      }
      if ( i == 0 ) {
	m_master = binnedLike;
      }
      toMerge.push_back( & (binnedLike->countsMap().data()) );
      energyBins.push_back(& (binnedLike->energies()) );
      size_t nPixComp = binnedLike->countsMap().pixels().size();
      size_t nEBinsComp = binnedLike->energies().size()-1;
      m_nPixelsByComp[i] = nPixComp;
      m_nEBinsByComp[i] = nEBinsComp;
      m_sizeByComp[i] = nPixComp*nEBinsComp;
      nPixTot += nPixComp;      
    }

    // Second Loop, merge the energy bins

    // These are the local indexes of the 
    std::vector<size_t> localIdx( summedLike.numComponents(), 0);
    std::vector<int> matches( summedLike.numComponents(), -1);
    double nextEnergy(0.);
    while ( nextEnergy >= 0 ) {      
      nextEnergy = findMinAndMatches( energyBins, localIdx, matches );
      if ( nextEnergy > 0 ) {
	m_energyBinLocal.push_back(matches);
	m_energiesMerged.push_back(nextEnergy);
      }
    }
    setDims(nPixTot,m_energiesMerged.size()-1);
    m_localData.resize(size());
    // Now fill the local data
    mergeVectors(toMerge,m_nPixelsByComp,m_energyBinLocal,m_localData);
  }


  void FitScanModelWrapper_Summed::extractModelFromSource(Source& aSrc,
							  std::vector<float>& model,							  
							  bool rescaleToNormOne) const {
    std::vector< const std::vector<float>* > toMerge;
    for ( size_t i(0); i < m_summedLike.numComponents(); i++ ) {
      LogLike* comp = m_summedLike.getComponent(i);
      BinnedLikelihood* binnedLike = dynamic_cast<BinnedLikelihood*>(comp);
      if ( binnedLike == 0 ) {
	throw std::runtime_error("FitScanModelWrapper can only be built from a sum of BinnedLikelihoods");
	continue;
      }
      size_t compSize = binnedLike->countsMap().data().size();
      std::vector<float>* modelVect = new std::vector<float>(compSize);
      FitUtils::extractModelFromSource(aSrc,*binnedLike,*modelVect,rescaleToNormOne);      
      toMerge.push_back( modelVect );
    }
    mergeVectors(toMerge,m_nPixelsByComp,m_energyBinLocal,model);
    for ( std::vector< const std::vector<float>* >::iterator itrDel = toMerge.begin();
	  itrDel != toMerge.end(); itrDel++ ) {
      std::vector<float>* toDel = const_cast<std::vector<float>*>(*itrDel);
      delete toDel;
    }
  }
  
  void FitScanModelWrapper_Summed::extractModels(const std::string& test_name,
						 std::vector<std::vector<float> >& templates,		       
						 std::vector<float>& fixed,
						 std::vector<float>& test_source_model,
						 std::vector<float>& refPars) const {

    static double tol = 1e-3;
    static int maxIter = 30;
    std::vector<FitScanCache*> caches;
    std::vector< std::vector< const std::vector<float>* > > templatesToMerge;
    std::vector< const std::vector<float>* > fixedToMerge;
    std::vector< const std::vector<float>* > testToMerge;

    std::vector< FitScanModelWrapper_Binned* > wrappers;

    for ( size_t i(0); i < m_summedLike.numComponents(); i++ ) {
      LogLike* comp = m_summedLike.getComponent(i);
      BinnedLikelihood* binnedLike = dynamic_cast<BinnedLikelihood*>(comp);
      if ( binnedLike == 0 ) {
	throw std::runtime_error("FitScanModelWrapper can only be built from a sum of BinnedLikelihoods");
	continue;
      }       
      FitScanModelWrapper_Binned* nbw = new FitScanModelWrapper_Binned(*binnedLike);
      wrappers.push_back(nbw);
      FitScanCache* cache = new FitScanCache(*nbw,test_name,tol,maxIter,false);
      fixedToMerge.push_back( &(cache->allFixed()) );
      testToMerge.push_back( &(cache->targetModel()) );
      refPars.resize(cache->refValues().size());
      if ( i == 0 ) {
	std::copy(cache->refValues().begin(),cache->refValues().end(),refPars.begin());
	templatesToMerge.resize( cache->allModels().size() );	 
      }
       
      for ( int j(0); j < cache->allModels().size(); j++ ) {
	const std::vector<float>* ptr = &(cache->allModels()[j]);
	templatesToMerge[j].push_back(ptr);
      }
    }

    mergeVectors(fixedToMerge,m_nPixelsByComp,m_energyBinLocal,fixed);
    mergeVectors(testToMerge,m_nPixelsByComp,m_energyBinLocal,test_source_model);

    templates.resize(templatesToMerge.size());
    for ( size_t itmp(0); itmp < templatesToMerge.size(); itmp++ ) {
      templates[itmp].resize(size());
      mergeVectors(templatesToMerge[itmp],m_nPixelsByComp,m_energyBinLocal,templates[itmp]);
    }

    for ( std::vector<FitScanModelWrapper_Binned*>::iterator itrd = wrappers.begin(); itrd != wrappers.end(); itrd++ ) {
      FitScanModelWrapper_Binned* toDel = *itrd;
      delete toDel;
    }
  
  }

  double FitScanModelWrapper_Summed::value() const {
    return m_summedLike.value();
  }

  void FitScanModelWrapper_Summed::addSource(Source* aSrc) {
    for ( size_t i(0); i < m_summedLike.numComponents(); i++ ) {
      LogLike* comp = m_summedLike.getComponent(i);
      comp->addSource(aSrc);
    }
  }

  void FitScanModelWrapper_Summed::syncParams() {
    m_summedLike.syncParams();
  }
    
  void FitScanModelWrapper_Summed::removeSource(const std::string& sourceName) {
    for ( size_t i(0); i < m_summedLike.numComponents(); i++ ) {
      LogLike* comp = m_summedLike.getComponent(i);
      BinnedLikelihood* binnedLike = dynamic_cast<BinnedLikelihood*>(comp);
      Source* delSrc = binnedLike->deleteSource(sourceName);
      if ( ! delSrc->fixedSpectrum() ) {
	binnedLike->eraseSourceMap(sourceName);
      }
      delete delSrc;
    }
  }

  int FitScanModelWrapper_Summed::writeFits_EnergyBins(const std::string& fitsFile, const evtbin::Binner* binner) const {
    m_master->countsMap().writeEbounds(fitsFile,binner);    
    return 0;
  }
  
  int FitScanModelWrapper_Summed::writeFits_GTIs(const std::string& fitsFile) const {
    m_master->countsMap().writeGti(fitsFile);
    return 0;
  } 

  const size_t FitScanModelWrapper_Summed::numComponents() const {
    return m_summedLike.numComponents();
  }
  
  const BinnedLikelihood* FitScanModelWrapper_Summed::getComponent(size_t idx) const {
    return dynamic_cast<const BinnedLikelihood*>(m_summedLike.getComponent(idx));
  }

  
  int FitScanModelWrapper_Summed::shiftTestSource(const std::vector<TestSourceModelCache*>& modelCaches,
						  const astro::SkyDir& newDir,
						  std::vector<float>& targetModel) const {
    if ( modelCaches.size() != m_summedLike.numComponents() ) {
      throw std::runtime_error("FitScanModelWrapper_Summed number of modelCaches should equal number of Likelihood components");
      return -1;
    }

    std::vector< std::vector<float> > targetModels(m_summedLike.numComponents());
    std::vector< const std::vector<float>* > targetModelsPtrs;
    for ( size_t i(0); i < modelCaches.size(); i++ ) {
      size_t compSize = m_sizeByComp[i];
      targetModels[i].resize(compSize);
      std::vector<float>* newVect = &(targetModels[i]);
      modelCaches[i]->translateMap(newDir,*newVect);
      targetModelsPtrs.push_back(newVect);
    }

    mergeVectors(targetModelsPtrs,m_nPixelsByComp,m_energyBinLocal,targetModel);

    return 0;      
  }
  
  void FitScanModelWrapper_Summed::set_klims(size_t kmin, size_t kmax) {
    for ( size_t i(0); i < m_summedLike.numComponents(); i++ ) {
      LogLike* comp = m_summedLike.getComponent(i);
      BinnedLikelihood* binnedLike = dynamic_cast<BinnedLikelihood*>(comp);
      int kmin_local = m_energyBinLocal[kmin][i];
      int kmax_local = m_energyBinLocal[kmin][i];
      if ( kmin_local < 0 || kmax_local < 0 ) {
	binnedLike->set_klims(0,0);
      } else {
	binnedLike->set_klims(kmin_local,kmax_local);
      }
    }
  }
  
  FitScanCache::FitScanCache(FitScanModelWrapper& modelWrapper,
			     const std::string& testSourceName,
			     double tol, int maxIter,
			     bool useReduced) 
    :m_modelWrapper(modelWrapper),
     m_testSourceName(testSourceName),
     m_tol(tol),
     m_maxIter(maxIter),
     m_nebins(m_modelWrapper.nEBins()),
     m_npix(m_modelWrapper.nPix()),
     m_data(modelWrapper.data()),
     m_allFixed(m_modelWrapper.size()),
     m_targetModel(m_modelWrapper.size()),
     m_useReduced(useReduced),
     m_loglike_ref(m_modelWrapper.value()),
     m_currentFixed(m_modelWrapper.size()),
     m_prior_test(0),
     m_prior_bkg(0),
     m_currentBestModel(m_modelWrapper.size()),
     m_currentTestSourceIndex(-1),
     m_currentLogLike(0.),
     m_currentEDM(0.),
     m_currentEnergyBin(-1),
     m_firstBin(0),
     m_lastBin(0){
    // Extract the reference values for everything
    m_modelWrapper.extractModels(m_testSourceName,
				 m_allModels,
				 m_allFixed,
				 m_targetModel,
				 m_refValues);
    // Set up the baseline fit.   
    std::vector<bool> freeSources(m_allModels.size(),true);
    std::vector<float> parScales(m_allModels.size(),1.);
    setEnergyBin(-1);

    if ( m_useReduced ) {
      reduceModels();
    }

    refactorModel(freeSources,parScales,false);  
  }
  
 
  FitScanCache::~FitScanCache() {
    delete m_prior_test;
    delete m_prior_bkg;
  }

  void FitScanCache::refactorModel(const std::vector<bool>& freeSources, 
				   const std::vector<float>& pars_scales,
				   bool include_test) {
    
    // Get the ptr fo the test source model, if requested
    const std::vector<float>* test_source_ptr = include_test ? 
      ( m_useReduced ? &m_targetRedModel : &m_targetModel ) : 0;
    
    // This is to remapping the fit parameters
    std::vector<float> pars_out;
    
    // Refactor the models, 
    //  1) adds fixed models to the m_currentFixed vector
    //  2) pushes free models onto the m_currentModels vector
    //  3) pushes free paramters to the pars_out vector
    if ( m_useReduced ) {						
      m_currentFixed.resize(m_dataRed.size());
      m_currentBestModel.resize(m_dataRed.size());
      FitUtils::refactorModels(m_allRedModels,m_allRedFixed,pars_scales,freeSources,
			       test_source_ptr,
			       m_currentModels,m_currentFixed,pars_out);
    } else {
      FitUtils::refactorModels(m_allModels,m_allFixed,pars_scales,freeSources,
			       test_source_ptr,
			       m_currentModels,m_currentFixed,pars_out);      
    }
    
    // Set the initial and current parameters and covariences for the current fit
    size_t nfree = m_currentModels.size();
    m_currentTestSourceIndex = include_test ? nfree - 1 : -1;
    FitUtils::Vector_Stl_to_Hep(pars_out,m_initPars);
    FitUtils::Vector_Stl_to_Hep(pars_out,m_currentPars);
    m_currentLogLike = 0.;
    m_currentEDM = 0.;
    
    m_currentCov = CLHEP::HepSymMatrix(nfree);
    m_currentGrad = CLHEP::HepVector(nfree);
    m_currentSourceIndices.clear();
    m_currentRefValues.clear();
    for ( size_t i(0); i < freeSources.size(); i++ ) {
      if ( freeSources[i] ) {
	m_currentSourceIndices.push_back(i);
	m_currentRefValues.push_back(m_refValues[i]);
      }
    }  

    delete m_prior_test;
    delete m_prior_bkg;
    m_prior_test = 0;
    m_prior_bkg = 0;
  }
  
  void FitScanCache::getParScales(std::vector<float>& pars_scales) {
    // resize the vector and fill it with ones
    pars_scales.resize(nBkgModel());
    FitUtils::setVectorValue(1.0,pars_scales.begin(),pars_scales.end());
    // loop on the vector of current indices and copy the values 
    // to the correct location in pars_scales
    size_t idx(0);
    for ( std::vector<int>::const_iterator itr = m_currentSourceIndices.begin();
	  itr != m_currentSourceIndices.end(); itr++, idx++ ) {
      pars_scales[*itr] = m_currentPars[idx];
    }
  }
  
  void FitScanCache::setEnergyBin(int energyBin) {
    // set the start and stop indices for the current energy bin
    m_currentEnergyBin = energyBin;
    if ( energyBin < 0 ) {
      // Do all bins.  Setting m_lastBin = 0 will do this
      m_firstBin = 0;
      m_lastBin = 0;
    } else {
      // Loop from m_npix*energyBin to m_npix*(energyBin+1)
      if ( m_useReduced ) {
	m_firstBin = energyBin == 0 ? 0 : m_energyBinStopIdxs[energyBin-1];
	m_lastBin = m_energyBinStopIdxs[energyBin];
      } else {
	m_firstBin = m_npix*energyBin;
	m_lastBin = m_firstBin+m_npix;
      }
    }
  }
  
  void FitScanCache::setTestSource(Source& aSrc) {
    // First remove the current version of the source
    removeTestSourceFromCurrent();
    // Extract the predicted counts from the SourceModel object
    // setting the last arguement to true sets the parameters scale to 1.0
    m_modelWrapper.extractModelFromSource(aSrc,m_targetModel,true);
    if ( m_useReduced ) {
      FitUtils::sparsifyModel(m_nonZeroBins,m_targetModel,m_targetRedModel);    
    }    
    // Add the new version of the source to the model
    addTestSourceToCurrent(0.0);
  } 
  
  int FitScanCache::shiftTestSource(const std::vector<TestSourceModelCache*>& modelCaches,
				    const astro::SkyDir& newDir) {
    // First remove the current version of the source
    removeTestSourceFromCurrent();
    
    int status = m_modelWrapper.shiftTestSource(modelCaches,newDir,m_targetModel);
    if ( status != 0 ) {
      // FIXME, do we throw an exception here?
      return status;
    }
    if ( m_useReduced ) {
      FitUtils::sparsifyModel(m_nonZeroBins,m_targetModel,m_targetRedModel);    
    }        
    // Add the new version of the source to the model
    addTestSourceToCurrent(0.0);
    return 0;
  }
  

  void FitScanCache::addTestSourceToCurrent(double initNorm) {
    // If the test source is already in the model, do nothing
    if ( m_currentTestSourceIndex >= 0 ) return;
    
    // Add the target source model to the vector of models
    if ( m_useReduced ) {      
      m_currentModels.push_back(&m_targetRedModel);
    } else {
      m_currentModels.push_back(&m_targetModel);
    }
    
    // Push the initNorm value to the back of the vector of parameters
    std::vector<float> pars;
    FitUtils::Vector_Hep_to_Stl(m_initPars,pars);
    pars.push_back(initNorm);
    FitUtils::Vector_Stl_to_Hep(pars,m_initPars);
    FitUtils::Vector_Stl_to_Hep(pars,m_currentPars);
    
    FitUtils::Vector_Hep_to_Stl(m_currentGrad,pars);
    pars.push_back(0.);
    FitUtils::Vector_Stl_to_Hep(pars,m_currentGrad);   
    
    // Copy the covariance matrix upper block elements
    m_currentTestSourceIndex = pars.size() -1;
    CLHEP::HepSymMatrix newCov(pars.size());
    newCov.sub(1,m_currentCov);
    m_currentCov = newCov;
  }
  
  void FitScanCache::removeTestSourceFromCurrent() {
    // If the test source not in the model, do nothing
    if ( m_currentTestSourceIndex < 0 ) return;
    m_currentTestSourceIndex = -1;
    // Remove the last element from the model list
    m_currentModels.pop_back();
    int nPars = m_currentModels.size();
    if ( nPars > 0 ) {
      // copy the first nPars elements of the parameter vectors
      // and the upper covariance matrix upper block elements
      m_initPars = m_initPars.sub(1,nPars);
      m_currentPars = m_currentPars.sub(1,nPars);
      m_currentCov = m_currentCov.sub(1,nPars);
      m_currentGrad = m_currentGrad.sub(1,nPars);
    } else {
      // No free parameters, make null versions of the 
      // parameter vectors and covariance matirx
      m_initPars = CLHEP::HepVector();
      m_currentPars = CLHEP::HepVector();
      m_currentCov = CLHEP::HepSymMatrix();
      m_currentGrad = CLHEP::HepVector();
    }
  }
  

  /* Set the prior */
  void FitScanCache::buildPriorsFromExternal(const CLHEP::HepVector& centralVals,
					     const CLHEP::HepSymMatrix& covariance,
					     const std::vector<bool>& constrainPars) {
    if ( m_prior_test ) {
      m_prior_test->update(centralVals,covariance,constrainPars,true);
    } else {
      m_prior_test = new FitScanMVPrior(centralVals,covariance,constrainPars,true);
    }
    if ( m_prior_bkg ) {
      m_prior_bkg->update(centralVals,covariance,constrainPars,false);
    } else {
      m_prior_bkg = new FitScanMVPrior(centralVals,covariance,constrainPars,false);
    }

  }
  
  /* Set the prior from the current fit*/
  void FitScanCache::buildPriorsFromCurrent(const std::vector<bool>& constrainPars,
					    double covScaleFactor) {
    if ( m_currentTestSourceIndex < 0 ) {
      CLHEP::HepSymMatrix scaledCov = covScaleFactor*m_currentCov;
      buildPriorsFromExternal(m_currentPars,scaledCov,constrainPars);    
      return;
    }
    CLHEP::HepVector redPars = m_currentPars.sub(1,nBkgModel());
    CLHEP::HepSymMatrix scaledCov = covScaleFactor*m_currentCov.sub(1,nBkgModel());
    buildPriorsFromExternal(redPars,scaledCov,constrainPars); 
  }
  

  int FitScanCache::fitCurrent(bool usePrior, int verbose) {
    // This just passes the cached data along to the FitUtils function
    // and latches the output into the internal cache
    FitScanMVPrior* prior(0);
    if ( usePrior ) {
      prior = m_currentTestSourceIndex >= 0 ? m_prior_test : m_prior_bkg;
    }
    int status = FitUtils::fitNorms_newton(m_useReduced ? m_dataRed : m_data,
					   m_initPars,
					   m_currentModels,
					   m_currentFixed,
					   prior,
					   m_tol,m_maxIter,
					   m_currentPars,
					   m_currentCov,
					   m_currentGrad,
					   m_currentBestModel,
					   m_currentEDM,
					   m_currentLogLike,
					   m_firstBin,
					   m_lastBin,
					   verbose);
    return status;
  }
  
  
  /* Calculate the log-likelihood for the currently cached values */
  int FitScanCache::calculateLoglikeCurrent(double& logLike) {
    
    FitUtils::sumModel(m_currentPars,m_currentModels,m_currentFixed,m_currentBestModel,
		       m_firstBin,m_lastBin);
        
    std::vector<float>::const_iterator data_start = m_useReduced ? m_dataRed.begin() + m_firstBin : m_data.begin() + m_firstBin;
    std::vector<float>::const_iterator data_end =  m_useReduced ? 
      ( m_lastBin == 0 ? m_dataRed.end() : m_dataRed.begin() + m_lastBin ) :
      ( m_lastBin == 0 ? m_data.end() : m_data.begin() + m_lastBin );
    std::vector<float>::const_iterator model_start = m_currentBestModel.begin() + m_firstBin;
    std::vector<float>::const_iterator model_end = m_lastBin == 0 ? m_currentBestModel.end() : m_currentBestModel.begin() + m_lastBin;
    
    logLike = FitUtils::negativeLogLikePoisson(data_start,data_end,model_start,model_end);
    return 0;
  }
  
  int FitScanCache::scanNormalization(int nnorm, double normSigma,
				      double posErr, double negErr,
				      std::vector<double>& norms,				      
				      std::vector<double>& logLikes) {
    
    // If there is no current test source we can't scan its normalization
    if ( m_currentTestSourceIndex < 0 ) {
      return -1;
    } 
    
    // Get the scan range and the scan step
    double norm_mle = m_currentPars[m_currentTestSourceIndex];
    double scan_val = std::max(1.0e-10,norm_mle - (normSigma*negErr));
    double scan_max = std::max(2.5e-01,norm_mle + (normSigma*posErr));
    double lin_step = (scan_max - scan_val) / float(nnorm-1);

    // Resize the output vectors
    norms.resize(nnorm);
    logLikes.resize(nnorm);
    
    // Remove test source from the models we are fitting
    std::vector<const std::vector<float>* > models_temp(m_currentModels);
    models_temp.pop_back();
    // Space for the sum of the fixed models, this will include the test source
    std::vector<float> fixed_temp(m_currentFixed.size());
    
    // Check to see if there are any other free sources
    bool do_profile = models_temp.size() > 0.;
    
    // Make space for the output
    CLHEP::HepVector init_pars;
    CLHEP::HepVector pars_temp;
    CLHEP::HepVector grad_temp;
    CLHEP::HepSymMatrix covs_temp;
    std::vector<float> model_temp(m_currentFixed.size());
    double edm_temp(0.);    
    
    // if we are doing profile fitting, copy over the initial pars
    if ( do_profile ) {
      init_pars = m_currentPars.sub(1,models_temp.size());
    }
    
    // Loop on the normalizations
    for ( int is(0); is < nnorm; is++ ) {
      norms[is] = scan_val;
      // Add the test model to the fixed models (with the correct normalization factor)
      FitUtils::vectorAdd(m_currentFixed.begin()+m_firstBin,m_currentFixed.begin()+m_lastBin,
			  m_targetModel.begin()+m_firstBin,m_targetModel.begin()+m_lastBin,
			  fixed_temp.begin()+m_firstBin,fixed_temp.begin()+m_lastBin,
			  1.,scan_val);
      if ( do_profile ) {
	// Fit the other sources
	int status = FitUtils::fitNorms_newton((m_useReduced ? m_dataRed : m_data),
					       init_pars,models_temp,fixed_temp,
					       m_prior_bkg,
					       m_tol,m_maxIter,
					       pars_temp,covs_temp,grad_temp,
					       model_temp,
					       edm_temp,logLikes[is],
					       m_firstBin,m_lastBin);
	if ( status ) {	    
	  std::cout << "Failed profile fit on energy bin " << is << ".  Status: " << status << std::endl;
	  return status;
	} 
      } else {
	// No need to do the fit, just build the total model and get the log-likelihood
	FitUtils::sumModel(init_pars,models_temp,fixed_temp,model_temp,m_firstBin,m_lastBin);
	std::vector<float>::const_iterator itrDataBeg = m_useReduced ? m_dataRed.begin()+m_firstBin : m_data.begin()+m_firstBin;
	std::vector<float>::const_iterator itrDataEnd = m_useReduced ? m_dataRed.begin()+m_lastBin : m_data.begin()+m_lastBin;
	logLikes[is] = FitUtils::negativeLogLikePoisson(itrDataBeg,itrDataEnd,
							model_temp.begin()+m_firstBin,model_temp.begin()+m_lastBin);   
      }
      // step the scan value
      scan_val += lin_step;
    }
    return 0;
  }  
  
  int FitScanCache::signalUncertainty_quad(double deltaLogLike,
					   double& posErr,
					   double& negErr) {
    // If there is no current test source we can't return its errors
    if ( m_currentTestSourceIndex < 0 ) {
      posErr = 0.;
      negErr = 0.;
      return -1;
    } 
    
    // This is just solving the quadratic equation
    double inv_curve = m_currentCov[m_currentTestSourceIndex][m_currentTestSourceIndex];
    double a_val = 1./inv_curve;
    double b_val = m_currentGrad[m_currentTestSourceIndex];
    double det = b_val*b_val;
    det += 4.*a_val*deltaLogLike;
    
    if ( det < 0 ) {
      posErr = 0.;
      negErr = 0.;      
      return -2;
    }
    
    double sqrt_det = std::sqrt(det);    
    posErr = (-b_val + sqrt_det)/(2.*a_val);
    // Don't let the negative error be greater than the current value....
    negErr = std::min((+b_val + sqrt_det)/(2.*a_val),m_currentPars[m_currentTestSourceIndex]);
    return 0;
  }
  

  void FitScanCache::reduceModels() {					
    if ( ! m_useReduced ) {
      return;
    }
    FitUtils::extractNonZeroBins(m_data,m_npix,m_nonZeroBins,m_dataRed,m_energyBinStopIdxs);
    int nred = m_dataRed.size();

    // for debugging
    if ( false ) {
      std::cout << "Reduced data has " << nred << " bins." << std::endl;    
      for ( std::vector<int>::const_iterator itrIdx = m_energyBinStopIdxs.begin();
	    itrIdx != m_energyBinStopIdxs.end(); itrIdx++ ) {
	std::cout << *itrIdx << ' ';
      }
      std::cout << std::endl;
    }

    for ( std::vector<std::vector<float> >::const_iterator itr = m_allModels.begin();
	  itr != m_allModels.end(); itr++ ) {
      std::vector<float> redModel;
      FitUtils::sparsifyModel(m_nonZeroBins,*itr,redModel);
      if ( redModel.size() != nred ) {
	throw std::runtime_error("Mistmatch between reduced data vector size and reduced model vector size.");
      }
      m_allRedModels.push_back(redModel);
    }
    FitUtils::sparsifyModel(m_nonZeroBins,m_allFixed,m_allRedFixed);
    FitUtils::sparsifyModel(m_nonZeroBins,m_targetModel,m_targetRedModel);
    if ( m_allRedFixed.size() != nred ) {
      throw std::runtime_error("Mistmatch between reduced data vector size and fixed component vector size.");
    }
    if ( m_targetRedModel.size() != nred ) {
      throw std::runtime_error("Mistmatch between reduced data vector size and target component vector size.");
    }
  }



  evtbin::Binner* FitScanner::buildEnergyBinner(const std::vector<double>& energies) {
    std::vector<evtbin::Binner::Interval> energy_intervals;
    for (unsigned int i = 0; i < energies.size()-1; i++) {
      energy_intervals.push_back(evtbin::Binner::Interval(energies[i], 
							  energies[i+1]));
    }
    evtbin::Binner* binner = new evtbin::OrderedBinner(energy_intervals,
						       "photon energy");
    return binner;
  }
  
  astro::SkyProj* FitScanner::buildSkyProj(const std::string &projName,
					   const astro::SkyDir& dir,
					   double pixSize,
					   int nPix, bool galactic) {
    double crpix[] = {nPix/2. + 0.5, nPix/2. + 0.5};
    double xref = galactic ? dir.l() : dir.ra();
    double yref = galactic ? dir.b() : dir.dec();
    double crval[] = {xref, yref};
    double cdelt[] = {-pixSize, pixSize};
    astro::SkyProj* skyProj = new astro::SkyProj(projName, crpix, crval, cdelt, 0, galactic);
    return skyProj;
  }

  // C'tor from WCS grid of directions
  FitScanner::FitScanner(BinnedLikelihood& binnedLike,
			 optimizers::Optimizer& optimizer,
			 const astro::SkyProj& proj,
			 int nx, int ny)
    :m_modelWrapper(new FitScanModelWrapper_Binned(binnedLike)),
     m_opt(&optimizer),
     m_proj(&proj),
     m_testSourceDir(),
     m_dir1_binner(new evtbin::LinearBinner(0,nx,1.,"pix_x")),
     m_dir2_binner(new evtbin::LinearBinner(0,ny,1.,"pix_y")),
     m_energy_binner(0),
     m_norm_binner(0),
     m_testSource(0),
     m_testSourceName("TestSource"),
     m_cache(0),
     m_testSourceCaches(1,0),
     m_verbose_null(0),
     m_verbose_bb(0),
     m_verbose_scan(0),
     m_writeTestImages(false),
     m_useReduced(true){
        
    // Build the energy binned from the energies in the BinnedLikelihood
    m_energy_binner = buildEnergyBinner(m_modelWrapper->energies());  
  }

  // C'tor from WCS grid of directions
  FitScanner::FitScanner(SummedLikelihood& summedLike,
			 optimizers::Optimizer& optimizer,
			 const astro::SkyProj& proj,
			 int nx, int ny)
    :m_modelWrapper(new FitScanModelWrapper_Summed(summedLike)),
     m_opt(&optimizer),
     m_proj(&proj),
     m_testSourceDir(),
     m_dir1_binner(new evtbin::LinearBinner(0,nx,1.,"pix_x")),
     m_dir2_binner(new evtbin::LinearBinner(0,ny,1.,"pix_y")),
     m_energy_binner(0),
     m_norm_binner(0),
     m_testSource(0),
     m_testSourceName("TestSource"),
     m_cache(0),
     m_testSourceCaches(summedLike.numComponents(),0),
     m_verbose_null(0),
     m_verbose_bb(0),
     m_verbose_scan(0),
     m_writeTestImages(false),
     m_useReduced(true){
        
    // Build the energy binned from the energies in the BinnedLikelihood
    m_energy_binner = buildEnergyBinner(m_modelWrapper->energies());  
  }

  

  
  // C'tor from HEALPix region set of directions
  /* 
  FitScanner::FitScanner(LogLike& logLike,
			 optimizers::Optimizer& optimizer,
			 const astro::HealpixProj& proj,
			 const std::string& region)
    :m_logLike(&logLike),
     m_opt(&optimizer),
     m_proj(&proj),
     m_testSourceDir(),
     m_dir1_binner(new evtbin::HealpixBinner(proj.healpix().Nside(), 
					     proj.healpix().Scheme(),
					     SET_NSIDE,
					     proj.isGalactic(),
					     region,
					     "HEALPIX")),
     m_dir2_binner(0),
     m_energy_binner(0),
     m_norm_binner(0),
     m_testSource(0),
     m_testSourceName("TestSource"),
     m_cache(0),
     m_testSourceCache(0){
    
    // Make sure we are doing a binned fit
    BinnedLikelihood* binnedLike = dynamic_cast<BinnedLikelihood*>(m_logLike);
    if ( binnedLike == 0 ) {
      throw std::runtime_error("FitScanner requires using Binned likelihood fitting.");
    }
    
    // Build the energy binned from the energies in the BinnedLikelihood
    m_energy_binner = buildEnergyBinner(binnedLike->energies());
  }
  */

  // D'tor, does cleanup
  FitScanner::~FitScanner() throw() {
    //if ( m_testSource != 0 ) {
    //  removeTestSourceFromModel();
    //}
    delete m_testSource;
    delete m_dir1_binner;
    delete m_dir2_binner;
    delete m_energy_binner;
    delete m_norm_binner;
    for ( std::vector< std::pair< std::string,std::pair<HistND*,std::string> > >::iterator itr = m_scanData.begin();
	  itr != m_scanData.end(); itr++ ) {
      delete itr->second.first;
    }
    m_scanData.clear();
    delete m_cache;
    deleteTestModelCaches();
    delete m_modelWrapper;
  }
  
  
  /* Build a TS map.
     This scans over the directions and calculates the Test Statistics w.r.t. the null 
     hypothesis for each */  
  int FitScanner::run_tsmap(double tol, int tolType, int maxIter, bool remakeTestSource) {
    return run_tscube(false,0,0.0,tol,tolType,maxIter,remakeTestSource);
  }
  
  /* Build a TS cube.
     For each point in a TS Map this also calculate the spectrum as a function of energy
     and can also scan over the normalization
  */  
  int FitScanner::run_tscube(bool doSED, int nNorm, double normSigma, 
			     double covScale,
			     double tol, int maxIter, int tolType, 
			     bool remakeTestSource,
			     int ST_scan_level) {

    if ( true ) {
      std::cout << "nNorm =     " << nNorm << std::endl;
      std::cout << "normSigma = " << normSigma << std::endl;
      std::cout << "covScale  = " << covScale << std::endl;
      std::cout << "tol       = " << tol << std::endl;
      std::cout << "maxIter   = " << maxIter << std::endl;
      std::cout << "tolType   = " << tolType << std::endl;
      std::cout << "remakeSrc = " << (remakeTestSource ? 'T' : 'F') << std::endl;
      std::cout << "st_scan   = " << ST_scan_level << std::endl;
    }

    // How much fitting do we do with ScienceTools fitter
    const bool baseline_st( ST_scan_level > 0 );
    const bool broadband_st( ST_scan_level > 1 );
    const bool sed_st( ST_scan_level > 2 );
    const bool normscan_st( ST_scan_level > 3 );
        
    int status(0);
    double loglike_null_st = m_modelWrapper->value();
    double loglike_null = m_modelWrapper->value();

    // If requested, do the baseline fit with the ScienceTools fitter
    // Note that this will re-fit the source spectra, if so directed by 
    // the input xml file
    if ( baseline_st ) {
      std::cout << "Doing baseline fit with standard fitter." << std::endl;
      status = baselineFit(tol,tolType);
      if ( status != 0 ) {
	std::cout << "Baseline fit failed with status " << status << std::endl;
	// return status;	
      }
      loglike_null_st = m_modelWrapper->value();
      std::cout << "Did baseline fit.  Likelihood before: " << loglike_null << ", after: " << loglike_null_st << std::endl;
    }
    
    // Build the cache object, deleting the old version if needed
    delete m_cache;
    m_cache = new FitScanCache(*m_modelWrapper,m_testSourceName,tol,maxIter,useReduced());        
    m_cache->calculateLoglikeCurrent(loglike_null);
    std::cout << "Got null likelihood " << loglike_null << std::endl;    

    // Do the baseline fit and latch the results.   
    status = m_cache->fitCurrent(false,verbose_null());
    if ( status != 0 ) return status;	
    
    loglike_null = m_cache->currentLogLike();
    std::cout << "Redid baseline fit with linear fitter.  Likelihood: " << loglike_null << std::endl;
    
    // Get the number of bins in the grid of the region
    // This allows for:
    //   1) standard WCS-grid: m_dir1_binner and m_dir2_binner both exist
    //   2) HEALpix pixelization: m_dir1_binner exists, dir2_binner does not
    long nxbins = m_dir1_binner->getNumBins();
    long nybins = m_dir2_binner ? m_dir2_binner->getNumBins() : 1;    
    
    // Make the binner for the normalization scan, if requested
    if ( nNorm > 0 ) {
      m_norm_binner = new evtbin::LinearBinner(0,nNorm,1.,"norm");
    }
    
    // Check to see if we should also do the normalization scan using the 
    // ScienceTools minimizer
    int nNorm_st = normscan_st ? nNorm : 0;

    // These are the names of the various histograms we store
    
    // 2D (or HEALPix) histograms (MAPS)
    static const std::string tsmap_name("TSMAP");
    static const std::string norm_map_name("N_MAP");
    static const std::string posErr_map_name("ERRP_MAP");
    static const std::string negErr_map_name("ERRN_MAP");

    // 3D (or HEALPix,energy) histograms (CUBES)    
    static const std::string tscube_name("TSCUBE");
    static const std::string norm_cube_name("N_CUBE");
    static const std::string posErr_cube_name("ERRPCUBE");
    static const std::string negErr_cube_name("ERRNCUBE");
    static const std::string nll_cube_name("NLL_CUBE");

    // 4D (or HEALPIX,energy,normalization) histograms (SCANS)
    static const std::string norm_name("NORMSCAN");
    static const std::string delta_ll_name("NLL_SCAN");
    

    // Build the varius histograms
    // Note that buildHist can return a null pointer
    
    // Histograms for Newton's Method fitting
    HistND* ts_map = buildHist(tsmap_name,true,false,false);
    HistND* norm_map = buildHist(norm_map_name,true,false,false);
    HistND* posErr_map = buildHist(posErr_map_name,true,false,false);
    HistND* negErr_map = buildHist(negErr_map_name,true,false,false);

    HistND* ts_cube = doSED ? buildHist(tscube_name,true,true,false) : 0;
    HistND* norm_cube = doSED ? buildHist(norm_cube_name,true,true,false) : 0;
    HistND* posErr_cube = doSED ? buildHist(posErr_cube_name,true,true,false) : 0;
    HistND* negErr_cube = doSED ? buildHist(negErr_cube_name,true,true,false) : 0;
    HistND* nll_cube = doSED ? buildHist(nll_cube_name,true,true,false) : 0;
    HistND* norm_vals = nNorm > 0 ? buildHist(norm_name,true,true,true) : 0;
    HistND* delta_ll_vals = nNorm > 0 ? buildHist(delta_ll_name,true,true,true) : 0;
    
    // Histograms for ScienceTools fitting
    HistND* ts_map_st = broadband_st ? buildHist(tsmap_name+"_ST",true,false,false) : 0;
    HistND* ts_cube_st = ( sed_st && doSED ) ? buildHist(tscube_name+"_ST",true,true,false) : 0;
    HistND* norm_cube_st = ( sed_st && doSED ) ? buildHist(norm_cube_name+"_ST",true,true,false) : 0;
    HistND* nll_cube_st = ( sed_st && doSED ) ? buildHist(nll_cube_name+"_ST",true,true,false) : 0;
    HistND* norm_vals_st =  nNorm_st > 0 ? buildHist(norm_name+"_ST",true,true,true) : 0;
    HistND* delta_ll_vals_st = nNorm_st > 0 ? buildHist(delta_ll_name+"_ST",true,true,true) : 0;
    
    
    // We want all the free sources to be refit at each grid point (for the broadband fits)
    std::vector<bool> freeSources(m_cache->nBkgModel(),true);
    // We use the current (null) fit to get the initial point for each grid point
    std::vector<float> parScales(m_cache->nBkgModel(),1.0);
    m_cache->getParScales(parScales);
    
    // Count the number of failed fits
    int nfailed_bb(0);
    int nfailed_bb_newton(0);
    int nfailed_scan(0);
    int nfailed_scan_newton(0);
    
    // We store the output by pixel, so these are useful
    long ipix(0);
    int npix = nPixels();
    
    // Note the loop order, outer loop is on Y, this matches HistND structure
    for ( long iy(0); iy < nybins; iy++ ) {      
      for ( long ix(0); ix < nxbins; ix++, ipix++ ) {
	
	// Set the test source direction from the grid
	status = setTestSourceDir(ix,iy);
	if ( status != 0 ) {
	  throw std::runtime_error("Failed to set test source direction.");
	  return -1;
	}

	// Add the test source to the SourceModel if needed.
	// This also recomputes the test source image	
	if ( broadband_st || remakeTestSource ) {
	  status = addTestSourceToModel();
	  if ( status != 0 ) {
	    throw std::runtime_error("Failed to add test source to model.");
	    return -1;
	  }	  
	}
	
	// Do the broadband fit with the ScienceTools, if requested
	if ( broadband_st ) {
	  status = fitTestSourceBroadband(tol,tolType);
	  if ( status != 0 ) {
	    // count the number of failed fits
	    nfailed_bb++;
	    continue;
	  }
	  
	  // get the TS value and copy to the output histogram
	  double loglike_bb_st = m_modelWrapper->value();
	  double tsval_st = 2*(loglike_bb_st - loglike_null_st);	
	  if ( ts_map_st ) {
	    ts_map_st->setBinDirect(ipix,tsval_st);
	  }
	}

	// This resets the cache to the null fit 
	m_cache->refactorModel(freeSources,parScales,false);

	// Add the current test source to the null fit
	if ( remakeTestSource ) {
	  // This version uses the SourceMap recomputed by the SourceModel
	  m_cache->setTestSource(*m_testSource);
	} else {
	  // This version just shifts the image by some number of pixels
	  m_cache->shiftTestSource(m_testSourceCaches,m_testSourceDir);
	  if ( writeTestImages() ) {
	    char buffer[255];
	    static const std::string testImages("TestImages.fits");
	    for ( size_t iComp(0); iComp < m_testSourceCaches.size(); iComp++ ) {
	      sprintf(buffer,"I_%d_%d_%d",iy,ix,iComp);
	      m_testSourceCaches[iComp]->writeTestSourceToFitsImage(testImages,buffer);
	    }
	  }
	}
	
	// Set the cache to do a broadband fit
	m_cache->setEnergyBin(-1);	
	status = m_cache->fitCurrent(false,verbose_bb());
	if ( status != 0 ) {
	  nfailed_bb_newton++;
	  continue;
	}
	
	// get the TS value and copy to the output histogram
	double tsval_newton = 2*(m_cache->currentLogLike() - loglike_null);
	double normVal = m_cache->currentPars()[m_cache->testSourceIndex()];
	double posErr(0.);
	double negErr(0.);

	m_cache->signalUncertainty_quad(0.5,posErr,negErr);
	
	if ( ts_map ) ts_map->setBinDirect(ipix,tsval_newton);
	if ( norm_map ) norm_map->setBinDirect(ipix,normVal);
	if ( posErr_map ) posErr_map->setBinDirect(ipix,posErr);
	if ( negErr_map ) negErr_map->setBinDirect(ipix,negErr);

	// if we are not doing the SED, we can move the next grid location
	if ( ! doSED ) continue;

	// Space for the output of the SED scan	
	std::vector<double> norm_mles;
	std::vector<double> pos_errs;
	std::vector<double> neg_errs;
	std::vector<double> logLike_mles; 
	std::vector<std::vector<double> > logLikes;
	std::vector<std::vector<double> > norms;

	// Space for the output of the SED (with ScienceTools fitting)
	std::vector<double> norm_mles_st;
	std::vector<double> logLike_mles_st;
	std::vector<std::vector<double> > logLikes_st;
	std::vector<std::vector<double> > norms_st;

	if ( sed_st ) {
	  Likelihood::ScanUtils::sed_binned(*m_modelWrapper,
					    m_testSourceName,
					    *m_opt,
					    tol,tolType,
					    nNorm_st,
					    norms_st,
					    norm_mles_st,
					    logLike_mles_st,
					    logLikes_st);
	}

	status = sed_binned_newton(nNorm,normSigma,covScale,
				   norm_mles,pos_errs,neg_errs,
				   logLike_mles,
				   norms,logLikes);
	
	if ( status != 0 ) {
	  nfailed_scan_newton++;
	  if ( status < 0 ) { 
	    std::cout << 'z' << std::flush;
	    continue;
	  } else if ( status < 10 ) {
	    std::cout << status << std::flush;
	  } else {
	    std::cout << 'x' << std::flush;
	  }
	} else {
	  std::cout << '.' << std::flush;
	}

	// This block copies the SED scan data to the output histograms

	// This is the pixel index
	int idx_sed = ipix;
	// This is the stride from one normalization set to the next
	int step_norm = nPixels() * nEBins();
	// Loop on energy bins
	for ( int iE(0); iE < nEBins(); iE++, idx_sed += nPixels() ) {
	  // Fill the histograms that use pixel and energy bin
	  double ts_val_bin = 2.*(logLike_mles[iE] - logLikes[iE][0]);
	  if ( ts_cube ) ts_cube->setBinDirect(idx_sed,ts_val_bin);
	  if ( norm_cube ) norm_cube->setBinDirect(idx_sed,norm_mles[iE]);
	  if ( posErr_cube ) posErr_cube->setBinDirect(idx_sed,pos_errs[iE]);
	  if ( negErr_cube ) negErr_cube->setBinDirect(idx_sed,neg_errs[iE]);
	  if ( nll_cube ) nll_cube->setBinDirect(idx_sed,logLike_mles[iE]);
	  if ( sed_st ) {
	    double ts_val_bin_st = 2.*(logLike_mles_st[iE] - logLikes_st[iE][0]);
	    if ( ts_cube_st ) ts_cube_st->setBinDirect(idx_sed,ts_val_bin_st);
	    if ( norm_cube_st ) norm_cube_st->setBinDirect(idx_sed,norm_mles_st[iE]);
	    if ( nll_cube_st ) nll_cube_st->setBinDirect(idx_sed,logLike_mles_st[iE]);
	  }
	  // This is the index for the pixel,energy bin
	  int idx_norm = idx_sed;
	  // Loop on normalization scan points
  	  for ( int iN(0); iN < nNorm; iN++, idx_norm += step_norm ) {
	    // Fill the histograms that use pixel, energy bin and normalization
	    double deltaLogLike = logLikes[iE][iN] - logLike_mles[iE] ;
 	    if ( norm_vals ) norm_vals->setBinDirect(idx_norm,norms[iE][iN]);
	    if ( delta_ll_vals ) delta_ll_vals->setBinDirect(idx_norm,deltaLogLike);
	    if ( normscan_st ) {
	      double deltaLogLike_st = logLikes_st[iE][iN] - logLike_mles_st[iE];
	      if ( norm_vals_st ) norm_vals_st->setBinDirect(idx_norm,norms_st[iE][iN]);
	      if ( delta_ll_vals_st ) delta_ll_vals_st->setBinDirect(idx_norm,deltaLogLike_st);
	    }
	  }	  
	}
	if ( broadband_st || remakeTestSource ) {
	  removeTestSourceFromModel();
	}
	m_cache->removeTestSourceFromCurrent();
      }
      std::cout << std::endl;      
    }

    // Warn about failed fits
    if ( nfailed_bb > 0 ) {
      std::cout << "There were " << nfailed_bb << " failed broadband fits with the standard fitter." << std::endl;
    }
    if ( nfailed_bb_newton > 0 ) {
      std::cout << "There were " << nfailed_bb_newton << " failed broadband fits with the Newton's method fitter." << std::endl;
    }
    if ( nfailed_scan > 0 ) {
      std::cout << "There were " << nfailed_scan << " failed normalization scans with the standard fitter." << std::endl;
    }
    if ( nfailed_scan_newton > 0 ) {
      std::cout << "There were " << nfailed_bb_newton << " failed normalization scans with the Newton's method fitter." << std::endl;
    }
    return 0;
  }


  /* Write the stored data to a FITS file */
  int FitScanner::writeFitsFile(const std::string& fitsFile,
				const std::string& creator,
				const std::string& fits_template,
				bool copyGTIs) const {

    int status(0);

    // Make the fits file
    // FIXME, should we check to see if it already exists?
    tip::IFileSvc::instance().createFile(fitsFile, fits_template);
    
    // Is this a WCS-based or HEALPix based scan
    bool image_based = m_dir2_binner != 0;    

    // Names of the columns we will write to the results table
    std::vector<std::pair<std::string,std::pair<HistND*,std::string> > > colData;

    static const std::string TSMAP_NAME("TSMAP");

    // Loop over the histogams and pull out the ones that 
    // can be represented as images
    // I.e., WCS-based maps or cubes
    // Push the other histograms onto the list of columns 
    // to save in the results table
    for ( std::vector< std::pair< std::string,std::pair<HistND*,std::string> > >::const_iterator itrHists = m_scanData.begin();
	  itrHists != m_scanData.end(); itrHists++ ) {
      if ( image_based && 
	   itrHists->second.second.find("NORM") == std::string::npos ) {
	if ( itrHists->first == TSMAP_NAME ) {
	  // Write the TSMAP to the primary
	  status = writeFitsImage(fitsFile,"",*(itrHists->second.first));
	} else {
	  status = writeFitsImage(fitsFile,itrHists->first,*(itrHists->second.first));
	}
	if ( status ) {
	  throw std::runtime_error("Failed to write fits image");
	  return status;
	}
      } else {
	if ( itrHists->first != TSMAP_NAME ) {
	  colData.push_back(*itrHists);
	}
      }
    }
	      
    // Ok, now save the remaining histograms to a single table
    status = writeFitsTable_byPixel(fitsFile,"SCANDATA",colData);
    if ( status ) {
      throw std::runtime_error("Failed to write fits table for scan data");
      return status;
    }

    // Save the energy bin edges
    // FIXME we could write flux conversion parameters here
    status = writeFits_EnergyBins(fitsFile);
    if ( status ) {
      throw std::runtime_error("Failed to write energy bins");
      return status;
    }

    // Write the GTIs, if requested
    if ( ! copyGTIs) return 0;
    status = writeFits_GTIs(fitsFile);
    if ( status ) {
      throw std::runtime_error("Failed to write GTIs");
      return status;
    }
    
    return 0;
    }
  
  // The number of pixels we are scanning over
  int FitScanner::nPixels() const {
    int retVal(1);
    retVal *= m_dir1_binner->getNumBins();
    // HEALPix only uses a single index
    if ( m_dir2_binner ) {
      retVal *= m_dir2_binner->getNumBins();
    }
    return retVal;
  }
  
  // The number of energy bins we are scanning over
  int FitScanner::nEBins() const {
    return m_energy_binner ? m_energy_binner->getNumBins() : 0;
  }
  
  // The number of normalization values we are scanning over
  int FitScanner::nNorms() const {
    return m_norm_binner ? m_norm_binner->getNumBins() : 0;
  }  
  
  /* Use a powerlaw point source */ 
  int FitScanner::setPowerlawPointTestSource(optimizers::FunctionFactory& factory, double index) {
    std::map<std::string, Source *>::const_iterator itrFind = m_modelWrapper->getMasterComponent().sources().find(m_testSourceName);
    if ( itrFind != m_modelWrapper->getMasterComponent().sources().end() ) {
      removeTestSourceFromModel();      
    }
    // delete the old test source
    delete m_testSource;
    m_testSource = 0;
    m_testSource = new Likelihood::PointSource();
    m_testSource->setName(m_testSourceName);
    optimizers::Function * pl = factory.create("PowerLaw");
    double parValues[] = {1., -1.*index, 100.};
    std::vector<double> pars(parValues, parValues + 3);
    pl->setParamValues(pars);
    optimizers::Parameter indexParam = pl->getParam("Index");
    indexParam.setBounds(-3.5, -1.);
    pl->setParam(indexParam);
    optimizers::Parameter prefactorParam = pl->getParam("Prefactor");
    prefactorParam.setBounds(1e-10, 1e3);
    prefactorParam.setScale(1e-9);
    pl->setParam(prefactorParam);
    m_testSource->setSpectrum(pl);   
    int status = buildTestModelCache();
    if ( status ) {
      throw std::runtime_error("Failed FitScanner::buildTestModelCache().");
      return -1;
    }
    return 0;    
  }


  /* Use a source from the model */
  int FitScanner::setTestSourceByName(const std::string& sourceName) {
    
    // Protect against the case that the sourceName is the same as the m_testSourceName
    std::map<std::string, Source *>::const_iterator itrFind;
    if ( m_testSourceName != sourceName ) {
      // Remove the old source, if it exists
      itrFind = m_modelWrapper->getMasterComponent().sources().find(m_testSourceName);
      if ( itrFind != m_modelWrapper->getMasterComponent().sources().end() ) {
	removeTestSourceFromModel();      
      }
      // Delete the old test soruce
      delete m_testSource;

      // Check to see if the new source in is the model 
      m_testSourceName = sourceName;
    }

    itrFind = m_modelWrapper->getMasterComponent().sources().find(m_testSourceName);
    if ( itrFind == m_modelWrapper->getMasterComponent().sources().end() ) {
      throw std::runtime_error("FitScanner::setTestSourceByName could not find test source in model");
      return -1;
    }
    m_testSource = itrFind->second;
    int status = buildTestModelCache();
    if ( status ) {
      throw std::runtime_error("Failed FitScanner::buildTestModelCache().");
      return -1;
    }
    return 0;        
  }


  /* This adds the test source to the source model */
  int FitScanner::addTestSourceToModel() {
    PointSource* ptrSrc = dynamic_cast<PointSource*>(m_testSource);
    if ( ptrSrc == 0 ) {
      throw std::runtime_error("Test source must be a point source, for now...");
      return -1;
    }
    ptrSrc->setDir(m_testSourceDir.ra(), m_testSourceDir.dec(),
		   false, false);
    Likelihood::ScanUtils::freezeSourceParams(*m_testSource);         
    optimizers::Parameter& normPar = m_testSource->spectrum().normPar();
    normPar.setValue( 1e-10 );
    normPar.setFree(true); 

    m_modelWrapper->addSource(ptrSrc);
    m_modelWrapper->syncParams();
    return 0;
  }

  /* This removes the test source from the source model */
  void FitScanner::removeTestSourceFromModel() {
    m_modelWrapper->removeSource(m_testSourceName);    
    return;
  }

  /* Set the direction of the test source, based on the loop parameters */
  int FitScanner::setTestSourceDir(int ix, int iy) {
    // No projection means we are not looping over direction
    if ( m_proj == 0 ) { 
      if ( ix != 0 ||
	   iy != 0 ) {
	
      }            
      return 0;
    }
    // Get the currect pixel indices
    double xpix = m_dir1_binner->getInterval(ix).midpoint();
    double ypix = m_dir2_binner ? m_dir2_binner->getInterval(iy).midpoint() : 0.;
    // astro::SkyDir doesn't have assignment operator, so we do this
    m_testSourceDir() = astro::SkyDir(xpix,ypix,*m_proj,false)();   
    return 0;
  }

  /* This does the baseline fit
     i.e., the fit without the test source */
  int FitScanner::baselineFit(double tol, int tolType) {
    int status = m_opt->find_min(false,tol,tolType);
    // FIXME, we should save the fit here, and come back to it
    // m_logLike->saveCurrentFit();
    return status;
  } 


  /* This does the baseline fit with Newton's Method,
     for the normalization parameters only */
  int FitScanner::baselineFit_Newton(double tol,int maxIter) {
    static const bool useReduced(true);
    if ( m_cache == 0 ) {
      // Build the cache
      m_cache = new FitScanCache(*m_modelWrapper,m_testSourceName,tol,maxIter,useReduced); 
    } else {
      // Reset the cache (all sources free, no test source)
      std::vector<bool> freePars(m_cache->nBkgModel(),true);
      std::vector<float> pars_scales(m_cache->nBkgModel(),1.0);
      m_cache->refactorModel(freePars,pars_scales,false);
    }
    int ok = m_cache->fitCurrent(false);
    return ok;
  }


  /* This does the broadband fit
     i.e., the fit with the source across the entire energy range */  
  int FitScanner::fitTestSourceBroadband(double tol, int tolType) {
    // Just pass it along to the optimizer
    int status = m_opt->find_min_only(false,tol,tolType); 
    return status;
  }


  /* This does the sed fitting with Newton's method,                                                       
     for the normalization paramters only */  
  int FitScanner::sed_binned_newton(int nnorm, double normSigma, 
				    double constrainScale,
				    std::vector<double>& norm_mles,
				    std::vector<double>& pos_errs,
				    std::vector<double>& neg_errs,
				    std::vector<double>& logLike_mles,
				    std::vector<std::vector<double> >& norms,
				    std::vector<std::vector<double> >& logLikes) {


    static const bool redoFailedVerbose(false);

    const double errorLevel = 0.5*normSigma*normSigma;

    // We can't do the fitting without a FitScanCache
    if ( m_cache == 0 ) return -1;

    // first we fix everything except the signal component to their current values 
    std::vector<float> par_scales;
    m_cache->getParScales(par_scales);
    bool usePrior(false);

    if ( constrainScale < 0 ) {
      std::vector<bool> freeSources(m_cache->nBkgModel(),false);
      m_cache->refactorModel(freeSources,par_scales,true);
    } else {
      usePrior = true;
      std::vector<bool> constrainPars(m_cache->nBkgModel(),true);
      m_cache->buildPriorsFromCurrent(constrainPars,constrainScale);
    }
    
    // Latch the index of the test source
    int test_idx = m_cache->testSourceIndex();

    // Allocate the output vectors
    norm_mles.resize(m_cache->nebins());
    logLike_mles.resize(m_cache->nebins());
    pos_errs.resize(m_cache->nebins());
    neg_errs.resize(m_cache->nebins());
    norms.resize(m_cache->nebins());    
    logLikes.resize(m_cache->nebins());
        
    // This is to keep track of failed fits.
    // Usually they just have to do with problem
    // with the matrix inversion, so we might not want to crash
    int nfailed(0);

    // Loop on the energy bins
    for ( size_t i(0); i < m_cache->nebins(); i++ ) {
      m_cache->setEnergyBin(i);
      int status = m_cache->fitCurrent(usePrior);
      if ( status ) {
	// if the fit failed, fill the output vectors, and move on.
	// for debugging, redo failed fits with verbose on
	if ( redoFailedVerbose ) {
	  m_cache->fitCurrent(4);
	}
	nfailed++;
	logLike_mles[i] = 0.;
	norm_mles[i] = -1.;
	pos_errs[i] = 0.;
	neg_errs[i] = 0.;
	norms[i].resize(nnorm,0.);
	logLikes[i].resize(nnorm,0.);
	continue;
      }
      // latch the information for the output vectors
      logLike_mles[i] = m_cache->currentLogLike();
      norm_mles[i] = m_cache->currentPars()[test_idx];
      m_cache->signalUncertainty_quad(0.5,pos_errs[i],neg_errs[i]);
      double negLim(0.);
      double posLim(0.);
      
      m_cache->signalUncertainty_quad(errorLevel,posLim,negLim);
      m_cache->scanNormalization(nnorm,1.0,posLim,negLim,norms[i],logLikes[i]);
      if ( false ) {
	std::cout << "Done scan " << i << ' ' << norm_mles[i] 
		  << ' ' << norms[i].back() << ' ' << ( logLikes[i].back() - logLike_mles[i] ) << std::endl;
      }
    }

    // Reset the cache to do broadband fitting
    m_cache->setEnergyBin(-1);
    return nfailed;
  } 

  /* Build and cache an image of the test source */
  int FitScanner::buildTestModelCache() {
    deleteTestModelCaches();
    if ( m_testSource == 0 ) {
      throw std::runtime_error("FitScanner needs a test source to build a test source cache.");
      return -1;
    }
    PointSource* ptrSrc = dynamic_cast<PointSource*>(m_testSource);
    if ( ptrSrc == 0 ) {      
      throw std::runtime_error("Test source must be a point source, for now...");
      return -1;
    }
    
    for ( size_t iComp(0); iComp < m_modelWrapper->numComponents(); iComp++ ) {      
      const BinnedLikelihood* binnedLike = m_modelWrapper->getComponent(iComp);      
      if ( binnedLike == 0 ) {
	throw std::runtime_error("FitScanner requires using Binned likelihood fitting.");
	return -1;
      }
      // astro::SkyDir doesn't have assignment operator, so we do this
      m_testSourceDir() = binnedLike->countsMap().refDir()();
      
      // Now set the direction in the source.   
      ptrSrc->setDir(m_testSourceDir.ra(), m_testSourceDir.dec(),
		     false, false);
      // Build the cache with the new test source
      m_testSourceCaches[iComp] = new TestSourceModelCache(*binnedLike,*ptrSrc);
    }
    return 0;
  }

  /* Build an n-dimensional histogram based on the loop parameters */
  HistND* FitScanner::buildHist(const std::string& name, 
				bool do_pix,
				bool do_energy,
				bool do_norm) {

    // We use this to keep track of the axes of the histograms by name
    std::string dimString;
    // The binners used to define the axes
    std::vector<evtbin::Binner*> binners;
    if ( do_pix ) {
      binners.push_back(const_cast<evtbin::Binner*>(m_dir1_binner));
      if ( m_dir2_binner ) {
	binners.push_back(const_cast<evtbin::Binner*>(m_dir2_binner));
      }
      dimString += "PIX:";
    }
    if ( do_energy ) {
      if ( m_energy_binner == 0 ) {
	throw std::runtime_error("No Energy binner in FitScanner::buildHist");
	return 0;
      }
      binners.push_back(const_cast<evtbin::Binner*>(m_energy_binner));
      dimString += "ENERGY:";
    }
    if ( do_norm ) {
      if ( m_norm_binner == 0 ) {
	throw std::runtime_error("No normalization binner in FitScanner::buildHist");
	return 0;
      }
      binners.push_back(const_cast<evtbin::Binner*>(m_norm_binner));
      dimString += "NORM:";
    }
    // Nothing was actually requested, return a null hist
    if ( binners.size() == 0 ) {
      return 0;
    }
    // Ok, build the histogram from the binners
    HistND* hist = new HistND(binners);
    std::pair<HistND*,std::string> pair(hist,dimString);
    m_scanData.push_back(std::pair<std::string,std::pair<HistND*,std::string> >(name,pair) );
    return hist;
  }

  /* Write a histogram as a FITS image */
  int FitScanner::writeFitsImage(const std::string& fitsFile,
				 const std::string& extName,
				 const HistND& hist) const {
    
    // Get the image axes
    long nxbins = m_dir1_binner->getNumBins();
    long nybins = m_dir2_binner ? m_dir2_binner->getNumBins() : 0;
    if ( nxbins == 0 || nybins == 0 ) {
      throw std::runtime_error("Histogram is not an image or series of images.");
      return -1;
    }

    // Get the data
    const std::vector<const evtbin::Binner*>& binners = hist.getBinners();
    // Make sure it is 2 or 3D (map or cube)
    if ( binners.size() < 2 || binners.size() > 3 ) {
      throw std::runtime_error("Wrong number of histogram axes");
      return -1;
    }
    
    if ( binners[0]->getNumBins() != nxbins ) {
      throw std::runtime_error("X-axis size does not match.");
      return -1;
    }
    
    if ( binners[1]->getNumBins() != nybins ) {
      throw std::runtime_error("Y-axis size does not match.");
      return -1;
    }

    // Get the axes size
    std::vector<long> naxes;
    naxes.push_back(nxbins);
    naxes.push_back(nybins);
    
    if ( binners.size() == 3 ) {
      naxes.push_back( binners[2]->getNumBins() );
    }
   
    // Add an image to the file
    if ( ! extName.empty() ) {
      tip::IFileSvc::instance().appendImage(fitsFile,extName,naxes);
    }
    tip::Image* image(tip::IFileSvc::instance().editImage(fitsFile, extName));

    // Set the image dimension and keywords
    image->setImageDimensions(naxes);
    tip::Header & header(image->getHeader());
    // FIXME, astro::ProjBase::setKeywords() should be a const function
    astro::ProjBase* nc_proj = const_cast<astro::ProjBase*>(m_proj);
    if ( nc_proj ) {
      nc_proj->setKeywords(header);
    }
    
    if (  binners.size() == 3 ) {
      header.setKeyword("CRPIX3",1);
      header.setKeyword("CRVAL3",binners[2]->getInterval(0).begin());
      header.setKeyword("CDELT3",binners[2]->getInterval(0).width());
      header.setKeyword("CTYPE3","photon energy");
    }
    
    // actually set the data
    image->set(hist.data());
    // deleting the image is the thing that actually writes the data to the file
    delete image;

    return 0;
  }
  

  /* Write a histogram as a FITS table */
  int FitScanner::writeFitsTable_byPixel(const std::string& fitsFile,
					 const std::string& extName,
					 const std::vector<std::pair<std::string,std::pair<HistND*,std::string> > >& colData) const {


    // We need to know the size of the table
    int ncol = colData.size();
    int npix = nPixels();

    // The table is already defined in the template, so we can just edit it
    tip::Table *table = tip::IFileSvc::instance().editTable(fitsFile, extName);  
    // Resize the table to have as many records as there are pixels
    table->setNumRecords(npix);

    tip::Header & header(table->getHeader());
    // FIXME, astro::ProjBase::setKeywords() should be a const function
    astro::ProjBase* nc_proj = const_cast<astro::ProjBase*>(m_proj);
    if ( nc_proj ) {
      nc_proj->setKeywords(header);
    }      
    
    // Loop over the columns
    for ( int icol(0); icol < ncol; icol++ ) {
      const std::string& colName = colData[icol].first;
      HistND* hist = colData[icol].second.first;
      const std::string& dimString =  colData[icol].second.second;
      const std::vector<float>& histData = hist->data();
      // We need to know the shape of the data in each column
      int nDataTotal = histData.size();
      int nData = nDataTotal / npix;      
      char buffer[255];
      sprintf(buffer,"%iE",nData);
      std::string formatString(buffer);
      table->appendField(colName,formatString); 
      std::string pixelDimString;
      bool has_dim = convertDimString(dimString,pixelDimString);
      if ( has_dim ) {
	setDimKeyword(header,icol,pixelDimString);
      }
      // Now add the column to the table, and copy the data into the column
      tip::IColumn* col = table->getColumn(table->getFieldIndex(colName));
      std::vector<float> writeValue(nData);
      for ( int ipix(0); ipix < npix; ipix++ ) {
	int idx(ipix);
	for ( int idata(0); idata < nData; idata++, idx+=npix ) {
	  writeValue[idata] = histData[idx];
	}
	col->set(ipix,writeValue);
      } 
    }

    delete table;
    return 0;
  }

  /* write the energy bins */
  int FitScanner::writeFits_EnergyBins(const std::string& fitsFile) const {
    if ( m_energy_binner == 0 ) {
      return 0;
    }
    return m_modelWrapper->writeFits_EnergyBins(fitsFile,m_energy_binner);
  }
  
  /* write the Good time intervals */
  int FitScanner::writeFits_GTIs(const std::string& fitsFile) const {  
    return m_modelWrapper->writeFits_GTIs(fitsFile);
  }

  /* Convert the dimension string to the format expected by FITS */
  bool FitScanner::convertDimString(const std::string& inString,
				    std::string& outString, 
				    bool do_pix,
				    bool do_energy,
				    bool do_norm) const {
    std::vector<std::string> tokens;
    outString.clear();
    static const std::string delims(":");
    facilities::Util::stringTokenize(inString,delims,tokens);
    outString = "(";
    int nd(0);
    for ( std::vector<std::string>::const_iterator itr = tokens.begin();
	  itr != tokens.end(); itr++ ) {
      char buffer[255];
      if ( do_pix && *itr == "PIX" ) {
	sprintf(buffer,"%i,",nPixels());
	outString += buffer;
	nd++;
      } else if ( do_energy && *itr == "ENERGY" ) {
	sprintf(buffer,"%i,",nEBins());
	outString += buffer;
	nd++;
      } else if ( do_norm && *itr == "NORM" ) {
	sprintf(buffer,"%i,",nNorms());
	outString += buffer;
	nd++;
      }
    }
    outString += ")";
    return nd > 1;
  }

  /* set the TDIM keyword */
  void FitScanner::setDimKeyword(tip::Header& header,
				 int icol,
				 const std::string& dimString) const {    
    char buf[255];
    sprintf(buf,"TDIM%i",icol+1);
    header.setKeyword(buf,dimString);
  }

  
  void FitScanner::deleteTestModelCaches() {
    for ( std::vector<TestSourceModelCache*>::iterator itrDel = m_testSourceCaches.begin();
	  itrDel != m_testSourceCaches.end(); itrDel++ ) {
      TestSourceModelCache* toDel = *itrDel;
      delete toDel;
      *itrDel = 0;
    }
  }


} // namespace Likelihood
