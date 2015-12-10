/**
 * @file Convolve.cxx
 * @brief Functions to perform convolutions of HEALPix maps
 * @author E. Charles
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/src/ConvolveHealpix.cxx,v 1.2 2015/03/03 06:00:00 echarles Exp $
 */


#include "Likelihood/ConvolveHealpix.h"
#include <set>
#include <cmath>
#include "alm_healpix_tools.h"
#include "alm.h"

namespace Likelihood {

  namespace ConvolveHealpix {

    void fillMapWithPSF_pole(const MeanPsf & psf, 
			     const double& energy, 
			     Healpix_Map<float>& hp) {
      // This takes advantage of the RING scheme, we just loop on the rings of the 
      // map until we reach the boundtry of the PSF
      // 
      // This is copied from gardian/Psf.cc
      hp.fill(0.);
      int nRing = 1;
      int startPix, nPixRing;
      double theta(0);
      bool shifted(false);
      int order = hp.Order();
      if ( order < 0 ) {
	throw std::runtime_error("ConvolveHealpix::fillMapWithPSF_pole"
				 "only works with standard nside = 2^order maps");
	return;
      }

      // This is the limit we fill out to
      double th  = astro::degToRad(1000.*psfMinPixSize(energy));
      // This is the fraction of the PSF in that limit
      double fraction = psf.integral( astro::radToDeg(th),energy );
      while ( theta < th ){
	hp.get_ring_info2(nRing, startPix, nPixRing, theta, shifted);
	//Difference in order to get the finer structure of the psf
	const int dOrder = std::max(std::min(13-order, 6 - 2*(nRing-1)),1);
	const int nPix = (1<<dOrder)*(1<<dOrder);
	//Create a finer healpix grid to get a better representation for the psf
	Healpix_Base hpf(order+dOrder, NEST);
	//Loop over the relevant pixels in the finer grid and set their value with
	//interpolation
	int sp = hp.ring2nest(startPix);
	double sum = 0;
	for (int i = sp*nPix; i < (sp+1)*nPix; ++i){
	  double cdist = astro::radToDeg(hpf.pix2ang(i).theta);
	  //Do linear interpolation although it isn't linear.  Finer binning in the
	  //input psf should solve that.
	  sum += std::max( psf(energy,cdist ), 0.0);
	}
	sum /= nPix;
	//Add the value to all the pixels in a ring
	for (int i = startPix; i < startPix+nPixRing; ++i){
	  hp[i] = sum;
	}
	++nRing;
      }
      //Integrate over the psf to assure proper normalization
      double integral_per_pixel = 1./(ASTRO_4PI*hp.average());
      //std::cout<<"PSF integral = "<<integral<<std::endl;
      //Scale the map with the integral/fraction
      for (int i=0; i<hp.Npix(); ++i){
	hp[i] *= integral_per_pixel*fraction;
      }
    }

    
    void fillMapWithPSF_refDir(const MeanPsf& psf, 
			       const double& energy, 
			       const astro::SkyDir& refDir,
			       bool use_lb,
			       Healpix_Map<float>& hp){
      // This does not take advantage of the ring scheme
      // but rather explicitly gets the location of all the sub-pixels for the integration
      // 
      // This is copied from gardian/Psf.cc

      int order = hp.Order();
      if ( order < 0 ) {
	throw std::runtime_error("ConvolveHealpix::fillMapWithPSF_refDir"
				 "only works with standard nside = 2^order maps");
	return;
      }
    
      // make a pointing object
      double refTheta = astro::degToRad ( astro::latToTheta_Deg( use_lb ? refDir.b() : refDir.dec() ) );
      double refPhi = astro::degToRad( use_lb ? refDir.l() : refDir.ra() );
      const pointing center(refTheta,refPhi);
      // This is the limit we fill out to
      double th = astro::degToRad(500*psfMinPixSize(energy));
      // This is the fraction of the PSF in that limit
      double fraction = psf.integral( astro::radToDeg(th),energy );

      //Find the pixels associated with the psf function
      std::vector<int> listpix;
      hp.query_disc_inclusive(center, th, listpix);

      //Use a finer grid in the neighboring pixels
      int cp = hp.ang2pix(center);
      fix_arr<int,8> nb;
      hp.neighbors(cp, nb);

      //Create a set of the pixel numbers
      std::set<int> nbset(&nb[0], &nb[0]+8);
      //Remove -1 from the set
      nbset.erase(-1);
      //Set for the total pixels
      std::set<int> listset(listpix.begin(),listpix.end());
      
      //Exclude the center pixel and its neighbors from the total
      //pixels.  We also only use the intersection of the neighbor
      //pixels and the total pixels, in case we don't want to include
      //all of the neighbors. (For a very narrow psf or small
      //fraction).
      listset.erase(cp);
      std::set<int> listReduced, nbReduced;
      std::set_difference(listset.begin(), listset.end(), nbset.begin(), nbset.end(), 
			  std::insert_iterator<std::set<int> >(listReduced,listReduced.begin()));
      std::set_intersection(listset.begin(), listset.end(), nbset.begin(), nbset.end(), 
			    std::insert_iterator<std::set<int> >(nbReduced,nbReduced.begin()));
      
      //Build the psf map
      vec3 cvec = center.to_vec3();
      std::set<int>::iterator it = listReduced.begin();
      for ( ; it != listReduced.end() ; ++it){
	// This will break if we go to extremely high resolution
	double cdist = astro::radToDeg(acos(dotprod(cvec, hp.pix2ang(*it).to_vec3())));
	hp[*it] = std::max( psf(energy,cdist), 0.0 );
      }
      //Neighbors are next, increase the order by 3
      it = nbReduced.begin();
      if (hp.Scheme() == RING) {
	for ( ; it != nbReduced.end(); ++it) {
	  const int dOrder = std::min(13-order, 3);
	  const int nPix = (1<<dOrder)*(1<<dOrder);
	  Healpix_Base hpf(order+dOrder, NEST);
	  const int sp = hp.ring2nest(*it);
	  double sum = 0;
	  for (int i = sp*nPix; i < (sp+1)*nPix; ++i) {
	    double cdist = astro::radToDeg(acos(dotprod(cvec, hpf.pix2ang(i).to_vec3())));
	    sum += std::max( psf(energy,cdist), 0.0 );
	  }
	  hp[*it] = sum/nPix;
	}
      } else {
	for ( ; it != nbReduced.end(); ++it) {
	  const int dOrder = std::min(13-order, 3);
	  const int nPix = (1<<dOrder)*(1<<dOrder);
	  Healpix_Base hpf(order+dOrder, NEST);
	  const int sp = *it;
	  double sum = 0;
	  for (int i = sp*nPix; i < (sp+1)*nPix; ++i) {
	    double cdist = astro::radToDeg(acos(dotprod(cvec, hpf.pix2ang(i).to_vec3())));
	    sum += std::max( psf(energy,cdist), 0.0 );
	  }
	  hp[*it] = sum/nPix;
	}
      }

      //And the center pixel, increase the order by 6
      const int dOrder = std::min(13-order, 6);
      const int nPix = (1<<dOrder)*(1<<dOrder);
      Healpix_Base hpf(order+dOrder, NEST);
      int sp;
      if (hp.Scheme() == RING)
	sp = hp.ring2nest(cp);
      else
	sp = cp;
      double sum = 0;
      for (int i = sp*nPix; i < (sp+1)*nPix; ++i) {
	double cdist = astro::radToDeg(acos(dotprod(cvec, hpf.pix2ang(i).to_vec3())));
	// std::cout << i << ' ' << cdist << std::endl;
	sum += std::max( psf(energy,cdist), 0.0 );
      }
      hp[cp] = sum/nPix;
      
      //Integrate over the psf to assure proper normalization
      double integral_per_pixel = 1./(hp.average()*hp.Npix());
      //std::cout<<"PSF integral = "<<integral<<std::endl;
      //Scale the map with the integral/fraction
      for (int i=0; i<hp.Npix(); ++i){
	hp[i] *= integral_per_pixel*fraction;
      }
    } 

    void convolve(const Healpix_Map<float>& inMap,
		  const Healpix_Map<float>& psf,
		  Healpix_Map<float>& outMap) {
      //
      // This is copied from gardian/Psf.cc
      //
      
      const int lmax = 3*outMap.Nside()-1; //3 seems to be the magic number

      //Average of the map, needed for scaling the convolution
      double avgBefore = inMap.average();

      //Convert the psf to Healpix, use a higher order to get a better
      //convolution

      Alm<xcomplex<float> > psfalm(lmax,0);  //Symmetric so only m = 0 needed.
      Alm<xcomplex<float> > mapalm(lmax,lmax);  //Not symmetric so m <= l needed.

      int niter = 3;
      map2alm_iter(psf, psfalm, niter);  //Using iterations requires us to scale the map to preserve counts
    

      //Harmonic transform map to alm
      map2alm_iter(inMap, mapalm, niter);

      //Multiply the alms
      for (int l = 0; l<=lmax; ++l) {
	float factor = sqrt(ASTRO_4PI/(2*l+1));
	for (int m = 0; m<=l; ++m) {
	  mapalm(l,m) *= factor*psfalm(l,0);
	}
      }
      
      //Convert alms back to the map
      alm2map(mapalm, outMap);

      //Scale the map, so number of counts are equal before and after
      double scale = avgBefore/outMap.average();
      for (int j = 0; j < outMap.Npix(); ++j){
	outMap[j] *= scale;
      }
      //Remove values from the map that are smaller in absolute value than the largest negative number
      float min, max;
      outMap.minmax(min,max);
      if (min < 0) {
	for (int j = 0; j < outMap.Npix(); ++j) {
	  if (outMap[j] < fabs(min)) {
	    outMap[j] = 0;
	  }
	} 
      }
    }

    double psfMinPixSize(double energy) {
      // This is the point beyond which it doesn't really
      // make sense to make smaller pixels
      static double min_angle(0.05*M_PI/180.);
      static double ebreak(1e4);
      if (energy > ebreak) {
	/// PSF containment flattens above ebreak.
	return min_angle;
      }
      /// Empirical power-law scaling.
      return min_angle*std::pow(energy/ebreak, -0.8);      
    }


  } // namespace ConvolveHealpix
 
} // namespace Likelihood
