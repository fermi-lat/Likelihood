/** 
 * @file HistND.h
 * @brief N-dimensional histogram.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/Likelihood/HistND.h,v 1.3 2015/03/03 05:59:56 echarles Exp $
 */

#ifndef Likelihood_HistND_h
#define Likelihood_HistND_h

#include <algorithm>
#include <vector>

#include "evtbin/Hist.h"

namespace evtbin {
   class Binner;
}

namespace Likelihood {

/** 
 * @class HistND
 * @brief N dimensional histogram.
 *
 */
class HistND : public evtbin::Hist {
public:
    
   /// @brief Create a three dimensional histogram which uses the given 
   /// binner objects to determine the indices.
   /// @param binners vector of references to binners.  The size of this 
   /// vector sets the dimensionality of the histogram.
   HistND(const std::vector<evtbin::Binner *> & binners);

   HistND(const HistND & rhs);
   
   virtual ~HistND() throw();
     
   /// @brief Increment the bin appropriate for the given value.  This
   /// is generic for N-dimensional histograms.
   /// @param value Vector giving the value being binned. The vector
   /// must have at least as many values as the dimensionality of the
   /// histogram.
   virtual void fillBin(const std::vector<double> & values, 
                        double weight = 1.);

   inline void setBinDirect(long ibin,double weight = 1.) { m_data[ibin] = weight; }

   /// @return The unrolled index of the point in N-dimensional space
   /// corresponding to the specified coordinates.  If the point is
   /// outside the bounded volume (minus an border region in pixels),
   /// return -1.
   /// @param values N-dimensional coordinates (x, y, z,...)
   /// @param border_size Border region to exclude in units of pixels.
   long binIndex(const std::vector<double> & values, long border_size=0) const;

   long binIndex(const std::vector<long> & ivalues) const;
   
   double operator[](long index) const;

   const std::vector<float> & data() const {return m_data;}
   
   void setData(const std::vector<float> & data) {m_data = data;}
   void setData(const std::vector<double> & data) {
      m_data.resize(data.size(), 0);
      std::copy(data.begin(), data.end(), m_data.begin());
   }

   /// set the data along a slice of the histogram
   /// 
   /// @param idim is the dimension to step along
   /// @param ivalues are index values for each dimemsion of the slice (ivalues[idim] is ignored)
   /// @param data is the input data    
   void setSlice(unsigned int idim,
		 const std::vector<unsigned int>& ivalues,
		 const std::vector<float>& data);
 
   /// set the data along a slice of the histogram
   /// 
   /// @param idim is the dimension to step along
   /// @param ivalues are index values for each dimemsion of the slice (ivalues[idim] is ignored)
   /// @param data is the input data    
   void setSlice(unsigned int idim,
		 const std::vector<unsigned int>& ivalues,
		 const std::vector<double>& data);


   /// get the data along a slice of the histogram
   /// 
   /// @param idim is the dimension to step along
   /// @param ivalues are index values for each dimemsion of the slice (ivalues[idim] is ignored)
   /// @param data is filled with the data from that slice
   void getSlice(unsigned int idim,
                 const std::vector<unsigned int>& ivalues,
                 std::vector<float>& data);

   /// get the data along a slice of the histogram
   /// 
   /// @param idim is the dimension to step along
   /// @param ivalues are index values for each dimemsion of the slice (ivalues[idim] is ignored)
   /// @param data is filled with the data from that slice
   void getSlice(unsigned int idim,
                 const std::vector<unsigned int>& ivalues,
                 std::vector<double>& data);
   
   /// Sum the values in a dimension of the histogram
   /// 
   /// @param idim is the dimension to sum over
   /// The sum runs from [firstBin,lastBin)
   HistND * sumRange(unsigned int idim, unsigned int firstBin, unsigned int lastBin) const;

   HistND * clone() const {return new HistND(*this);}


private:

   std::vector<float> m_data;
     
   unsigned int m_ndims;
     
   std::vector<unsigned int> m_strides;

};

} // namespace Likelihood

#endif // Likelihood_HistND_h
