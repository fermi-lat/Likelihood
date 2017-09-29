/**
 * @file SparseVector.h
 * @brief Implmentation of a sparse vector (similar to an std::map, but with less functionality, and using less memory)
 * @author E. Charles
 *

 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SparseVector.h,v 1.3 2016/11/01 23:49:20 echarles Exp $
 */

#ifndef Likelihood_SparseVector_h
#define Likelihood_SparseVector_h

#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace Likelihood {
  
  template <typename T>
  class SparseVector {

  public:

    /* This gets used in std::lower_bound */
    static bool key_compare_lower(const std::pair<size_t,T>& a,
				  size_t b){
      return a.first < b;
    }

    /* This gets used in std::uppwer_bound */
    static bool key_compare_upper(size_t a, 
				  const std::pair<size_t,T>& b){
      return a < b.first;
    }


  public:

    /* There are all of the types defined for a std::vector<T> */
    typedef size_t key_type;
    typedef T value_type;
    typedef T* pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef size_t size_type;
    typedef int difference_type;


    /* There are all of the types defined for a std::vector< std::pair<size_t,T> > */
    typedef typename std::vector< std::pair<size_t,T> >::iterator iterator;
    typedef typename std::vector< std::pair<size_t,T> >::const_iterator const_iterator;
    typedef typename std::vector< std::pair<size_t,T> >::reverse_iterator reverse_iterator;
    typedef typename std::vector< std::pair<size_t,T> >::const_reverse_iterator const_reverse_iterator;

    /* These are here for convienice */
    typedef typename std::pair<size_t,T> container_type;
    typedef typename std::pair<size_t,T>* container_pointer_type;
    typedef typename std::pair<size_t,T>& container_reference_type;
    typedef const typename std::pair<size_t,T>& container_const_reference_type;
       
  public:

    /* Null c'tor, creates an empty vector*/
    SparseVector(T null = 0);   

    /* Create a vector with n elements, all set to null */    
    SparseVector(size_type n, T null = 0);   
   
    /* Copies in a vector, removing any null elements */    
    SparseVector(const std::vector<T>& vect, T null = 0);

    /* Copy c'tor */    
    SparseVector(const SparseVector<T>& other);

    /* Assignment operator */
    SparseVector<T>& operator=(const SparseVector<T>&);

    /* -------------------- Functions from the std::vector<T> interaface ------------------ */

    /* Standard stl iterators. 
     Note: these only loop on non-zero values */
    inline iterator begin() { return m_data.begin(); }
    inline iterator end() { return m_data.end(); }
    inline const_iterator begin() const { return m_data.begin(); }
    inline const_iterator end()  const { return m_data.end(); } 

    /* Standard stl reverse interators 
       Note: these only loop on non-zero elements */
    inline reverse_iterator rbegin() { return m_data.rbegin(); }
    inline reverse_iterator rend() { return m_data.rend(); }
    inline const_reverse_iterator rbegin() const { return m_data.rbegin(); }
    inline const_reverse_iterator rend()  const { return m_data.rend(); }


    /* Size of the vector, including null elements */
    inline size_type size() const { return m_size; }

    /* Maximum possible size of the underlying data structure */
    inline size_type max_size() const { return m_data.max_size(); }
    /* Currently allocated size of the underlying data structure */
    inline size_type capacity() const { return m_data.capacity(); }
 
    /* Test if size() == 0 */
    inline bool empty() const { return m_size == 0; }
    
    /* First and last non-zero elements */
    container_reference_type front() { return m_data.front(); }
    container_const_reference_type front() const { return m_data.front(); }
    container_reference_type back() { return m_data.back(); }
    container_const_reference_type back() const { return m_data.back(); }

    /* Clear the underlying data struture and set the size to zero */
    void clear() { m_data.clear(); m_size = 0; }

    /* Clear the underlying data struture, but leave the size as is */
    void clear_data() { m_data.clear(); }
    

    /* -------------------- Functions from the std::vector< std::pair< size_t, T> > interaface ------------------ */

    /* Find a particular element */
    iterator find(size_type k);
    const_iterator find(size_type k) const;

    /* Find the first element whose key is not less than k */
    iterator lower_bound(size_type k);
    const_iterator lower_bound(size_type k) const;

    /* Find the first element whose key is greater than k */
    iterator upper_bound(size_type k);
    const_iterator upper_bound(size_type k) const;

    /* Return a particular element, or the null value if the element is not present */
    const_reference operator[](size_type n) const;    
    void push_back(const T& val);
    void swap(SparseVector<T>& other);
    void resize(size_type n);

    /* ------------------- Implementation specific functions ----------------- */

    /* Return the null value */
    inline const_reference null() const { return m_null; }

    /* Test if a particular value is null */
    inline bool test_value(const T& val) {  return val != m_null; }

    /* Clear out the non-null values */
    inline void set_to_null() { m_data.clear(); }

    /* Free the memory */
    void free_memory() {
      clear();
      std::vector< std::pair<size_t,T> > nullV;
      m_data.swap(nullV);
    }

    /* Return a const reference to the non-null elements */
    inline const std::vector<std::pair<size_t,T> >& non_null() const { return m_data; }
    
    /* Fill this from a vector */
    void fill_from_vect(const std::vector<T>& vect);

    /* Fill this from a map */
    void fill_from_map(const std::map<size_t,T>& map);

    /* Fill this from key and value vectors */
    void fill_from_key_and_value(const std::vector<size_t>& keys,
				 const std::vector<T>& values);
    

    /* Fill a vector */
    void fill_vect(std::vector<T>& vect) const;
    void fill_vect(std::vector<T>& vect, size_t first, size_t last) const;

    /* Fill a map */
    void fill_map(std::map<size_t,T>& vect) const;
    void fill_map(std::map<size_t,T>& vect, size_t first, size_t last) const;

    /* Fill key and value vectors */
    void fill_key_and_value(std::vector<size_t>& keys,
			    std::vector<T>& values) const;
    void fill_key_and_value(std::vector<size_t>& keys,
			    std::vector<T>& values, 
			    size_t first, size_t last) const;
    

    /* Add to a vector.  This does not add the non-null elements. */
    void add_to_vect(std::vector<T>& vect) const;
  
    /* Subtract from a vector.  This does not add the non-null elements.  */
    void subtract_from_vect(std::vector<T>& vect) const;
  
  private:

    /* The full size of the array, include null elements */
    size_type m_size;
    /* The value for null elements, typically zero */
    value_type m_null;
    /* The indices and vlaues for non-null elements */
    std::vector< std::pair<size_t,T> > m_data;    

  };  


  template <typename T>
  SparseVector<T>::SparseVector(T null)
    :m_size(0),
     m_null(null){
  }
  
  template <typename T>
  SparseVector<T>::SparseVector(size_type n, T null)
    :m_size(n),
     m_null(null){
  }
 
  template <typename T>
  SparseVector<T>::SparseVector(const std::vector<T>& vect, T null)
    :m_size(vect.size()),
     m_null(null){
    fill_from_vect(vect);
  }

  template <typename T>
  SparseVector<T>::SparseVector(const SparseVector<T>& other)
    :m_size(other.size()),
     m_null(other.m_null),
     m_data(other.m_data){
  }

  template <typename T>
  SparseVector<T>&  SparseVector<T>::operator=(const SparseVector<T>& other) {
    if ( this != &other) {
      m_size = other.m_size;
      m_null = other.m_null;
      m_data = other.m_data;
    }
    return *this;
  }

  template <typename T>
  typename SparseVector<T>::const_reference SparseVector<T>::operator[](size_type n) const {
    const_iterator it = find(n);
    return it != end() ? it->second : m_null;
  }
  
  template <typename T>
  void SparseVector<T>::push_back(const T& val) {
    if ( test_value(val) ) {
      m_data.push_back(std::make_pair<size_t,T>(m_size,val));
    }
    m_size += 1;
  }
  
  template <typename T>  
  void SparseVector<T>::swap(SparseVector<T>& other) {
    if ( this == &other ) return;
    size_t tt = other.m_size;
    other.m_size = m_size;
    m_size = tt;
    T nn = other.m_null;
    other.m_null = m_null;
    m_null = nn;
    m_data.swap(other.m_data);
  }

  template <typename T>
  void SparseVector<T>::resize(size_type n) {
    if ( n >= m_size ) {
      m_size = n;
      return;
    }
    iterator itrDel = lower_bound(n);
    m_data.erase(itrDel,m_data.end());
  }

  template <typename T>
  typename SparseVector<T>::const_iterator SparseVector<T>::find(size_type k) const {
    const_iterator it = lower_bound(k);
    if ( it == end() ) return end();
    if ( it->first == k ) return it;
    return end();
  }

  template <typename T>
  typename SparseVector<T>::iterator SparseVector<T>::find(size_type k) {
    iterator it = lower_bound(k);
    if ( it == end() ) return end();
    if ( it->first == k ) return it;
    return end();    
  }

  template <typename T>
  typename SparseVector<T>::const_iterator SparseVector<T>::lower_bound(size_type k) const {
    return std::lower_bound(begin(),end(),k,key_compare_lower);
  }

  template <typename T>
  typename SparseVector<T>::iterator SparseVector<T>::lower_bound(size_type k) {
    return std::lower_bound(begin(),end(),k,key_compare_lower);
  }
  
  template <typename T>
  typename SparseVector<T>::const_iterator SparseVector<T>::upper_bound(size_type k) const {
    return std::upper_bound(begin(),end(),k,key_compare_upper);
  }

  template <typename T>
  typename SparseVector<T>::iterator SparseVector<T>::upper_bound(size_type k) {
    return std::upper_bound(begin(),end(),k,key_compare_upper);
  }
  
  template <typename T>
  void SparseVector<T>::fill_from_vect(const std::vector<T>& vect) {
    if ( m_size != vect.size() ) {
      throw std::runtime_error("SparseVector<T>::fill_from_vect, input size does not match");
    }
    m_data.clear();
    size_t idx(0);
    for ( typename std::vector<T>::const_iterator itr = vect.begin(); itr!= vect.end(); itr++, idx++ ) {
      if ( test_value(*itr) ) {
	m_data.push_back(std::make_pair<size_t,T>(idx,*itr));
      }
    }
  }
  
  template <typename T>
  void SparseVector<T>::fill_from_map(const std::map<size_t,T>& map) {
    m_data.clear();
    for ( typename std::map<size_t,T>::const_iterator itr = map.begin(); itr!= map.end(); itr++ ) {
      m_data.push_back(std::make_pair<size_t,T>(itr->first,itr->second));
    }    
  }
  
  template <typename T>
  void SparseVector<T>::fill_from_key_and_value(const std::vector<size_t>& keys,
						const std::vector<T>& values) {
    m_data.clear();
    if ( keys.size() != values.size() ) {
      throw std::runtime_error("SparseVector<T>::fill_from_key_and_value: key size does not equal value size");
    }
    std::vector<size_t>::const_iterator itrKey = keys.begin();
    typename std::vector<T>::const_iterator itrVal = values.begin();
    
    for ( typename std::vector<T>::const_iterator itrVal = values.begin(); itrVal != values.end(); itrKey++, itrVal++ ) {
      m_data.push_back(std::make_pair<size_t,T>(*itrKey,*itrVal));
    }     
  }
  
  template <typename T>
  void SparseVector<T>::fill_vect(std::vector<T>& vect) const {
    vect.resize(size());
    for ( const_iterator itr = begin(); itr != end(); itr++ ) {
      vect[itr->first] = itr->second;
    }
  }

  template <typename T>
  void SparseVector<T>::fill_vect(std::vector<T>& vect, size_t first, size_t last) const {
    vect.resize(last-first);
    for ( const_iterator itr = lower_bound(first); itr < upper_bound(last) ; itr++ ) {
      vect[itr->first - first] = itr->second;
    }
  }

  template <typename T>
  void SparseVector<T>::fill_map(std::map<size_t,T>& map) const {
    map.clear();
    for ( const_iterator itr = begin(); itr != end(); itr++ ) {
      map[itr->first] = itr->second;
    }    

  }
  
  template <typename T>
  void SparseVector<T>::fill_map(std::map<size_t,T>& map, size_t first, size_t last) const {
    map.clear();
    for ( const_iterator itr = lower_bound(first); itr < upper_bound(last) ; itr++ ) {
      map[itr->first - first] = itr->second;
    }    
  }
  
  template <typename T>
  void SparseVector<T>::fill_key_and_value(std::vector<size_t>& keys,
					   std::vector<T>& values) const {
    keys.resize(m_data.size());
    values.resize(m_data.size(),m_null);
    size_t idx(0);   
    for ( const_iterator itr = begin(); itr != end(); itr++, idx++ ) {
      keys[idx] = itr->first;
      values[idx] = itr->second;
    }
 
  }

  template <typename T>
  void SparseVector<T>::fill_key_and_value(std::vector<size_t>& keys,
					   std::vector<T>& values, 
					   size_t first, size_t last) const {
    const_iterator itr_begin = lower_bound(first);
    const_iterator itr_end = upper_bound(last);
    keys.resize(itr_end-itr_begin);
    values.resize(itr_end-itr_begin,m_null);
    size_t idx(0);   
    for ( const_iterator itr = itr_begin; itr != itr_end; itr++, idx++ ) {
      keys[idx-first] = itr->first;
      values[idx-first] = itr->second;
    }
  }
  
  template <typename T>
  void SparseVector<T>::add_to_vect(std::vector<T>& vect) const {
    for ( const_iterator itr = begin(); itr != end(); itr++ ) {
      vect[itr->first] += itr->second;
    }    
  }

  template <typename T>
  void SparseVector<T>::subtract_from_vect(std::vector<T>& vect) const {
    for ( const_iterator itr = begin(); itr != end(); itr++ ) {
      vect[itr->first] -= itr->second;
    }    
  }

} // namespace Likelihood

#endif // Likelihood_FitUtils_h
