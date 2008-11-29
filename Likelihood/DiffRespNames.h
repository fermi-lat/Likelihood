/**
 * @file DiffRespNames.h
 * @brief Class to keep track of diffuse response column names.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_DiffRespNames_h
#define Likelihood_DiffRespNames_h

namespace Likelihood {

/**
 * @class DiffRespNameError
 * @brief Exception class for handling error from DiffRespNames
 */

class DiffRespNameError : public std::exception {

public :

   DiffRespNameError() {}

   DiffRespNameError(std::string errorString) : m_what(errorString) {}

   virtual ~DiffRespNameError() throw() {}

   virtual const char * what() const throw() {return m_what.c_str();}

protected:

   std::string m_what;

};

/**
 * @class DiffRespNames
 *
 */

class DiffRespNames {

public:

   DiffRespNames(const std::string & rootName="DIFRSP");

   size_t size() const;

   const std::vector<std::string> & colnames() const;

   const std::string & operator[](size_t indx) const;

   const std::string & operator[](const std::string & key) const;

   void addColumn(const std::string & diffRspName);

   std::string key(const std::string & colname) const;

private:

   std::string m_rootName;

   std::vector<std::string> m_colnames;

};

} // namespace Likelihood

#endif // Likelihood_DiffRespNames_h

