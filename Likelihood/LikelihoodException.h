#ifndef LIKELIHOOD_EXCEPTION_H
#define LIKELIHOOD_EXCEPTION_H
#include <exception>
#include <string>

namespace Likelihood {
  class LikelihoodException: public std::exception {
  public:
    LikelihoodException(std::string errorString, int code=0) : 
      m_what(errorString), m_code(code) 
      {}
    virtual const char *what() const {return m_what.c_str();}
    virtual const int code() const {return m_code;}
  private:
    std::string m_what;
    int m_code;
  };
}
#endif //LIKELIHOOD_EXCEPTION_H
