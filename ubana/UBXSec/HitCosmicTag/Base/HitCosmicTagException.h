/**
 * \file HitCosmicTagException.h
 *
 * \ingroup HitCosmicTag
 * 
 * \brief Class def header for exception classes in HitCosmicTag package
 *
 * @author Marco Del Tutto
 */

/** \addtogroup HitCosmicTag

    @{*/
#ifndef HITCOSMICTAG_EXCEPTION_H
#define HITCOSMICTAG_EXCEPTION_H

#include <iostream>
#include <exception>

namespace cosmictag {
  /**
     \class HitCosmicTagException
     Generic (base) exception class for HitCosmicTag package
  */
  class HitCosmicTagException : public std::exception{

  public:
    /// Default ctor w/ error message input (optional)
    HitCosmicTagException(const std::string& msg="");

    virtual ~HitCosmicTagException() throw(){};

    virtual const char* what() const throw();

  private:
    /// Error message
    std::string _msg;
  };

}
#endif
/** @} */ // end of doxygen group 

