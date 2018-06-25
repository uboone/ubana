#ifndef HITCOSMICTAG_EREXCEPTION_CXX
#define HITCOSMICTAG_EREXCEPTION_CXX

#include "HitCosmicTagException.h"

namespace cosmictag {

  HitCosmicTagException::HitCosmicTagException(const std::string& msg)
    : std::exception()
  {
      _msg = "\033[93m EXCEPTION \033[00m\033[95m";
      _msg += msg;
      _msg += "\033[00m\n";
  }

  const char* HitCosmicTagException::what() const throw()
  { return _msg.c_str(); }

}
#endif
