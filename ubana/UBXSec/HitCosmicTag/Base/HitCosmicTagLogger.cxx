#ifndef __HITCOSMICTAGLOGGER_CXX__
#define __HITCOSMICTAGLOGGER_CXX__

#include "HitCosmicTagLogger.h"

namespace cosmictag {

  std::map<std::string,logger> *logger::_logger_m = nullptr;

  msg::Level_t logger::_level_default = msg::kNORMAL;
  
  std::ostream& logger::send(const msg::Level_t level) const
  {
    (*_ostrm)  << msg::kStringPrefix[level].c_str()
	       << " ";
    return (*_ostrm);
  }
  
  std::ostream& logger::send(const msg::Level_t level,
			     const std::string& function ) const
  {
    auto& strm(send(level));
    strm << "<" << _name << "::" << function.c_str() << "> ";
    return strm;
  }
  
  std::ostream& logger::send(const msg::Level_t level,
			     const std::string& function,
			     const unsigned int line_num ) const
  {
    auto& strm(send(level));
    strm << "<" << _name << "::" << function.c_str() << "::L" << line_num << "> ";
    return strm;
  }
  
  std::ostream& logger::send(const msg::Level_t level,
			     const std::string& function,
			     const unsigned int line_num,
			     const std::string& file_name) const
  {
    auto& strm(send(level,function));
    // FIXME temporary operation to fetch file name from full path
    strm << file_name.substr(file_name.rfind("/")+1,file_name.size()).c_str() << "::L" << line_num << " ";
    return strm;
  }
}
#endif
