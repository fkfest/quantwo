#include "utilities.h"
inline void error(std::string what, std::string where)
{
  std::cerr << "  ERROR: " << what << std::endl;
  if (where.size() > 0)
    std::cerr << "  in " << where << std::endl;
  exit(1);
}
inline void say(std::string what, std::string where)
{
  _xout0(what << std::endl);
  if (where.size() > 0)
    _xout0("  in " << where << std::endl);
}
inline std::string exepath(){
  std::string path;
  char buff[1024];
#ifdef __linux__
// Linux
  ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
  if (len != -1) {
    buff[len] = '\0';
    path = std::string(buff);
#elif defined TARGET_OS_MAC
// Mac
  int  bufsize = sizeof(buff);
  if (_NSGetExecutablePath(buff, &bufsize) == 0){
    path = std::string(buff);
#else
#error "unknown platform"
  if (false){
#endif
  } else {
     /* handle error condition */
     error("Problem by detection of exe-path");
  }
  return DirName(path);
}

template <class T>
inline
bool str2num(T& t, const std::string& s,
             std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

template <class T>
inline
std::string num2str(const T& t,
             std::ios_base& (*f)(std::ios_base&))
{
  std::ostringstream oss;
  oss << f << t;
  return oss.str();
}


