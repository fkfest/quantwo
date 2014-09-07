#include "utilities.h"
#ifdef _WIN32
#include <windows.h>
#endif
void error(std::string what, std::string where)
{
  std::cerr << "  ERROR: " << what << std::endl;
  if (where.size() > 0)
    std::cerr << "  in " << where << std::endl;
  exit(1);
}
void say(std::string what, std::string where)
{
  _xout0(what << std::endl);
  if (where.size() > 0)
    _xout0("  in " << where << std::endl);
}
std::string exepath(){
  std::string path;
  char buff[1024];
#if defined __linux__ || defined __CYGWIN__
// Linux
  ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
  if (len != -1) {
    buff[len] = '\0';
    path = std::string(buff);
#elif defined __MACH__
// Mac
  uint  bufsize = sizeof(buff);
  if (_NSGetExecutablePath(buff, &bufsize) == 0){
    path = std::string(buff);
#elif defined _WIN32
// windows
  HMODULE hModule = GetModuleHandleW(NULL);
  if (hModule != NULL){
    GetModuleFileName(hModule,buff, sizeof(buff));
    path = std::string(buff);
    std::size_t bslash = path.rfind("\\");
    if (bslash != std::string::npos )
      path = path.substr(0,bslash);
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

std::size_t curlyfind(const std::string& str, const std::string& what, std::size_t ipos)
{
  std::size_t ipos1, ires(ipos), icur(0);
  int level = 0;
  bool found(false);
  for ( ipos1 = ipos; ipos1 < str.size() && !found; ++ipos1 ){
    if ( level == 0 ){
      // search at the current level only!
      if ( str[ipos1] == what[icur] ){
        if ( icur == 0 ) ires = ipos1;
        ++icur;
        if ( icur == what.size() ) found = true;
      } else if ( icur > 0 ){
        icur = 0;
        ipos1 = ires;
        continue;
      }
    }
    if ( str[ipos1] == '{' ) ++level;
    if ( str[ipos1] == '}' ) --level;
  }
  if (!found)
    ires = std::string::npos;
  return ires;
}
