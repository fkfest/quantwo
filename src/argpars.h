#ifndef ARGPARS_H
#define ARGPARS_H

#include <string>
#include <cstring>
#include <vector>

/*!
 * Argument list parser
 */
class ArgPars {
public:
  ArgPars() : _argc(0), _curarg(0), _nextlet(0) {}
  ArgPars(int argc, char **argv) : _argc(argc), _argv(argv), 
      _options(argc,false), _curarg(0), _nextlet(0) {}
  // next option
  bool nextoption(std::string & opt);
  // next argument after the current option. Using shift one can get other arguments
  bool optarg(std::string& arg, int shift = 0);
  // mark the next argument as a part of options
  void markasoption( int shift = 0 ) {
    unsigned int pos = _curarg + 1 + shift;
    if ( pos < _options.size() ) _options[pos] = true;
  }
  // unhandled arguments left
  bool nextremaining(std::string& arg);
private:
  // arguments
  int _argc;
  char **_argv;
  // status: is it an option?
  std::vector<bool> _options;
  // current argument
  int _curarg;
  // current letter in multi-option
  unsigned int _nextlet;
};

bool ArgPars::nextoption(std::string & opt) {
  if ( _curarg >= _argc ) {
    opt.clear();
    return false;
  }
  if ( _nextlet > 0 && _nextlet < strlen(_argv[_curarg])) {
    // next letter
    opt = std::string(1,_argv[_curarg][_nextlet]);
    ++_nextlet;
    return true;
  }
  
  // next argument
  ++_curarg;
  _nextlet = 0;
  if ( _curarg == _argc ) {
    opt.clear();
    return false;
  }
  
  if ( !_options[_curarg] && _argv[_curarg][0] == '-' ) {
    // get options
    _options[_curarg] = true;
    if ( strlen(_argv[_curarg]) > 2 && _argv[_curarg][1] == '-' ){
      // long option "--word"
      opt = std::string(&_argv[_curarg][1]);
      return true;
    } else {
      // short options "-wo"
      ++_nextlet;
    }    
  }
  return nextoption(opt);
}
bool ArgPars::optarg(std::string& arg, int shift) {
  int pos = _curarg + 1 + shift;
  if ( pos < _argc ) {
    arg = std::string(_argv[pos]);
    return true;
  } else {
    arg.clear();
    return false;
  }
}

bool ArgPars::nextremaining(std::string& arg) {
  for ( int iarg = 1; iarg < _argc; ++iarg ) {
    if ( !_options[iarg] ) {
      _options[iarg] = true;
      arg = std::string(_argv[iarg]);
      return true;
    }
  }
  arg.clear();
  return false;
}

#endif
