#ifndef ARGPARS_H
#define ARGPARS_H

#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>

//! Main functionality in class ArgPars
namespace ArgParser {
/*!
 * Option (can be used to generate help print)
 */
class ArgOpt {
public:
  ArgOpt(const std::string& desc, const std::vector<std::string>& opt ) : _opt(opt), _desc(desc) {}
  ArgOpt(const std::string& desc, const std::string& opt1 ) : _desc(desc) { _opt.push_back(opt1);}
  ArgOpt(const std::string& desc, const std::string& opt1, const std::string& opt2 ) : _desc(desc)
  { _opt.push_back(opt1);_opt.push_back(opt2);}
  ArgOpt(const std::string& desc, const std::string& opt1, const std::string& opt2, const std::string& opt3 ) : _desc(desc)
  { _opt.push_back(opt1);_opt.push_back(opt2);_opt.push_back(opt3);}
  ArgOpt(const std::string& desc, const std::string& opt1, const std::string& opt2, const std::string& opt3, const std::string& opt4  ) : _desc(desc)
  { _opt.push_back(opt1);_opt.push_back(opt2);_opt.push_back(opt3);_opt.push_back(opt4);}
  ArgOpt(const std::string& desc, const std::string& opt1, const std::string& opt2, const std::string& opt3, const std::string& opt4, const std::string& opt5 ) : _desc(desc)
  { _opt.push_back(opt1);_opt.push_back(opt2);_opt.push_back(opt3);_opt.push_back(opt4);_opt.push_back(opt5);}
  // check if opt is equal to one of the option-keywords
  bool equal(const std::string& opt) const {
    for ( unsigned i = 0; i < _opt.size(); ++i )
      if ( opt == _opt[i] )
        return true;
    return false;
  }
  void print(std::ostream & o, int optwidth) const {
    std::stringstream str;
    for ( std::vector<std::string>::const_iterator it=_opt.begin(); it != _opt.end(); ++it ){
      str << "-" << *it << " ";
    }
    o << std::left << std::setw(optwidth) << str.str() << std::setw(0);
    o << _desc;
  }
private:
  // option
  std::vector< std::string > _opt;
  // description
  std::string _desc;
};

/*!
 * All options
 */
class ArgOpts : public std::vector<ArgOpt> {
public:
  ArgOpts() : std::vector<ArgOpt>() {}
  // add new option
  ArgOpts::iterator add(const ArgOpt& aopt) { return insert(end(),aopt); }
  // add and check option
  bool check(const std::string& opt, const ArgOpt& aopt){ return add(aopt)->equal(opt);}
  // print help
  // with usage information and description. The width reserved for options can be changed using optwidth
  void printhelp(std::ostream & o, const std::string& usage = "", const std::string& desc = "", int optwidth = 16) const {
    if (!usage.empty()) o << "Usage: " << usage << std::endl << std::endl;
    if (!desc.empty()) o << desc << std::endl;
    if (!empty()) {
      o << "Options:" << std::endl;
      for ( const_iterator it = begin(); it != end(); ++it ) {
        it->print(o,optwidth);
        o << std::endl;
      }
    }
  }
};

/*!
 * Argument list parser
 */
class ArgPars {
public:
  ArgPars() : _argc(0), _curarg(0), _nextlet(0) {}
  ArgPars(int argc, char **argv) : _argc(argc), _argv(argv),
      _options(argc,false), _curarg(0), _nextlet(0) {}
  // next option, if clear is true: clear stored option info
  bool nextoption(std::string & opt, bool clear = true);
  // next option, if clear is true: clear stored option info
  bool nextoption(bool clear = true);
  // add and check option
  bool check(const ArgOpt& aopt) { return _argopts.check(_opt,aopt); }
  // get current option
  std::string get_current_option() const { return _opt; }
  // next argument after the current option. Using shift one can get other arguments
  // if tolerate_minus = true: arguments starting with "-" are allowed
  bool optarg(std::string& arg, int shift = 0, bool tolerate_minus = false);
  // mark the next argument as a part of options
  void markasoption( int shift = 0 ) {
    unsigned int pos = _curarg + 1 + shift;
    if ( pos < _options.size() ) _options[pos] = true;
  }
  // unhandled arguments left
  bool nextremaining(std::string& arg);
  // print help
  // with usage information and description. The width reserved for options can be changed using optwidth
  void printhelp(std::ostream & o, const std::string& usage = "", const std::string& desc = "", int optwidth = 16) const {
          _argopts.printhelp(o,usage,desc,optwidth);
        }
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
  // current option
  std::string _opt;
  // list of possible options
  ArgOpts _argopts;
};


bool ArgPars::nextoption(bool clear) {
  if (clear) _argopts.clear();
  if ( _curarg >= _argc ) {
    _opt.clear();
    return false;
  }
  if ( _nextlet > 0 && _nextlet < strlen(_argv[_curarg])) {
    // next letter
    _opt = std::string(1,_argv[_curarg][_nextlet]);
    ++_nextlet;
    return true;
  }

  // next argument
  ++_curarg;
  _nextlet = 0;
  if ( _curarg == _argc ) {
    _opt.clear();
    return false;
  }

  if ( !_options[_curarg] && _argv[_curarg][0] == '-' ) {
    // get options
    _options[_curarg] = true;
    if ( strlen(_argv[_curarg]) > 2 && _argv[_curarg][1] == '-' ){
      // long option "--word"
      _opt = std::string(&_argv[_curarg][1]);
      return true;
    } else {
      // short options "-wo"
      ++_nextlet;
    }
  }
  return nextoption();
}
bool ArgPars::nextoption(std::string & opt, bool clear) {
  bool ret = nextoption(clear);
  opt = _opt;
  return ret;
}
bool ArgPars::optarg(std::string& arg, int shift, bool tolerate_minus) {
  int pos = _curarg + 1 + shift;
  if ( pos >= _argc || ( !tolerate_minus && _argv[pos][0] == '-' ) ) {
    arg.clear();
    return false;
  } else {
    arg = std::string(_argv[pos]);
    return true;
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

}
#endif
