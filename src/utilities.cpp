inline void error(std::string what, std::string where)
{
  std::cerr << "  ERROR: " << what << std::endl;
  if (where!="")
    std::cerr << "  in " << where << std::endl;
  exit(1);
}
inline void say(std::string what, std::string where)
{
  std::cout << what << std::endl;
  if (where!="")
    std::cout << "  in " << where << std::endl;
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

