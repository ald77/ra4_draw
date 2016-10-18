#ifndef H_CONFIG_PARSER
#define H_CONFIG_PARSER

#include <string>
#include <map>
#include <sstream>
#include <limits>

class ConfigParser{
public:
  ConfigParser();

  ConfigParser & Load(const std::string &file_path,
                      const std::string &option_set);
  void Save(const std::string &file_path,
            const std::string &option_set) const;

  bool HaveOpt(const std::string &option) const;

  template<typename T = std::string>
  T GetOpt(const std::string &option) const;

  template<typename T = std::string>
  ConfigParser & SetOpt(const std::string &option, const T& value);

  const std::map<std::string, std::string> & Options() const;

private:
  std::map<std::string, std::string> options_;

  void WriteToStream(std::ofstream &f, const std::string &option_set) const;
};

std::ostream & operator<<(std::ostream &stream, const ConfigParser &cp);


template<typename T>
T ConfigParser::GetOpt(const std::string &option) const{
  std::istringstream iss(options_.at(option));
  T val;
  iss >> val;
  return val;
}

template<typename T>
ConfigParser & ConfigParser::SetOpt(const std::string &option, const T& value){
  std::ostringstream oss;
  oss.precision(std::numeric_limits<T>::digits10);
  oss << value;
  options_[option] = oss.str();
  return *this;
}

#endif
