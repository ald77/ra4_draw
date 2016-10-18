#include "core/config_parser.hpp"

#include <fstream>

#include <sys/stat.h>
#include <unistd.h>

#include "core/utilities.hpp"

using namespace std;

ConfigParser::ConfigParser():
  options_(){
}

ConfigParser & ConfigParser::Load(const string &file_path,
                                  const string &option_set){
  ifstream file(file_path.c_str());
  string line, opt_set, opt_name, opt_value;
  while(getline(file, line)){
    line = Strip(line);

    //Comment line
    if(line.size()==0 || line.front() == '#') continue;
    
    //New option set
    if(line.front() == '[' && line.back() == ']'){
      opt_set = Strip(line.substr(1, line.size()-2));
      continue;
    }

    //Option
    if(opt_set != option_set) continue;
    auto eq_pos = line.find('=');
    opt_name = Strip(line.substr(0, eq_pos));
    opt_value = Strip(line.substr(eq_pos+1));
    options_[opt_name] = opt_value;
  }
  return *this;
}

void ConfigParser::Save(const string &file_path,
                        const string &option_set) const{
  bool file_exists = FileExists(file_path);
  if(file_exists){
    string out_path = MakeTemp(file_path);
    ifstream in_file(file_path.c_str());
    ofstream out_file(out_path.c_str());
    string line, opt_set, opt_name, opt_value;
    while(getline(in_file, line)){
      string clean_line = Strip(line);

      //Blank line
      if(clean_line.size()==0) continue;
      
      //Comment line
      if(clean_line.front() == '#'){
        out_file << line << endl;
        continue;
      }

      //New option set
      if(clean_line.front() == '[' && clean_line.back() == ']'){
        opt_set = Strip(clean_line.substr(1, line.size()-2));
        if(opt_set != option_set) out_file << endl << line << endl;
        continue;
      }

      if(opt_set != option_set){
        out_file << line << endl;
        continue;
      }
    }
    WriteToStream(out_file, option_set);
    out_file.close();
    execute("mv "+out_path+" "+file_path);
  }else{
    ofstream out_file(file_path.c_str());
    WriteToStream(out_file, option_set);
  }
}

bool ConfigParser::HaveOpt(const std::string &option) const{
  return options_.find(option) != options_.cend();
}

const map<string, string> & ConfigParser::Options() const{
  return options_;
}

void ConfigParser::WriteToStream(ofstream &f, const string &option_set) const{
  f<< endl;
  f << '[' << option_set << ']' << endl;
  for(const auto &opt: options_){
    f << "  " << opt.first << " = " << opt.second << endl;
  }
}

ostream & operator<<(ostream &stream, const ConfigParser &cp){
  const auto & options = cp.Options();
  string body = accumulate(options.cbegin(), options.cend(), string(),
                           [](const string &result, const pair<string, string> &opt){
                             string to_add = opt.first+": "+opt.second;
                             return (result.empty() ? string() : result+", ")+to_add;
                           });
  stream << '{' << body << '}';
  return stream;
}

template<>
std::string ConfigParser::GetOpt<std::string>(const std::string &option) const{
  return options_.at(option);
}
