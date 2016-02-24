#include "generate_baby.hpp"

#define ERROR(x) throw std::runtime_error(string("Error in file ")+__FILE__+" at line "+to_string(__LINE__)+" (in "+__func__+"): "+x);
#define DBG(x) std::cout << "In " << __FILE__ << " at line " << __LINE__ << " (in function " << __func__ << "): " << x << std::endl;

#include <cctype>

#include <string>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

int main(int argc, char *argv[]){
  set<string> files;
  for(int argi = 1; argi < argc; ++argi){
    files.insert(argv[argi]);
  }

  set<Variable> vars = GetVariables(files);
  WriteBaseHeader(vars, files);
  WriteBaseSource(vars);
  for(const auto &file: files){
    WriteSpecializedHeader(vars, file);
    WriteSpecializedSource(vars, file);
  }
}

SimpleVariable::SimpleVariable(const string &type, const string &name):
  type_(type),
  name_(name){
  }

Variable::Variable(const string &name):
  name_(name),
  type_map_(){
}

SimpleVariable Variable::GetEntry(const string &baby_type) const{
  return SimpleVariable(baby_type, Type());
}

const std::string & Variable::Name() const{
  return name_;
}

std::string & Variable::Name(){
  return name_;
}

string Variable::Type() const{
  if(ImplementInBase() || VirtualInBase()){
    for(const auto &var: type_map_){
      if(var.second != "") return var.second;
    }
  }
  return "";
}

string Variable::Type(const string &baby_type) const{
  if(HasEntry(baby_type)){
    return type_map_.at(baby_type);
  }else{
    return "";
  }
}

string Variable::DecoratedType() const{
  string type = Type();
  if(type.find("std::") != string::npos){
    type = type+"*";
  }
  return type;
}

string Variable::DecoratedType(const string &baby_type) const{
  string type = Type(baby_type);
  if(type.find("std::") != string::npos){
    type = type+"*";
  }
  return type;
}

bool Variable::HasEntry(const string &baby_type) const{
  return type_map_.find(baby_type) != type_map_.cend();
}

void Variable::SetEntry(const string &baby_type,
                        const string &type){
  type_map_[baby_type] = type;
}

bool Variable::ImplementInBase() const{
  set<string> type_set = GetTypeSet();
  return type_set.size()==1 && *(type_set.cbegin()) != "";
}

bool Variable::VirtualInBase() const{
  set<string> type_set = GetTypeSet();
  if(type_set.size() != 2) return false;
  auto iter = type_set.cbegin();
  string first_type = *iter;
  ++iter;
  string second_type = *iter;
  return (first_type != "" && second_type == "")
    || (first_type == "" && second_type != "");
}

bool Variable::ImplementIn(const std::string &baby_type) const{
  return VirtualInBase() && Type(baby_type)!="";
}

bool Variable::NotInBase() const{
  set<string> type_set = GetTypeSet();
  if(type_set.size() == 2){
    auto iter = type_set.cbegin();
    string first_type = *iter;
    ++iter;
    string second_type = *iter;
    return (first_type != second_type)
      && (first_type != "")
      && (second_type != "");
  }else{
    return type_set.size() >= 3;
  }
}

bool Variable::EverythingIn(const string &baby_type) const{
  return NotInBase() && Type(baby_type)!="";
}

bool Variable::operator<(const Variable &other) const{
  return name_ < other.name_;
}

set<string> Variable::GetTypeSet() const{
  set<string> types;
  for(const auto &type: type_map_){
    types.insert(type.second);
  }
  return types;
}

set<Variable> GetVariables(const set<string> &files){
  vector<Variable> vars;
  for(const auto &file: files){
    ifstream ifs("txt/variables/"+file);
    for(string line; std::getline(ifs, line); ){
      if(IsComment(line)) continue;
      SimpleVariable simple_var = GetVariable(line);
      Variable var(simple_var.name_);
      vector<Variable>::iterator this_var = vars.end();
      for(auto iter = vars.begin();
          iter != vars.end();
          ++iter){
        if(iter->Name() == var.Name()) this_var = iter;
      }
      if(this_var != vars.end()){
        this_var->SetEntry(file, simple_var.type_);
      }else{
        var.SetEntry(file, simple_var.type_);
        vars.push_back(var);
      }
    }
  }

  for(auto &var: vars){
    for(const auto &file: files){
      if(!var.HasEntry(file)) var.SetEntry(file, "");
    }
  }

  return {vars.cbegin(), vars.cend()};
}

bool IsComment(const string &line){
  if(line.size() <= 2) return true;
  for(const auto &letter: line){
    if(letter == ' ') continue;
    if(isalpha(letter) || letter == '_'){
      return false;
    }else{
      return true;
    }
  }
  return true;
}

SimpleVariable GetVariable(string line){
  RemoveExtraSpaces(line);
  auto loc = line.rfind(';');
  if(loc != string::npos) line = line.substr(0,loc);
  loc = line.rfind(' ');
  if(loc == string::npos){
    ERROR("Could not separate type and variable in "+line);
    return SimpleVariable("", line);
  }else{
    return SimpleVariable(line.substr(0,loc), line.substr(loc+1));
  }
}

void RemoveExtraSpaces(string &line){
  size_t loc = line.find_first_not_of(" ");
  if(loc != string::npos){
    line = line.substr(loc);
  }
  while((loc = line.find("  ", loc)) != string::npos){
    line.replace(loc, 2, " ");
    ++loc;
  }
  loc = line.find_last_not_of(" ");
  if(loc != string::npos){
    line = line.substr(0, loc+1);
  }
}

string ToUpper(string x){
  for(size_t i = 0; i < x.size(); ++i){
    x.at(i) = toupper(x.at(i));
  }
  return x;
}

string ToLower(string x){
  for(size_t i = 0; i < x.size(); ++i){
    x.at(i) = tolower(x.at(i));
  }
  return x;
}

void WriteBaseHeader(const set<Variable> &vars,
                     const set<string> &types){
  ofstream file("inc/baby.hpp");
  file << "#ifndef H_BABY\n";
  file << "#define H_BABY\n\n";

  file << "#include <vector>\n";
  file << "#include <set>\n";
  file << "#include <memory>\n";
  file << "#include <string>\n\n";

  file << "#include \"TChain.h\"\n\n";

  file << "class Baby{\n";
  file << "public:\n";
  file << "  explicit Baby(const std::set<std::string> &file_names);\n";
  file << "  Baby(Baby &&) = default;\n";
  file << "  Baby& operator=(Baby &&) = default;\n";
  file << "  virtual ~Baby() = default;\n\n";

  file << "  long GetEntries() const;\n";
  file << "  virtual void GetEntry(long entry);\n\n";

  for(const auto &var: vars){
    if(var.ImplementInBase()){
      file << "  "
           << var.DecoratedType() << " const & "
           << var.Name() << "() const;\n";
    }else if(var.VirtualInBase()){
      file << "  virtual "
           << var.DecoratedType() << " const & "
           << var.Name() << "() const = 0;\n";
    }
  }
  file << "\n";

  file << "protected:\n";
  file << "  virtual void Initialize();\n\n";

  file << "  std::unique_ptr<TChain> chain_;\n";
  file << "  long entry_;\n\n";

  file << "private:\n";
  file << "  Baby() = delete;\n";
  file << "  Baby(const Baby &) = delete;\n";
  file << "  Baby& operator=(const Baby &) = delete;\n\n";

  file << "  std::set<std::string> file_names_;\n";
  file << "  mutable long total_entries_;\n";
  file << "  mutable bool cached_total_entries_;\n\n";

  for(const auto &var: vars){
    if(!var.ImplementInBase()) continue;
    file << "  "
         << var.DecoratedType() << " "
         << var.Name() << "_;\n";
    file << "  TBranch *b_" << var.Name() << "_;\n";
    file << "  mutable bool c_" << var.Name() << "_;\n";
  }
  file << "};\n\n";

  for(const auto &type: types){
    file << "#include \"baby_" << type << ".hpp\"\n";
  }
  file << '\n';
  
  file << "#endif" << endl;
  file.close();
}

void WriteBaseSource(const set<Variable> &vars){
  ofstream file("src/baby.cpp");
  file << "#include \"baby.hpp\"\n\n";

  file << "#include <mutex>\n";
  file << "#include <stdexcept>\n\n";

  file << "#include \"utilities.hpp\"\n\n";

  file << "using namespace std;\n\n";

  file << "Baby::Baby(const set<string> &file_names):\n";
  file << "  chain_(nullptr),\n";
  file << "  file_names_(file_names),\n";
  file << "  total_entries_(0),\n";
  auto last_base = vars.cbegin();
  bool found_in_base = false;
  for(auto iter = vars.cbegin(); iter != vars.cend(); ++iter){
    if(iter->ImplementInBase()){
      last_base = iter;
      found_in_base = true;
    }
  }
  if(vars.size() == 0 || !found_in_base){
    file << "  cached_total_entries_(false){\n";
  }else{
    file << "  cached_total_entries_(false),\n";
    for(auto var = vars.cbegin(); var != last_base; ++var){
      if(!var->ImplementInBase()) continue;
      file << "  " << var->Name() << "_{},\n";
      file << "  b_" << var->Name() << "_(nullptr),\n";
      file << "  c_" << var->Name() << "_(false),\n";
    }
    file << "  " << last_base->Name() << "_{},\n";
    file << "  b_" << last_base->Name() << "_(nullptr),\n";
    file << "  c_" << last_base->Name() << "_(false){\n";
  }
  file << "  lock_guard<mutex> lock(Multithreading::root_mutex);\n";
  file << "  chain_ = unique_ptr<TChain>(new TChain(\"tree\"));\n";
  file << "  for(const auto &file: file_names){\n";
  file << "    chain_->Add((file+\"/tree\").c_str());\n";
  file << "  }\n";
  file << "}\n";

  file << "long Baby::GetEntries() const{\n";
  file << "  if(!cached_total_entries_){\n";
  file << "    cached_total_entries_ = true;\n";
  file << "    lock_guard<mutex> lock(Multithreading::root_mutex);\n";
  file << "    total_entries_ = chain_->GetEntries();\n";
  file << "  }\n";
  file << "  return total_entries_;\n";
  file << "}\n\n";

  file << "void Baby::GetEntry(long entry){\n";
  for(const auto &var: vars){
    if(!var.ImplementInBase()) continue;
    file << "  c_" << var.Name() << "_ = false;\n";
  }
  file << "  lock_guard<mutex> lock(Multithreading::root_mutex);\n";
  file << "  entry_ = chain_->LoadTree(entry);\n";
  file << "}\n\n";

  file << "void Baby::Initialize(){\n";
  file << "  lock_guard<mutex> lock(Multithreading::root_mutex);\n";
  file << "  chain_->SetMakeClass(1);\n";
  for(const auto &var: vars){
    if(!var.ImplementInBase()) continue;
    file << "  chain_->SetBranchAddress(\"" << var.Name() << "\", &" << var.Name() << "_, &b_" << var.Name() << "_);\n";
  }
  file << "}\n";

  for(const auto &var: vars){
    if(!var.ImplementInBase()) continue;
    file << var.DecoratedType() << " const & Baby::" << var.Name() << "() const{\n";
    file << "  if(!c_" << var.Name() << "_ && b_" << var.Name() << "_){\n";
    file << "    b_" << var.Name() << "_->GetEntry(entry_);\n";
    file << "    c_" << var.Name() << "_ = true;\n";
    file << "  }\n";
    file << "  return " << var.Name() << "_;\n";
    file << "}\n\n";
  }
  file << flush;
  file.close();
}

void WriteSpecializedHeader(const set<Variable> &vars, const string &type){
  ofstream file("inc/baby_"+type+".hpp");
  file << "#ifndef H_BABY_" << ToUpper(type) << "\n";
  file << "#define H_BABY_" << ToUpper(type) << "\n\n";

  file << "#include \"baby.hpp\"\n\n";

  file << "class Baby_" << type << ": virtual public Baby{\n";
  file << "public:\n";
  file << "  explicit Baby_" << type << "(const std::set<std::string> &file_names);\n";
  file << "  virtual ~Baby_" << type << "() = default;\n\n";

  file << "  virtual void GetEntry(long entry);\n\n";

  for(const auto &var: vars){
    if(var.VirtualInBase()){
      if(var.ImplementIn(type)){
        file << "  virtual " << var.DecoratedType() << " const & " << var.Name() << "() const;\n";
      }else{
        file << "  __attribute__((noreturn)) virtual " << var.DecoratedType()
             << " const & " << var.Name() << "() const;\n";
      }
    }else if(var.EverythingIn(type)){
      file << "  " << var.DecoratedType(type) << " const & " << var.Name() << "() const;\n";
    }
  }
  file << "\n";

  file << "private:\n";
  file << "  Baby_" << type << "() = delete;\n";
  file << "  Baby_" << type << "(const Baby_" << type << " &) = delete;\n";
  file << "  Baby_" << type << "& operator=(const Baby_" << type << " &) = delete;\n";
  file << "  Baby_" << type << "(Baby_" << type << " &&) = delete;\n";
  file << "  Baby_" << type << "& operator=(Baby_" << type << " &&) = delete;\n";

  file << "  virtual void Initialize();\n\n";

  for(const auto &var: vars){
    if(var.ImplementIn(type) || var.EverythingIn(type)){
      file << "  " << var.DecoratedType(type) << " " << var.Name() << "_;\n";
      file << "  TBranch *b_" << var.Name() << "_;\n";
      file << "  mutable bool c_" << var.Name() << "_;\n";
    }
  }
  file << "};\n\n";

  file << "#endif" << endl;
  file.close();
}

void WriteSpecializedSource(const set<Variable> &vars, const string &type){
  ofstream file("src/baby_"+type+".cpp");
  file << "#include \"baby_" << type << ".hpp\"\n\n";

  file << "#include \"utilities.hpp\"\n\n";
  
  file << "using namespace std;\n\n";

  file << "Baby_" << type << "::Baby_" << type << "(const set<string> &file_names):\n";
  int implemented_here = 0;
  for(const auto & var: vars){
    if(var.ImplementIn(type) || var.EverythingIn(type)) ++implemented_here;
  }
  if(implemented_here == 0){
    file << "  Baby(file_names){\n";
  }else{
    file << "  Baby(file_names),\n";
    set<Variable>::const_iterator last = vars.cend();
    for(auto var = vars.cbegin(); var != vars.cend(); ++var){
      if(var->ImplementIn(type) || var->EverythingIn(type)){
        last = var;
      }
    }
    for(auto var = vars.cbegin(); var != vars.cend(); ++var){
      if(var->ImplementIn(type) || var->EverythingIn(type)){
        file << "  " << var->Name() << "_{},\n";
        file << "  b_" << var->Name() << "_(nullptr),\n";
        if(var != last){
          file << "  c_" << var->Name() << "_(false),\n";
        }else{
          file << "  c_" << var->Name() << "_(false){\n";
        }
      }
    }
  }
  file << "  Initialize();\n";
  file << "}\n\n";

  file << "void Baby_" << type << "::GetEntry(long entry){\n";
  for(const auto &var: vars){
    if(var.ImplementIn(type) || var.EverythingIn(type)){
      file << "  c_" << var.Name() << "_ = false;\n";
    }
  }
  file << "  Baby::GetEntry(entry);\n";
  file << "}\n\n";

  file << "void Baby_" << type << "::Initialize(){\n";
  file << "  Baby::Initialize();\n";
  if(vars.size() > 0){
    file << "  lock_guard<mutex> lock(Multithreading::root_mutex);\n";
  }
  for(const auto &var: vars){
    if(var.ImplementIn(type) || var.EverythingIn(type)){
      file << "  chain_->SetBranchAddress(\"" << var.Name() << "\", &"
           << var.Name() << "_, &b_" << var.Name() << "_);\n";
    }
  }
  file << "}\n";

  for(const auto &var: vars){
    if(var.ImplementIn(type) || var.EverythingIn(type)){
      file << var.DecoratedType(type) << " const & Baby_" << type << "::" << var.Name() << "() const{\n";
      file << "  if(!c_" << var.Name() << "_ && b_" << var.Name() << "_){\n";
      file << "    b_" << var.Name() << "_->GetEntry(entry_);\n";
      file << "    c_" << var.Name() << "_ = true;\n";
      file << "  }\n";
      file << "  return " << var.Name() << "_;\n";
      file << "}\n\n";
    }else if(var.VirtualInBase()){
      file << var.DecoratedType() << " const & Baby_" << type << "::" << var.Name() << "()  const{\n";
      file << "  ERROR(\"" << var.DecoratedType() << ' ' << var.Name()
           << " not available in babies of type " << type << ".\")\n";
      file << "}\n\n";
    }
  }

  file << flush;
  file.close();
}
