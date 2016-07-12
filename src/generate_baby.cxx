/*! \class SimpleVariable

  \brief Pairs a type and name of a variable accessible via Baby

  Used only to generate Baby classes, but not in subsequent analysis code.
*/

/*! \class Variable

  \brief A variable to be accessible in Baby classes.

  Stores information about which concrete classes in which to implement the
  variable.

  Used only to generate Baby classes, but not in subsequent analysis code.
*/
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

/*!\brief Standard constructor

  \param[in] type Type of variable (e.g., int, float, etc.)

  \param[in] name Name of variable (e.g., ht, met, etc.)
*/
SimpleVariable::SimpleVariable(const string &type, const string &name):
  type_(type),
  name_(name){
  }

/*!\brief Standard constructor

  \param[in] name Name of variable (e.g., ht, met, etc.)
*/
Variable::Variable(const string &name):
  name_(name),
  type_map_(){
}

/*!\brief Get variable name

  \return Variable name
*/
const std::string & Variable::Name() const{
  return name_;
}

/*!\brief Get variable name

  \return Variable name
*/
std::string & Variable::Name(){
  return name_;
}

/*!\brief Get one type for variable across all Baby classes

  \return If type is the same in all derived Baby classes, return that
  type. Otherwise, returns empty string.

  \see Variable::Type()
*/
string Variable::Type() const{
  if(ImplementInBase() || VirtualInBase()){
    for(const auto &var: type_map_){
      if(var.second != "") return var.second;
    }
  }
  return "";
}

/*!\brief Get type of variable for a derived Baby class

  \param[in] baby_type Name of derived Baby class (e.g., basic, full) for which
  to get type

  \return Type of variable as used in the derived Baby class

  \see Variable::DecoratedType(const std::string &baby_type) const
*/
string Variable::Type(const string &baby_type) const{
  if(HasEntry(baby_type)){
    return type_map_.at(baby_type);
  }else{
    return "";
  }
}

/*!\brief Get type with "std::" and "*" if needed

  \return If type is the same in all Baby classes, returns that type with
  "std::" and "*" added if needed. Otherwise, returns empty string.

  \see Variable::Type()
*/
string Variable::DecoratedType() const{
  string type = Type();
  if(type.find("std::") != string::npos){
    type = type+"*";
  }
  return type;
}

/*!\brief Get type with "std::" and "*" if needed

  \param[in] baby_type Name of derived Baby class (e.g., basic, full) for which
  to get type

  \return Type of variable as used in the derived Baby class, with "std::" and
  "*" if needed

  \see Variable::Type(const std::string &baby_type) const
*/
string Variable::DecoratedType(const string &baby_type) const{
  string type = Type(baby_type);
  if(type.find("std::") != string::npos){
    type = type+"*";
  }
  return type;
}

/*!\brief Check if a derived Baby class has this variable

  \param[in] baby_type Name of derived Baby class (e.g., basic, full) to check

  \return True if and only if the derived Baby type contains this variable
*/
bool Variable::HasEntry(const string &baby_type) const{
  return type_map_.find(baby_type) != type_map_.cend();
}

/*!\brief Set the type of this variable for use in a derived Baby class

  \param[in] baby_type Specifies the derived Baby class (e.g., basic, full) for
  which to set variable type

  \param[in] type Variable type to set (e.g., int, float)
*/
void Variable::SetEntry(const string &baby_type,
                        const string &type){
  type_map_[baby_type] = type;
}

/*!\brief Check if variable can be implemented directly in Baby

  Implementation in Baby requires that the variable be present in all derived
  Baby classes with the same type

  \return True if variable accessor can be implemented in Baby
*/
bool Variable::ImplementInBase() const{
  set<string> type_set = GetTypeSet();
  return type_set.size()==1 && *(type_set.cbegin()) != "";
}

/*!\brief Check if variable should have pure virtual method in Baby

  Pure virtual accessor is needed if and only if the variable has exactly one
  type across all derived Baby classes, but does is only present in a subset of
  the classes.

  \return True if variable should have virtual accessor in Baby
*/
bool Variable::VirtualInBase() const{
  set<string> type_set = GetTypeSet();

  //type_set should contain one real type and an empty type indicating absence
  //from some Baby classes
  if(type_set.size() != 2) return false;
  auto iter = type_set.cbegin();
  string first_type = *iter;
  ++iter;
  string second_type = *iter;
  return (first_type != "" && second_type == "")
    || (first_type == "" && second_type != "");
}

/*!\brief Check if variable needs accessor implemented in a derived Baby class

  \param[in] baby_type Type of derived Baby to check

  \return True if variable is virtual in Baby, and implemented in the specified
  derived class
*/
bool Variable::ImplementIn(const std::string &baby_type) const{
  return VirtualInBase() && Type(baby_type)!="";
}

/*!\brief Check if variable is completely absent in base Baby

  This occus if the variable has different types in different derived Baby
  classes

  \return True if variable is only present in derived Baby classes and not in
  Baby
*/
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

/*!\brief Check if variable needs declaration and definition in a derived class

  \param[in] baby_type Specifies derived Baby class to check

  \return True if variable exists in specified derived class, but does not have
  virtual function in base Baby
*/
bool Variable::EverythingIn(const string &baby_type) const{
  return NotInBase() && Type(baby_type)!="";
}

/*!\brief Comparison operator allows storing Variable in set

  Sorts alphabetically by name

  \param[in] other Variable to compare to

  \return True if *this is alphabetically before other
*/
bool Variable::operator<(const Variable &other) const{
  return name_ < other.name_;
}

/*!\brief Get all types (int, float, etc.) used across derived Baby classes

  May include an empty string if variable is absent from some derived Baby
  classes

  \return List of types used in all derived Baby classes
*/
set<string> Variable::GetTypeSet() const{
  set<string> types;
  for(const auto &type: type_map_){
    types.insert(type.second);
  }
  return types;
}

/*!/brief Reads files to get variable names and associated types

  \param[in] files Set of files to be read. Typically everything in
  txt/variables

  \return All variables with associated types used for each derived Baby class
*/
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

/*!\brief Check if line in variable list file is a comment (blank)

  Not very robust. Currently, a line is a "comment" if is is less than 3
  characters long or contains neither letters nor underscores. Should add
  support for "#" starting a comment line and more robust checking.

  \return True if line meets above definition of a comment
*/
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

/*!\brief Extracts variable type and name from line from variable text files

  \param[in] line The line to be parsed

  \return SimpleVariable with name and type set
*/
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

/*!\brief Removes leading, trailing, and double spaced from a text line

  \param[in,out] line Line of text from which spaces are to be removed
*/
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

/*!\brief Replaces all lower case letters with upper case equivalent

  \param[in] x String to capitalize

  \return Copy of x with all lower case letters replaced with upper case
  equivalent
*/
string ToUpper(string x){
  for(size_t i = 0; i < x.size(); ++i){
    x.at(i) = toupper(x.at(i));
  }
  return x;
}

/*!\brief Replaces all upper case letters with lower case equivalent

  \param[in] x String to convert to all lower case

  \return Copy of x with all upper case letters replaced with lower case
  equivalent
*/
string ToLower(string x){
  for(size_t i = 0; i < x.size(); ++i){
    x.at(i) = tolower(x.at(i));
  }
  return x;
}

/*!\brief Writes inc/baby.hpp

  \param[in] vars All variables for all Baby classes, with type information

  \param[in] types Names of derived Baby classes (basic, full, etc.)
*/
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

  file << "class NamedFunc;\n\n";

  file << "class Baby{\n";
  file << "public:\n";
  file << "  explicit Baby(const std::set<std::string> &file_names);\n";
  file << "  Baby(Baby &&) = default;\n";
  file << "  Baby& operator=(Baby &&) = default;\n";
  file << "  virtual ~Baby() = default;\n\n";

  file << "  long GetEntries() const;\n";
  file << "  virtual void GetEntry(long entry);\n\n";

  file << "  const std::set<std::string> & FileNames() const;\n\n";

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

  file << "  const std::unique_ptr<TChain> & GetTree() const;\n\n";

  file << "  static NamedFunc GetFunction(const std::string &var_name);\n\n";

  file << "protected:\n";
  file << "  virtual void Initialize();\n\n";

  file << "  std::unique_ptr<TChain> chain_;//!<Chain to load variables from\n";
  file << "  long entry_;//!<Current entry\n\n";

  file << "private:\n";
  file << "  Baby() = delete;\n";
  file << "  Baby(const Baby &) = delete;\n";
  file << "  Baby& operator=(const Baby &) = delete;\n\n";

  file << "  std::set<std::string> file_names_;//!<Files loaded into TChain\n";
  file << "  mutable long total_entries_;//!<Cached number of events in TChain\n";
  file << "  mutable bool cached_total_entries_;//!<Flag if cached event count up to date\n\n";

  for(const auto &var: vars){
    if(!var.ImplementInBase()) continue;
    file << "  "
         << var.DecoratedType() << " "
         << var.Name() << "_;//!<Cached value of " << var.Name() << '\n';
    file << "  TBranch *b_" << var.Name() << "_;//!<Branch from which "
         << var.Name() << " is read\n";
    file << "  mutable bool c_" << var.Name() << "_;//!<Flag if cached "
         << var.Name() << " up to date\n";
  }
  file << "};\n\n";

  for(const auto &type: types){
    file << "#include \"baby_" << type << ".hpp\"\n";
  }
  file << '\n';

  file << "#endif" << endl;
  file.close();
}

/*!\brief Writes src/baby.cpp

  \param[in] vars All variables for all Baby classes, with type information
*/
void WriteBaseSource(const set<Variable> &vars){
  ofstream file("src/baby.cpp");
  file << "/*! \\class Baby\n\n";

  file << "  \\brief Abstract base class for access to ntuple variables\n\n";

  file << "  Loads variables on demand and caches for fast repeated use within an event.\n\n";

  file << "  A derived class is used for each known ntuple format. Variables and functions\n";
  file << "  are kept in this base class whenever possible, and placed in the derived classes\n";
  file << "  only when necessary. All variables that are identical across derived Baby\n";
  file << "  classes are implemented fully in this base class. Others are implemented with\n";
  file << "  dummy virtual accessor functions and internal member variables in the base\n";
  file << "  class, with derived classes providing a real implementation of the accessor if\n";
  file << "  the variable is defined in the corresponding ntuple format. If the variable has\n";
  file << "  inconsistent types across the ntuple formats, then each derived class must\n";
  file << "  provide all necessary accessors and internal variables; in such a case, access\n";
  file << "  through this abstract base is not possible.\n";
  file << "*/\n";

  file << "#include \"baby.hpp\"\n\n";

  file << "#include <mutex>\n";
  file << "#include <type_traits>\n";
  file << "#include <utility>\n";
  file << "#include <stdexcept>\n\n";

  file << "#include \"named_func.hpp\"\n";
  file << "#include \"utilities.hpp\"\n\n";

  file << "using namespace std;\n\n";

  file << "namespace{\n";
  file << "  using ScalarType = NamedFunc::ScalarType;\n";
  file << "  using VectorType = NamedFunc::VectorType;\n";
  file << "  using ScalarFunc = NamedFunc::ScalarFunc;\n";
  file << "  using VectorFunc = NamedFunc::VectorFunc;\n\n";

  file << "  /*!\\brief Get dummy NamedFunc in case of substitution failure\n\n";

  file << "    \\param[in] name Name of function/variable\n\n";

  file << "    \\return Dummy NamedFunc that always returns 0\n";
  file << "  */\n";
  file << "  template<typename T>\n";
  file << "    NamedFunc GetFunction(T,\n";
  file << "                          const string &name){\n";
  file << "    DBG(\"Could not find appropriate type for \\\"\" << name << \".\\\"\");\n";
  file << "    return NamedFunc(name, [](const Baby &){return 0.;}, false);\n";
  file << "  }\n\n";

  file << "  /*!\\brief Get NamedFunc for a function returning a scalar\n\n";

  file << "    \\param[in] baby_func Member function pointer to variable accessor\n\n";

  file << "    \\param[in] name Name of function/variable\n\n";

  file << "    \\return NamedFunc that returns appropriate scalar\n";
  file << "  */\n";
  file << "  template<typename T>\n";
  file << "    NamedFunc GetFunction(T const &(Baby::*baby_func)() const,\n";
  file << "                          const string &name){\n";
  file << "    return NamedFunc(name,\n";
  file << "                     [baby_func](const Baby &b){\n";
  file << "                       return ScalarType((b.*baby_func)());\n";
  file << "                     });\n";
  file << "  }\n\n";

  file << "  /*!\\brief Get NamedFunc for a function returning a vector\n\n";

  file << "    \\param[in] baby_func Member function pointer to variable accessor\n\n";

  file << "    \\param[in] name Name of function/variable\n\n";

  file << "    \\return NamedFunc that returns appropriate vectorr\n";
  file << "  */\n";
  file << "  template<typename T>\n";
  file << "    NamedFunc GetFunction(vector<T>* const &(Baby::*baby_func)() const,\n";
  file << "                          const string &name){\n";
  file << "    return NamedFunc(name,\n";
  file << "                     [baby_func](const Baby &b){\n";
  file << "                       const auto &raw = (b.*baby_func)();\n";
  file << "                       return VectorType(raw->cbegin(), raw->cend());\n";
  file << "                     });\n";
  file << "  }\n\n";

  bool have_vector_double = false;
  for(auto var = vars.cbegin(); var != vars.cend() && !have_vector_double; ++var){
    if(var->Type().find("vector<double>") != string::npos){
      have_vector_double = true;
    }
  }

  if(have_vector_double){
    file << "    template<>\n";
    file << "      NamedFunc GetFunction<vector<double>* const &(Baby::*)() const>(vector<double>* const &(Baby::*baby_func)() const,\n";
    file << "                                                                      const string &name){\n";
    file << "      return NamedFunc(name, [baby_func](const Baby &b){return *((b.*baby_func)());}, true);\n";
    file << "  }\n";
  }
  file << "}\n\n";

  file << "/*!\\brief Standard constructor\n\n";

  file << "  \\param[in] file_names ntuple files to read from\n";
  file << "*/\n";
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
  file << "    chain_->Add((file).c_str());\n";
  //file << "    chain_->Add((file+\"/tree\").c_str());\n";
  file << "  }\n";
  file << "}\n";

  file << "/*!\\brief Get number of entries in TChain and cache it\n\n";

  file << "  \\return Number of entries in TChain\n";
  file << "*/\n";
  file << "long Baby::GetEntries() const{\n";
  file << "  if(!cached_total_entries_){\n";
  file << "    cached_total_entries_ = true;\n";
  file << "    lock_guard<mutex> lock(Multithreading::root_mutex);\n";
  file << "    total_entries_ = chain_->GetEntries();\n";
  file << "  }\n";
  file << "  return total_entries_;\n";
  file << "}\n\n";

  file << "/*!\\brief Change current entry\n\n";

  file << "  \\param[in] entry Entry number to load\n";
  file << "*/\n";
  file << "void Baby::GetEntry(long entry){\n";
  for(const auto &var: vars){
    if(!var.ImplementInBase()) continue;
    file << "  c_" << var.Name() << "_ = false;\n";
  }
  file << "  lock_guard<mutex> lock(Multithreading::root_mutex);\n";
  file << "  entry_ = chain_->LoadTree(entry);\n";
  file << "}\n\n";

  file << "const std::set<std::string> & Baby::FileNames() const{\n";
  file << "  return file_names_;\n";
  file << "}\n\n";

  file << "/*! \\brief Get underlying TChain for this Baby\n\n";

  file << "  \\return Pointer to underlying TChain\n";
  file << "*/\n";
  file << "const unique_ptr<TChain> & Baby::GetTree() const{\n";
  file << "  return chain_;\n";
  file << "}\n\n";

  file << "/*! \\brief Get a NamedFunc accessing specified variable\n\n";

  file << "  \\return NamedFunc which returns specified variable from a Baby\n";
  file << "*/\n";
  file << "NamedFunc Baby::GetFunction(const std::string &var_name){\n";
  if(vars.size() != 0){
    file << "  if(var_name == \"" << vars.cbegin()->Name() << "\"){\n";
    file << "    return ::GetFunction(&Baby::" << vars.cbegin()->Name() << ", \"" << vars.cbegin()->Name() << "\");\n";
    for(auto var = ++vars.cbegin(); var != vars.cend(); ++var){
      file << "  }else if(var_name == \"" << var->Name() << "\"){\n";
      file << "    return ::GetFunction(&Baby::" << var->Name() << ", \"" << var->Name() << "\");\n";
    }
    file << "  }else{\n";
    file << "    DBG(\"Function lookup failed for \\\"\" << var_name << \".\\\"\");\n";
    file << "    return NamedFunc(var_name,\n";
    file << "                     [](const Baby &){\n";
    file << "                       return 0.;\n";
    file << "                     });\n";
    file << "  }\n";
  }else{
    file << "  DBG(\"No variables defined in Baby.\");\n";
    file << "  return NamedFunc(var_name,\n";
    file << "                   [](const Baby &){\n";
    file << "                     return 0.;\n";
    file << "                   });\n";
  }
  file << "}\n\n";

  file << "/*! \\brief Setup all branches\n";
  file << "*/\n";
  file << "void Baby::Initialize(){\n";
  file << "  lock_guard<mutex> lock(Multithreading::root_mutex);\n";
  file << "  chain_->SetMakeClass(1);\n";
  for(const auto &var: vars){
    if(!var.ImplementInBase()) continue;
    file << "  chain_->SetBranchAddress(\"" << var.Name() << "\", &" << var.Name() << "_, &b_" << var.Name() << "_);\n";
  }
  file << "}\n\n";

  for(const auto &var: vars){
    if(!var.ImplementInBase()) continue;
    file << "/*! \\brief Get " << var.Name() << " for current event and cache it\n\n";

    file << "  \\return " << var.Name() << " for current event\n";
    file << "*/\n";
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

/*!\brief Writes a derived Baby header file

  \param[in] vars All variables for all Baby classes, with type information

  \param[in] type Name of derived Baby class (basic, full, etc.)
*/
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
      file << "  " << var.DecoratedType(type) << " "
           << var.Name() << "_;//!<Cached value of " << var.Name() << '\n';
      file << "  TBranch *b_" << var.Name() << "_;\n//!<Branch from which "
           << var.Name() << " is read\n";
      file << "  mutable bool c_" << var.Name() << "_;//!<Flag if cached "
           << var.Name() << " up to date\n";
    }
  }
  file << "};\n\n";

  file << "#endif" << endl;
  file.close();
}

/*!\brief Writes a derived Baby source file

  \param[in] vars All variables for all Baby classes, with type information

  \param[in] type Name of derived Baby class (basic, full, etc.)
*/

void WriteSpecializedSource(const set<Variable> &vars, const string &type){
  ofstream file("src/baby_"+type+".cpp");
  file << "/*! \\class Baby_" << type << "\n\n";

  file << "  \\brief Derived class to access variables in " << type << " format ntuples\n\n";

  file << "  For variables not shared by all ntuple formats, the abstract base class Baby\n";
  file << "  cannot implement functions to get them, and derived classes must do the work.\n";
  file << "  This class implements getter functions for variables in the " << type << " format\n";
  file << "  ntuples, and dummy getters that throw an error for any variable not in " << type;
  file << "  format ntuples.\n";
  file << "*/\n";
  file << "#include \"baby_" << type << ".hpp\"\n\n";

  file << "#include \"utilities.hpp\"\n\n";

  file << "using namespace std;\n\n";

  file << "/*!\\brief Standard constructor\n\n";

  file << "  \\param[in] file_names ntuple files to read from\n";
  file << "*/\n";
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

  file << "/*!\\brief Change current entry\n\n";

  file << "  \\param[in] entry Entry number to load\n";
  file << "*/\n";
  file << "void Baby_" << type << "::GetEntry(long entry){\n";
  for(const auto &var: vars){
    if(var.ImplementIn(type) || var.EverythingIn(type)){
      file << "  c_" << var.Name() << "_ = false;\n";
    }
  }
  file << "  Baby::GetEntry(entry);\n";
  file << "}\n\n";

  file << "/*! \\brief Setup all branches\n";
  file << "*/\n";
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
      file << "/*!\\brief Get " << var.Name() << " for current event and cache it\n\n";

      file << "  \\return " << var.Name() << " for current event\n";
      file << "*/\n";
      file << var.DecoratedType(type) << " const & Baby_" << type << "::" << var.Name() << "() const{\n";
      file << "  if(!c_" << var.Name() << "_ && b_" << var.Name() << "_){\n";
      file << "    b_" << var.Name() << "_->GetEntry(entry_);\n";
      file << "    c_" << var.Name() << "_ = true;\n";
      file << "  }\n";
      file << "  return " << var.Name() << "_;\n";
      file << "}\n\n";
    }else if(var.VirtualInBase()){
      file << "/*!\\brief Dummy getter for " << var.Name() << ". Throws error\n\n";

      file << "  \\return Never returns. Throws error.\n";
      file << "*/\n";
      file << var.DecoratedType() << " const & Baby_" << type << "::" << var.Name() << "()  const{\n";
      file << "  ERROR(\"" << var.DecoratedType() << ' ' << var.Name()
           << " not available in babies of type " << type << ".\");\n";
      file << "}\n\n";
    }
  }

  file << flush;
  file.close();
}
