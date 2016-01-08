#ifndef H_GENERATE_BABY
#define H_GENERATE_BABY

#include <string>
#include <map>
#include <set>

struct SimpleVariable{
  SimpleVariable(const std::string &type, const std::string &name);
  std::string type_;
  std::string name_;
};

class Variable{
public:
  Variable(const std::string &name);

  SimpleVariable GetEntry(const std::string &baby_type) const;

  const std::string & Name() const;
  std::string & Name();

  std::string Type() const;
  std::string Type(const std::string &baby_type) const;
  std::string DecoratedType() const;
  std::string DecoratedType(const std::string &baby_type) const;

  bool HasEntry(const std::string &baby_type) const;
  void SetEntry(const std::string &baby_type,
                const std::string &type);

  bool ImplementInBase() const;
  bool VirtualInBase() const;
  bool ImplementIn(const std::string &baby_type) const;
  bool NotInBase() const;
  bool EverythingIn(const std::string &baby_type) const;

  bool operator<(const Variable& other) const;

private:
  std::set<std::string> GetTypeSet() const;

  std::string name_;
  std::map<std::string, std::string> type_map_;
};

std::set<Variable> GetVariables(const std::set<std::string> &files);

bool IsComment(const std::string &line);

SimpleVariable GetVariable(std::string line);

void RemoveExtraSpaces(std::string &line);

void WriteBaseHeader(const std::set<Variable> &vars,
                     const std::set<std::string> &types);

void WriteBaseSource(const std::set<Variable> &vars);

void WriteSpecializedHeader(const std::set<Variable> &vars,
                            const std::string &type);

void WriteSpecializedSource(const std::set<Variable> &vars,
                            const std::string &type);

void WriteMergedHeader(const std::set<Variable> &vars,
                       const std::set<std::string> &types);

#endif
