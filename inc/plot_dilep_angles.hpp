#ifndef H_DILEP_ANGLES
#define H_DILEP_ANGLES

#include <cstddef>

#include <memory>
#include <string>
#include <set>

#include "baby.hpp"
#include "process.hpp"
#include "named_func.hpp"

template<typename T>
std::shared_ptr<Process> Proc(const std::string process_name, Process::Type type,
			      int color, const std::set<std::string> &files, const std::string &cut = "1"){
  return std::make_shared<Process>(process_name, type, color,
				   std::unique_ptr<Baby>(new T(files)),
				   cut);
}

bool IsGoodMuon(const Baby &b, std::size_t imu);
bool IsGoodElectron(const Baby &b, std::size_t iel);
bool IsGoodTrack(const Baby &b, std::size_t itk);
bool IsGoodJet(const Baby &b, std::size_t ijet);
void GetAngles(const Baby &b, double &phi1, double &phi2);
NamedFunc::ScalarType MinDeltaPhiLepMet(const Baby &b);
NamedFunc::ScalarType MaxDeltaPhiLepMet(const Baby &b);
NamedFunc::ScalarType MinDeltaPhiLepJet(const Baby &b);
NamedFunc::ScalarType MaxDeltaPhiLepJet(const Baby &b);
NamedFunc::ScalarType MinDeltaPhiMetJet(const Baby &b);
NamedFunc::ScalarType MaxDeltaPhiMetJet(const Baby &b);

#endif
