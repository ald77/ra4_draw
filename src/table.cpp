#include "table.hpp"

#include <fstream>

#include "RooStats/NumberCountingUtils.h"

#include "utilities.hpp"

using namespace std;

Table::TableColumn::TableColumn(const Table &table,
				const shared_ptr<Process> &process):
  FigureComponent(table, process),
  sumw_(table.rows_.size(), 0.),
  sumw2_(table.rows_.size(), 0.){
}

void Table::TableColumn::RecordEvent(const Baby &baby,
				     const NamedFunc &process_cut){
  const Table& table = static_cast<const Table&>(figure_);
  for(size_t irow = 0; irow < table.rows_.size(); ++irow){
    const TableRow& row = table.rows_.at(irow);
    if(!row.is_data_row_) continue;
    NamedFunc cut = process_cut && row.cut_;
    const NamedFunc &wgt = row.weight_;

    bool have_vector = false;
    size_t min_vec_size = 0;

    bool cut_is_scalar = cut.IsScalar();
    NamedFunc::VectorType cut_vector;
    if(cut_is_scalar){
      if(!cut.GetScalar(baby)) return;
    }else{
      cut_vector = cut.GetVector(baby);
      if(!have_vector || cut_vector.size() < min_vec_size){
	have_vector = true;
	min_vec_size = cut_vector.size();
      }
    }

    bool wgt_is_scalar = wgt.IsScalar();
    NamedFunc::ScalarType wgt_scalar = 0.;
    NamedFunc::VectorType wgt_vector;
    if(wgt_is_scalar){
      wgt_scalar = wgt.GetScalar(baby);
    }else{
      wgt_vector = wgt.GetVector(baby);
      if(!have_vector || wgt_vector.size() < min_vec_size){
	have_vector = true;
	min_vec_size = wgt_vector.size();
      }
    }

    if(!have_vector){
      sumw_.at(irow) += wgt_scalar;
      sumw2_.at(irow) += wgt_scalar*wgt_scalar;
    }else{
      for(size_t iobject = 0; iobject < min_vec_size; ++iobject){
	NamedFunc::ScalarType this_cut = cut_is_scalar ? true : cut_vector.at(iobject);
	if(!this_cut) continue;
	NamedFunc::ScalarType this_wgt = wgt_is_scalar ? wgt_scalar : wgt_vector.at(iobject);
	sumw_.at(irow) += this_wgt;
	sumw2_.at(irow) += this_wgt*this_wgt;
      }
    }
  }
}

Table::Table(const string &name,
             const vector<TableRow> &rows,
             const vector<shared_ptr<Process> > &processes):
  Figure(),
  name_(name),
  rows_(rows),
  backgrounds_(),
  signals_(),
  datas_(){
  for(const auto &process: processes){
    switch(process->type_){
    case Process::Type::data:
      datas_.emplace_back(new TableColumn(*this, process));
      break;
    case Process::Type::background:
      backgrounds_.emplace_back(new TableColumn(*this, process));
      break;
    case Process::Type::signal:
      signals_.emplace_back(new TableColumn(*this, process));
      break;
    default:
      break;
    }
  }
}

void Table::Print(double luminosity){
  std::ofstream file("tables/"+name_+".tex");
  PrintHeader(file);
  for(size_t i = 0; i < rows_.size(); ++i){
    PrintRow(file, i, luminosity);
  }
  PrintFooter(file);
}

set<shared_ptr<Process> > Table::GetProcesses() const{
  set<shared_ptr<Process> > processes;
  for(const auto &proc: backgrounds_){
    processes.insert(proc->process_);
  }
  for(const auto &proc: signals_){
    processes.insert(proc->process_);
  }
  for(const auto &proc: datas_){
    processes.insert(proc->process_);
  }
  return processes;
}

Figure::FigureComponent * Table::GetComponent(const shared_ptr<Process> &process){
  const auto &component_list = GetComponentList(process);
  for(const auto &component: component_list){
    if(component->process_ == process){
      return component.get();
    }
  }
  DBG("Could not find histogram for process "+process->name_+".");
  return nullptr;
}

const vector<unique_ptr<Table::TableColumn> >& Table::GetComponentList(const shared_ptr<Process> &process){
  switch(process->type_){
  case Process::Type::data:
    return datas_;
  case Process::Type::background:
    return backgrounds_;
  case Process::Type::signal:
    return signals_;
  default:
    ERROR("Did not understand process type "+to_string(static_cast<long>(process->type_))+".");
    return backgrounds_;
  }
}
void Table::PrintHeader(ofstream &file) const{
  file << "  \\documentclass[10pt,oneside]{report}\n";
  file << "  \\usepackage{graphicx,xspace,amssymb,amsmath,colordvi,colortbl, verbatim,multicol}\n";
  file << "  \\usepackage{multirow, rotating}\n\n";

  file << "  \\linespread{1.3}\n";
  file << "  \\addtolength{\\oddsidemargin}{-1.8in}\n";
  file << "  \\addtolength{\\textheight}{1.4in}\n";
  file << "  \\thispagestyle{empty}\n\n";

  file << "  \\begin{document}\n\n";

  file << "  \\begin{sidewaystable}[tbp!]\n";
  file << "    \\centering\n";
  file << "    \\begin{tabular}{ l";
  
  if(backgrounds_.size() > 1){
    file << " | ";
    for(size_t i = 0; i < backgrounds_.size(); ++i){
      file << 'r';
    }
    file << " | r";
  }else if(backgrounds_.size() == 1){
    file << "r";
  }
  
  if(datas_.size() > 1){
    file << " | ";
    for(size_t i = 0; i < datas_.size(); ++i){
      file << 'r';
    }
    file << " | r";
  }else if(datas_.size() == 1){
    file << "r";
  }
  
  for(size_t i = 0; i < signals_.size(); ++i){
    file << " | rr";
  }

  file << " }\n";
  file << "    \\hline\\hline\n";
  file << "    Cut";

  if(backgrounds_.size() > 1){
    for(size_t i = 0; i < backgrounds_.size(); ++i){
      file << " & " << backgrounds_.at(i)->process_->name_;
    }
    file << " & SM Bkg.";
  }else if(backgrounds_.size() == 1){
    file << " & " << backgrounds_.front()->process_->name_;
  }
  
  if(datas_.size() > 1){ 
    file << " | ";
    for(size_t i = 0; i < datas_.size(); ++i){
      file << " & " << datas_.at(i)->process_->name_;
    }
    file << " & Data Tot.";
  }else if(datas_.size() == 1){
    file << " & " << datas_.front()->process_->name_;
  }
  
  for(size_t i = 0; i < signals_.size(); ++i){
    file << " & " << signals_.at(i)->process_->name_ << " & Z_{\\text{Bi}}";
  }

  file << "\\\\\n";
}

void Table::PrintRow(ofstream &file, size_t irow, double luminosity) const{
  const TableRow& row = rows_.at(irow);
  if(row.lines_before_ > 0){
    for(size_t i = 0; i < row.lines_before_; ++i){
      file << "\\hline";
    }
    file << "\n";
  }
  file << row.label_;

  if(row.is_data_row_){
    if(backgrounds_.size() > 1){
      for(size_t i = 0; i < backgrounds_.size(); ++i){
        file << " & " << luminosity*backgrounds_.at(i)->sumw_.at(irow);
      }
      file << " & " << luminosity*GetYield(backgrounds_, irow) << "\\pm" << luminosity*GetError(backgrounds_, irow);
    }else if(backgrounds_.size() == 1){
      file << " & " << luminosity*GetYield(backgrounds_, irow) << "\\pm" << luminosity*GetError(backgrounds_, irow);
    }

    if(datas_.size() > 1){
      for(size_t i = 0; i < datas_.size(); ++i){
        file << " & " << luminosity*datas_.at(i)->sumw_.at(irow);
      }
      file << " & " << luminosity*GetYield(datas_, irow);
    }else if(datas_.size() == 1){
      file << " & " << luminosity*GetYield(datas_, irow);
    }

    for(size_t i = 0; i < signals_.size(); ++i){
      file << " & " << luminosity*signals_.at(i)->sumw_.at(irow) << " & "
           << RooStats::NumberCountingUtils::BinomialExpZ(luminosity*signals_.at(i)->sumw_.at(irow),
                                                          luminosity*GetYield(backgrounds_, irow),
                                                          hypot(GetError(backgrounds_, irow)/GetYield(backgrounds_, irow), 0.3));
    }
  }else{
    file << "\\multicolumn{1}{c}{" << row.label_ << "}";
  }

  file << "\\\\\n";

  if(row.lines_after_ > 0){
    for(size_t i = 0; i < row.lines_after_; ++i){
      file << "\\hline";
    }
    file << "\n";
  }
}

void Table::PrintFooter(ofstream &file) const{
  file << "    Cut";

  if(backgrounds_.size() > 1){
    for(size_t i = 0; i < backgrounds_.size(); ++i){
      file << " & " << backgrounds_.at(i)->process_->name_;
    }
    file << " & SM Bkg.";
  }else if(backgrounds_.size() == 1){
    file << " & " << backgrounds_.front()->process_->name_;
  }
  
  if(datas_.size() > 1){ 
    file << " | ";
    for(size_t i = 0; i < datas_.size(); ++i){
      file << " & " << datas_.at(i)->process_->name_;
    }
    file << " & Data Tot.";
  }else if(datas_.size() == 1){
    file << " & " << datas_.front()->process_->name_;
  }
  
  for(size_t i = 0; i < signals_.size(); ++i){
    file << " & " << signals_.at(i)->process_->name_ << " & Z_{\\text{Bi}}";
  }

  file << "\\\\\n";

  file << "\\end{sidewaystable}\n";
  file << "\\end{document}\n";
}

size_t Table::NumColumns() const{
  return 1
    + (backgrounds_.size() <= 1 ? backgrounds_.size() : backgrounds_.size()+1)
    + (signals_.size() <= 1 ? signals_.size() : signals_.size()+1)
    + 2*signals_.size();
}

double Table::GetYield(const vector<unique_ptr<TableColumn> > &columns,
                       size_t irow){
  double yield = 0.;
  for(const auto &column: columns){
    yield += column->sumw_.at(irow);
  }
  return yield;
}

double Table::GetError(const vector<unique_ptr<TableColumn> > &columns,
                       size_t irow){
  double error = 0.;
  for(const auto &column: columns){
    error += column->sumw2_.at(irow);
  }
  return error;
}
