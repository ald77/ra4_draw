#ifndef H_TABLE
#define H_TABLE

#include <memory>
#include <vector>
#include <string>
#include <fstream>

#include "core/figure.hpp"
#include "core/table_row.hpp"
#include "core/process.hpp"
#include "core/gamma_params.hpp"
#include "core/plot_opt.hpp"

class Table final: public Figure{
public:
  class TableColumn final: public Figure::FigureComponent{
  public:
    TableColumn(const Table &table,
		const std::shared_ptr<Process> &process);
    ~TableColumn() = default;

    void RecordEvent(const Baby &baby) final;

    std::vector<double> sumw_, sumw2_;

  private:
    TableColumn() = delete;
    TableColumn(const TableColumn &) = delete;
    TableColumn& operator=(const TableColumn &) = delete;
    TableColumn(TableColumn &&) = delete;
    TableColumn& operator=(TableColumn &&) = delete;

    std::vector<NamedFunc> proc_and_table_cut_;
    NamedFunc::VectorType cut_vector_, wgt_vector_, val_vector_;
  };

  Table(const std::string &name,
	const std::vector<TableRow> &rows,
	const std::vector<std::shared_ptr<Process> > &processes,
	bool do_zbi=true,
	bool print_table=true,
	bool print_pie=false,
	bool print_titlepie=true);
  Table(Table &&) = default;
  Table& operator=(Table &&) = default;
  ~Table() = default;

  void Print(double luminosity,
             const std::string &subdir) final;
  
  std::vector<GammaParams> Yield(const Process *process, double luminosity) const;
  std::vector<GammaParams> BackgroundYield(double luminosity) const;
  std::vector<GammaParams> DataYield() const;
  
  std::set<const Process*> GetProcesses() const final;

  FigureComponent * GetComponent(const Process *process) final;
  
  std::string name_;
  std::vector<TableRow> rows_;
  bool do_zbi_;
  bool print_table_;
  bool print_pie_;
  bool print_titlepie_;
  std::vector<PlotOpt> plot_options_;//!<Styles with which to draw pie chart

private:
  std::vector<std::unique_ptr<TableColumn> > backgrounds_;//!<Background components of the figure
  std::vector<std::unique_ptr<TableColumn> > signals_;//!<Signal components of the figure
  std::vector<std::unique_ptr<TableColumn> > datas_;//!<Data components of the figure

  Table(const Table &) = delete;
  Table& operator=(const Table &) = delete;
  Table() = delete;

  const std::vector<std::unique_ptr<TableColumn> >& GetComponentList(const Process *process) const;

  void PrintHeader(std::ofstream &file, double luminosity) const;
  void PrintRow(std::ofstream &file, std::size_t irow, double luminosity) const;
  void PrintPie(std::size_t irow, double luminosity) const;
  void PrintFooter(std::ofstream &file) const;

  std::size_t NumColumns() const;

  static double GetYield(const std::vector<std::unique_ptr<TableColumn> > &columns,
                         std::size_t irow);
  static double GetError(const std::vector<std::unique_ptr<TableColumn> > &columns,
                         std::size_t irow);
};

#endif
