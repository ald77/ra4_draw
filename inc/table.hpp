#ifndef H_TABLE
#define H_TABLE

#include <memory>
#include <vector>
#include <string>
#include <fstream>

#include "figure.hpp"
#include "table_row.hpp"
#include "process.hpp"

class Table final: public Figure{
public:
  class TableColumn final: public Figure::FigureComponent{
  public:
    TableColumn(const Table &table,
		const std::shared_ptr<Process> &process);
    TableColumn(const TableColumn &) = default;
    TableColumn& operator=(const TableColumn &) = default;
    TableColumn(TableColumn &&) = default;
    TableColumn& operator=(TableColumn &&) = default;
    ~TableColumn() = default;

    void RecordEvent(const Baby &baby,
                     const NamedFunc &process_cut) final;

    std::vector<double> sumw_, sumw2_;

  private:
    TableColumn() = delete;
  };

  Table(const std::string &name,
	const std::vector<TableRow> &rows,
	const std::vector<std::shared_ptr<Process> > &processes);
  Table(const Table &) = default;
  Table& operator=(const Table &) = default;
  Table(Table &&) = default;
  Table& operator=(Table &&) = default;
  ~Table() = default;

  void Print(double luminosity) final;

  std::set<std::shared_ptr<Process> > GetProcesses() const final;

  FigureComponent * GetComponent(const std::shared_ptr<Process> &process) final;

  std::string name_;
  std::vector<TableRow> rows_;
private:
  std::vector<std::unique_ptr<TableColumn> > backgrounds_;//!<Background components of the figure
  std::vector<std::unique_ptr<TableColumn> > signals_;//!<Signal components of the figure
  std::vector<std::unique_ptr<TableColumn> > datas_;//!<Data components of the figure

  Table() = delete;

  const std::vector<std::unique_ptr<TableColumn> >& GetComponentList(const std::shared_ptr<Process> &process);

  void PrintHeader(std::ofstream &file) const;
  void PrintRow(std::ofstream &file, std::size_t irow, double luminosity) const;
  void PrintFooter(std::ofstream &file) const;

  std::size_t NumColumns() const;

  static double GetYield(const std::vector<std::unique_ptr<TableColumn> > &columns,
                         std::size_t irow);
  static double GetError(const std::vector<std::unique_ptr<TableColumn> > &columns,
                         std::size_t irow);
};

#endif
