#ifndef H_TABLE_ROW
#define H_TABLE_ROW

#include <string>

#include "named_func.hpp"

class TableRow{
public:
  explicit TableRow(const std::string &label,
                    std::size_t lines_before = 1,
                    std::size_t lines_after = 1);

  TableRow(const std::string &label,
           const NamedFunc &cut,
           std::size_t lines_before = 0,
           std::size_t line_after = 0,
           const NamedFunc &weight = "weight");
  TableRow(const TableRow &) = default;
  TableRow& operator=(const TableRow &) = default;
  TableRow(TableRow &&) = default;
  TableRow& operator=(TableRow &&) = default;
  ~TableRow() = default;

  std::string label_;
  NamedFunc cut_, weight_;
  std::size_t lines_before_, lines_after_;
  bool is_data_row_;

private:
  TableRow() = delete;
};

#endif
