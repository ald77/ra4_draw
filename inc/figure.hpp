#ifndef H_FIGURE
#define H_FIGURE

#include <memory>
#include <vector>

#include "process.hpp"
#include "baby.hpp"
#include "named_func.hpp"

class Figure{
public:
  class FigureComponent{
  public:
    FigureComponent(const Figure &figure,
                    const std::shared_ptr<Process> &process);
    FigureComponent(const FigureComponent &) = default;
    FigureComponent& operator=(const FigureComponent &) = default;
    FigureComponent(FigureComponent &&) = default;
    FigureComponent& operator=(FigureComponent &&) = default;
    virtual ~FigureComponent() = default;

    virtual void RecordEvent(const Baby &baby) = 0;

    const Figure& figure_;//!<Reference to figure containing this component
    std::shared_ptr<Process> process_;//!<Process associated to this part of the figure
    
  private:
    FigureComponent() = delete;
  };

  Figure() = default;
  Figure(const Figure &) = default;
  Figure& operator=(const Figure &) = default;
  Figure(Figure &&) = default;
  Figure& operator=(Figure &&) = default;
  virtual ~Figure() = default;

  virtual void Print(double luminosity,
                     const std::string &subdir) = 0;

  virtual std::set<std::shared_ptr<Process> > GetProcesses() const = 0;

  virtual FigureComponent * GetComponent(const std::shared_ptr<Process> &process) = 0;
  bool print_figure_;
};

#endif
