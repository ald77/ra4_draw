#ifndef H_FIGURE
#define H_FIGURE

#include <memory>
#include <vector>
#include <mutex>

#include "core/process.hpp"
#include "core/baby.hpp"
#include "core/named_func.hpp"

class Figure{
public:
  class FigureComponent{
  public:
    FigureComponent(const Figure &figure,
                    const std::shared_ptr<Process> &process);
    virtual ~FigureComponent() = default;

    virtual void RecordEvent(const Baby &baby) = 0;

    const Figure& figure_;//!<Reference to figure containing this component
    std::shared_ptr<Process> process_;//!<Process associated to this part of the figure
    std::mutex mutex_;

  private:
    FigureComponent() = delete;
    FigureComponent(const FigureComponent &) = delete;
    FigureComponent& operator=(const FigureComponent &) = delete;
    FigureComponent(FigureComponent &&) = delete;
    FigureComponent& operator=(FigureComponent &&) = delete;
  };

  Figure() = default;
  Figure(const Figure &) = default;
  Figure& operator=(const Figure &) = default;
  Figure(Figure &&) = default;
  Figure& operator=(Figure &&) = default;
  virtual ~Figure() = default;

  virtual void Print(double luminosity,
                     const std::string &subdir) = 0;

  virtual std::set<const Process*> GetProcesses() const = 0;

  virtual FigureComponent * GetComponent(const Process *process) = 0;
};

#endif
