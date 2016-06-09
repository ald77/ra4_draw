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

    virtual void RecordEvent(const Baby &baby,
                             const NamedFunc &process_cut) = 0;

    const Figure& figure_;//!<Reference to figure containing this component
    std::shared_ptr<Process> process_;//!<Process associated to this part of the figure
    
  private:
    FigureComponent() = delete;
  };

  Figure(const std::vector<std::shared_ptr<Process> > &processes);
  Figure(const Figure &) = default;
  Figure& operator=(const Figure &) = default;
  Figure(Figure &&) = default;
  Figure& operator=(Figure &&) = default;
  virtual ~Figure() = default;

  virtual void Print(double luminosity) = 0;

  std::set<std::shared_ptr<Process> > GetProcesses() const;

  FigureComponent * GetComponent(const std::shared_ptr<Process> &process);

protected:
  std::vector<std::unique_ptr<FigureComponent> > backgrounds_;//!<Background components of the figure
  std::vector<std::unique_ptr<FigureComponent> > signals_;//!<Signal components of the figure
  std::vector<std::unique_ptr<FigureComponent> > datas_;//!<Data components of the figure

private:
  Figure() = delete;

  std::vector<std::unique_ptr<FigureComponent> >& GetComponentList(const std::shared_ptr<Process> &process);
};

#endif
