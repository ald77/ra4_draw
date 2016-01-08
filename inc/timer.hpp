#ifndef H_TIMER
#define H_TIMER

#include <ctime>

class Timer{
public:
  explicit Timer(unsigned long num_its = 0, double auto_print = 0.);
  Timer(const Timer &) = default;
  Timer & operator=(const Timer &) = default;
  Timer(Timer &&) = default;
  Timer & operator=(Timer &&) = default;
  ~Timer() = default;  

  void SetAutoPrint(double auto_print);
  void SetNumIterations(unsigned long num_its);

  void Start();
  void Iterate();

  double GetRemainingTime() const;
  void PrintRemainingTime() const;

private:
  time_t start_time_;
  mutable time_t last_print_;
  unsigned long num_its_, cur_its_;
  double auto_print_;
};

#endif
