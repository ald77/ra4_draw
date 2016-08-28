#ifndef H_TIMER
#define H_TIMER

#include <chrono>
#include <ostream>
#include <mutex>
#include <string>

class Timer{
public:
  using Clock = std::chrono::high_resolution_clock;
  using TimeType = Clock::time_point;

  explicit Timer(std::size_t num_iterations = 0,
                 double auto_print = -1.,
                 bool erase_lines = false);
  Timer(std::size_t num_iterations,
        std::chrono::duration<double> auto_print,
        bool erase_lines = false);
  explicit Timer(const std::string &label,
                 std::size_t num_iterations = 0,
                 double auto_print = -1.,
                 bool erase_lines = false);
  Timer(const std::string &label,
        std::size_t num_iterations,
        std::chrono::duration<double> auto_print,
        bool erase_lines = false);
  Timer(const Timer &) = default;
  Timer & operator=(const Timer &) = default;
  Timer(Timer &&) = default;
  Timer & operator=(Timer &&) = default;
  ~Timer() = default;

  void Iterate();

  void Restart();
  void Restart(std::size_t num_iterations);

  std::chrono::duration<double> ElapsedTime() const;
  std::chrono::duration<double> RemainingTime() const;

  friend std::ostream & operator<<(std::ostream &stream, const Timer &timer);

  std::size_t Iteration() const;
  Timer & Iteration(size_t iteration);

  std::size_t NumIterations() const;
  Timer & NumIterations(std::size_t num_iterations);

  std::chrono::duration<double> AutoPrintTime() const;
  Timer & AutoPrintTime(double auto_print);
  Timer & AutoPrintTime(std::chrono::duration<double> auto_print);

  const std::string & Label() const;
  Timer & Label(const std::string &label);

private:
  std::string label_;
  TimeType start_time_;
  mutable TimeType last_print_;
  std::size_t iteration_, num_iterations_;
  std::chrono::duration<double> auto_print_;
  bool erase_lines_;
  static std::mutex mutex_;
};

#endif
