#include "timer.hpp"

#include <cmath>

#include <iostream>
#include <iomanip>

using namespace std;

namespace{
  void PrintTime(ostream &stream, double seconds){
    int seconds_round = round(seconds);
    int hours = seconds_round/3600;
    seconds_round -= 3600*hours;
    int minutes = seconds_round/60;
    seconds_round -= 60*minutes;
    auto fill = stream.fill();
    stream.fill('0');
    if(hours){
      stream << hours << ':' << setw(2) << minutes << ':' << setw(2) << seconds_round;
    }else{
      stream << minutes << ':' << setw(2) << seconds_round;
    }
    stream.fill(fill);
  }
}

mutex Timer::mutex_{};

Timer::Timer(size_t num_iterations,
	     double auto_print,
	     bool erase_lines):
  start_time_(Clock::now()),
  last_print_(start_time_),
  iteration_(static_cast<size_t>(-1)),
  num_iterations_(num_iterations),
  auto_print_(auto_print),
  erase_lines_(erase_lines){
}

Timer::Timer(size_t num_iterations,
	     chrono::duration<double> auto_print,
	     bool erase_lines):
  start_time_(Clock::now()),
  last_print_(start_time_),
  iteration_(static_cast<size_t>(-1)),
  num_iterations_(num_iterations),
  auto_print_(auto_print),
  erase_lines_(erase_lines){
}

void Timer::Iterate(){
  ++iteration_;
  if(auto_print_.count() >= 0.
     && chrono::duration<double>(Clock::now() - last_print_) >= auto_print_){
    if(erase_lines_) clog << "\r\33[2K" << *this;
    else clog << *this << '\n';
  }
}

void Timer::Restart(){
  iteration_ = -1;
  start_time_ = Clock::now();
  last_print_ = start_time_;
}

void Timer::Restart(size_t num_iterations){
  num_iterations_ = num_iterations;
  iteration_ = -1;
  start_time_ = Clock::now();
  last_print_ = start_time_;
}

chrono::duration<double> Timer::ElapsedTime() const{
  return chrono::duration<double>(Clock::now() - start_time_);
}

chrono::duration<double> Timer::RemainingTime() const{
  if(iteration_ == static_cast<size_t>(-1) || iteration_ == 0){
    return static_cast<chrono::duration<double> >(0.);
  }else{
    return ElapsedTime()*(num_iterations_ - iteration_)/iteration_;
  }
}

ostream & operator<<(ostream &stream, const Timer &timer){
  double elapsed = timer.ElapsedTime().count();
  double remaining = timer.RemainingTime().count();
  auto eta_raw = chrono::system_clock::now()
    +static_cast<chrono::seconds>(static_cast<chrono::seconds::rep>(remaining));
  auto eta_time = chrono::system_clock::to_time_t(eta_raw);
  string eta = ctime(&eta_time);
  if(eta.size() && eta.back() == '\n') eta.pop_back();
  {
    lock_guard<mutex> lock(Timer::mutex_);
    stream << timer.Iteration() << '/' << timer.NumIterations() << " in ";
    PrintTime(stream, elapsed);
    auto precision = stream.precision();
    stream.precision(3);
    if(timer.NumIterations()){
      stream << " (" << 100.*timer.Iteration()/timer.NumIterations() << "%, ";
    }else{
      stream << " (" << 100. << "%, ";
    }
    if(elapsed > 0.){
      stream << timer.Iteration()/elapsed << " Hz). ";
    }else{
      stream << 0. << " Hz). ";
    }
    stream.precision(precision);
    PrintTime(stream, remaining);
    stream << " left. ETA: " << eta;
  }
  timer.last_print_ = Timer::Clock::now();
  return stream;
}

size_t Timer::Iteration() const{
  return iteration_;
}

Timer & Timer::Iteration(size_t iteration){
  iteration_ = iteration;
  return *this;
}

size_t Timer::NumIterations() const{
  return num_iterations_;
}

Timer & Timer::NumIterations(size_t num_iterations){
  num_iterations_ = num_iterations;
  return *this;
}

chrono::duration<double> Timer::AutoPrintTime() const{
  return auto_print_;
}

Timer & Timer::AutoPrintTime(double auto_print){
  auto_print_ = static_cast<chrono::duration<double> >(auto_print);
  return *this;
}

Timer & Timer::AutoPrintTime(chrono::duration<double> auto_print){
  auto_print_ = auto_print;
  return *this;
}
