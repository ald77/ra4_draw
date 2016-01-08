#include "timer.hpp"
#include <ctime>
#include <cstdio>
#include <cmath>

Timer::Timer(unsigned long num_its, double auto_print):
  start_time_(0),
  last_print_(0),
  num_its_(num_its),
  cur_its_(0),
  auto_print_(auto_print){
  }

void Timer::SetAutoPrint(double auto_print){
  auto_print_ = auto_print;
}

void Timer::SetNumIterations(unsigned long num_its){
  num_its_=num_its;
}

void Timer::Start(){
  cur_its_=0;
  time(&start_time_);
  last_print_ = start_time_;
}

void Timer::Iterate(){
  ++cur_its_;
  if(auto_print_>0){
    time_t cur_time;
    time(&cur_time);
    if(difftime(cur_time, last_print_)>auto_print_){
      PrintRemainingTime();
    }
  }
}

double Timer::GetRemainingTime() const{
  if(cur_its_==0){
    return 0.0;
  }else{
    time_t cur_time;
    time(&cur_time);
    return (difftime(cur_time,start_time_)*(num_its_-cur_its_))/cur_its_;
  }
}

void Timer::PrintRemainingTime() const{
  int secs(static_cast<int>(floor(GetRemainingTime()+0.5)));
  time_t endtime;
  time(&endtime);
  endtime+=secs;
  const short hours(secs/3600);
  secs-=3600*hours;
  const short minutes(secs/60);
  secs-=60*minutes;
  printf("Iteration %10ld of %10ld. %4hd:",cur_its_,num_its_,hours);
  if(minutes<10){
    printf("0");
  }
  printf("%hd:",minutes);
  if(secs<10){
    printf("0");
  }
  printf("%d remaining. Expected finish: %s",secs, ctime(&endtime));
  fflush(stdout);
  time(&last_print_);
}
