#include "thread_pool.hpp"

using namespace std;

ThreadPool::ThreadPool():
  tasks_(),
  threads_(),
  stop_thread_now_(),
  stop_at_empty_(),
  mutex_(),
  cv_(){
  size_t num_threads = thread::hardware_concurrency();
  if(num_threads > 2){
    --num_threads;
  }else{
    num_threads = 1;
  }
  Resize(num_threads);
}

ThreadPool::ThreadPool(size_t num_threads):
  tasks_(),
  threads_(),
  stop_thread_now_(),
  stop_at_empty_(),
  mutex_(),
  cv_(){
  Resize(num_threads);
}

ThreadPool::~ThreadPool(){
  stop_at_empty_ = true;
  {
    lock_guard<mutex> lock(mutex_);
    cv_.notify_all();
  }

  for(size_t ithread = 0; ithread < Size(); ++ithread){
    if(threads_.at(ithread)->joinable()){
      threads_.at(ithread)->join();
    }
  }
}

size_t ThreadPool::Size() const{
  return threads_.size();
}

void ThreadPool::Resize(size_t num_threads){
  size_t old_num_threads = Size();
  if(num_threads > old_num_threads){
    threads_.resize(num_threads);
    stop_thread_now_.resize(num_threads);

    for(size_t ithread = old_num_threads; ithread < num_threads; ++ithread){
      threads_.at(ithread).reset(new thread(&ThreadPool::DoTasksFromQueue, this, ithread));
    }
  }else if(num_threads < old_num_threads){
    for(size_t ithread = num_threads; ithread < old_num_threads; ++ithread){
      *stop_thread_now_.at(ithread) = true;
      threads_.at(ithread)->detach();
    }

    {
      lock_guard<mutex> lock(mutex_);
      cv_.notify_all();
    }

    threads_.resize(num_threads);
    stop_thread_now_.resize(num_threads);
  }
}

void ThreadPool::DoTasksFromQueue(size_t ithread){
  unique_ptr<function<void()> > task = tasks_.Pop();
  while(true){
    while(task != nullptr){
      (*task)();
      if(stop_thread_now_.at(ithread)) return;
      task = tasks_.Pop();
    }

    unique_lock<mutex> lock(mutex_);
    cv_.wait(lock, bind(&ThreadPool::ReadyToAct, this, ithread, ref(task)));
    if(task == nullptr) return;
  }
}

bool ThreadPool::ReadyToAct(size_t ithread, unique_ptr<function<void()> > &task){
  if(stop_thread_now_.at(ithread)){
    return true;
  }else{
    task = tasks_.Pop();
    return stop_at_empty_ || task != nullptr;
  }
}

void ThreadPool::Queue::Push(FuncPtr &func){
  lock_guard<mutex> lock(mutex_);
  queue_.push(unique_ptr<function<void()> >());
  queue_.back() = move(func);
}

ThreadPool::Queue::FuncPtr ThreadPool::Queue::Pop(){
  lock_guard<mutex> lock(mutex_);
  if(queue_.empty()){
    return FuncPtr();
  }else{
    Queue::FuncPtr func = move(queue_.front());
    queue_.pop();
    return func;
  }
}
