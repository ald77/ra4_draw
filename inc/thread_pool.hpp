#ifndef H_THREAD_POOL
#define H_THREAD_POOL

#include <algorithm>
#include <thread>
#include <future>
#include <memory>
#include <functional>
#include <queue>
#include <mutex>
#include <vector>
#include <atomic>
#include <condition_variable>

class ThreadPool{
public:
  ThreadPool();
  explicit ThreadPool(std::size_t num_threads);
  ~ThreadPool();

  std::size_t Size() const;
  void Resize(size_t num_threads);
  
  template<typename FuncType, typename...ArgTypes>
  auto Push(FuncType &&func, ArgTypes&&... args) -> std::future<decltype(func(args...))>;

private:
  ThreadPool(const ThreadPool &) = delete;
  ThreadPool& operator=(const ThreadPool &) = delete;
  ThreadPool(ThreadPool &&) = delete;
  ThreadPool& operator=(ThreadPool &&) = delete;

  void DoTasksFromQueue(size_t ithread);
  bool ReadyToAct(size_t ithread, std::unique_ptr<std::function<void()> > &task);

  class Queue{
  public:
    using FuncPtr = std::unique_ptr<std::function<void()> >;

    Queue() = default;
    ~Queue() = default;

    void Push(FuncPtr &func);
    FuncPtr Pop();

  private:
    Queue(const Queue &) = delete;
    Queue& operator=(const Queue &) = delete;
    Queue(Queue &&) = delete;
    Queue& operator=(Queue &&) = delete;

    std::queue<FuncPtr> queue_;
    std::mutex mutex_;
  };

  Queue tasks_;
  std::vector<std::unique_ptr<std::thread> > threads_;
  std::vector<std::shared_ptr<std::atomic<bool> > > stop_thread_now_;
  std::atomic<bool> stop_at_empty_;

  std::mutex mutex_;
  std::condition_variable cv_;
};

template<typename FuncType, typename...ArgTypes>
auto ThreadPool::Push(FuncType &&func, ArgTypes&&... args) -> std::future<decltype(func(args...))>{
  auto task =  std::make_shared<std::packaged_task<decltype(func(args...))()> >(std::bind(std::forward<FuncType>(func), std::forward<ArgTypes>(args)...));
  std::unique_ptr<std::function<void()> > pkg_func(new std::function<void()>([task](){(*task)();}));
  tasks_.Push(pkg_func);
  std::lock_guard<std::mutex> lock(mutex_);
  cv_.notify_one();
  return task->get_future();
}

#endif
