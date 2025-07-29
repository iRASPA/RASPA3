module;

#ifdef USE_LEGACY_HEADERS
#include <atomic>
#include <chrono>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <deque>
#include <exception>
#include <functional>
#include <future>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <queue>
#include <semaphore>
#include <stop_token>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>
#endif

#if __has_include(<omp.h>)
#include <omp.h>
#endif

export module threadpool;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

// https://github.com/DeveloperPaul123/thread-pool
namespace ThreadPool
{
/**
 * @brief Simple concept for the Lockable and Basic Lockable types as defined by the C++
 * standard.
 * @details See https://en.cppreference.com/w/cpp/named_req/Lockable and
 * https://en.cppreference.com/w/cpp/named_req/BasicLockable for details.
 */
template <typename Lock>
concept is_lockable = requires(Lock &&lock) {
  lock.lock();
  lock.unlock();
  { lock.try_lock() } -> std::convertible_to<bool>;
};

template <typename T, typename Lock = std::mutex>
  requires is_lockable<Lock>
class thread_safe_queue
{
 public:
  using value_type = T;
  using size_type = typename std::deque<T>::size_type;

  thread_safe_queue() = default;

  void push_back(T &&value)
  {
    std::scoped_lock lock(mutex_);
    data_.push_back(std::forward<T>(value));
  }

  void push_front(T &&value)
  {
    std::scoped_lock lock(mutex_);
    data_.push_front(std::forward<T>(value));
  }

  [[nodiscard]] bool empty() const
  {
    std::scoped_lock lock(mutex_);
    return data_.empty();
  }

  [[nodiscard]] std::optional<T> pop_front()
  {
    std::scoped_lock lock(mutex_);
    if (data_.empty()) return std::nullopt;

    auto front = std::move(data_.front());
    data_.pop_front();
    return front;
  }

  [[nodiscard]] std::optional<T> pop_back()
  {
    std::scoped_lock lock(mutex_);
    if (data_.empty()) return std::nullopt;

    auto back = std::move(data_.back());
    data_.pop_back();
    return back;
  }

  [[nodiscard]] std::optional<T> steal()
  {
    std::scoped_lock lock(mutex_);
    if (data_.empty()) return std::nullopt;

    auto back = std::move(data_.back());
    data_.pop_back();
    return back;
  }

  void rotate_to_front(const T &item)
  {
    std::scoped_lock lock(mutex_);
    auto iter = std::find(data_.begin(), data_.end(), item);

    if (iter != data_.end())
    {
      std::ignore = data_.erase(iter);
    }

    data_.push_front(item);
  }

  [[nodiscard]] std::optional<T> copy_front_and_rotate_to_back()
  {
    std::scoped_lock lock(mutex_);

    if (data_.empty()) return std::nullopt;

    auto front = data_.front();
    data_.pop_front();

    data_.push_back(front);

    return front;
  }

 private:
  std::deque<T> data_{};
  mutable Lock mutex_{};
};
}  // namespace ThreadPool

export namespace ThreadPool
{
namespace details
{

#ifdef __cpp_lib_move_only_function
using default_function_type = std::move_only_function<void()>;
#else
using default_function_type = std::function<void()>;
#endif
}  // namespace details

enum class ThreadingType : std::size_t
{
  Serial = 0,
  ThreadPool = 1,
  OpenMP = 2,
  GPU_Offload = 3
};

template <typename FunctionType = details::default_function_type, typename ThreadType = std::jthread>
  requires std::invocable<FunctionType> && std::is_same_v<void, std::invoke_result_t<FunctionType>>
class ThreadPool
{
 public:
  // singleton added
  inline static ThreadPool<details::default_function_type, std::jthread> &instance()
  {
    static auto singleton = ThreadPool();
    return singleton;
  }

  inline ThreadingType getThreadingType() { return threadingType; }

  inline std::size_t getThreadCount() { return number_of_threads; }

  template <typename InitializationFunction = std::function<void(std::size_t)>>
    requires std::invocable<InitializationFunction, std::size_t> &&
             std::is_same_v<void, std::invoke_result_t<InitializationFunction, std::size_t>>
  void init(const std::size_t nthreads, ThreadingType threading_type)
  {
    if (threading_type == ThreadingType::ThreadPool)
    {
      number_of_threads = nthreads > 0 ? nthreads - 1 : std::thread::hardware_concurrency() - 1;
      // number_of_threads =  nthreads > 0 ? nthreads:  std::thread::hardware_concurrency();
      threadingType = threading_type;
      tasks_ = std::deque<task_item>(number_of_threads);
      InitializationFunction init = [](std::size_t) {};
      std::size_t current_id = 0;
      for (std::size_t i = 0; i < number_of_threads; ++i)
      {
        priority_queue_.push_back(std::size_t(current_id));
        try
        {
          threads_.emplace_back(
              [&, id = current_id, init](const std::stop_token &stop_tok)
              {
                init(id);
                do
                {
                  // wait until signaled
                  tasks_[id].signal.acquire();

                  do
                  {
                    // invoke the task
                    while (auto task = tasks_[id].tasks.pop_front())
                    {
                      try
                      {
                        unassigned_tasks_.fetch_sub(1, std::memory_order_release);
                        std::invoke(std::move(task.value()));
                        completed_tasks_.fetch_sub(1, std::memory_order_release);
                      }
                      catch (...)
                      {
                      }
                    }

                    // try to steal a task
                    for (std::size_t j = 1; j < tasks_.size(); ++j)
                    {
                      const std::size_t index = (id + j) % tasks_.size();
                      if (auto task = tasks_[index].tasks.steal())
                      {
                        // steal a task
                        unassigned_tasks_.fetch_sub(1, std::memory_order_release);
                        std::invoke(std::move(task.value()));
                        completed_tasks_.fetch_sub(1, std::memory_order_release);
                        // stop stealing once we have invoked a stolen task
                        break;
                      }
                    }
                    // check if there are any unassigned tasks before rotating to the
                    // front and waiting for more work
                  } while (unassigned_tasks_.load(std::memory_order_acquire) > 0);

                  priority_queue_.rotate_to_front(id);
                  // check if all tasks are completed and release the barrier (binary
                  // semaphore)
                  if (completed_tasks_.load(std::memory_order_acquire) == 0)
                  {
                    threads_done_.release();
                  }

                } while (!stop_tok.stop_requested());
              });
          // increment the thread id
          ++current_id;
        }
        catch (...)
        {
          // catch all

          // remove one item from the tasks
          tasks_.pop_back();

          // remove our thread from the priority queue
          std::ignore = priority_queue_.pop_back();
        }
      }
    }
  }

  ~ThreadPool()
  {
    wait_for_tasks();

    // stop all threads
    for (std::size_t i = 0; i < threads_.size(); ++i)
    {
      threads_[i].request_stop();
      tasks_[i].signal.release();
      threads_[i].join();
    }
  }

  /// thread pool is non-copyable
  ThreadPool(const ThreadPool &) = delete;
  ThreadPool &operator=(const ThreadPool &) = delete;

  /**
   * @brief Enqueue a task into the thread pool that returns a result.
   * @details Note that task execution begins once the task is enqueued.
   * @tparam Function An invokable type.
   * @tparam Args Argument parameter pack
   * @tparam ReturnType The return type of the Function
   * @param f The callable function
   * @param args The parameters that will be passed (copied) to the function.
   * @return A std::future<ReturnType> that can be used to retrieve the returned value.
   */
  template <typename Function, typename... Args, typename ReturnType = std::invoke_result_t<Function &&, Args &&...>>
    requires std::invocable<Function, Args...>
  [[nodiscard]] std::future<ReturnType> enqueue(Function f, Args... args)
  {
#ifdef __cpp_lib_move_only_function
    // we can do this in C++23 because we now have support for move only functions
    std::promise<ReturnType> promise;
    auto future = promise.get_future();
    auto task = [func = std::move(f), ... largs = std::move(args), promise = std::move(promise)]() mutable
    {
      try
      {
        if constexpr (std::is_same_v<ReturnType, void>)
        {
          func(largs...);
          promise.set_value();
        }
        else
        {
          promise.set_value(func(largs...));
        }
      }
      catch (...)
      {
        promise.set_exception(std::current_exception());
      }
    };
    enqueue_task(std::move(task));
    return future;
#else
    /*
     * use shared promise here so that we don't break the promise later (until C++23)
     *
     * with C++23 we can do the following:
     *
     * std::promise<ReturnType> promise;
     * auto future = promise.get_future();
     * auto task = [func = std::move(f), ...largs = std::move(args),
                      promise = std::move(promise)]() mutable {...};
     */
    auto shared_promise = std::make_shared<std::promise<ReturnType>>();
    auto task = [func = std::move(f), ... largs = std::move(args), promise = shared_promise]()
    {
      try
      {
        if constexpr (std::is_same_v<ReturnType, void>)
        {
          func(largs...);
          promise->set_value();
        }
        else
        {
          promise->set_value(func(largs...));
        }
      }
      catch (...)
      {
        promise->set_exception(std::current_exception());
      }
    };

    // get the future before enqueuing the task
    auto future = shared_promise->get_future();
    // enqueue the task
    enqueue_task(std::move(task));
    return future;
#endif
  }

  /**
   * @brief Enqueue a task to be executed in the thread pool that returns void.
   * @tparam Function An invokable type.
   * @tparam Args Argument parameter pack for Function
   * @param func The callable to be executed
   * @param args Arguments that will be passed to the function.
   */
  template <typename Function, typename... Args>
    requires std::invocable<Function, Args...> && std::is_same_v<void, std::invoke_result_t<Function &&, Args &&...>>
  void enqueue_detach(Function &&func, Args &&...args)
  {
    enqueue_task(std::move(
        [f = std::forward<Function>(func), ... largs = std::forward<Args>(args)]() mutable -> decltype(auto)
        {
          // suppress exceptions
          try
          {
            std::invoke(f, largs...);
          }
          catch (...)
          {
          }
        }));
  }

  [[nodiscard]] auto size() const { return threads_.size(); }

  /**
   * @brief Wait for all tasks to finish.
   * @details This function will block until all tasks have been completed.
   */
  void wait_for_tasks()
  {
    if (completed_tasks_.load(std::memory_order_acquire) > 0)
    {
      // wait for all tasks to finish
      threads_done_.acquire();
    }
  }

 private:
  ThreadPool() : tasks_() {};

  std::size_t number_of_threads;
  ThreadingType threadingType;

  template <typename Function>
  void enqueue_task(Function &&f)
  {
    auto i_opt = priority_queue_.copy_front_and_rotate_to_back();
    if (!i_opt.has_value())
    {
      // would only be a problem if there are zero threads
      return;
    }
    auto i = *(i_opt);
    unassigned_tasks_.fetch_add(1, std::memory_order_relaxed);
    completed_tasks_.fetch_add(1, std::memory_order_relaxed);
    tasks_[i].tasks.push_back(std::forward<Function>(f));
    tasks_[i].signal.release();
  }

  struct task_item
  {
    thread_safe_queue<FunctionType> tasks{};
    std::binary_semaphore signal{0};
  };

  std::vector<ThreadType> threads_;
  std::deque<task_item> tasks_;
  thread_safe_queue<std::size_t> priority_queue_;
  std::atomic_int_fast64_t unassigned_tasks_{}, completed_tasks_{};
  std::binary_semaphore threads_done_{0};
};

/**
 * @example mandelbrot/source/main.cpp
 * Example showing how to use thread pool with tasks that return a value. Outputs a PPM image of
 * a mandelbrot.
 */
}  // namespace ThreadPool
