module;

#include <omp.h>

export module threadpool;

import <atomic>;
import <chrono>;
import <exception>;
import <functional>;
import <future>;
import <iostream>;
import <memory>;
import <mutex>;
import <queue>;
import <deque>;
import <thread>;
import <type_traits>;
import <utility>;
import <optional>;
import <semaphore>;

import threading;



// https://github.com/DeveloperPaul123/thread-pool

template <typename T>
class thread_safe_queue 
{
  public:
    using value_type = T;
    using size_type = typename std::deque<T>::size_type;

    thread_safe_queue() = default;

    void push(T&& value) {
        std::lock_guard lock(mutex_);
        data_.push_back(std::forward<T>(value));
    }

    [[nodiscard]] bool empty() const {
        std::lock_guard lock(mutex_);
        return data_.empty();
    }

    [[nodiscard]] std::optional<T> pop() {
        std::lock_guard lock(mutex_);
        if (data_.empty()) return std::nullopt;

        auto front = std::move(data_.front());
        data_.pop_front();
        return front;
    }

    [[nodiscard]] std::optional<T> steal() {
        std::lock_guard lock(mutex_);
        if (data_.empty()) return std::nullopt;

        auto back = std::move(data_.back());
        data_.pop_back();
        return back;
    }

  private:
    using mutex_type = std::mutex;
    std::deque<T> data_{};
    mutable mutex_type mutex_{};
};


export class ThreadPool 
{
  public:
    enum class ThreadingType : size_t
    {
      Serial = 0,
      ThreadPool = 1,
      OpenMP = 2,
      GPU_Offload = 3
    };

    explicit ThreadPool(const unsigned int &number_of_threads = std::thread::hardware_concurrency(), ThreadingType type = ThreadingType::ThreadPool) 
        : threadingType(type), tasks_(number_of_threads) 
    {
      if(type == ThreadingType::OpenMP)
      {
        omp_set_dynamic(0);     // Explicitly disable dynamic teams
        omp_set_num_threads(static_cast<int>(number_of_threads + 1));
      }
      else 
      {
        for (std::size_t i = 0; i < number_of_threads; ++i) 
        {
          try 
          {
            threads_.emplace_back([&, id = i](const std::stop_token &stop_tok) 
            {
              do 
              {
                // wait until signaled
                tasks_[id].signal.acquire();

                do 
                {
                  // invoke the task
                  while (auto task = tasks_[id].tasks.pop()) 
                  {
                    try 
                    {
                      pending_tasks_.fetch_sub(1, std::memory_order_release);
                      std::invoke(std::move(task.value()));
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
                      pending_tasks_.fetch_sub(1, std::memory_order_release);
                      std::invoke(std::move(task.value()));
                      // stop stealing once we have invoked a stolen task
                      break;
                    }
                  }
                } while (pending_tasks_.load(std::memory_order_acquire) > 0);
              } while (!stop_tok.stop_requested());
            });
          } 
          catch (...) 
          {
            // catch all
          }
        }
      }
    }

    static ThreadPool &createPool(size_t count = 1, ThreadingType type = ThreadingType::ThreadPool)
    {
      static ThreadPool instance(static_cast<unsigned int>(count > 0 ? count - 1 :  std::thread::hardware_concurrency()), type);
      return instance;
    }

    static ThreadPool &instance()
    {
      return createPool();
    }

    size_t getThreadCount() {return tasks_.size();}

    ThreadingType threadingType{ThreadingType::OpenMP};

    ~ThreadPool() 
    {
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
     * @tparam Function An invokable type.
     * @tparam Args Argument parameter pack
     * @tparam ReturnType The return type of the Function
     * @param f The callable function
     * @param args The parameters that will be passed (copied) to the function.
     * @return A std::future<ReturnType> that can be used to retrieve the returned value.
     */
    template <typename Function, typename... Args,
              typename ReturnType = std::invoke_result_t<Function &&, Args &&...>>
    requires std::invocable<Function, Args...>
    [[nodiscard]] std::future<ReturnType> enqueue(Function f, Args... args) 
    {
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
        auto task = [func = std::move(f), ... largs = std::move(args),
                     promise = shared_promise]() {
            try {
                promise->set_value(func(largs...));
            } catch (...) {
                promise->set_exception(std::current_exception());
            }
        };

        // get the future before enqueuing the task
        auto future = shared_promise->get_future();
        // enqueue the task
        enqueue_task(std::move(task));
        return future;
    }

    /**
     * @brief Enqueue a task to be executed in the thread pool that returns void.
     * @tparam Function An invokable type.
     * @tparam Args Argument parameter pack for Function
     * @param func The callable to be executed
     * @param args Arguments that will be passed to the function.
     */
    template <typename Function, typename... Args>
    requires std::invocable<Function, Args...> &&
        std::is_same_v<void, std::invoke_result_t<Function &&, Args &&...>>
    void enqueue_detach(Function &&func, Args &&...args) {
        enqueue_task(
            std::move([f = std::forward<Function>(func),
                       ... largs = std::forward<Args>(args)]() mutable -> decltype(auto) {
                // suppress exceptions
                try {
                    std::invoke(f, largs...);
                } catch (...) {
                }
            }));
    }

  private:
    template <typename Function>
    void enqueue_task(Function &&f) {
        const std::size_t i = count_++ % tasks_.size();
        pending_tasks_.fetch_add(1, std::memory_order_relaxed);
        tasks_[i].tasks.push(std::forward<Function>(f));
        tasks_[i].signal.release();
    }

    struct task_item {
        thread_safe_queue<std::function<void()>> tasks{};
        std::binary_semaphore signal{0};
    };

    std::vector<std::jthread> threads_;
    std::deque<task_item> tasks_;
    std::size_t count_{};
    std::atomic_int_fast64_t pending_tasks_{};
};
