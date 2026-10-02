module;

export module spatial_decomposition_worker_team;

import std;

/**
 * \brief A persistent team of worker threads with a shared barrier, for the spatial-decomposition MD engine.
 *
 * The team owns `size() - 1` background `std::jthread`s; the calling thread takes part as member 0. `run(task)`
 * executes `task(member)` once on every member and returns when all of them are done; inside the task the members
 * separate their phases with `sync()`, a `std::barrier` over the whole team. Each member keeps the same index for
 * the lifetime of the team, so per-thread state (sub-domain, neighbour lists, mesh buffer) stays attached to the
 * same core and cache from step to step.
 *
 * A team of one has no background threads: `run` calls the task inline and `sync` is a no-op, so a serial run goes
 * through exactly the same code path as a parallel one.
 *
 * Exceptions thrown inside a task are caught on the throwing member and rethrown by `run` on the caller after all
 * members have finished; a member that has failed still arrives at every subsequent barrier of the task (see
 * `phase`), so the others cannot dead-lock on it.
 */
export class WorkerTeam
{
 public:
  explicit WorkerTeam(std::size_t numberOfMembers)
      : members(std::max<std::size_t>(1, numberOfMembers)), barrier(static_cast<std::ptrdiff_t>(members))
  {
    for (std::size_t member = 1; member < members; ++member)
    {
      workers.emplace_back([this, member](std::stop_token stopToken) { workerLoop(member, stopToken); });
    }
  }

  ~WorkerTeam()
  {
    {
      std::scoped_lock lock(mutex);
      stopping = true;
    }
    wakeUp.notify_all();
    for (std::jthread& worker : workers) worker.request_stop();
  }

  WorkerTeam(const WorkerTeam&) = delete;
  WorkerTeam& operator=(const WorkerTeam&) = delete;

  std::size_t size() const { return members; }

  /// Shared counter for handing out the work items of a task dynamically among the members: member 0 calls
  /// resetWork() before the task, and every member claims the first index of its next `count` items with
  /// claimWork(count) until the result reaches the number of items.
  void resetWork() { workCounter.store(0, std::memory_order_relaxed); }
  std::size_t claimWork(std::size_t count) { return workCounter.fetch_add(count, std::memory_order_relaxed); }

  /// Executes task(member) on every member (the caller is member 0) and waits for all of them.
  void run(std::function<void(std::size_t)> task)
  {
    if (members == 1)
    {
      currentTask = std::move(task);
      execute(0);
      currentTask = nullptr;
      if (failed())
      {
        std::exception_ptr error;
        {
          std::scoped_lock errorLock(errorMutex);
          error = std::exchange(firstError, nullptr);
          errorCount.store(0, std::memory_order_release);
        }
        std::rethrow_exception(error);
      }
      return;
    }

    {
      std::scoped_lock lock(mutex);
      currentTask = std::move(task);
      ++generation;
      pending = members - 1;
    }
    wakeUp.notify_all();

    execute(0);

    std::unique_lock lock(mutex);
    finished.wait(lock, [this] { return pending == 0; });
    currentTask = nullptr;

    if (failed())
    {
      std::exception_ptr error;
      {
        std::scoped_lock errorLock(errorMutex);
        error = std::exchange(firstError, nullptr);
        errorCount.store(0, std::memory_order_release);
      }
      lock.unlock();
      std::rethrow_exception(error);
    }
  }

  /// Barrier over all members; a no-op for a team of one.
  void sync()
  {
    if (members > 1) barrier.arrive_and_wait();
  }

  /**
   * \brief Runs one phase of a task on this member and then synchronizes with the team.
   *
   * An exception in the body is recorded (the first one wins) and the phase still arrives at the barrier;
   * `failed()` lets later phases skip their work once any member has failed.
   */
  template <typename Body>
  void phase(Body&& body)
  {
    if (!failed())
    {
      try
      {
        body();
      }
      catch (...)
      {
        recordError(std::current_exception());
      }
    }
    sync();
  }

  bool failed() const { return errorCount.load(std::memory_order_acquire) != 0; }

 private:
  void workerLoop(std::size_t member, std::stop_token stopToken)
  {
    std::uint64_t seen = 0;
    while (true)
    {
      {
        std::unique_lock lock(mutex);
        wakeUp.wait(lock, [&] { return stopping || generation != seen || stopToken.stop_requested(); });
        if (stopping || stopToken.stop_requested()) return;
        seen = generation;
      }
      execute(member);
      bool last = false;
      {
        std::scoped_lock lock(mutex);
        last = (--pending == 0);
      }
      if (last) finished.notify_one();
    }
  }

  void execute(std::size_t member)
  {
    try
    {
      currentTask(member);
    }
    catch (...)
    {
      recordError(std::current_exception());
    }
  }

  void recordError(std::exception_ptr error)
  {
    std::scoped_lock lock(errorMutex);
    if (!firstError) firstError = error;
    errorCount.fetch_add(1, std::memory_order_release);
  }

  std::size_t members;
  std::barrier<> barrier;
  std::vector<std::jthread> workers;

  std::mutex mutex;
  std::condition_variable wakeUp;
  std::condition_variable finished;
  std::function<void(std::size_t)> currentTask;
  std::uint64_t generation{0};
  std::size_t pending{0};
  bool stopping{false};

  std::mutex errorMutex;
  std::exception_ptr firstError;
  std::atomic<std::size_t> workCounter{0};
  std::atomic<std::size_t> errorCount{0};
};
