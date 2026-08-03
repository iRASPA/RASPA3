module;

module structure_parallel;

import std;

import threadpool;

namespace
{

// Set for as long as this thread is inside `forEachIndex` or `forEachBlock`. A loop reached from inside one of
// these is over work that is already being done in parallel at a coarser grain, and nesting a second lot of
// threads under the first would ask for more of them than were allowed and slow the first lot down competing
// for them.
thread_local bool insideParallelLoop = false;

ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>& threadPool()
{
  return ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance();
}

// The first exception any lane let out, kept until they have all stopped. There is no way to recall a thread
// that is already running, so the ones that have not failed run to the end of their own work and the failure is
// raised in the calling thread afterwards.
class FirstError
{
 public:
  void capture()
  {
    const std::scoped_lock lock(mutex_);
    if (!error_) error_ = std::current_exception();
  }

  void rethrowIfAny() const
  {
    if (error_) std::rethrow_exception(error_);
  }

 private:
  std::exception_ptr error_{};
  std::mutex mutex_{};
};

// How many lanes are actually run, given what the caller asked for: one where there is no pool to run them on,
// which is also the route a serial run takes.
std::size_t runnableLanes(std::size_t workers)
{
  auto& pool = threadPool();
  const std::size_t runnable =
      pool.getThreadingType() == ThreadPool::ThreadingType::ThreadPool ? pool.getThreadCount() + 1 : 1uz;
  return laneCount(workers, runnable);
}

// `lane(worker)` for every worker below `lanes`, the calling thread taking lane 0 and the pool the rest, and
// whatever a lane let out rethrown once they have all stopped.
void runLanes(std::size_t lanes, FirstError& error, const std::function<void(std::size_t worker)>& lane)
{
  auto guarded = [&](std::size_t worker)
  {
    const bool outermost = !insideParallelLoop;
    insideParallelLoop = true;
    try
    {
      lane(worker);
    }
    catch (...)
    {
      error.capture();
    }
    if (outermost) insideParallelLoop = false;
  };

  std::vector<std::future<void>> pending(lanes - 1);
  for (std::size_t worker = 1; worker < lanes; ++worker) pending[worker - 1] = threadPool().enqueue(guarded, worker);
  guarded(0);

  for (std::future<void>& task : pending) task.get();
  error.rethrowIfAny();
}

}  // namespace

std::size_t workersAvailable()
{
  if (insideParallelLoop) return 1;

  auto& pool = threadPool();
  if (pool.getThreadingType() != ThreadPool::ThreadingType::ThreadPool) return 1;
  return pool.getThreadCount() + 1;
}

std::size_t laneCount(std::size_t requested, std::size_t available)
{
  return std::min(std::max<std::size_t>(1, requested), std::max<std::size_t>(1, available));
}

std::pair<std::size_t, std::size_t> blockOfLane(std::size_t count, std::size_t lanes, std::size_t worker)
{
  if (lanes == 0 || worker >= lanes) return {count, count};

  const std::size_t share = count / lanes;
  const std::size_t remainder = count % lanes;

  // The lanes that take one index more are the first `remainder` of them, so a lane's own start is its share of
  // the whole plus however many of those came before it.
  const std::size_t begin = worker * share + std::min(worker, remainder);
  const std::size_t end = begin + share + (worker < remainder ? 1 : 0);
  return {begin, end};
}

void forEachIndex(std::size_t count, std::size_t workers,
                  const std::function<void(std::size_t worker, std::size_t index)>& body)
{
  const std::size_t lanes = runnableLanes(workers);
  if (lanes == 1 || count <= 1)
  {
    forEachIndexOfLane(count, 1, 0, [&](std::size_t index) { body(0, index); });
    return;
  }

  // Which indices a lane takes is `forEachIndexOfLane`'s business and a matter of arithmetic, not of who
  // finishes first, which is what makes the run repeatable; interleaving them rather than dealing out blocks
  // is what keeps the workers busy, one sphere or one probe radius costing an order of magnitude more than
  // another and the expensive ones coming in runs.
  FirstError error;
  runLanes(lanes, error,
           [&](std::size_t worker)
           {
             forEachIndexOfLane(count, lanes, worker,
                                [&](std::size_t index)
                                {
                                  try
                                  {
                                    body(worker, index);
                                  }
                                  catch (...)
                                  {
                                    error.capture();
                                  }
                                });
           });
}

void forEachBlock(std::size_t count, std::size_t workers,
                  const std::function<void(std::size_t worker, std::size_t begin, std::size_t end)>& body)
{
  const std::size_t lanes = runnableLanes(workers);
  if (lanes == 1 || count <= 1)
  {
    body(0, 0, count);
    return;
  }

  FirstError error;
  runLanes(lanes, error,
           [&](std::size_t worker)
           {
             const auto [begin, end] = blockOfLane(count, lanes, worker);
             body(worker, begin, end);
           });
}
