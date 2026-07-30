module;

module exact_parallel;

import std;

import threadpool;

namespace
{

// Set for as long as this thread is inside `forEachIndex`. A loop reached from inside one of these is over
// work that is already being done in parallel at a coarser grain, and nesting a second lot of threads under
// the first would ask for more of them than were allowed and slow the first lot down competing for them.
thread_local bool insideParallelLoop = false;

}  // namespace

std::size_t workersAvailable()
{
  if (insideParallelLoop) return 1;

  auto& pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance();
  if (pool.getThreadingType() != ThreadPool::ThreadingType::ThreadPool) return 1;
  return pool.getThreadCount() + 1;
}

std::size_t laneCount(std::size_t requested, std::size_t available)
{
  return std::min(std::max<std::size_t>(1, requested), std::max<std::size_t>(1, available));
}

void forEachIndex(std::size_t count, std::size_t workers,
                  const std::function<void(std::size_t worker, std::size_t index)>& body)
{
  auto& pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance();

  // A caller that asks for lanes without a pool behind them gets the serial route.
  const std::size_t runnable =
      pool.getThreadingType() == ThreadPool::ThreadingType::ThreadPool ? pool.getThreadCount() + 1 : 1uz;
  const std::size_t lanes = laneCount(workers, runnable);
  if (lanes == 1 || count <= 1)
  {
    forEachIndexOfLane(count, 1, 0, [&](std::size_t index) { body(0, index); });
    return;
  }

  std::exception_ptr firstError;
  std::mutex errorMutex;

  // Which indices a lane takes is `forEachIndexOfLane`'s business and a matter of arithmetic, not of who
  // finishes first, which is what makes the run repeatable; interleaving them rather than dealing out blocks
  // is what keeps the workers busy, one sphere or one probe radius costing an order of magnitude more than
  // another and the expensive ones coming in runs.
  auto lane = [&](std::size_t worker)
  {
    const bool outermost = !insideParallelLoop;
    insideParallelLoop = true;
    forEachIndexOfLane(count, lanes, worker,
                       [&](std::size_t index)
                       {
                         try
                         {
                           body(worker, index);
                         }
                         catch (...)
                         {
                           const std::scoped_lock lock(errorMutex);
                           if (!firstError) firstError = std::current_exception();
                         }
                       });
    if (outermost) insideParallelLoop = false;
  };

  std::vector<std::future<void>> pending(lanes - 1);
  for (std::size_t worker = 1; worker < lanes; ++worker) pending[worker - 1] = pool.enqueue(lane, worker);
  lane(0);

  for (std::future<void>& task : pending) task.get();
  if (firstError) std::rethrow_exception(firstError);
}
