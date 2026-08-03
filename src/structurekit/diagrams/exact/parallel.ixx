module;

export module exact_parallel;

import std;

import threadpool;

// Running the independent parts of an exact analysis at the same time.
//
// The expensive loops of these analyses are over things that do not consult one another: one sphere's exposed
// surface is settled by that sphere and its neighbours, one probe diameter's geometry is built from the same
// fixed atoms as every other. So the loops are over independent work and the only question is how to hand out
// the indices and what to do about the answers.
//
// Two properties are wanted of the answer, and they are not the same property.
//
// It has to be the same every run. A structural analysis is a measurement of a fixed input, and a measurement
// that moves in the last digits between two runs of the same command is one nobody can check. So which worker
// takes which index is fixed by the index and the number of workers, and never by which worker happens to be
// free: `forEachIndex` deals index i to worker i modulo the count, and a caller that accumulates into
// per-worker partials and then reduces them in worker order gets an answer that depends on nothing but the
// number of workers.
//
// It cannot be the same as the serial answer, and this is worth being plain about rather than quiet about. A
// sum of floating-point numbers depends on the order they are added in, and a reduction over partials adds
// them in a different order from a single running total however carefully it is arranged. So an analysis run
// on several threads can differ from the same analysis run on one, in the last two or three digits of a
// quantity that is a sum of some hundreds of thousands of arcs. Serial is the default, and the serial route
// is the one that is unchanged; asking for threads is asking for the answer sooner and not for the same bits.
//
// Loops that write one slot per index rather than accumulate --- the pore-size sweep is one, each diameter
// filling its own row --- have no reduction and so no difference at all, on any number of threads.

// How many indices can be under way at once: the pool's helper threads and the calling thread, or one where
// there is no pool. One inside a loop that is already running here, so that a parallel loop reached from
// inside another one runs serially rather than asking for threads that are already busy.
export std::size_t workersAvailable();

// How many lanes a request for `requested` of them turns into when there are `available` to run them on: never
// none, and never more than either number. More lanes than there are threads would leave the extra ones queued
// behind the calling thread's own lane, which finishes them, but only after everything else.
export std::size_t laneCount(std::size_t requested, std::size_t available);

// The indices lane `worker` of `lanes` takes, in the order it takes them: `worker`, `worker + lanes`, and so on
// below `count`. This is the whole of the dealing-out. It is a function of its arguments and of nothing else,
// which is what makes a run repeatable and what lets the dealing-out be checked without any threads to run it
// on; `forEachIndex` is this loop, once per lane, with somewhere to run the lanes.
export template <typename Body>
void forEachIndexOfLane(std::size_t count, std::size_t lanes, std::size_t worker, Body&& body)
{
  if (lanes == 0) return;
  for (std::size_t index = worker; index < count; index += lanes) body(index);
}

// `body(worker, index)` for every index below `count`, on `workers` threads including the calling one.
// `worker` is below `workers` and no two calls with the same worker run at the same time, so a caller may
// give each worker scratch of its own and an accumulator of its own and share nothing else.
//
// An exception thrown by `body` comes back out of this call. On one worker it comes out at once and the
// indices behind it are not reached; on several, each worker carries on with its own and the first exception
// is rethrown once they have all stopped, there being no way to recall a thread that is already running.
export void forEachIndex(std::size_t count, std::size_t workers,
                         const std::function<void(std::size_t worker, std::size_t index)>& body);
