1) Splitting the computation in framework-molecule and molecule-molecule is _faster_
   The reason is probably that for molecule-molecule you need an additional check whether you
   compute the interaction with itself.

2) OpenMP is 20% slower for a serial run.

   900 Methane in 4x4x4 MFI, 100 cycles:

   threads   threadpool  openmp    serial
   --------------------------------------
   1         57.0 s      71.3 s    57.0 s
   2         33.3 s      38.5 s    
   4         21.2 s      21.6 s
   6         18.3 s      16.6 s
   8         16.1 s      13.8 s

3) possible OpenMP options:

   export OMP_NUM_THREADS=4
   export OMP_WAIT_POLICY=active
   export OMP_DYNAMIC=false
   export OP_PROC_BIND=true

4) Intel Threading Building Blocks:

   sudo apt install libtbb-dev
   link with: -ltbb
   
   https://www.intel.com/content/www/us/en/developer/articles/guide/get-started-with-parallel-stl.html
   
   use as:
   #include <execution>
   #include <algorithm>
   std::sort(std::execution::par_unseq, input.begin(), input.end());



5) dU/dlambda

For generic different lambdas:

if (!noCharges && rr < cutOffChargeSquared)
{
  double r = std::sqrt(rr);
  ForceFactor energyFactor = potentialCoulombGradient(forceField, scaleCoulombA * scaleCoulombB, r, chargeA, chargeB);

  energy += EnergyFactor(energyFactor.energy, (scaleCoulombA < 1.0 ? scaleCoulombB * energyFactor.dUdlambda : 0.0) +
                                              (scaleCoulombB < 1.0 ? scaleCoulombA * energyFactor.dUdlambda : 0.0));

  const double3 f = energyFactor.forceFactor * dr;

  it1->gradient += f;
  it2->gradient -= f;
}
