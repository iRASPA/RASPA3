#include <gtest/gtest.h>

import std;

import archive;
import double3;
import simd_quatd;
import atom;
import molecule;
import molecule_property_settings;
import property_end_to_end_acf;

namespace
{
EndToEndACFSettings acfSettings(std::size_t numberOfBlockElements, std::size_t sampleEvery,
                                std::optional<std::size_t> writeEvery)
{
  return EndToEndACFSettings{
      .sampleEvery = sampleEvery, .writeEvery = writeEvery, .numberOfBlockElements = numberOfBlockElements};
}

// Ornstein-Uhlenbeck process for the end-to-end vector: R(t + dt) = R(t) e^{-dt/tau} + sigma sqrt(1 - e^{-2 dt/tau}) xi
// with xi standard normal per Cartesian component, so that <R(0).R(t)> = 3 sigma^2 e^{-t/tau}.
struct OrnsteinUhlenbeck
{
  double tau;
  double sigma;
  double dt;
  std::mt19937_64 generator{12345};
  std::normal_distribution<double> normal{0.0, 1.0};

  double3 initial() { return sigma * double3(normal(generator), normal(generator), normal(generator)); }

  double3 step(const double3 &r)
  {
    const double decay = std::exp(-dt / tau);
    const double noise = sigma * std::sqrt(1.0 - decay * decay);
    return decay * r + noise * double3(normal(generator), normal(generator), normal(generator));
  }
};

template <typename Accessor>
std::optional<double> interpolate(const std::vector<EndToEndAutoCorrelationFunctionData> &data, double time,
                                  Accessor column)
{
  for (std::size_t i = 1; i < data.size(); ++i)
  {
    if (data[i].time >= time)
    {
      const double f = (time - data[i - 1].time) / (data[i].time - data[i - 1].time);
      return column(data[i - 1]) + f * (column(data[i]) - column(data[i - 1]));
    }
  }
  return std::nullopt;
}

std::optional<double> interpolateNormalized(const std::vector<EndToEndAutoCorrelationFunctionData> &data, double time)
{
  return interpolate(data, time, [](const EndToEndAutoCorrelationFunctionData &p) { return p.normalized; });
}

double endToEndSquaredACF(const EndToEndAutoCorrelationFunctionData &p) { return p.endToEndSquaredACF.value(); }
double radiusOfGyrationSquaredACF(const EndToEndAutoCorrelationFunctionData &p)
{
  return p.radiusOfGyrationSquaredACF.value();
}

// Rodrigues rotation of v about the unit axis by angle
double3 rotate(const double3 &v, const double3 &axis, double angle)
{
  const double c = std::cos(angle), s = std::sin(angle);
  return c * v + s * double3::cross(axis, v) + (1.0 - c) * double3::dot(axis, v) * axis;
}
}  // namespace

TEST(end_to_end_acf, lags_follow_the_order_n_blocking_scheme)
{
  // one component of one molecule, n = 4 block elements, sampled every cycle with a 0.5 ps step
  PropertyEndToEndAutoCorrelationFunction property({1}, {std::array<std::size_t, 2>{0, 1}},
                                                   {acfSettings(4, 1, std::nullopt)}, 1, 0.5);

  // constant vector: C(t) = |R|^2 at every lag
  const double3 r(1.0, 2.0, 2.0);
  for (std::size_t i = 0; i < 70; ++i) property.addSampleVectors(std::span<const double3>(&r, 1));

  const std::vector<EndToEndAutoCorrelationFunctionData> data = property.result(0);
  ASSERT_FALSE(data.empty());

  // in samples: block 0 holds lags 0, 1, 2, 3; block 1: 4, 8, 12; block 2: 16, 32, 48; block 3 (lag 64) needs
  // 128 samples before its first non-zero lag has data
  std::vector<double> expectedLags{0.0, 0.5, 1.0, 1.5, 2.0, 4.0, 6.0, 8.0, 16.0, 24.0};
  ASSERT_GE(data.size(), expectedLags.size());
  for (std::size_t i = 0; i < expectedLags.size(); ++i)
  {
    EXPECT_DOUBLE_EQ(data[i].time, expectedLags[i]);
    EXPECT_NEAR(data[i].acf, 9.0, 1e-12);
    EXPECT_NEAR(data[i].normalized, 1.0, 1e-12);
  }
  // the lags are strictly increasing, no duplicates across blocks
  for (std::size_t i = 1; i < data.size(); ++i) EXPECT_GT(data[i].time, data[i - 1].time);
  EXPECT_EQ(static_cast<std::size_t>(data[0].numberOfSamples), 70uz);
}

TEST(end_to_end_acf, recovers_the_relaxation_time_of_an_exponentially_decorrelating_vector)
{
  const std::size_t numberOfMolecules = 64;
  const double tau = 20.0;   // ps
  const double sigma = 3.0;  // Angstrom per component
  const double dt = 1.0;     // ps between samples (time step 0.1 ps, sampled every 10 cycles)

  PropertyEndToEndAutoCorrelationFunction property({numberOfMolecules}, {std::array<std::size_t, 2>{0, 1}},
                                                   {acfSettings(25, 10, std::nullopt)}, numberOfMolecules, 0.1);

  OrnsteinUhlenbeck process{tau, sigma, dt};
  std::vector<double3> vectors(numberOfMolecules);
  for (double3 &r : vectors) r = process.initial();

  const std::size_t numberOfSamples = 40000;
  for (std::size_t s = 0; s < numberOfSamples; ++s)
  {
    property.addSampleVectors(vectors);
    for (double3 &r : vectors) r = process.step(r);
  }

  const std::vector<EndToEndAutoCorrelationFunctionData> data = property.result(0);
  ASSERT_GT(data.size(), 10uz);

  // C(0) = <R^2> = 3 sigma^2
  EXPECT_NEAR(data[0].acf, 3.0 * sigma * sigma, 0.05 * 3.0 * sigma * sigma);

  // C(tau)/C(0) = 1/e, C(2 tau)/C(0) = 1/e^2
  const std::optional<double> atTau = interpolateNormalized(data, tau);
  const std::optional<double> atTwoTau = interpolateNormalized(data, 2.0 * tau);
  ASSERT_TRUE(atTau.has_value());
  ASSERT_TRUE(atTwoTau.has_value());
  EXPECT_NEAR(atTau.value(), std::exp(-1.0), 0.04);
  EXPECT_NEAR(atTwoTau.value(), std::exp(-2.0), 0.04);

  // the three tau_R estimates agree with the input relaxation time
  const EndToEndRelaxationTimes times = property.relaxationTimes(0);
  EXPECT_NEAR(times.meanSquaredEndToEnd, 3.0 * sigma * sigma, 0.05 * 3.0 * sigma * sigma);
  ASSERT_TRUE(times.oneOverE.has_value());
  EXPECT_NEAR(times.oneOverE.value(), tau, 0.12 * tau);
  ASSERT_TRUE(times.exponentialFit.has_value());
  EXPECT_NEAR(times.exponentialFit.value(), tau, 0.12 * tau);
  ASSERT_TRUE(times.integrated.has_value());
  EXPECT_NEAR(times.integrated.value(), tau, 0.20 * tau);
}

TEST(end_to_end_acf, end_to_end_vectors_are_taken_from_the_molecule_atoms)
{
  // component 0: two molecules of three atoms, ends (0, 2); component 1: one molecule without end-to-end atoms
  std::vector<std::size_t> numberOfMoleculesPerComponent{2, 1};
  std::vector<std::optional<std::array<std::size_t, 2>>> ends{std::array<std::size_t, 2>{0, 2}, std::nullopt};

  std::vector<Molecule> molecules;
  std::vector<Atom> atoms;
  auto addMolecule = [&](std::size_t componentId, std::vector<double3> positions)
  {
    Molecule molecule(double3(0.0, 0.0, 0.0), simd_quatd(0.0, 0.0, 0.0, 1.0), 1.0, componentId, positions.size());
    molecule.atomIndex = atoms.size();
    for (const double3 &p : positions)
    {
      Atom atom;
      atom.position = p;
      atoms.push_back(atom);
    }
    molecules.push_back(molecule);
  };
  addMolecule(0, {double3(0.0, 0.0, 0.0), double3(1.0, 0.0, 0.0), double3(3.0, 4.0, 0.0)});  // |R| = 5
  addMolecule(0, {double3(1.0, 1.0, 1.0), double3(0.0, 0.0, 0.0), double3(1.0, 1.0, 3.0)});  // |R| = 2
  addMolecule(1, {double3(0.0, 0.0, 0.0), double3(9.0, 9.0, 9.0)});

  // component 1 is not sampled (no settings)
  PropertyEndToEndAutoCorrelationFunction property(numberOfMoleculesPerComponent, ends,
                                                   {acfSettings(25, 5, std::nullopt), std::nullopt}, molecules.size(),
                                                   0.001);

  property.addSample(0, molecules, atoms);  // sampled (0 % 5 == 0)
  property.addSample(1, molecules, atoms);  // skipped
  property.addSample(5, molecules, atoms);  // sampled

  EXPECT_TRUE(property.hasData(0));
  EXPECT_FALSE(property.hasData(1));
  EXPECT_EQ(property.dataPerComponent[0].count, 2uz);
  EXPECT_EQ(property.dataPerComponent[1].count, 0uz);

  const std::vector<EndToEndAutoCorrelationFunctionData> data = property.result(0);
  ASSERT_GE(data.size(), 2uz);
  EXPECT_DOUBLE_EQ(data[0].time, 0.0);
  EXPECT_NEAR(data[0].acf, 0.5 * (25.0 + 4.0), 1e-12);  // average over the two molecules
  EXPECT_NEAR(data[1].acf, 0.5 * (25.0 + 4.0), 1e-12);  // static configuration: no decay
  EXPECT_DOUBLE_EQ(data[1].time, 5.0 * 0.001);
  EXPECT_TRUE(property.result(1).empty());

  // the squared radius of gyration (uniform weights): molecule 0 has its center at (4/3, 4/3, 0)
  EXPECT_NEAR(PropertyEndToEndAutoCorrelationFunction::radiusOfGyrationSquared(molecules[0], atoms), 46.0 / 9.0,
              1e-12);
  // the conformations do not change between the two samples (the molecules differ, but each one is frozen): the
  // scalars are static per molecule but differ between the molecules, so the fluctuation functions are flat at 1
  EXPECT_NEAR(data[0].endToEndSquaredACF.value(), 1.0, 1e-9);
  EXPECT_NEAR(data[1].endToEndSquaredACF.value(), 1.0, 1e-9);
  EXPECT_NEAR(data[0].radiusOfGyrationSquaredACF.value(), 1.0, 1e-9);
  EXPECT_NEAR(data[1].radiusOfGyrationSquaredACF.value(), 1.0, 1e-9);
  const EndToEndRelaxationTimes times = property.relaxationTimes(0);
  EXPECT_TRUE(times.endToEndSquared.available);
  EXPECT_NEAR(times.endToEndSquared.mean, 0.5 * (25.0 + 4.0), 1e-12);
  EXPECT_FALSE(times.endToEndSquared.times.oneOverE.has_value());  // never decays
  EXPECT_TRUE(times.endToEndSquared.times.integratedIsLowerBound);

  // a changed number of molecules is refused
  molecules.pop_back();
  EXPECT_THROW(property.addSample(10, molecules, atoms), std::runtime_error);
}

TEST(end_to_end_acf, rigid_tumbling_decays_the_vector_function_but_not_the_conformational_functions)
{
  // a closed "hairpin" of fixed shape that only rotates: <R(0).R(t)> decays with the rotational diffusion while
  // R^2 and Rg^2 are constants, so their fluctuation functions are not available (zero variance)
  const std::vector<double3> hairpin{double3(0.0, 0.0, 0.0),  double3(1.5, 0.0, 0.0), double3(3.0, 0.0, 0.0),
                                     double3(4.5, 0.0, 0.0),  double3(4.5, 1.5, 0.0), double3(4.5, 3.0, 0.0),
                                     double3(3.0, 3.0, 0.0),  double3(1.5, 3.0, 0.0), double3(0.0, 3.0, 0.0)};
  const std::size_t numberOfAtoms = hairpin.size();

  std::vector<Molecule> molecules;
  std::vector<Atom> atoms;
  Molecule molecule(double3(0.0, 0.0, 0.0), simd_quatd(0.0, 0.0, 0.0, 1.0), 1.0, 0, numberOfAtoms);
  molecule.atomIndex = 0;
  molecules.push_back(molecule);
  atoms.resize(numberOfAtoms);

  PropertyEndToEndAutoCorrelationFunction property({1}, {std::array<std::size_t, 2>{0, numberOfAtoms - 1}},
                                                   {acfSettings(10, 1, std::nullopt)}, 1, 1.0);

  // rotational random walk: each sample, rotate the whole molecule by a random small angle about a random axis
  std::mt19937_64 generator{987};
  std::normal_distribution<double> normal{0.0, 1.0};
  std::vector<double3> body = hairpin;
  const double expectedRgSquared = PropertyEndToEndAutoCorrelationFunction::radiusOfGyrationSquared(
      molecule,
      [&]
      {
        std::vector<Atom> a(numberOfAtoms);
        for (std::size_t i = 0; i < numberOfAtoms; ++i) a[i].position = hairpin[i];
        return a;
      }());
  for (std::size_t s = 0; s < 20000; ++s)
  {
    double3 axis(normal(generator), normal(generator), normal(generator));
    axis = axis.normalized();
    const double angle = 0.15 * normal(generator);
    double3 center(0.0, 0.0, 0.0);
    for (const double3 &p : body) center += p;
    center /= static_cast<double>(numberOfAtoms);
    for (double3 &p : body) p = center + rotate(p - center, axis, angle) + double3(0.01, 0.0, 0.0);  // plus drift
    for (std::size_t i = 0; i < numberOfAtoms; ++i) atoms[i].position = body[i];
    property.addSample(s, molecules, atoms);
  }

  const std::vector<EndToEndAutoCorrelationFunctionData> data = property.result(0);
  ASSERT_GT(data.size(), 10uz);
  EXPECT_NEAR(data[0].acf, 9.0, 1e-9);  // |R|^2 = 3^2 of the closed hairpin
  // the vector function has tumbled away
  const EndToEndRelaxationTimes times = property.relaxationTimes(0);
  ASSERT_TRUE(times.oneOverE.has_value());
  EXPECT_LT(times.oneOverE.value(), 500.0);
  EXPECT_LT(data.back().normalized, 0.3);

  // the shape has not changed: no conformational relaxation can be measured, the channels are flagged
  EXPECT_FALSE(times.endToEndSquared.available);
  EXPECT_FALSE(times.radiusOfGyrationSquared.available);
  EXPECT_NEAR(times.endToEndSquared.mean, 9.0, 1e-9);
  EXPECT_NEAR(times.radiusOfGyrationSquared.mean, expectedRgSquared, 1e-9);
  for (const EndToEndAutoCorrelationFunctionData &point : data)
  {
    EXPECT_FALSE(point.endToEndSquaredACF.has_value());
    EXPECT_FALSE(point.radiusOfGyrationSquaredACF.has_value());
  }
}

TEST(end_to_end_acf, conformational_functions_recover_the_relaxation_of_the_scalars)
{
  // Gaussian OU vector with relaxation time tau: the fluctuation function of R^2 = |R|^2 decays as e^{-2t/tau}
  // (Isserlis: <x^2(0) x^2(t)> - <x^2>^2 = 2 <x(0) x(t)>^2 per component), i.e. with tau/2. The "Rg^2" channel is
  // fed an independent scalar OU process with its own relaxation time tauG, whose fluctuation function is
  // e^{-t/tauG}.
  const std::size_t numberOfMolecules = 64;
  const double tau = 20.0;   // ps
  const double tauG = 8.0;   // ps
  const double sigma = 3.0;  // Angstrom per component
  const double dt = 1.0;     // ps between samples

  PropertyEndToEndAutoCorrelationFunction property({numberOfMolecules}, {std::array<std::size_t, 2>{0, 1}},
                                                   {acfSettings(25, 10, std::nullopt)}, numberOfMolecules, 0.1);

  OrnsteinUhlenbeck process{tau, sigma, dt};
  OrnsteinUhlenbeck scalarProcess{tauG, 1.0, dt};
  scalarProcess.generator.seed(777);
  std::vector<double3> vectors(numberOfMolecules);
  std::vector<double3> scalarCarrier(numberOfMolecules);  // only the x component is used
  std::vector<double> radiiOfGyrationSquared(numberOfMolecules);
  for (double3 &r : vectors) r = process.initial();
  for (double3 &g : scalarCarrier) g = scalarProcess.initial();

  const std::size_t numberOfSamples = 40000;
  for (std::size_t s = 0; s < numberOfSamples; ++s)
  {
    for (std::size_t m = 0; m < numberOfMolecules; ++m) radiiOfGyrationSquared[m] = 10.0 + scalarCarrier[m].x;
    property.addSampleVectors(vectors, radiiOfGyrationSquared);
    for (double3 &r : vectors) r = process.step(r);
    for (double3 &g : scalarCarrier) g = scalarProcess.step(g);
  }

  const std::vector<EndToEndAutoCorrelationFunctionData> data = property.result(0);
  ASSERT_GT(data.size(), 10uz);

  // normalized to one at zero lag
  for (const EndToEndAutoCorrelationFunctionData &point : data)
  {
    ASSERT_TRUE(point.endToEndSquaredACF.has_value());
    ASSERT_TRUE(point.radiusOfGyrationSquaredACF.has_value());
  }
  EXPECT_NEAR(data[0].endToEndSquaredACF.value(), 1.0, 1e-9);
  EXPECT_NEAR(data[0].radiusOfGyrationSquaredACF.value(), 1.0, 1e-9);

  // R^2: e^{-2t/tau}
  const std::optional<double> r2AtHalfTau = interpolate(data, 0.5 * tau, endToEndSquaredACF);
  const std::optional<double> r2AtTau = interpolate(data, tau, endToEndSquaredACF);
  ASSERT_TRUE(r2AtHalfTau.has_value());
  ASSERT_TRUE(r2AtTau.has_value());
  EXPECT_NEAR(r2AtHalfTau.value(), std::exp(-1.0), 0.05);
  EXPECT_NEAR(r2AtTau.value(), std::exp(-2.0), 0.05);

  // Rg^2 channel: e^{-t/tauG}
  const std::optional<double> gAtTauG = interpolate(data, tauG, radiusOfGyrationSquaredACF);
  ASSERT_TRUE(gAtTauG.has_value());
  EXPECT_NEAR(gAtTauG.value(), std::exp(-1.0), 0.05);

  const EndToEndRelaxationTimes times = property.relaxationTimes(0);
  ASSERT_TRUE(times.endToEndSquared.available);
  EXPECT_NEAR(times.endToEndSquared.mean, 3.0 * sigma * sigma, 0.05 * 3.0 * sigma * sigma);
  // Var(R^2) = 3 * 2 sigma^4
  EXPECT_NEAR(times.endToEndSquared.variance, 6.0 * std::pow(sigma, 4), 0.10 * 6.0 * std::pow(sigma, 4));
  ASSERT_TRUE(times.endToEndSquared.times.oneOverE.has_value());
  EXPECT_NEAR(times.endToEndSquared.times.oneOverE.value(), 0.5 * tau, 0.12 * 0.5 * tau);
  ASSERT_TRUE(times.endToEndSquared.times.exponentialFit.has_value());
  EXPECT_NEAR(times.endToEndSquared.times.exponentialFit.value(), 0.5 * tau, 0.15 * 0.5 * tau);

  ASSERT_TRUE(times.radiusOfGyrationSquared.available);
  EXPECT_NEAR(times.radiusOfGyrationSquared.mean, 10.0, 0.1);
  EXPECT_NEAR(times.radiusOfGyrationSquared.variance, 1.0, 0.1);
  ASSERT_TRUE(times.radiusOfGyrationSquared.times.oneOverE.has_value());
  EXPECT_NEAR(times.radiusOfGyrationSquared.times.oneOverE.value(), tauG, 0.12 * tauG);

  // the vector function is unaffected by the scalars: tau as before
  ASSERT_TRUE(times.oneOverE.has_value());
  EXPECT_NEAR(times.oneOverE.value(), tau, 0.12 * tau);

  // without Rg^2 input the Rg^2 channel is not available
  PropertyEndToEndAutoCorrelationFunction withoutRg({numberOfMolecules}, {std::array<std::size_t, 2>{0, 1}},
                                                    {acfSettings(25, 10, std::nullopt)}, numberOfMolecules, 0.1);
  for (std::size_t s = 0; s < 100; ++s)
  {
    withoutRg.addSampleVectors(vectors);
    for (double3 &r : vectors) r = process.step(r);
  }
  const std::vector<EndToEndAutoCorrelationFunctionData> noRg = withoutRg.result(0);
  ASSERT_FALSE(noRg.empty());
  EXPECT_FALSE(noRg[0].radiusOfGyrationSquaredACF.has_value());
  EXPECT_TRUE(noRg[0].endToEndSquaredACF.has_value());
  EXPECT_FALSE(withoutRg.relaxationTimes(0).radiusOfGyrationSquared.available);
  EXPECT_TRUE(withoutRg.relaxationTimes(0).endToEndSquared.available);
}

TEST(end_to_end_acf, binary_archive_round_trip_continues_the_accumulation)
{
  const std::size_t numberOfMolecules = 8;
  PropertyEndToEndAutoCorrelationFunction original({numberOfMolecules, 3},
                                                   {std::array<std::size_t, 2>{0, 1}, std::nullopt},
                                                   {acfSettings(5, 10, 5000), std::nullopt}, numberOfMolecules + 3, 0.002);

  // a trajectory of 137 + 61 samples, generated once
  OrnsteinUhlenbeck process{5.0, 2.0, 0.02};
  std::vector<std::vector<double3>> trajectory;
  {
    std::vector<double3> vectors(numberOfMolecules + 3);
    for (double3 &r : vectors) r = process.initial();
    for (std::size_t s = 0; s < 137 + 61; ++s)
    {
      trajectory.push_back(vectors);
      for (double3 &r : vectors) r = process.step(r);
    }
  }
  // a stand-in Rg^2 per molecule and sample, derived from the trajectory
  auto radiiOfGyrationSquared = [&](std::size_t s)
  {
    std::vector<double> rg(numberOfMolecules + 3);
    for (std::size_t m = 0; m < rg.size(); ++m) rg[m] = 0.25 * double3::dot(trajectory[s][m], trajectory[s][m]) + 1.0;
    return rg;
  };
  for (std::size_t s = 0; s < 137; ++s) original.addSampleVectors(trajectory[s], radiiOfGyrationSquared(s));

  const std::filesystem::path path = std::filesystem::temp_directory_path() / "raspa3_end_to_end_acf_round_trip.bin";
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << original;
  }
  PropertyEndToEndAutoCorrelationFunction restored;
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);

  EXPECT_EQ(restored.dataPerComponent[0].count, original.dataPerComponent[0].count);
  EXPECT_EQ(restored.dataPerComponent[0].numberOfBlocks, original.dataPerComponent[0].numberOfBlocks);
  EXPECT_EQ(restored.sampleEvery(0), 10uz);
  EXPECT_EQ(restored.writeEvery(0), std::optional<std::size_t>{5000});
  EXPECT_EQ(restored.numberOfBlockElements(0), 5uz);
  EXPECT_FALSE(restored.isSampled(1));
  EXPECT_FALSE(restored.hasData(1));

  // continuing both copies with the same samples gives bit-identical functions
  for (std::size_t s = 137; s < 137 + 61; ++s)
  {
    original.addSampleVectors(trajectory[s], radiiOfGyrationSquared(s));
    restored.addSampleVectors(trajectory[s], radiiOfGyrationSquared(s));
  }

  const std::vector<EndToEndAutoCorrelationFunctionData> a = original.result(0);
  const std::vector<EndToEndAutoCorrelationFunctionData> b = restored.result(0);
  ASSERT_EQ(a.size(), b.size());
  for (std::size_t i = 0; i < a.size(); ++i)
  {
    EXPECT_DOUBLE_EQ(a[i].time, b[i].time);
    EXPECT_DOUBLE_EQ(a[i].acf, b[i].acf);
    EXPECT_DOUBLE_EQ(a[i].numberOfSamples, b[i].numberOfSamples);
    ASSERT_TRUE(a[i].endToEndSquaredACF.has_value());
    ASSERT_TRUE(a[i].radiusOfGyrationSquaredACF.has_value());
    EXPECT_DOUBLE_EQ(a[i].endToEndSquaredACF.value(), b[i].endToEndSquaredACF.value());
    EXPECT_DOUBLE_EQ(a[i].radiusOfGyrationSquaredACF.value(), b[i].radiusOfGyrationSquaredACF.value());
  }
  const EndToEndRelaxationTimes ta = original.relaxationTimes(0), tb = restored.relaxationTimes(0);
  EXPECT_DOUBLE_EQ(ta.endToEndSquared.mean, tb.endToEndSquared.mean);
  EXPECT_DOUBLE_EQ(ta.endToEndSquared.variance, tb.endToEndSquared.variance);
  EXPECT_DOUBLE_EQ(ta.radiusOfGyrationSquared.mean, tb.radiusOfGyrationSquared.mean);
  EXPECT_DOUBLE_EQ(ta.radiusOfGyrationSquared.variance, tb.radiusOfGyrationSquared.variance);
}
