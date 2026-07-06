#include <gtest/gtest.h>

import std;

import janaf;
import nasa_polynomials;

// The janaf and nasa_polynomials modules use the same energy reference (dissociated atoms in
// their electronic ground state), so ln(q/V) of species present in both databases must agree.
// Small deviations remain because the databases are independent: Burcat incorporates updated
// (ATcT) enthalpies of formation (about 1 kJ/mol difference for NO), and the janaf module
// interpolates linearly between tabulated points.
TEST(nasa_polynomials, consistent_with_janaf_tables)
{
  const double tolerance = 0.15;  // in ln(q/V)
  const double energyTolerance = 1.5e3 / 8.314;  // deviations equivalent to at most 1.5 kJ/mol, in K

  for (std::string_view name : {"H2", "N2", "O2", "H2O", "NH3", "CO", "CO2", "CH4", "SO2", "H2S", "NO", "HCl"})
  {
    ASSERT_TRUE(NASAPolynomials::contains(name)) << name;
    ASSERT_TRUE(JANAF::contains(name)) << name;
    for (double temperature : {300.0, 600.0, 1000.0, 2000.0})
    {
      EXPECT_NEAR(NASAPolynomials::logPartitionFunction(name, temperature),
                  JANAF::logPartitionFunction(name, temperature),
                  std::max(tolerance, energyTolerance / temperature))
          << name << " at " << temperature << " K";
    }
  }
}

// Regression values for the propene-metathesis species of examples/advanced/8_mc_reaction_mc,
// which are not available in the JANAF tables.
TEST(nasa_polynomials, propene_metathesis_species)
{
  EXPECT_NEAR(NASAPolynomials::logPartitionFunction("propene", 450.0), 925.378, 0.01);
  EXPECT_NEAR(NASAPolynomials::logPartitionFunction("ethene", 450.0), 607.764, 0.01);
  EXPECT_NEAR(NASAPolynomials::logPartitionFunction("cis-2-butene", 450.0), 1240.033, 0.01);
  EXPECT_NEAR(NASAPolynomials::logPartitionFunction("trans-2-butene", 450.0), 1240.524, 0.01);
  EXPECT_NEAR(NASAPolynomials::logPartitionFunction("1-butene", 450.0), 1238.612, 0.01);
  EXPECT_NEAR(NASAPolynomials::logPartitionFunction("isobutene", 450.0), 1241.162, 0.01);
}

TEST(nasa_polynomials, name_matching)
{
  // exact database names
  EXPECT_TRUE(NASAPolynomials::contains("C3H6 propylene"));
  EXPECT_TRUE(NASAPolynomials::contains("C4H8,tr2-butene"));
  EXPECT_TRUE(NASAPolynomials::contains("C6H6 BENZENE"));

  // aliases, case-insensitive and ignoring separators
  EXPECT_TRUE(NASAPolynomials::contains("Trans-2-Butene"));
  EXPECT_TRUE(NASAPolynomials::contains("carbon dioxide"));
  EXPECT_TRUE(NASAPolynomials::contains("n-hexane"));
  EXPECT_TRUE(NASAPolynomials::contains("methanol"));
  EXPECT_TRUE(NASAPolynomials::contains("argon"));
  EXPECT_FALSE(NASAPolynomials::contains("unobtainium"));

  EXPECT_EQ(NASAPolynomials::logPartitionFunction("propene", 450.0),
            NASAPolynomials::logPartitionFunction("C3H6 propylene", 450.0));

  EXPECT_GT(NASAPolynomials::availableSpecies().size(), std::size_t{1900});

  // unique substring match on the descriptive part of database names
  EXPECT_TRUE(NASAPolynomials::contains("naphthalene"));
  EXPECT_TRUE(NASAPolynomials::contains("styrene"));
  EXPECT_TRUE(NASAPolynomials::contains("1-hexene"));
  EXPECT_FALSE(NASAPolynomials::contains("pentene"));  // ambiguous: 1-pentene, 2-pentene, ...
  EXPECT_THROW(std::ignore = NASAPolynomials::logPartitionFunction("pentene", 450.0), std::runtime_error);

  EXPECT_THROW(std::ignore = NASAPolynomials::logPartitionFunction("unobtainium", 450.0), std::runtime_error);
  EXPECT_THROW(std::ignore = NASAPolynomials::logPartitionFunction("propene", 100.0), std::runtime_error);
  EXPECT_THROW(std::ignore = NASAPolynomials::logPartitionFunction("propene", 7000.0), std::runtime_error);
}

// m-Xylene is a local supplement (absent from the Burcat database) and p-xylene is a local
// replacement (the Burcat record, based on Draeger & Scott 1981, is an entropy outlier flagged
// in the NIST WebBook: ~15 J/mol/K below the recommended value). Both use low ranges fitted to
// the recommended ideal-gas Cp (Draeger 1985 / Chao 1984/1986), anchored at experimental
// ΔfH°(298.15 K) and S°(298.15 K), with the Burcat high-range Cp coefficients.
TEST(nasa_polynomials, xylene_supplements)
{
  EXPECT_TRUE(NASAPolynomials::contains("m-xylene"));
  EXPECT_TRUE(NASAPolynomials::contains("p-xylene"));
  EXPECT_TRUE(NASAPolynomials::contains("1,3-dimethylbenzene"));  // substring of the database name
  EXPECT_TRUE(NASAPolynomials::contains("1,4-dimethylbenzene"));

  EXPECT_NEAR(NASAPolynomials::logPartitionFunction("m-xylene", 450.0), 2113.158, 0.01);
  EXPECT_NEAR(NASAPolynomials::logPartitionFunction("m-xylene", 600.0), 1595.393, 0.01);
  EXPECT_NEAR(NASAPolynomials::logPartitionFunction("p-xylene", 450.0), 2112.337, 0.01);
  EXPECT_NEAR(NASAPolynomials::logPartitionFunction("p-xylene", 600.0), 1594.624, 0.01);

  // G°(298.15) = ΔfH° - T S°: 17.2 - 298.15*0.35769 kJ/mol (meta), 17.9 - 298.15*0.3524 (para)
  EXPECT_NEAR(NASAPolynomials::standardGibbsEnergy("m-xylene", 298.15), -89.445, 0.01);
  EXPECT_NEAR(NASAPolynomials::standardGibbsEnergy("p-xylene", 298.15), -87.168, 0.01);

  // mutual consistency of the three isomers at 450 K, e.g. G°(m) - G°(o) is about
  // (17.2 - 19.0) - 450*(357.69 - 353.6)/1000 = -3.6 kJ/mol. The meta isomer is the most
  // stable, largely because its entropy is higher by ~R ln 2 (symmetry number 2 vs 4).
  const double gibbsOrtho = NASAPolynomials::standardGibbsEnergy("o-xylene", 450.0);
  const double gibbsMeta = NASAPolynomials::standardGibbsEnergy("m-xylene", 450.0);
  const double gibbsPara = NASAPolynomials::standardGibbsEnergy("p-xylene", 450.0);
  EXPECT_LT(gibbsMeta, gibbsOrtho);
  EXPECT_LT(gibbsMeta, gibbsPara);
  EXPECT_NEAR(gibbsMeta - gibbsOrtho, -3.6, 2.0);
  EXPECT_NEAR(gibbsMeta - gibbsPara, -3.1, 2.0);
}

// The standard Gibbs energy of the reference elements is -T*S°(T) + [H°(T) - H°(298.15 K)];
// spot-check N2 against the JANAF table (G° = -T*S° + [H-H298]: at 1000 K,
// S° = 228.170 J/mol/K and H-H298 = 21.463 kJ/mol -> G° = -206.707 kJ/mol).
TEST(nasa_polynomials, standard_gibbs_energy)
{
  EXPECT_NEAR(NASAPolynomials::standardGibbsEnergy("N2", 1000.0), -206.707, 0.05);
  EXPECT_NEAR(NASAPolynomials::standardGibbsEnergy("H2", 298.15), -298.15 * 130.680 / 1000.0, 0.05);
}
