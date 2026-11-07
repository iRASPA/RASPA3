#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import double3;
import transition_matrix;

// The diagnostic columns of the transition-matrix statistics file: the running sums of the forward and reverse
// transition log-probabilities along the detailed-balance walk, whose difference is ln Pi. They start at zero
// at the lowest visited state of the current rebuild, follow ln Pi flat across an unsupported link and outside
// the visited range, and are reset together with the statistics they are derived from.

namespace
{

class ScopedCurrentPath
{
 public:
  explicit ScopedCurrentPath(const std::filesystem::path& path) : original_(std::filesystem::current_path())
  {
    std::filesystem::create_directories(path);
    std::filesystem::current_path(path);
  }

  ~ScopedCurrentPath()
  {
    std::error_code ignored;
    const std::filesystem::path used = std::filesystem::current_path();
    std::filesystem::current_path(original_);
    std::filesystem::remove_all(used, ignored);
  }

 private:
  std::filesystem::path original_;
};

TransitionMatrix makeMatrix()
{
  TransitionMatrix matrix;
  matrix.doTMMC = true;
  matrix.useBias = true;
  matrix.useTMBias = true;
  matrix.minMacrostate = 0;
  matrix.maxMacrostate = 4;
  matrix.initialize();
  // states 1, 2, 3 visited; the links 1 -> 2 and 2 -> 3 have support in both directions
  matrix.cmatrix[1] = double3(0.0, 2.0, 3.0);
  matrix.cmatrix[2] = double3(1.0, 3.0, 2.0);
  matrix.cmatrix[3] = double3(2.0, 4.0, 0.0);
  matrix.histogram[1] = matrix.histogram[2] = matrix.histogram[3] = 1uz;
  return matrix;
}

double total(const double3& row) { return row.x + row.y + row.z; }

}  // namespace

TEST(TRANSITION_MATRIX_BOOKKEEPING, forward_and_reverse_sums_start_at_the_lowest_visited_state_and_give_lnpi)
{
  ScopedCurrentPath scratch(std::filesystem::temp_directory_path() /
                            std::format("raspa3_test_{}_transition_matrix", test_support::uniqueSuffix()));
  TransitionMatrix matrix = makeMatrix();
  matrix.recomputeLnPiAndBias();
  EXPECT_EQ(matrix.numberOfUpdates, 1uz);

  EXPECT_DOUBLE_EQ(matrix.forward_lnpi[1], 0.0);
  EXPECT_DOUBLE_EQ(matrix.reverse_lnpi[1], 0.0);
  for (std::size_t i = 1; i < 3; ++i)
  {
    const double forward = std::log(matrix.cmatrix[i].z / total(matrix.cmatrix[i]));
    const double reverse = std::log(matrix.cmatrix[i + 1].x / total(matrix.cmatrix[i + 1]));
    EXPECT_NEAR(matrix.forward_lnpi[i + 1] - matrix.forward_lnpi[i], forward, 1e-12) << "link " << i;
    EXPECT_NEAR(matrix.reverse_lnpi[i + 1] - matrix.reverse_lnpi[i], reverse, 1e-12) << "link " << i;
  }
  // the difference of the sums is ln Pi (up to its normalization)
  for (std::size_t i = 1; i <= 3; ++i)
  {
    EXPECT_NEAR(matrix.lnpi[i] - matrix.lnpi[1], (matrix.forward_lnpi[i] - matrix.reverse_lnpi[i]), 1e-12)
        << "state " << i;
  }
  // flat outside the visited range, like ln Pi
  EXPECT_DOUBLE_EQ(matrix.forward_lnpi[0], matrix.forward_lnpi[1]);
  EXPECT_DOUBLE_EQ(matrix.reverse_lnpi[0], matrix.reverse_lnpi[1]);
  EXPECT_DOUBLE_EQ(matrix.forward_lnpi[4], matrix.forward_lnpi[3]);
  EXPECT_DOUBLE_EQ(matrix.reverse_lnpi[4], matrix.reverse_lnpi[3]);

  // a later rebuild with another lowest visited state starts its sums afresh there
  matrix.histogram[1] = 0uz;
  matrix.recomputeLnPiAndBias();
  EXPECT_EQ(matrix.numberOfUpdates, 2uz);
  EXPECT_DOUBLE_EQ(matrix.forward_lnpi[2], 0.0);
  EXPECT_DOUBLE_EQ(matrix.reverse_lnpi[2], 0.0);
  EXPECT_DOUBLE_EQ(matrix.forward_lnpi[1], 0.0);
  EXPECT_DOUBLE_EQ(matrix.forward_lnpi[0], 0.0);
  EXPECT_NEAR(matrix.forward_lnpi[3], std::log(matrix.cmatrix[2].z / total(matrix.cmatrix[2])), 1e-12);
}

TEST(TRANSITION_MATRIX_BOOKKEEPING, clearing_the_statistics_resets_the_derived_columns_and_the_update_count)
{
  ScopedCurrentPath scratch(std::filesystem::temp_directory_path() /
                            std::format("raspa3_test_{}_transition_matrix_clear", test_support::uniqueSuffix()));
  TransitionMatrix matrix = makeMatrix();
  matrix.recomputeLnPiAndBias();
  ASSERT_NE(matrix.forward_lnpi[3], 0.0);
  ASSERT_EQ(matrix.numberOfUpdates, 1uz);

  // production start: the statistics go, the bias stays, the rebuild count starts over
  TransitionMatrix kept = matrix;
  kept.clearStatisticsKeepBias();
  EXPECT_EQ(kept.numberOfUpdates, 0uz);
  EXPECT_EQ(kept.bias, matrix.bias);

  // after initialization: everything derived from the statistics goes with them
  matrix.rezeroAfterInitialization = true;
  matrix.clearCMatrix();
  EXPECT_EQ(matrix.numberOfUpdates, 0uz);
  for (std::size_t i = 0; i < matrix.forward_lnpi.size(); ++i)
  {
    EXPECT_DOUBLE_EQ(matrix.forward_lnpi[i], 0.0) << "state " << i;
    EXPECT_DOUBLE_EQ(matrix.reverse_lnpi[i], 0.0) << "state " << i;
    EXPECT_DOUBLE_EQ(matrix.lnpi[i], 0.0) << "state " << i;
  }
}
