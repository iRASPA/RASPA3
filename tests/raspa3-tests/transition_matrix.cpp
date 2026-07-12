#include <gtest/gtest.h>

import std;

import double3;
import transition_matrix;

TEST(TRANSITION_MATRIX, bounds_are_checked_before_indexing)
{
  TransitionMatrix matrix;
  matrix.doTMMC = true;
  matrix.useBias = true;
  matrix.useTMBias = true;
  matrix.rejectOutOfBound = true;
  matrix.minMacrostate = 2;
  matrix.maxMacrostate = 4;
  matrix.initialize();

  matrix.bias = {0.0, std::log(2.0), std::log(8.0)};

  EXPECT_DOUBLE_EQ(matrix.biasFactor(3, 2), 2.0);
  EXPECT_DOUBLE_EQ(matrix.biasFactor(1, 2), 0.0);
  EXPECT_DOUBLE_EQ(matrix.biasFactor(5, 4), 0.0);

  matrix.updateMatrix(double3(1.0, 0.0, 0.0), 1);
  matrix.updateMatrix(double3(0.0, 0.0, 1.0), 5);
  for (const double3& row : matrix.cmatrix)
  {
    EXPECT_DOUBLE_EQ(row.x, 0.0);
    EXPECT_DOUBLE_EQ(row.y, 0.0);
    EXPECT_DOUBLE_EQ(row.z, 0.0);
  }
}

TEST(TRANSITION_MATRIX, collection_probabilities_are_clamped_and_reset)
{
  TransitionMatrix matrix;
  matrix.doTMMC = true;
  matrix.rezeroAfterInitialization = true;
  matrix.minMacrostate = 0;
  matrix.maxMacrostate = 2;
  matrix.initialize();

  matrix.updateMatrix(double3(-1.0, 0.25, 2.0), 1);
  EXPECT_DOUBLE_EQ(matrix.cmatrix[1].x, 0.0);
  EXPECT_DOUBLE_EQ(matrix.cmatrix[1].y, 0.25);
  EXPECT_DOUBLE_EQ(matrix.cmatrix[1].z, 1.0);

  matrix.updateHistogram(1);
  matrix.numberOfSteps = 7;
  matrix.clearCMatrix();

  EXPECT_EQ(matrix.numberOfSteps, 0);
  EXPECT_EQ(matrix.histogram, std::vector<std::size_t>({0, 0, 0}));
  for (const double3& row : matrix.cmatrix)
  {
    EXPECT_DOUBLE_EQ(row.x, 0.0);
    EXPECT_DOUBLE_EQ(row.y, 0.0);
    EXPECT_DOUBLE_EQ(row.z, 0.0);
  }
}

TEST(TRANSITION_MATRIX, rejects_an_inverted_macrostate_range)
{
  TransitionMatrix matrix;
  matrix.doTMMC = true;
  matrix.minMacrostate = 3;
  matrix.maxMacrostate = 2;

  EXPECT_THROW(matrix.initialize(), std::invalid_argument);
}
