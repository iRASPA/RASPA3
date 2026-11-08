#include <gtest/gtest.h>

import std;

import symmetric_matrix_layout;

TEST(symmetric_matrix_layout, Test_dynamic_extents)
{
  // the dynamic-extent constructor checks squareness with assert
  std::dextents<std::size_t, 2> extents(4, 4);
  layout_symmetric_matrix::mapping<std::dextents<std::size_t, 2>> mapping(extents);

  EXPECT_EQ(mapping.required_span_size(), 10uz);
  EXPECT_EQ(mapping(0, 0), 0uz);
  EXPECT_EQ(mapping(3, 3), 9uz);
  EXPECT_EQ(mapping(1, 3), mapping(3, 1));

  std::vector<double> storage(mapping.required_span_size());
  std::mdspan<double, std::dextents<std::size_t, 2>, layout_symmetric_matrix> matrix(storage.data(), extents);
  matrix[2, 0] = 5.0;
  EXPECT_EQ((matrix[0, 2]), 5.0);
}

TEST(symmetric_matrix_layout, Test_static_extents)
{
  layout_symmetric_matrix::mapping<std::extents<std::size_t, 3, 3>> mapping(std::extents<std::size_t, 3, 3>{});

  EXPECT_EQ(mapping.required_span_size(), 6uz);
  EXPECT_EQ(mapping(2, 1), mapping(1, 2));
}
