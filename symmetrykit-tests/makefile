CXX=clang++
GTEST_DIR=/usr/local/include/googletest/googletest/
GTEST_INCLUDE=/usr/local

CPPFLAGS += -isystem $(GTEST_INCLUDE)/include

MODULEPATHS = -fprebuilt-module-path=../mathkit -fprebuilt-module-path=../foundationkit -fprebuilt-module-path=../symmetrykit

CXXFLAGS += -std=c++2b -fmodules -fbuiltin-module-map -stdlib=libc++ -O3 -march=native -ffast-math -fopenmp -Wall -Wextra $(MODULEPATHS)


TESTS = test_symmetrykit

GTEST_HEADERS = $(GTEST_INCLUDE)/include/gtest/*.h \
                $(GTEST_INCLUDE)/include/gtest/internal/*.h

all : $(TESTS)

clean :
	rm -f $(TESTS) gtest.a gtest_main.a *.o

gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_INCLUDE) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_INCLUDE) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Builds a sample test.  A test should link with either gtest.a or
# gtest_main.a, depending on whether it defines its own main()
# function.

test_symmetrykit.o : main.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c main.cpp

transformation_matrix.o: $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c transformation_matrix.cpp

find_conventional_cell.o: $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c find_conventional_cell.cpp

find_smallest_primitive_cell.o: $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c find_smallest_primitive_cell.cpp

find_smallest_primitive_cell_no_partial_occupancies.o: $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c find_smallest_primitive_cell_no_partial_occupancies.cpp

find_pointgroup.o: $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c find_pointgroup.cpp

find_pointgroup_no_partial_occupancies.o: $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c find_pointgroup_no_partial_occupancies.cpp

find_spacegroup.o: $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c find_spacegroup.cpp

find_spacegroup_no_partial_occupancies.o: $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c find_spacegroup_no_partial_occupancies.cpp

test_symmetrykit : transformation_matrix.o find_conventional_cell.o find_conventional_cell.o find_smallest_primitive_cell.o find_smallest_primitive_cell_no_partial_occupancies.o find_pointgroup.o find_pointgroup_no_partial_occupancies.o find_pointgroup_no_partial_occupancies.o find_spacegroup.o find_spacegroup_no_partial_occupancies.o  main.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -L../mathkit -L../foundationkit -L../symmetrykit -lsymmetrykit -lfoundationkit -lmathkit -llapack -lpthread $^ -o $@
