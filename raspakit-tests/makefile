CXX=clang++
GTEST_DIR=/usr/local/include/googletest/googletest/
GTEST_INCLUDE=/usr/local
#GTEST_DIR=/opt/homebrew/include/googletest/googletest/
#GTEST_INCLUDE=/opt/homebrew

CPPFLAGS += -isystem $(GTEST_INCLUDE)/include

MODULEPATHS = -fprebuilt-module-path=../mathkit -fprebuilt-module-path=../foundationkit -fprebuilt-module-path=../symmetrykit -fprebuilt-module-path=../raspakit

CXXFLAGS += -std=c++2b -fmodules -fbuiltin-module-map -stdlib=libc++ -O3 -march=native -ffast-math -fopenmp -Wall -Wextra $(MODULEPATHS)


TESTS = test_raspakit

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

test_raspakit.o : main.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c main.cpp

static_energy.o: static_energy.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c static_energy.cpp

gradients.o: gradients.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c gradients.cpp

pbc.o: pbc.c $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c pbc.cpp

pressure.o: pressure.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c pressure.cpp

ewald.o: ewald.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c ewald.cpp

dudlambda.o: dudlambda.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c dudlambda.cpp

spacegroup.o: spacegroup.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c spacegroup.cpp

insertion_deletion.o: insertion_deletion.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c insertion_deletion.cpp

energy_decomposition.o: energy_decomposition.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c energy_decomposition.cpp

test_raspakit : static_energy.o gradients.o pressure.o ewald.o dudlambda.o spacegroup.o cbmc_interactions.o insertion_deletion.o energy_decomposition.o main.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -L../mathkit -L../foundationkit -L../symmetrykit -L../raspakit  -lraspakit -lsymmetrykit -lfoundationkit -lmathkit -llapack -lpthread $^ -o $@
