CXX=clang++
GTEST_DIR=/usr/local/include/googletest/googletest/
GTEST_INCLUDE=/usr/local

PYTHONLDFLAGS = $(shell python3-config --ldflags)

CPPFLAGS += -isystem $(GTEST_INCLUDE)/include

MODULEPATHS = -fprebuilt-module-path=../../src -fprebuilt-module-path=../../src/foundationkit -fprebuilt-module-path=../../src/mathkit -fprebuilt-module-path=../../src/symmetrykit -fprebuilt-module-path=../../src/raspakit -fprebuilt-module-path=../../src/raspakit/interactions -fprebuilt-module-path=../../src/raspakit/mc_moves -fprebuilt-module-path=../../src/raspakit/potentials -fprebuilt-module-path=../../src/raspakit/properties -fprebuilt-module-path=../../src/raspakit/cbmc

CXXFLAGS += -std=c++2b -fmodules -fbuiltin-module-map -stdlib=libc++ -O3 -march=native -ffast-math -fopenmp -Wall -Wextra $(MODULEPATHS)

TESTS = test_foundationkit

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

test_foundationkit.o : main.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c main.cpp

archive.o: $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c archive.cpp

test_foundationkit : archive.o main.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PYTHONLDFLAGS) ../../src/raspakit.so  -llapack -lpthread -lpython3.12 $^ -o $@
