CPP    := $(wildcard *.cpp)
PCM    := $(CPP:.cpp=.pcm)
OBJ    := $(CPP:.cpp=.o)
OBJPCM := $(CPP:.cpp=.pcm.o)

.PHONY: default
default: libfoundationkit.a libfoundationkit.so;

.SUFFIXES:
.SUFFIXES: .pcm .ixx *.cpp

.SECONDARY:

$(OBJ): %.o: %.cpp %.pcm
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(PCM): %.pcm: %.ixx
	$(CXX) $(CXXFLAGS) -x c++-module --precompile -o $@ $<
	$(CXX) $(CXXFLAGS) -x c++-module $< -c -o $@.o	

scanner.pcm: characterset.pcm
threadpool.pcm: threading.pcm

libfoundationkit.a: $(OBJ)
	ar -rs libfoundationkit.a $(OBJPCM) $(OBJ)
	ranlib -c libfoundationkit.a

libfoundationkit.so: $(OBJ)
	$(CXX) $(CXXFLAGS) -shared -undefined dynamic_lookup `python3 -m pybind11 --includes` `python3-config --ldflags` $(OBJPCM) $(OBJ)  -o $@

clean:
	rm -f *.pcm *.o *.a
