CPP    := $(wildcard *.cpp)
PCM    := $(CPP:.cpp=.pcm)
OBJ    := $(CPP:.cpp=.o)
OBJPCM := $(CPP:.cpp=.pcm.o)

$(info    CPP is $(CPP))
$(info    IXX is $(IXX))
$(info    PCM is $(PCM))
$(info    OBJ is $(OBJ))

.PHONY: default
default: libfoundationkit;

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

libfoundationkit: $(OBJ)
	ar -rs libfoundationkit.a $(OBJPCM) $(OBJ)
	ranlib -c libfoundationkit.a

clean:
	rm -f *.pcm *.o *.a
