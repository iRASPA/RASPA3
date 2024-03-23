CPP    := $(wildcard *.cpp)
PCM    := $(CPP:.cpp=.pcm)
OBJ    := $(CPP:.cpp=.o)
OBJPCM := $(CPP:.cpp=.pcm.o)

.PHONY: default
default: libmathkit.a libmathkit.so;

.SUFFIXES:
.SUFFIXES: .pcm .ixx *.cpp

.SECONDARY:

$(OBJ): %.o: %.cpp %.pcm
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(PCM): %.pcm: %.ixx
	$(CXX) $(CXXFLAGS) -x c++-module --precompile -o $@ $<
	$(CXX) $(CXXFLAGS) -x c++-module $< -c -o $@.o	

# dependencies
int3.pcm: hashcombine.pcm ring.pcm
float3.pcm: int3.pcm
double3.pcm: hashcombine.pcm int3.pcm bool3.pcm
randomnumbers.pcm: double3.pcm double3x3.pcm simd_quatd.pcm
float4.pcm: float3.pcm double3.pcm double4.pcm
int3x3.pcm: hashcombine.pcm int3.pcm double3.pcm
double3x3.pcm: double3.pcm int3x3.pcm simd_quatd.pcm
double4x3.pcm: double3.pcm double4.pcm
double4x4.pcm: double3.pcm double4.pcm simd_quatd.pcm double3x3.pcm
simd_quatd.pcm: double3.pcm
float3x3.pcm: float3.pcm simd_quatd.pcm
float3x4.pcm: float3.pcm float4.pcm
float4x4.pcm: float3.pcm float4.pcm double3x3.pcm double4x4.pcm simd_quatd.pcm float3x3.pcm
ringmatrix.pcm: int3x3.pcm
matrix.pcm: ringmatrix.pcm
symmetric_matrix_layout.pcm: mdspan.pcm


libmathkit.a: $(OBJ)
	ar -rs libmathkit.a $(OBJPCM) $(OBJ)
	ranlib -c libmathkit.a

libmathkit.so: $(OBJ)
	$(CXX) $(CXXFLAGS) -shared $(LINK_FLAGS) `python3 -m pybind11 --includes` `python3-config --ldflags` $(OBJPCM) $(OBJ) -o $@

clean:
	rm -f *.pcm *.o *.a *.so
