CPP    := $(wildcard *.cpp)
PCM    := $(CPP:.cpp=.pcm)
OBJ    := $(CPP:.cpp=.o)
OBJPCM := $(CPP:.cpp=.pcm.o)

.PHONY: default
default: libsymmetrykit.a libsymmetrykit.so;

.SUFFIXES:
.SUFFIXES: .pcm .ixx *.cpp

.SECONDARY:

$(OBJ): %.o: %.cpp %.pcm
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(PCM): %.pcm: %.ixx
	$(CXX) $(CXXFLAGS) -x c++-module --precompile -o $@ $<
	$(CXX) $(CXXFLAGS) -x c++-module $< -c -o $@.o	

# dependencies
skasymmetricatom.pcm: skatomcopy.pcm skelement.pcm
skcell.pcm: skboundingbox.pcm
skstructure.pcm: skcell.pcm skasymmetricatom.pcm
skparser.pcm: skstructure.pcm
skposcarparser.pcm: skparser.pcm 
skposcarlegacyparser.pcm: skparser.pcm skstructure.pcm skcell.pcm skasymmetricatom.pcm
sktransformationmatrix.pcm: skrotationmatrix.pcm
skintegerchangeofbasis.pcm: sktransformationmatrix.pcm
skseitzintegermatrix.pcm: skrotationmatrix.pcm skonethirdseitzmatrix.pcm
skintegersymmetryoperationset.pcm: skdefinitions.pcm sksymmetrycell.pcm skseitzintegermatrix.pcm
skseitzmatrix.pcm: skrotationmatrix.pcm skonethirdseitzmatrix.pcm
skpointsymmetryset.pcm: skrotationmatrix.pcm
sksymmetryoperationset.pcm: skdefinitions.pcm skpointgroup.pcm skrotationmatrix.pcm sktransformationmatrix.pcm skseitzmatrix.pcm sksymmetrycell.pcm
skpointgroup.pcm: skdefinitions.pcm skrotationaloccurancetable.pcm skrotationmatrix.pcm skpointsymmetryset.pcm sktransformationmatrix.pcm
skspacegroupsetting.pcm: skdefinitions.pcm skasymmetricunit.pcm sktransformationmatrix.pcm skintegersymmetryoperationset.pcm skseitzintegermatrix.pcm skpointgroup.pcm
skrotationalchangeofbasis.pcm: skrotationmatrix.pcm skseitzintegermatrix.pcm
skspacegroupdatabase.pcm: skspacegroupsetting.pcm
sksymmetrycell.pcm: skdefinitions.pcm skrotationmatrix.pcm sktransformationmatrix.pcm skpointsymmetryset.pcm skpointgroup.pcm
skspacegroup.pcm: skseitzmatrix.pcm skpointgroup.pcm skpointsymmetryset.pcm sksymmetryoperationset.pcm skspacegroupsetting.pcm skrotationalchangeofbasis.pcm skspacegroupdatabase.pcm sksymmetrycell.pcm

libsymmetrykit.a: $(OBJ)
	ar -rs libsymmetrykit.a $(OBJPCM) $(OBJ)
	ranlib -c libsymmetrykit.a

libsymmetrykit.so: $(OBJ)
	$(CXX) $(CXXFLAGS) -shared -undefined dynamic_lookup `python3 -m pybind11 --includes` `python3-config --ldflags` $(OBJPCM) $(OBJ) -o $@

clean:
	rm -f *.pcm *.o *.a
