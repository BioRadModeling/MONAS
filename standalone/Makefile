CXX = g++
CXXFLAGS = -std=c++11 -pthread

monas: main.o TsGetSurvivalRBEQualityFactor.o TsGSM2parallel.o TsSpecificEnergy.o TsLinealEnergy.o
	$(CXX) $(CXXFLAGS) $^ -o $@

TsGetSurvival.o: TsGetSurvivalRBEQualityFactor.cc TsGetSurvival.hh TsGSM2parallel.hh
	$(CXX) $(CXXFLAGS) -c TsGetSurvival.cc

TsGSM2parallel.o: TsGSM2parallel.cc TsGSM2parallel.hh TsSpecificEnergy.hh
	$(CXX) $(CXXFLAGS) -c TsGSM2parallel.cc

TsSpecificEnergy.o: TsSpecificEnergy.cc TsSpecificEnergy.hh
	$(CXX) $(CXXFLAGS) -c TsSpecificEnergy.cc

TsLinealEnergy.o: TsLinealEnergy.cc TsLinealEnergy.hh
	$(CXX) $(CXXFLAGS) -c TsLinealEnergy.cc

.PHONY = all clean
all: monas TsGetSurvivalRBEQualityFactor.o TsGSM2parallel.o TsSpecificEnergy.o TsLinealEnergy.o

clean:
	rm -rf *.o
	rm -rf monas 
