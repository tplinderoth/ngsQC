CXX = g++
CXXFLAGS = -O3 -Wall
BIN = vcfCleaner expandrf vcf2bimbam vcf2sweepfinder vcf2eigenstrat

all: $(BIN)

vcfCleaner: vcfCleaner.cpp vcfCleaner.h
	$(CXX) $(CXXFLAGS) vcfCleaner.cpp -o vcfCleaner -lz -lboost_iostreams
expandrf: expandrf.cpp
	$(CXX) $(CXXFLAGS) expandrf.cpp -o expandrf
vcf2bimbam: vcf2bimbam.cpp
	$(CXX) $(CXXFLAGS) vcf2bimbam.cpp -o vcf2bimbam
vcf2sweepfinder: vcf2sweepfinder.cpp
	$(CXX) $(CXXFLAGS) vcf2sweepfinder.cpp -o vcf2sweepfinder
vcf2eigenstrat: vcf2eigenstrat.cpp
	$(CXX) $(CXXFLAGS) vcf2eigenstrat.cpp -o vcf2eigenstrat

clean:
	rm -f $(BIN) *.o *.d

.PHONY: clean all
