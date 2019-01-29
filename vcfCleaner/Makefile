CXX = g++
CXXFLAGS = -O3 -Wall
BIN = vcfCleaner expandrf

all: $(BIN)

vcfCleaner: vcfCleaner.cpp vcfCleaner.h
	$(CXX) $(CXXFLAGS) vcfCleaner.cpp -o vcfCleaner
expandrf: expandrf.cpp
	$(CXX) $(CXXFLAGS) expandrf.cpp -o expandrf

clean:
	rm -f $(BIN) *.o *.d

.PHONY: clean all
