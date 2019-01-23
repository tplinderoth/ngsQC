CXX = g++
CXXFLAGS = -O3 -Wall
CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)
BIN = vcfCleaner

all: $(BIN)

.PHONY: clean

-include $(OBJ:.o=.d)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp > $*.d

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(BIN) *.o

clean:
	rm -f $(BIN) *.o *.d
