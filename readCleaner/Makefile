CXX = g++
CXXFLAGS = -O3 -Wall
CXXSRC = $(wildcard *.cpp)
LIBS = -lboost_system -lboost_filesystem
OBJ = $(CXXSRC:.cpp=.o)

all: readCleaner

.PHONY: clean

-include $(OBJ:.o=.d)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $*.cpp $(LIBS)
	$(CXX) -MM $(CXXFLAGS) $*.cpp $(LIBS) > $*.d

readCleaner: $(OBJ)
	$(CXX) $(CXXFLAGS) -o readCleaner *.o $(LIBS)

clean:
	rm -f readCleaner *.o *.d
