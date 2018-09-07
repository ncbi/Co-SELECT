CXX = g++ -O3
CPP_FILES := $(wildcard *.cpp)
CPP_FILES = \
  OptionParser.cpp \
  Output.cpp  \
  ReadData.cpp \
  SNP.cpp \
  ShapePrediction.cpp \
  TableCompiler.cpp \
  main.cpp \
  map.cpp \
  prediction.cpp \
  properties.cpp \
  utilities.cpp

OBJ_FILES = $(patsubst %.cpp,%.o,$(CPP_FILES))

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

dnashape : $(OBJ_FILES) 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@
