appname := gravlen

CXX := g++
INC := -I/home/epeng/Program/lsstsw/stack/Linux64/eigen/3.2.5/include -I/home/epeng/Program/phosim_core/source/cfitsio/include
CXXFLAGS := -Wall -O3 -g -std=c++11 $(INC)
LDFLAGS := -L. -L/home/epeng/Program/phosim_core/source/cfitsio/lib
LDLIBS := -lcfitsio -lgsl -lgslcblas -lm -lfortranstuff

SRCS := main.cpp Image.cpp commons.cpp Model.cpp gl_crit.cpp parafit.cpp
OBJS  := $(patsubst %.cpp, %.o, $(SRCS))

all: $(appname)

$(appname): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(appname) $(OBJS) $(LDLIBS)

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(OBJS)

dist-clean: clean
	rm -f *~ .depend

include .depend
