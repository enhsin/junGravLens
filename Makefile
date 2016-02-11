NAME=PC
NAME=LAPTOP
#ifeq ($(NAME), PC)
#INC=-I/usr/local/Cellar/boost/1.58.0/include -I/Users/cheng109/toberemoved/phosim/phosim_core/source/cfitsio/include -I/usr/local/include -I/Users/juncheng/work/phosim_core/source/cfitsio/include -I/usr/include/python2.7/
#LIB=-L/usr/local/Cellar/boost/1.58.0/lib -L/usr/local/lib -L/Users/cheng109/toberemoved/phosim/phosim_core/source/cfitsio/lib -L/Users/juncheng/work/phosim_core/source/cfitsio/lib -L/opt/local/lib
#INC_EIGEN=-I/Users/cheng109/eigen-eigen-c58038c56923/
#endif
ifeq ($(NAME), LAPTOP)
INC=-I/usr/local/Cellar/boost/1.58.0/include -I/Users/juncheng/work/cfitsio -I/usr/local/include -I/usr/local/Cellar/gsl/1.16/include
LIB=-L/usr/local/Cellar/boost/1.58.0/lib -L/Users/juncheng/work/cfitsio -L/opt/local/lib -L/usr/local/Cellar/gsl/1.16/lib
INC_EIGEN=-I/Users/juncheng/work/eigen-eigen-b30b87236a1b
endif
CC=g++ #clang++
CFLAGS=-Wall -O3
LDFLAGS=-lcfitsio -lgsl 
# -larmadillo -lboost_iostreams -lboost_system #-fopenmp 
SRCS=main.cpp Image.cpp commons.cpp Model.cpp
#OMP_NUM_THREADS=8


#all: main.o Image.o commons.o Model.o
#	g++ $(INC) $(LIB) $(CFLAG) $(LDFLAGS) commons.cpp main.cpp Image.cpp Model.cpp -o junGL 
#	./junGL

all: #main.cpp Image.o commons.o Model.o 
	
	$(CC) $(CFLAGS) $(INC) $(INC_EIGEN) $(LIB) $(LDFLAGS) commons.cpp main.cpp Image.cpp Model.cpp gl_crit.cpp parafit.cpp -o junGL #-mmacosx-version-min=10.8
	#valgrind --tool=memcheck --leak-check=full --verbose --log-file=memcheck.log --track-origins=yes ./junGL
	./junGL

	
	